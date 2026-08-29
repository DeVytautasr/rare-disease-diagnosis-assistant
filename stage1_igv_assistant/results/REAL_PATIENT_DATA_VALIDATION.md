# Real Patient Data Validation — Stage 1 Breakpoint Tools

**Date:** 2026-08-29
**Data:** two aligned WGS BAMs, ~39 GB and ~41.6 GB, referred to throughout as
SAMPLE_A and SAMPLE_B. Real sample identifiers are deliberately absent from
this document; the mapping is held outside the repository.

---

## Scope, and what this document does not claim

**All tested positions were arbitrary code-exercise coordinates, chosen to
exercise code paths — not candidate breakpoints.** They were selected for
ordinary coverage and for region class (centromeric, gene-dense, contig edge,
and so on), never for signal. Several were rejected and replaced precisely
*because* they showed anomalous coverage.

**No clinical interpretation is implied by anything in this document.** Where
an evidence score, a depth ratio, a discordant count, or a gene annotation
appears below, it is reported as the output of a function under test. None of
it is a finding about either sample, about any variant, or about any
individual. The X/Y depth ratios in Task 3 are reported as a property of the
coverage distribution relevant to a scoring threshold, and carry no
phenotypic or sex inference.

This was a validation run. **No code was changed.** No duplicate filtering was
added; that decision was deliberately deferred until after these results.

### Method

Tools were called directly from `stage1_igv_assistant.tools.bam_tools` via
`python3` against the real BAMs. Every call was wrapped so that a failure
produced a traceback rather than aborting the matrix. Scratch files were
written to absolute `/tmp` paths. Reference is hs38DH: 3,366 contigs, 26
primary, 261 `*_alt`, 42 `*_random`, 2,512 `chrUn_*`, and 525 `HLA-*` contigs
which carry **no** `chr` prefix. Both BAMs are coordinate-sorted Illumina,
`bwa mem -M`, and pass `samtools quickcheck -v`.

---

## Headline

**The machinery is sound. The interpretation surface is not.**

Across 70 primary calls in Task 1, 16 probe groups in Task 2, and the full
Task 3 edge matrix, there were **zero exceptions on real coordinates**, zero
malformed returns, zero missing keys, no memory pressure (peak RSS 68 MB, 0
major page faults), and no BAM-backed call slower than 0.096 s. Contig
resolution works in both naming conventions. Zero-coverage handling is exactly
correct. SA-tag parsing is exact against an independent parser.

Every serious defect found is in what the numbers *mean*, not in whether the
code runs. The most severe is that **the headline evidence score can go up
when data quality goes down** — reaching a maximum of 100/100 "strong" in
windows where essentially all reads are unmappable.

---

## Findings, ranked by how badly they could mislead

### 1. Losing data raises the score — up to 100/100 "strong" on unusable data

`summarize_breakpoint_evidence` excludes unassessable layers from **both the
numerator and the denominator** of `evidence_score` (documented at
`bam_tools.py:274-275`). Separately, `get_read_depth_profile` takes **no
`min_mapq` parameter at all**, while the discordant, soft-clip and split
layers all filter at `min_mapq=20`.

In a pericentromeric window where nearly every read is MAPQ 0, the three
filtered layers see nothing and drop out; the depth layer, which never
filtered, still scores. The denominator falls to 1, and a 25/25 depth score
normalises to 100.

Observed at chr1:122,040,000, chr1:122,240,000 and chr1:122,260,000 in **both**
BAMs, verbatim:

> Evidence strength: strong. Score: 100.0/100 normalised over 1 applicable
> layer(s) (1/1 showing signal); raw score over all 4 layers was 25.0/100.

Of 28 (locus, BAM) cells meeting the all-MAPQ-0 criterion, **16 scored ≥ 70
("strong")**.

The MCP server instructions direct the assistant to report `evidence_score`
and `evidence_strength` as the primary finding, on the grounds that the raw
score is "artificially capped by layers this data could never produce". That
reasoning is sound for a genuine technology limit — a Novoalign BAM with no SA
tags — but it is exactly inverted for a *quality* limit, where the layers
could have produced signal and were silenced by bad mapping. The normalisation
cannot currently tell those two cases apart.

The same mechanism fires from a parameter alone: `window_bp=0` at a clean 32x
MAPQ-60 locus returns `SCORE=60.0 "moderate"` from `raw=15.0`.

### 2. "No reads in window" is reported when reads exist but were MAPQ-filtered

The reason string accompanying finding 1 is itself wrong. At an HLA contig the
stats layer reports 4 (SAMPLE_A) and 7 (SAMPLE_B) reads while the other layers
report `"no reads in window"`. In the chr1 pericentromere the same string is
emitted at windows holding 61, 256, 304, 402 and 452 reads.

A reader is told "no data here". The truth is "all coverage here is
multi-mapping" — a completely different statement, and the one that would have
warned them off finding 1.

### 3. `verdict` and `is_balanced` can contradict each other

In `check_reciprocal_breakpoint`, every `verdict` branch requires
`primary_disc >= 5` (`bam_tools.py:1954-1958`), while `is_balanced` is
`primary_disc >= 3 and reciprocal_disc >= 3` (line 1978). Anything in the 3–4
band on both sides returns:

```
verdict     : INSUFFICIENT EVIDENCE at both positions
is_balanced : True
```

Reproduced on real data. Two callers reading the same return object get
opposite answers, and both fields look authoritative.

Separately, that `else` branch message asserts "at both positions" even when
one side has signal — it fired on a case with 3 discordant pairs and 1
back-pointing read at the reciprocal position.

### 4. The depth layer's false-positive rate at ~31x is not small

Across 120 pseudo-random autosomal loci per BAM (fixed seed, 117 scored):

| Measure | SAMPLE_A | SAMPLE_B |
|---|---|---|
| depth layer scored > 0 | 30/117 (26%) | 30/117 (26%) |
| any `evidence_score > 0` | 69% | 79% |
| median `depth_ratio_min_to_mean` | 0.739 | 0.753 |

The sentence *"Read depth drops to N% of the window mean near the queried
position … consistent with a possible deletion"* was emitted at **26% of
ordinary loci chosen at random**.

The 0.7 scoring boundary sits *inside* the noise distribution rather than
outside it. Its empirical calibration came from 300x GIAB data, where the same
statistic has roughly 3x less relative scatter. On chrX, at roughly half
autosomal depth, the depth layer fired at 56–58% of loci.

Task 1's independent 42-locus control grid agrees: **32/42 (76%) of arbitrary
ordinary-coverage positions returned `evidence_strength: "weak"`** rather than
`"none"`, because both the discordant and split layers award 7.5/25 for *any*
non-zero fraction (`bam_tools.py:1411-1412`, `1484-1485`). At 30x WGS a ±500 bp
window holds ~250 reads, and one read with an interchromosomal mate is
ordinary background.

The module's own provenance comment already states these cutoffs are heuristic
and were "not checked against true-negative/normal-coverage regions for a
false-positive rate". This run supplies that missing number.

The soft-clip layer, scored on `max_clips_at_position` rather than a fraction,
was clean at **0/42**. That layer's redesign is holding up on real data.

### 5. Windows overrunning a contig end are not clamped, and the full width is reported

`get_bam_stats_at_locus(chrM, 16000, 17500)` returns the same 11,840 reads as
`16000-16569`, but `mean_depth` 1,011.06 instead of 2,665.36 — exactly
×569/1500. The denominator includes 931 bases that do not exist, and the `end`
field echoes the requested 17,500. At `16000-26569` the reported depth is
143.49.

`get_read_depth_profile` emits bins entirely past the contig end: 9 of 20 at
chrM:16,400, 20 of 20 at chr21:47,709,983, and 9,918 of 10,000 for
`chrM 0-2,000,000`, where mean depth reads 21.45 against a true ~2,589.

This manufactures a maximum depth score. At chrM:16,400 in both BAMs:
`evidence_score 65.0`, `moderate`, `4/4 layers`, `depth_score 25.0`, with
*"Read depth drops to 0% of the window mean … consistent with a possible
deletion."* The 0.0 bins are the ones past base 16,569. `dip_is_at_focus` is
True *because* the query is near the edge, so the localisation gate that exists
to suppress off-position dips reinforces this artefact instead.

Related: `consensus_clip_position` was returned as **16,569 on a 16,569 bp
contig** (last valid 0-based base is 16,568), with 1,113 supporting reads.

### 6. Fabricated contig names for HLA partners

`_normalize_chrom` prepends `chr` unconditionally. Every contig class in this
reference is already prefixed *except* the 525 `HLA-*` contigs, so:

```
HLA-A*01:01:01:01  ->  chrHLA-A*01:01:01:01   (in header: True -> False)
```

In an MHC window, **16 of 44 partner keys named contigs that do not exist in
the header**, carrying 102 SA entries. No error, no warning. Broad scan volume:
3,428 and 3,746 HLA-targeting SA entries across 151 and 213 distinct HLA
contigs.

`example_partner_loci` compounds it: the display string is
`f"{first_rname}:{first_pos}"`, and HLA contig names contain colons, yielding
`'HLA-DRB1*12:17:4992'`. Only `rsplit(":", 1)` can recover that. The same
partner therefore appears under two different broken spellings in one result.

### 7. `min_mapq` never reaches the SA record; strand is discarded

`get_split_reads` reads only `fields[0]` and `fields[1]` of each SA entry.
Strand, CIGAR, mapQ and NM are parsed past and dropped.

`min_mapq` filters the **primary** alignment only. In the chr1 window, 144/619
(SAMPLE_A) and 224/316 (SAMPLE_B) SA records have mapQ 0 — **71% for
SAMPLE_B** — and those partner contigs are counted at full weight regardless of
the caller's `min_mapq`. Evidence that looks filtered is not.

Strand is dropped even though 12.3% and 10.1% of SA entries flip strand
relative to the primary. Orientation is the field distinguishing an
inversion-type junction from a direct one.

### 8. Alt contigs are treated as separate chromosomes from their own primary

On `chr15_KI270905v1_alt` the largest single mate partner is chr15 itself,
counted as a discordant inter-chromosomal pair. At MAPQ 60 with no filtering
to suppress it, this yields discordant fractions of 0.172 and 0.223 and
composite scores of 32.5 "weak" and **40.0 "moderate"**.

At chr6:32,578,000 (MHC), SAMPLE_B returns split fraction 0.63 with partners on
chr6 alt contigs → 25/25 split score, composite **47.5 "moderate"**, and the
observation text reads *"scattered across 4 chromosomes — no dominant
partner"* for what is one chromosome plus its own alt haplotypes.

On an HLA contig with `min_mapq=0` the discordant fraction is **1.000** (all
mates on chr6) — the top tier, generated entirely by routine alt-aware mapping.

### 9. `split_reads` and `partner_chromosomes` are in different units

`split += 1` once per read; `partner_chroms[...] += 1` once per SA *entry*. So
`sum(partner_chromosomes.values()) != split_reads` (619 vs 603; 316 vs 303).
Multi-segment reads are common — 3,807 and 3,282 reads with more than one SA
entry, some with four. A read with two segments on the same contig
double-counts as two pieces of evidence for that partner.

The docstring says `{chr_name: count}` without saying count of what.

### 10. Two entry points to the depth layer disagree

`summarize_breakpoint_evidence` hardcodes its own depth geometry
(`bam_tools.py:1358-1361`), ignoring `window_bp`:

```python
depth_profile = get_read_depth_profile(
    bam_path, chromosome, max(0, position - 2000), position + 2000,
    window_size=200, focus_position=position
)
```

The documented tool order has the operator call `get_read_depth_profile`
standalone first. At a ±5000/500 bp geometry the two verdicts **conflict at 4
of 10 positions** — standalone `likely_deletion: False`, `summarize` `True`
with `depth_score: 15.0`. `dip_is_at_focus` conflicts at 3 of 10, in both
directions.

`depth_ratio_min_to_mean` is strongly window-size dependent and the 0.7
threshold was calibrated at `window_size=200`. Holding the span fixed at ±2000:
ratio 0.523 at ws=100, 0.585 at ws=200, 0.830 at ws=500, 0.868 at ws=1000 —
crossing the threshold between 200 and 500. Any threshold quoted for this layer
is meaningless without its window size, and nothing enforces or warns about the
coupling.

Relatedly, `likely_deletion` is essentially always True for wide scans: True at
4 kb/20 bins, 4 kb/80, 20 kb/100, 100 kb/500 and 400 kb/2,000 bins at one
ordinary locus, and True for whole-chrM profiles at both `window_size=100` and
`window_size=1`.

### 11. A stale unit label in the sentence most likely to be quoted

`bam_tools.py:1556`:

```python
f"vs mean {depth_profile['summary']['mean_depth']} reads/window) — "
```

Since the 2026-08-11 FIX-1, `min_depth` and `mean_depth` are true **per-base**
depth. The function's own docstring says so: *"Computes true mean per-base read
depth in sliding windows"*. This is the last `reads/window` in the code, and it
sits in the one sentence an assistant is most likely to paste verbatim into a
report. The adjacent off-position branch prints the same quantities with no
unit at all, so the two branches are also inconsistent with each other.

### 12. A per-read SA coordinate is embedded in prose output

`bam_tools.py:1498` interpolates a concrete partner locus taken from a single
read's SA tag into the `supporting_observations` sentence:

```python
example_suffix = f" (e.g. {examples[0]})" if examples else ""
```

`example_partner_loci` as a structured field is trivial to strip before a
report is produced. Inside a sentence it is not. On patient data this is a
per-read position of exactly the kind that should not leave the analysis, and
it arrives pre-embedded in the field a reporting assistant is most likely to
copy wholesale. Flagged as a data-handling consideration for anything built on
top of these tools.

### 13. `clinical_note` asserts a breakpoint from a coordinate lookup

`bam_tools.py:1842-1843` returns, for **any** queried coordinate:

```python
"clinical_note": f"Breakpoint directly disrupts {len(gene_list)} gene(s)."
```

It emitted this at arbitrary control positions where there is no breakpoint and
no evidence of one. The field name and the verb "disrupts" both assert
something a coordinate lookup cannot know.

### 14. Out-of-range coordinates return clean zeros rather than an error

`stats(chr21, 46,709,983, ...)` past the contig end returns
`reads=0, mean_depth 0.0` with the requested span echoed and no out-of-range
signal. pysam did not raise for any out-of-range fetch tested, so
`_fetch_or_error`'s `invalid_region` path never fired; the only `invalid_region`
seen came from an unknown contig name. **A typo'd coordinate is
indistinguishable from a genuine coverage gap.**

Degenerate parameters likewise return plausible "no data" rather than errors:
`window_bp=0` → "no reads in window" at a clean 32x locus; `window_size` larger
than the span → a single bin with ratio 1.000, which can never show a dip — a
structural false negative. Float coordinates are accepted silently.

---

## Crashes

Three, all parameter-domain — none occurred at any real coordinate. All three
reproduced independently on synthetic data:

| Input | Exception | Site |
|---|---|---|
| `get_read_depth_profile(..., window_size=0)` | `ValueError: range() arg 3 must not be zero` | `bam_tools.py:999` |
| `get_read_depth_profile(..., window_size=-100)` | `IndexError: list index out of range` | `bam_tools.py:1008` |
| `get_bam_stats_at_locus(..., start="...", end="...")` | `TypeError: '<' not supported between instances of 'str' and 'int'` | `bam_tools.py:343` |

The `window_size=-100` case is data-dependent: over a zero-coverage region the
same call returns a clean `assessable: false` result, so it only crashes where
reads exist. The string-coordinate case bypasses `_validate_range`'s contract —
every other malformed input returns a structured error dict.

---

## What worked

Worth recording as explicitly as the defects.

- **Zero-coverage handling is exactly correct.** No division by zero, no 0/0.
  Fractions are `None`, `assessable` is `false` with a reason, `evidence_score`
  is `None`, `evidence_score_raw` is 0, `signal_layers` is `0/0`, strength is
  `NOT ASSESSABLE`. Zero coverage is never reported as a negative result.
- **SA-tag parsing is exact.** Against an independent regex parser over
  chr1:1,000,000-1,100,000: 603/603 and 303/303 reads, 0 malformed entries,
  `partner_chromosomes` matching key-for-key. No parse errors across ~470k SA
  entries scanned.
- **The split-read layer does not depend on supplementary flags.** `bwa mem -M`
  produces **zero** `0x800` reads — confirmed at chromosome scale over ~9-10 M
  reads per chromosome per BAM — and `get_split_reads` still returns non-zero
  throughout, because it keys on the SA tag of the primary record.
  Incidentally, `SA_any − SA_on_primary` equals the secondary count exactly in
  all ten region-samples, so the `is_secondary` skip is doing real work;
  without it the counts would roughly double.
- **Partner ordering is deterministic** — count-descending, ties by file order,
  stable across repeated calls.
- **Contig resolution works in both conventions**, and on `_alt`, `_random`,
  `chrUn_` and exact-match `HLA-*` names. `get_bam_stats_at_locus(bam, "1", …)`
  returns identical numbers to `"chr1"`.
- **Error paths return structured dicts**, not exceptions: unknown contig →
  `invalid_region` with a header sample; inverted range and negative start →
  `invalid_parameters`; an invalid layer name → `invalid_parameters` listing
  valid values.
- **`detect_applicable_layers` is correct on both BAMs** — all four layers
  applicable, from paired reads and SA tags present, matching the header.
  Because all four are applicable, `evidence_score == evidence_score_raw` at
  every Task 1 position; the normalisation is a no-op on this data.

### Performance

No BAM-backed call exceeded 0.1 s. Peak RSS 68 MB with 0 major page faults —
file size is irrelevant to memory, since htslib seeks via the BAI and reads
only the needed BGZF blocks.

| Call | Median | Max |
|---|---|---|
| `get_bam_stats_at_locus` | 0.018 s | 0.019 s |
| `count_discordant_pairs` | 0.016 s | 0.020 s |
| `count_soft_clipped_reads` | 0.016 s | 0.018 s |
| `get_split_reads` | 0.017 s | 0.017 s |
| `get_read_depth_profile` (10 kb, 500 bp bins) | 0.035 s | 0.044 s |
| `summarize_breakpoint_evidence` | 0.090 s | 0.096 s |
| `detect_applicable_layers` | 0.022 s | 0.025 s |
| `get_gene_at_locus` (network) | 1.276 s | 3.584 s |

One scale-dependent behaviour a small fixture cannot show: opening the BAM and
loading the 10.7 MB, 3,366-contig `.bai` costs **14.2 ms**, while the actual
500 bp fetch on an open handle costs **0.7 ms**. Roughly 95% of every call is
index reloading. Each tool opens and closes independently, and
`summarize_breakpoint_evidence` performs **five separate opens**, which is why
it costs 0.090 s against ~0.017 s for a single-layer call. This is not a
correctness problem and is fine for tens or hundreds of loci, but a genome-wide
scan would be dominated by index reloading rather than by reading data.

### Network

All `get_gene_at_locus` calls succeeded — no timeouts, no HTTP errors. Latency
0.171–3.584 s, two orders of magnitude slower than anything touching the BAMs
and the only source of variance in the run. Two observability notes rather than
faults: the function returns no attempt count on success and its backoff starts
at 1 s, so a slow call including a silent retry is indistinguishable from a
single slow request; and one returned gene had no `external_name`, so
`g.get("external_name", "unknown")` put the literal string `"unknown"` into
`gene_name` rather than falling back to the Ensembl ID that was present.

---

## Status

Reported only. **No code was changed**, no fixes were applied, and no duplicate
filtering was added. The duplicate-filtering decision, and any response to the
findings above, follow from these results rather than preceding them.
