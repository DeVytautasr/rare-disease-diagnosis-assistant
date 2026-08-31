# Real Patient Data Validation — Stage 1 Breakpoint Tools

> **Which validation document is this?** This one covers the **real patient
> BAMs**. For the public GIAB HG002 / PacBio-Illumina validation see
> `GIAB_PUBLIC_DATA_VALIDATION.md` (renamed from `REAL_DATA_VALIDATION.md`
> on 2026-08-30 to keep the two apart). The two datasets reach different
> conclusions about threshold calibration; do not read a figure from one as
> if it came from the other.

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

> **Annotation added 2026-08-30.** The findings below are unchanged. Nothing
> measured in the validation run has been edited, re-run in place, or removed —
> those numbers are the evidence base the decisions were made from, and a
> decision cannot be checked against evidence that was rewritten to agree with
> it. What follows is a status column recording what was done about each
> finding, plus a re-measurement section at the end.

| # | Status | Commit |
|---:|---|---|
| 1 | FIXED | `b37b347` |
| 2 | FIXED | `b37b347` |
| 3 | FIXED | `89512e5` |
| 4 | PART FIXED, PART DELIBERATELY NOT | `8364b23, 59281e3` |
| 5 | FIXED | `c0dcfbc` |
| 6 | FIXED | `016fc06` |
| 7 | FIXED | `b7cc665` |
| 8 | FIXED, WITH A DOCUMENTED CAVEAT | `e29fae1` |
| 9 | FIXED | `92d54e3` |
| 10 | FIXED | `f4b8e27` |
| 11 | FIXED | `99bbb71` |
| 12 | FIXED | `b1ea12c` |
| 13 | FIXED | `7da77c7` |
| 14 | FIXED | `ffac4b6, c0dcfbc` |
| crashes | FIXED | `ffac4b6` |

All fixes carry a regression test in `stage1_igv_assistant/tests/`, named for
the condition that exposed them.

> Two further findings, **15 and 16**, were added on 2026-08-31. They did not
> come from this validation run — see the section at the end of this document.

### 1. Losing data raises the score — up to 100/100 "strong" on unusable data
> **Status: FIXED** — `b37b347`. Tri-state assessment plus a `LOW_MAPQ_QUALITY_GATE = 0.4` gate: a window where >40% of reads fail the MAPQ filter now returns `evidence_score: None` / `QUALITY-LIMITED` instead of a normalised score.


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
> **Status: FIXED** — `b37b347`. `assessable` is now false ONLY when the window held no reads at all; reads-present-but-filtered is a distinct `quality_limited` state with its own reason.


The reason string accompanying finding 1 is itself wrong. At an HLA contig the
stats layer reports 4 (SAMPLE_A) and 7 (SAMPLE_B) reads while the other layers
report `"no reads in window"`. In the chr1 pericentromere the same string is
emitted at windows holding 61, 256, 304, 402 and 452 reads.

A reader is told "no data here". The truth is "all coverage here is
multi-mapping" — a completely different statement, and the one that would have
warned them off finding 1.

### 3. `verdict` and `is_balanced` can contradict each other
> **Status: FIXED** — `89512e5`. `verdict` and `is_balanced` are derived from one decision rather than computed twice.


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
> **Status: PART FIXED, PART DELIBERATELY NOT** — `8364b23, 59281e3`. First half fixed: `MIN_ABSOLUTE_SUPPORT = 3` gates the bottom tier of the discordant and split layers. Second half deliberately unchanged: the 0.7 depth threshold was NOT moved, because moving a calibration constant on the strength of one cohort of two would trade a measured false-positive rate for an unmeasured false-negative one. The 26% figure is recorded in the threshold-provenance comment instead. See the re-measurement below.


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
> **Status: FIXED** — `c0dcfbc`. Windows are clamped to the contig length, `clamped_to_contig` and `contig_length` are returned, and a start past the end is an `out_of_range` error rather than clean zeros.


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
> **Status: FIXED** — `016fc06`. `_canonical_chrom` resolves against the BAM header in both directions and returns the name unchanged when it cannot resolve. No contig name is invented.


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
> **Status: FIXED** — `b7cc665`. SA entries are filtered by their own mapQ, `sa_entries_below_min_mapq` reports the drops, strand is kept as `partner_strand_concordant`/`partner_strand_flipped`, and `get_split_reads`'s `min_mapq` default moved from 0 to 20 to match its neighbours. See the re-measurement below — the drop is large.


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
> **Status: FIXED, WITH A DOCUMENTED CAVEAT** — `e29fae1`. `_primary_contig` maps the well-defined `chrN_*` pattern to `chrN` for the discordance decision only; the exact contig is still reported, and same-primary mates are counted in `same_primary_alt_mates`. The 525 `HLA-*` contigs are NOT mapped to chr6 — their names encode no such link, so mapping them would mean hardcoding a guess. Known caveat, asserted in the test so it stays visible.


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
> **Status: FIXED** — `92d54e3`. Partner counts are deduplicated per read, and `sa_entries_total` preserves the raw entry count. A read naming two distinct partners still contributes to each — the docstring now says so.


`split += 1` once per read; `partner_chroms[...] += 1` once per SA *entry*. So
`sum(partner_chromosomes.values()) != split_reads` (619 vs 603; 316 vs 303).
Multi-segment reads are common — 3,807 and 3,282 reads with more than one SA
entry, some with four. A read with two segments on the same contig
double-counts as two pieces of evidence for that partner.

The docstring says `{chr_name: count}` without saying count of what.

### 10. Two entry points to the depth layer disagree
> **Status: FIXED** — `f4b8e27`. `depth_window_bp` and `depth_window_size` are parameters, defaulted to the previous 2000/200 so behaviour is unchanged, and echoed in the result so a caller can reproduce the ratio. The window-size coupling is stated with the threshold constant.


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
> **Status: FIXED** — `99bbb71`. The sentence names per-base depth, matching what the tool computes.


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
> **Status: FIXED** — `b1ea12c`. The per-read SA coordinate is out of the prose; `example_partner_loci` carries it as a structured field.


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
> **Status: FIXED** — `7da77c7`. `clinical_note` states what the coordinate lookup established, not that a breakpoint disrupts a gene.


`bam_tools.py:1842-1843` returns, for **any** queried coordinate:

```python
"clinical_note": f"Breakpoint directly disrupts {len(gene_list)} gene(s)."
```

It emitted this at arbitrary control positions where there is no breakpoint and
no evidence of one. The field name and the verb "disrupts" both assert
something a coordinate lookup cannot know.

### 14. Out-of-range coordinates return clean zeros rather than an error
> **Status: FIXED** — `ffac4b6, c0dcfbc`. Out-of-range coordinates return a structured `out_of_range` error.


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

### Addendum, 2026-08-30: post-fix re-measurement

The paragraph above describes the validation run itself and stays true of it.
The fixes came afterwards, in the commits named in the status column. This
section records what changed when the fixed code was re-run against the same
two BAMs. It is additive: nothing above was re-run in place or edited.

Method: three git worktrees of this repository at `ffac4b6` (before any of
this group of fixes), `8364b23` (plus `MIN_ABSOLUTE_SUPPORT`) and `e29fae1`
(plus findings 7 and 8), imported in turn so that the BAM, the coordinates and
the parameters are identical and the only variable is the code.

#### The 42-locus control grid (finding 4, first half)

The original grid's coordinates were not recorded and its scratch files are
gone, so the identical positions could not be recovered. A fresh grid was
drawn by a rule fixed in advance — seed 20260830, uniform over chr1–chr22
weighted by contig length, 5 Mb excluded at each contig end, 150 candidates
drawn, then accepted in draw order where the ±500 bp read count falls within
0.6–1.4× the sample's own median non-zero window depth, until 42 are held. No
position was chosen for the answer it gave. Both samples selected the same 42
coordinates, as the fixed seed requires; their accepted depth bands differ
(145–337 and 154–358 reads, medians 241 and 256).

`evidence_strength` over the 42 positions:

| | SAMPLE_A before | SAMPLE_A after | SAMPLE_B before | SAMPLE_B after |
|---|---:|---:|---:|---:|
| `none` | 14 (33%) | **30 (71%)** | 14 (33%) | **25 (60%)** |
| `weak` | 28 (67%) | **12 (29%)** | 28 (67%) | **17 (40%)** |
| `moderate` / `strong` | 0 | 0 | 0 | 0 |

Layer firing rate (score > 0), before → after:

| Layer | SAMPLE_A | SAMPLE_B |
|---|---|---|
| discordant | 25/42 (60%) → **7/42 (17%)** | 25/42 (60%) → **9/42 (21%)** |
| split | 2/42 (5%) → **0/42** | 2/42 (5%) → **0/42** |
| soft-clip | 1/42 (2%) → 1/42 (2%) | 2/42 (5%) → 2/42 (5%) |
| depth | 0/42 → 0/42 | 0/42 → 0/42 |

Mean `evidence_score` fell from 8.04 to 4.46 (SAMPLE_A) and 9.46 to 6.25
(SAMPLE_B); the median fell from 7.50 to 0.00 in both. The verdict changed at
16/42 and 11/42 positions, in every case from `weak` to `none`, never the
reverse. The `minsup` and `post` columns are identical at all 42 rows in both
samples — findings 7 and 8 do not move these ordinary autosomal control loci,
which is what should be expected of them.

Two honest qualifications. First, this is a **different grid** from the one
that produced the 32/42 (76%) figure above, so the two "before" numbers —
28/42 (67%) here, 32/42 (76%) there — are not the same measurement. They agree
in size and direction, which is the most that can be claimed. Second, the
soft-clip and depth layers are unchanged by this fix, so their rates are a
control on the method rather than a result: they did not move, and they were
not expected to.

The bottom tier no longer fires on background. It has not been shown to still
fire on real events — this grid contains no known breakpoint, so it can only
measure false positives. That asymmetry is why the depth threshold was left
alone.

#### Finding 7: SA-record MAPQ filtering

At chr1:1,000,000–1,100,000, the window used for the SA-parsing check above:

| Call | SAMPLE_A before → after | SAMPLE_B before → after |
|---|---|---|
| explicit `min_mapq=20` | 293 → **207** (−29.4%) | 101 → **40** (−60.4%) |
| explicit `min_mapq=0` | 603 → 603 (unchanged) | 303 → 303 (unchanged) |
| default, no `min_mapq` | 603 → **207** (−65.7%) | 303 → **40** (−86.8%) |

SA entries dropped for failing their own mapQ: 91/298 and 62/102. The
`min_mapq=0` row is the control — the parameter is real in both directions,
and a caller who asks for no filtering still gets none. The default row moves
most because two changes compound there: the SA filter, and the default itself
moving from 0 to 20.

Strand, previously discarded entirely, at `min_mapq=0`: 447 concordant /
172 flipped (SAMPLE_A) and 144 / 172 (SAMPLE_B). SAMPLE_B's SA entries in this
window are close to evenly split between orientations.

#### Finding 8: alt contigs

At chr6:32,578,000 (MHC), defaults:

| | SAMPLE_A before → after | SAMPLE_B before → after |
|---|---|---|
| `evidence_score` | 15.0 → **7.5** | 47.5 → **22.5** |
| `evidence_strength` | weak → weak | **moderate → weak** |
| discordant pairs | 16/296 → 15/296 | 24/219 → **11/219** |
| `same_primary_alt_mates` | — → 1 | — → **13** |
| split reads | 10 → 0 | 80 → **2** |

SAMPLE_B's split partners before the fix were `chr6_GL000256v2_alt` (47),
`chr6_GL000253v2_alt` (15), `chr6_GL000252v2_alt` (10), `chr6_GL000254v2_alt`
(8) — the observation text called that *"scattered across 4 chromosomes — no
dominant partner"* for what is chr6 and three of its own alt haplotypes. That
sentence is gone. Most of the split collapse here is finding 7's MAPQ filter
rather than finding 8; the two ship together and were measured together.

On `chr15_KI270905v1_alt` the result is weaker than expected, and it is
recorded as measured rather than as hoped. The contig midpoint has **zero
coverage** in both samples — `NOT ASSESSABLE`, not a negative control. Scanning
all 517 windows at 10 kb stride found exactly **one** window with ≥100 reads,
at position 1,840,500, in both samples. There:

| | SAMPLE_A | SAMPLE_B |
|---|---|---|
| discordant before → after (MAPQ 60) | 47/147 → 43/147 | 39/175 → 38/175 |
| fraction | 0.320 → 0.293 | 0.223 → 0.217 |
| `same_primary_alt_mates` | 4 | 1 |
| `evidence_score` | 40.0 → **40.0** | 40.0 → **40.0** |

SAMPLE_B's 0.223 matches the figure reported above exactly, so this is very
likely the same window the original run used for that sample; SAMPLE_A's 0.320
does not match the 0.172 reported above, so its original position was a
different one and could not be recovered.

**The fix does not rescue this locus, and the finding above overstates the
cause.** The top mate partners at this window are chr14, chr7 and chr13 — not
chr15. Only 4 and 1 mates were on chr15 or its own scaffolds. The 40.0
"moderate" score survives the fix because it is driven by mates on genuinely
different chromosomes, which is what mismapping onto a 5 Mb alt haplotype
produces. Same-primary alt pairing was real but small here. Finding 8 is
correct that alt contigs were miscounted; it is not correct that this was what
made the alt contig score "moderate". Alt-contig loci remain a known weak spot,
now for a reason the fix was never going to address.


---

## Findings 15 and 16, 2026-08-31 — found by an LLM session, not by this run

> **Provenance.** Findings 1–14 above came from the 2026-08-29 validation run,
> in which tools were called directly from Python. These two did not. They were
> found on 2026-08-31 by the first LLM session to use the MCP tools after the
> fourteen fixes, against the same two BAMs
> (`LLM_SESSION_5_PATIENT_DATA_qwen.md`). They are numbered into this document
> because they are defects in the same tool surface, reported against the code
> state `e29fae1` that the fourteen fixes produced — but they were **not**
> measured by the run above, and the ranking there does not include them.

| # | Status | Fixed in | Present in `e29fae1` | Cause |
|---:|---|---|---|---|
| 15 | FIXED | `96b5e25` | yes | latent before the fixes, **amplified** by `8364b23` |
| 16 | FIXED | `96b5e25` | yes | introduced by `8364b23` |

Both carry regression tests in `tests/test_subthreshold_observations.py`
(35 assertions).

### 15. A score of zero is reported as an empty window, directly beneath the counts that contradict it

`supporting_observations` ended with a blanket denial:

> No discordant pairs, soft-clipping, split reads, or depth changes detected
> near this position.

emitted whenever `evidence_score == 0` and the window held any reads. The
session saw it printed in the same list as, and immediately after, a sentence
reporting one discordant pair with a named mate chromosome. Two contradicting
sentences, with nothing to tell a reader which is authoritative — and the
denial is the one a model summarising the list is most likely to carry into a
report, because it is the last and the most general.

**The obvious explanation is wrong, and the grid says so.** It is natural to
blame `MIN_ABSOLUTE_SUPPORT` (`8364b23`): once 1–2 supporting reads score
nothing, a zero score stops implying an empty window. But the contradiction
fires on **7/42 and 8/42** control loci in the *pre-fix* code as well, because
two other counted-but-unscored cases already existed — soft clips below the
pileup threshold, and depth dips ruled off-position. `MIN_ABSOLUTE_SUPPORT`
took it from 7 and 8 to **22/42 and 18/42**, more than half of all loci in
`SAMPLE_A`. It tripled a defect it did not create.

That distinction determined the fix. Gating the denial on `MIN_ABSOLUTE_SUPPORT`
specifically would have left the soft-clip and depth cases firing. It is
instead gated on the **actual absence of all four layers**, so any future
layer-level gate inherits the correct behaviour. Where reads exist but sit
under a threshold they are named, with the bar they missed:

> No layer reached its scoring threshold. Support is present but sub-threshold:
> 1 discordant pair(s), below the 3-read minimum for this 500bp window.

Reads discarded by `min_mapq` count as presence for this purpose too: "we
filtered them out" is a different claim from "there was nothing there", and the
per-layer quality note already sat one line above the denial that contradicted
it.

Measured after the fix: **0/42 in both samples, at both window widths.** The
blanket sentence still fires verbatim where it is true, which is asserted
separately so the fix cannot degenerate into never denying anything.

*Caveat on the count.* The detector flags a locus when the denial co-occurs
with a non-zero discordant, split or soft-clip count. It does not count
off-position depth dips, so 7/42 and 8/42 are **lower bounds** on the pre-fix
rate.

### 16. The minimum-support threshold is not window-invariant, so the verdict moves with an argument

`MIN_ABSOLUTE_SUPPORT = 3` was deliberately an absolute count rather than a
fraction — a fraction cannot distinguish 1-in-250 from 1-in-4. But a bare count
is not invariant to the window it is counted in. The number of background reads
in a ±`window_bp` window scales with `window_bp`, so 3 is a strict filter at
500 bp and a loose one at 1500 bp.

The session hit this directly: one position read `none` with 1 discordant pair
at the default 500 bp and `weak` with 4 pairs at `window_bp=1000`. The verdict
changed because an argument changed, not because the evidence did — and
`window_bp` is a parameter the model chooses.

**Decision: scale the threshold, *and* report the window.** Not either/or.

Scaling is the substantive fix, and the shape is derivable rather than guessed:
at uniform coverage the expected background count is proportional to the number
of reads in the window, which is proportional to its width. The threshold is
now stated as 3 per `MIN_ABSOLUTE_SUPPORT_WINDOW_BP = 500` and scaled to the
window in use, floored at 2 so that no window is narrow enough to let one read
plus its duplicate score. But scaling alone would have *hidden* the coupling
rather than fixed it — the mistake the depth threshold was explicitly not
allowed to make — so `window_bp` and `min_supporting_reads` are now returned
fields and appear in `interpretation_template`, exactly as `depth_window_bp`
and `depth_window_size` already were, and for the same reason: a verdict quoted
without its window is not a reproducible statement.

#### The 42-locus grid, re-run at two window widths

Same selection rule as the 2026-08-30 addendum (seed 20260830, 5 Mb excluded at
each contig end, 150 candidates, accepted in draw order within 0.6–1.4× the
sample's own median candidate depth, until 42 held). The acceptance band always
uses ±500 bp, so **both window widths are evaluated at the same 42
coordinates**; the only variables are the code state and `window_bp`.

Three code states: `before` = `ffac4b6` (no minimum at all), `fixed3` =
`e29fae1` (a bare count of 3), `scaled` = `96b5e25`.

The `scaled` measurement was in fact taken against `52421ba`, the pre-amend
form of that commit. Its message was later corrected — the pre-fix
contradiction rate below disproved a causal claim it made — and two docstrings
were reworded in the same amend. The two trees are identical under an AST
comparison with docstrings stripped, so no measured number is affected; it is
noted because the grid was run before the SHA it is attributed to existed.

`evidence_strength`, count of `weak` out of 42 — the false-positive measure,
since this grid contains no known breakpoint:

| | SAMPLE_A @500 | SAMPLE_A @1000 | SAMPLE_B @500 | SAMPLE_B @1000 |
|---|---:|---:|---:|---:|
| `before` | 26 (62%) | 34 (81%) | 25 (60%) | 34 (81%) |
| `fixed3` | 11 (26%) | 18 (43%) | 15 (36%) | 21 (50%) |
| `scaled` | 11 (26%) | **7 (17%)** | 15 (36%) | **11 (26%)** |

**The improvement does not survive the window change under a bare count.**
Widening to 1000 bp gives back roughly half the gain — `SAMPLE_A` 26% → 43%,
`SAMPLE_B` 36% → 50%. Under scaling it survives and strengthens.

Discordant-layer firing rate (component score > 0), out of 42:

| | SAMPLE_A @500 | SAMPLE_A @1000 | SAMPLE_B @500 | SAMPLE_B @1000 |
|---|---:|---:|---:|---:|
| `before` | 23 | 32 | 22 | 33 |
| `fixed3` | 7 | **14** | 8 | **17** |
| `scaled` | 7 | 2 | 8 | 2 |

The `fixed3` row is the measurement that justifies the linear form: **the
firing rate doubles exactly when the window doubles** — 7→14 and 8→17 — which
is what proportional background predicts and what a window-invariant threshold
must cancel. This is no longer only an argument.

Mean `evidence_score` over the 39 scoreable loci (3 are `QUALITY-LIMITED` in
every state and every window, and are excluded):

| | SAMPLE_A @500 | SAMPLE_A @1000 | SAMPLE_B @500 | SAMPLE_B @1000 |
|---|---:|---:|---:|---:|
| `before` | 7.88 | 9.62 | 9.23 | 11.35 |
| `fixed3` | 4.42 | 5.77 | 6.15 | 7.88 |
| `scaled` | 4.42 | 3.46 | 6.15 | 5.00 |

Finding 15's contradiction rate over the same grid, out of 42:

| | SAMPLE_A @500 | SAMPLE_A @1000 | SAMPLE_B @500 | SAMPLE_B @1000 |
|---|---:|---:|---:|---:|
| `before` | 7 | 1 | 8 | 2 |
| `fixed3` | 22 | 17 | 18 | 15 |
| `scaled` | **0** | **0** | **0** | **0** |

The `before` row falls at the wider window only because more loci score
non-zero there and never reach the denial branch at all — it is not an
improvement.

**At the default width the scaling is exactly a no-op.** Every *scoring*
column of `fixed3@500` and `scaled@500` is identical, in both samples: the same
11 and 15 `weak` verdicts, the same 7 and 8 discordant firings, the same means
and medians. Every figure in the 2026-08-30 addendum still stands unaltered.
The one column that does move at 500 bp is the contradiction count (22→0,
18→0), which is finding 15's fix, not finding 16's — the two ship in one commit
and must not be read as one effect.

#### Two honest qualifications

**This is a third grid, not the addendum's.** The selection rule was re-derived
from its written description rather than from the original script, which is
gone. It draws a slightly different sample: median candidate depth 242 and 258
here against 241 and 256 there, and this run reports 3 `QUALITY-LIMITED` loci
per sample which the addendum's table did not separate out. The `before@500`
`weak` rates — 62% and 60% here, 67% and 67% there — agree in size and
direction. That is the most that can be claimed, and it is the same caveat the
addendum itself carried against the original 76% figure.

**Linear scaling over-corrects; it is conservative, not invariant.** True
verdict-invariance would hold the firing rate roughly constant across widths.
It does not: `scaled@1000` fires on 2/42 in both samples, against 7 and 8 at
the default. The reason is structural — a threshold that grows like the mean
outruns a count distribution whose spread grows like the root of the mean, so
the same nominal cutoff sits further into the tail at the wider window. A
variance-aware rule (scaling like mean + k·√mean) would be closer to invariant,
and was **not** adopted, because choosing *k* would be a second uncalibrated
judgement stacked on the first.

The cost of over-correcting is a loss of sensitivity at wide windows, and this
grid **cannot price it**: it contains no true positive, so it can only measure
false positives, and every number above is a false-positive number. That is the
same asymmetry that stopped the 0.7 depth threshold from being moved. The
difference here is that the alternative is not "leave it alone" but a bare
count that is demonstrably wrong in the *unsafe* direction — it inflates
evidence as the window widens. Erring conservative at non-default widths, with
the default untouched and the coupling now printed alongside every verdict, is
the position taken. It is a judgement call and is recorded as one.
