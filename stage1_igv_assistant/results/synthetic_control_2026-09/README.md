# Synthetic positive control, rebuilt 2026-09-25 — stopped at the IMP01 gate

Two dated entries follow: the pre-registered gate, which failed and stopped the
run, and — added later the same day as a separate entry — the revised gate, after
which the run continued through delly. The first entry is kept as written.

A new experiment, not a reproduction of Phase 6: new random draws, an
ALT-aware index, and partly new breakpoints. The Phase 6 generator was never
committed and was lost in the reinstall; the rebuild follows the approved
Phase 6 design and lives in `scripts/synthetic/`, committed before each step
ran. Large data (BAM, FASTQ, ART SAM) stays under `~/public_data/sim/`.

## Status

**Stopped at the pre-registered gate on IMP01.** Criterion (a) failed:

| Criterion | Required | IMP01 | Phase 6 |
|---|---|---|---|
| (a) truly spanning reads with an SA tag to the partner chromosome | ≥ 70% | **18 / 33 = 54.5%** | 50 / 61 = 82% |
| (b) SA tags pointing at a wrong partner | 0 | 0 | 0 |
| (c) depth change at chr20:200,000 / chr21:14,100,000 | within ±10% | −7.3% / +1.6% | −0.8% / +5.3% |

Nothing downstream of the gate was run. One implant BAM exists (IMP01); the
other eleven were prepared and aligned but not built, delly was not run, and
the demo bundle was not rebuilt.

What the per-read record (`IMP01_gate.json`) shows: every read that carries a
partner SA tag overlaps the far side of its junction by at least 31 bp; every
read without one overlaps it by at most 30 bp and carries no SA tag at all.
The committed gate defines a truly spanning read as one with at least one base
on each side of the junction. The Phase 6 record that set the 70% threshold
does not say how Phase 6 defined it.

## Files

| File | What it holds |
|---|---|
| `run_record.json` | status, code, tool versions (as invoked), every seed, class thresholds, the recorded-coordinate checks, the gate result, per-implant compensation counts |
| `background_measurements.json` | base-quality bins, fragment length (mean 442.64, SD 104.58), isolated-mismatch error rate (0.002104) |
| `art_profiles.json` | the five 150 bp ART profiles measured the same way; HS25 closest (0.001905) |
| `scan_candidates.json` | the 978 scanned candidate breakends and their statistics |
| `selection.json` | class thresholds, pool sizes, the 12 implants with their sources and seeds |
| `pooled_alignment.json` | the single bwa call for all twelve implants |
| `IMP01_gate.json` | the gate, with every truly spanning read |

## 2026-09-25, later the same day — revised criterion (a); continued through delly

A new entry. The pre-registered result above stands as recorded, in
`IMP01_gate.json` and `run_record.json`.

**Why criterion (a) changed** (the user's decision, taken after seeing the
result). The 70% threshold came from the original IMP01 gate (50/61 = 82%),
whose definition of a truly spanning read was lost with the generator. Under the
definition kept here — at least one base on each side of the junction — a 150 bp
read crossing a junction at a uniformly random position has a shorter side of
30 bp or more, the least that can reach bwa mem's `-T 30`, at 91 of 149
positions: 61.1%. Counting the errors ART actually placed, 239 of the 390 truly
spanning reads across all twelve implants (61.3%) have a shorter side that still
reaches score 30 with a 19 bp exact seed (`-k 19`). The user also quotes
183/304 = 60.2% from the original experiment; no surviving record holds that
figure, so it could not be checked.

**Replacement criterion, defined after the result and recorded as such:**

| Criterion | Result |
|---|---|
| (a1) every truly spanning read without a partner SA tag is explained by `-T 30` / `-k 19` | 15 such reads, **0 unexplained**: 14 have a shorter side of 4–28 bp; the one with 30 bp carries an ART mismatch (best score 25, longest exact match 17 bp) |
| (a2) observed SA fraction against expectation — not pass/fail | 18/33 = 54.5%, against 61.1% geometric (P(≤ 18) = 0.27), 61.3% with ART's errors, 60.2% quoted |
| (b), unchanged | 0 SA tags to a wrong partner |
| (c), unchanged | depth −7.3% at chr20:200,000 and +1.6% at chr21:14,100,000 (without duplicates −4.6% and +3.7%) |

Revised gate: **PASS** (`IMP01_gate_revised_2026-09-25.json`, with every read).

**Dose.** From the pre-implant depth, the fragment length and the read length,
IMP01 should receive about 50.0 junction fragments and show about 17 truly
spanning reads per junction; it received 52 fragments and shows 15 at J20 and
18 at J21. The original gate's 61 spanning reads would need about 91 fragments
at this depth. The revised gate record holds the same table for all twelve
implants.

**Then continued.** The other eleven implants were built and delly v2.6.0
(`sr -h 1`, the unmodified exclude template) was run on all twelve and on the
untouched background. Every build and every delly job exited 0; every BCF names
one sample, NA12878; peak RSS 315,948–317,532 kB and wall time 0:43–1:10 per run.
Records per BCF: background 916, IMP01–IMP08 918, IMP09–IMP12 916 — reported, not
interpreted. The evidence chain and the recovery analysis have not been run.

The implants are aligned as their NA12878 background was (bwa `-Y`, no `-M`); the
patient BAMs were aligned with `-M` and without `-Y`. This is recorded as a
limitation in the run record.

The demo bundle was rebuilt from the new IMP01 and IMP10 (1.59 MB, `demo_bundle/`,
not committed). Its `DEMO.md` still carries the old implants' hard-coded claims —
a score of 47.5/100 `moderate` for DEMO_CLEAN, WITHHELD for DEMO_REPEAT — which
have not been checked against the new implants.

| New file | What it holds |
|---|---|
| `IMP01_gate_revised_2026-09-25.json` | the revised criterion, per read, with (a2), (b), (c) and the dose table |
| `implants_ground_truth.json` | all twelve implants: both breakpoints, each junction's orientation and BND notation, class, pre-implant depth and low_mapq_fraction at both ends, compensation counts, every seed |
| `run_record_2026-09-25_continued.json` | status, versions, seeds, thresholds, both gates, delly exit codes, samples, records, peak RSS, limitations |

## 2026-09-25, evening — the evidence chain and the recovery analysis

A new entry; nothing above is edited. Everything here is under
`analysis_2026-09-25/`, produced by scripts committed before they ran
(`scripts/synthetic/evidence_chain.py`, `adjudicate_reads.py`, `delly_ladder.py`,
`realign_noalt.py`), every tool call made through the MCP dispatch the interface
uses and recorded with its return. A new experiment: the figures quoted from the
first run, whose records were lost, are set beside these only for comparison.

**The chain.** Each of the thirteen BCFs was loaded through the bridge, filtered
by the standard funnel (PASS, PE ≥ 3, SR ≥ 1, primary contigs, outside the exclude
template; no svtype step, as in the first run) and the evidence tools were run at
both breakends of every survivor: 5,398 calls, no dispatch error, no error return.
An implanted junction counts as detected when a survivor matches it in orientation
with both breakends within 500 bp.

| Figure | Rebuild | First run (quoted) |
|---|---|---|
| Junctions detected | **16 of 24** | 14 of 24 |
| Translocations detected | **8 of 12**, each with both junctions | 7 of 12 |
| By class | clean 8/8 · repeat-adjacent **8/8** · low-mappability 0/8 | 8/8 · 6/8 · 0/8 |
| Where the losses happened | all 8 at discovery: delly wrote no record for them | all at discovery |
| Background funnel | 894 → 513 → 155 → 27 (primary 27, unmasked 27) | 894 → 513 → 155 → 27 |
| Per implant BAM | 29 survivors where detected, 27 where not; the 27 non-implanted survivors are the background's 27, identical by candidate_id, in all twelve | each detected implant added exactly 2 |
| Localisation of the 16 | 6 at 0 bp, 10 at 1 bp (larger offset of the two ends) | 12/14 at 0 bp, 2/14 at 2 bp |
| PE × SR grid | 16/24 in every cell; non-implanted survivors per BAM 156 / 156 / 155 / 118 at SR ≥ 0 and 28 / 28 / 27 / 25 at SR ≥ 1 (PE ≥ 1 / 2 / 3 / 5); SR ≥ 2 identical to SR ≥ 1 | 14/24 in every cell; SR ≥ 1 cut 156 → 28; SR ≥ 2 identical |

**Class robustness.** Under the stricter upper-decile rule (a dinucleotide run of
8 bp counts as repeat) IMP01 becomes repeat-adjacent: clean 6/6 and repeat-adjacent
10/10. **The per-class denominators change; the per-class rates (100 %, 100 %,
0 %) do not.**

**The reused coordinates.** IMP01 detected (as before); IMP10 missed (as before);
**IMP06 detected, where the first run missed it.** Realigned without the `.alt` —
all twelve implants' reads pooled exactly as the original, with a control alignment
that reproduced the committed one record for record — IMP06 loses both junctions:
delly reports neither and the BAM's survivors fall back to 27. Its junction pairs
with both mates at MAPQ ≥ 20 go from 37 to 0. IMP01 and IMP10, realigned as
controls, change in no read and no call. **ALT-aware alignment, not the random
draw, explains the difference** (`noalt/`). The first run's index is not recorded;
this shows the rebuild would have missed IMP06 without the `.alt`, which is the
condition that reproduces the first run's result.

**The scoring ceiling** (32 breakends of the 16 detected junctions): evidence
scores 32.5–55.0 (30 moderate, 2 weak); the observed discordant fraction at most
0.164 against the 0.5 of the top tier; the tool's own attainable ceiling 57.5 at
every breakend, below the strong band's 70; strong reachable nowhere. The depth
layer contributed 15 points at 7 breakends of IMP04, IMP05 and IMP07, where no
copy-number change was built in. First run: 40.0–55.0 all moderate, maximum
0.175, 57.5 at IMP01 chr20, spurious depth at IMP02 and IMP04.

**Rescue.** The four missed implants (IMP09–IMP12, all low-mappability) queried at
their true coordinates in fresh processes: every return reads provenance
`caller_supplied` (13 of 13 calls each), and both controls hold. The strongest,
IMP09 chr20:33,700,000, shows 10 discordant pairs, 7 soft-clipped reads (at most 1
at any position) and 9 split reads and scores 22.5 `weak`; IMP11 chr20:31,100,000
shows 11, 4 and 41 but is `QUALITY-LIMITED`; the rest are `QUALITY-LIMITED` with
0–1 discordant pairs. The first run's rescue example, IMP06 (44, 20 and 14), is
detected by the caller in this rebuild.

**Per read, against ART's ground truth** (`adjudication/`). The enumeration
reproduces every tool count exactly at all 24 breakpoints. Of 390 truly spanning
reads, the split_reads tool (its MCP default, min_mapq 0) counts 154 genuinely
(39.5 %); delly's SR sums to 137 (35.1 %). The tool's raw count is closer to the
truth than delly's SR in 10 of 12 implants (2 ties). Its false-positive rate by
class is 11.1 % / 20.4 % / 87.7 % (clean / repeat-adjacent / low-mappability),
r = 0.668 with the low_mapq_fraction; 41 of its counted reads truly span the
junction but carry only a mapQ-0 SA entry naming a wrong chromosome — the IMP01
gate's "no wrong partner" held at IMP01 only — and counting those as true gives
0 % / 5.8 % / 60.0 % and r = 0.298. At the low-mappability class the false
positives are concentrated at IMP11 chr20:31,100,000 (40 of 41). 149 of the 390
reads carry no SA tag; the soft-clip layer catches 77 of them. For the 72 it
misses, the first reason in the layer's own filter order: 19 lie outside both
breakpoint windows, 19 have MAPQ under 20, 31 have a shorter side under 10 bp and
3 end in a clip under 10 bp despite a longer overhang. 11 SA-less reads have a
27–29 bp shorter side; as at the gate, none of them — nor any other SA-less read
under 30 bp — could reach -T 30 with a 19 bp seed through chance matches.
Raising split_reads' min_mapq from 0 to 20 would lose 30 genuine reads (154 → 124)
and remove 85 of 86 false positives, and would zero the layer at all 8
low-mappability breakpoints. The default is unchanged. First run: 60.5 % against
54.4 %; FP 0 % / 3.9 % / 56.9 % driven by one locus, r = +0.142; 121 of 304 without
SA, 114 caught.

**Caller versus tool at the background's 27 survivors.** Where both breakends have a
low_mapq_fraction under 0.01 the tool's count (the larger of the two breakends)
differs from delly's SR by −5 to +9 reads; where a breakend reaches 0.4 or more the
divergence reaches 100-fold (1,700 against 17 at chr20:34.23 Mb, low_mapq_fraction
0.905); r(|difference|, low_mapq_fraction) = 0.689 (`figures.json`,
`caller_vs_tool_background`).

**The delly MAPQ ladder** (`ladder/`), on the four missed implants and the
background, twenty jobs at once, all exit 0; the default setting rerun is identical
record for record to the committed BCFs. At the defaults and at -q 1 -r 5 none of
the 8 missed junctions is called; at -q 0 (-r 5 or -r 0) one is (IMP11 J20,
`LowQual`) and it does not survive the funnel. On the background the raw BND count
rises 13 → 81 (-q 1 -r 5), 93 (-q 0 -r 5), 85 (-q 0 -r 0) and the survivors 27 → 27,
28, 28. First run: at -q 0, 2 of 10 missed junctions called and 1 survived; raw BND
about 7×; survivors 27 → 28.

**Checks that failed and were fixed, kept on record.** The no-ALT alignment's first
verification counted lines naming ALT contigs, so bwa's "read 0 ALT contigs" failed
it (`noalt/align_first_check_defective.json`); the check now reads the number.

| New file (under `analysis_2026-09-25/`) | What it holds |
|---|---|
| `chain/<label>.json.gz` | per BCF: the bridge load, every funnel and grid list, the evidence tools at both breakends of every survivor — every call and return |
| `rescue/IMPxx.json.gz` | the evidence tools at each implant's true breakpoints, with the provenance controls |
| `figures.json`, `figures.log.txt` | the figures above, with the per-implant matching |
| `adjudication/IMPxx.json.gz`, `summary.json` | per read, against ART's ground truth |
| `ladder/run.json`, `analysis.json` | the twenty ladder jobs and what they call |
| `noalt/` | the no-ALT realignment, its control, the three rebuilt BAMs' records, delly and the comparison |
