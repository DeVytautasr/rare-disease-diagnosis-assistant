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
