# Synthetic positive control, rebuilt 2026-09-25 — stopped at the IMP01 gate

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
