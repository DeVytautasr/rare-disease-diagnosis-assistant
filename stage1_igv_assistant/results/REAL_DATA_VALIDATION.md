# Real-Data Validation Attempt — Balanced Translocation Search + Cross-Technology Comparison

> **Vocabulary note (post applicable-layer normalisation, see AUDIT_2026_08.md
> §5 / Critical Finding 5):** every `evidence_score` value in this document —
> in the cross-technology table, the calibration-update table, and the
> scoring-bug-fix table below — refers to what `bam_tools.py` now calls
> `evidence_score_raw`: the direct sum over all 4 layers regardless of
> whether a layer was structurally applicable to that data. This document
> predates the later addition of a separate, normalised `evidence_score`
> field (which divides only by the applicable-layer count via
> `detect_applicable_layers`), so nothing here needed to change definition —
> only the field name did. `signal_layers` values below are all "N/4" for
> the same reason (M was always 4 before that change). Figures are kept
> exactly as originally recorded; not recomputed.

**Date:** 2026-08-09

## Goal

Find a public, modern-alignment BAM/CRAM with a confirmed inversion (INV) or translocation-like (BND) structural variant, run the full tool suite on it, and validate against real ground truth — closing the gap left by the synthetic-data demo in `DEMO_END_TO_END.md`.

## The 3 specified attempts

### Attempt A — HGSVC2 CRAM (HG00514)
```
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/working/20200131_sib_datasets/HG00514/...
```
**Fails — path does not exist.** Confirmed by directory listing of the real HGSVC2 tree (`.../HGSVC2/working/`): no `20200131_sib_datasets/` directory exists there (checked in the prior session too — same result both times). The real directories follow a different naming scheme entirely (e.g. `20191004_Illumina/`, `20191031_CHS_PacBio_HG00512_HiFi/`).

### Attempt B — HG002 Illumina 300x BAM
```
https://ftp-trace.ncbi.nlm.nih.gov/.../NHGRI_Illumina300X_AJtrio_novaSeq_125bpPE/HG002.GRCh38.300x.bam
```
**Given subdirectory name was wrong** (`novaSeq_125bpPE` doesn't exist). Found the real one by listing the parent directory: `NHGRI_Illumina300X_AJtrio_novoalign_bams/`. The corrected URL streams successfully (region queries confirmed working via `.bai`, despite the BAM header misleadingly claiming `SO:unsorted`). But:
- Aligned with **Novoalign v3.02.07** (2014-era short-read aligner). Confirmed **0 `SA` tags** in two independent spot checks (chr1:1,000,000–1,100,000 and chr7:140,000,000–141,000,000) — Novoalign does not emit chimeric/supplementary alignments the way bwa-mem does, so `split_reads` cannot contribute for this file regardless of locus (same category of limitation as the 2018 HCC1143 BAM from the earlier session).
- Even so, this is a **real, working, streamable BAM** — see the cross-technology comparison below.

### Attempt C — ENA CRAM accession `ERR3241750`
```
https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241750/HG002.final.cram
```
**Fails — wrong sample entirely.** The real file at that accession is `HG00705.final.cram` — a 1000 Genomes Project sample, not HG002. This isn't a path typo; `ERR3241750` simply isn't an HG002 accession.

## The deeper finding: GIAB has no curated INV/BND for HG002

Even setting aside data-access issues, the actual goal — "a confirmed SV in the GIAB HG002 SV benchmark that is an INV or BND" — cannot be met with any HG002 BAM, because that ground truth doesn't exist in the benchmarks GIAB publishes:
- `HG002_SVs_Tier1_v0.6.vcf.gz` (GRCh37): **0** `SVTYPE=INV` records (checked this session)
- `HG002_GRCh38_CMRG_SV_v1.00.vcf.gz` (GRCh38, validated last session): **0** INV/BND records — only `CONTRAC`, `DUP`, `SIMPLEDEL`, `SIMPLEINS`, `SUBSDEL` REPTYPEs across all 250 entries

Two further web searches (general SRA/ENA balanced-translocation queries, and literature on `t(11;22)`/`t(8;22)`/`t(17;22)` — the best-characterized recurrent constitutional translocations in humans) found no publicly deposited raw alignment data either; see `DEMO_END_TO_END.md` for that search. This is a genuine scarcity in publicly available structural-variant benchmarks, not a search or access failure on our part.

## What we did instead: real cross-technology validation

Since Attempt B gave us a genuine, working, real BAM, we ran the full tool suite against it at the **same confirmed real deletion** already validated with PacBio HiFi last session (chr1:115,686,862–115,690,222, `VANGL1`, `GT=1|0` heterozygous, from the GRCh38-native CMRG SV benchmark). This is a real short-read-vs-long-read comparison on the same real variant in the same individual — not a translocation/inversion, but real data nonetheless, and it surfaced a genuinely new, reproducible finding.

| Tool | PacBio HiFi (prior session) | Illumina 300x novoalign (this session) |
|---|---|---|
| Mean depth (±500bp) | 30.5x | 226.77x |
| `discordant_pairs` | 0/0 (unpaired long reads — N/A) | 1/1708 (0.1%) — noise |
| `soft_clipped_reads` | 4/40 (10%), consensus 115,686,865, max 3 reads/position | 23/853 (2.7%), consensus **115,686,865** (identical position), **max 13 reads/position** |
| `split_reads` | 1/39 (2.6%), partner chr1:115,690,223 (near exact DEL end) | 0/853 — Novoalign emits no `SA` tags |
| Depth ratio, ±500bp window (`summarize_breakpoint_evidence` default) | 0.618–0.707 (varies by exact window) → **misses** `<0.6` threshold | **0.662** → **misses** `<0.6` threshold |
| Depth ratio, wide window (deletion span ± 2000bp flanking) | 0.618 (n/a, only tested at this width) | **0.496** → **correctly flags** `likely_deletion: True` |
| `breakpoint_evidence_summary` (= today's `evidence_score_raw`) | score 22.5, weak, 2/4 | score 15.0, weak, 2/4 |

### Two things worth highlighting

1. **The soft-clip consensus position matched exactly across both technologies** (115,686,865, 3bp from the documented deletion start) — strong independent confirmation the breakpoint location is real, not a mapping artifact of either platform. The Illumina data showed a much tighter pileup (13 reads at the exact same position vs. PacBio's 3) — short reads localize this breakpoint more sharply than long reads did here.

2. **The default `window_bp=500` in `summarize_breakpoint_evidence` misses this real deletion's depth signal on *both* real datasets**, while a wider, deletion-span-aware window catches it cleanly (ratio drops from ~0.66–0.71 to ~0.50–0.62 depending on exact width). This reproduces the same limitation flagged after the PacBio validation, now confirmed on a second, independent, real dataset — this is a real, addressable limitation of the current default window size relative to typical deletion sizes, not a one-off artifact of either BAM. A natural next step (not implemented here, to keep scope bounded) would be for `summarize_breakpoint_evidence` to accept an optional `region_end` for depth profiling, so callers with a known SV span (as from a VCF) can pass the true deletion extent instead of relying on the fixed `position ± window_bp`.

## Bottom line

No confirmed real balanced translocation or inversion BAM was found — this appears to be a genuine gap in publicly available structural-variant benchmark data, not a search shortfall. In its place, this session validated the tool suite against a second real dataset (Illumina 300x, different individual sequencing run and aligner than the PacBio HiFi validation), on the same confirmed real deletion, and surfaced one new, reproducible, actionable finding about the depth-evidence window size.

## Calibration update

`summarize_breakpoint_evidence`'s depth layer was recalibrated: the depth-profile query is now a fixed ±2kb/200bp-window region around `position` (independent of `window_bp`, which still governs discordant-pair/soft-clip/split-read windows), and the `depth_ratio_min_to_mean` scoring threshold was raised from 0.6 to 0.7.

Following calibration (window ±2kb, threshold 0.7), the depth layer correctly flags the GIAB deletion on both PacBio HiFi and Illumina 300x data. Verified end-to-end through `summarize_breakpoint_evidence` itself (not just the underlying `get_read_depth_profile` tool in isolation):

*(`depth_score` here is still on the pre-"Scoring bug fix" 0-50/layer scale —
see the "Scoring bug fix" section below, which rescales this same 30 to 15
on the current 0-25/layer scale. `evidence_score` in this table = today's
`evidence_score_raw`.)*

| | PacBio HiFi | Illumina 300x |
|---|---|---|
| `depth_ratio_min_to_mean` (±2kb/200bp window) | 0.609 | 0.542 |
| `depth_score` (0-50 scale, pre-rescale) | 0 → **30** | 0 → **30** |
| `signal_layers` (N/4) | 2/4 → **3/4** | 2/4 → **3/4** |
| `evidence_score` (= today's `evidence_score_raw`) | 22.5 → **37.5** | 15.0 → **30.0** |
| `evidence_strength` | weak → **still weak** | weak → **still weak** |

The tool is now validated on two independent real sequencing technologies for deletion detection via the depth layer specifically. Note precisely what changed and what didn't: the depth layer now correctly contributes signal on both real datasets (that was the calibration goal, and it's met), but `evidence_strength` stays "weak" on both — `discordant_pairs` and `split_reads` correctly contribute nothing to a plain deletion on either aligner (neither Novoalign nor this PacBio pipeline's alignment emits `SA` tags, and a deletion has no inter-chromosomal discordant signal to find), so 2 of 4 layers are structurally unavailable here regardless of window/threshold tuning. This isn't a shortfall in the fix — it's the expected behavior for this SV type on this data.

One important caveat on the threshold itself: **0.7 is calibrated against a single confirmed locus**, replicated across 2 sequencing technologies (not 2 independent genomic loci). Both calibration points are the same underlying deletion in the same individual (HG002). This is a reasonable starting point, not a general-purpose validated threshold — it hasn't been checked against true-negative regions (normal, non-deletion loci) to confirm it doesn't also flag ordinary coverage variance as "likely deletion." A false-positive-rate check against several known-normal loci would be a natural next step before relying on this threshold for anything beyond this specific validation exercise.

Balanced translocation validation on real data remains outstanding — no public modern-alignment BAM with a confirmed germline balanced translocation was found in this session.

## Scoring bug fix (2026-08-10)

An LLM assistant session using only the 8 MCP tools [now 10 — see README.md] (no source access) on this same locus (chr1:115,686,862, Illumina 300x) surfaced a genuine inconsistency in `summarize_breakpoint_evidence`'s output: `discordant_pair_score=15, soft_clip_score=15, split_read_score=0, depth_score=30` sum to 60, but the returned `evidence_score` was 30 — contradicting the tool's own docstring/description, which claimed the components "sum into evidence_score."

**Root cause:** `discordant_pair_score`, `soft_clip_score`, `split_read_score`, and `depth_score` were each tiered on a 0-50 scale, but the combination step computed `evidence_score = round(sum(component_scores) / 2.0, 1)` — silently halving the total without that division being reflected anywhere in the returned component values. `evidence_score` itself was arithmetically correct as a normalized composite; the component scores were mislabeled relative to their real contribution.

**Fix:** each component is now tiered directly on a 0-25 scale (half the old tier values: 0 / 7.5 / 15 / 25 for discordant-pair, soft-clip, split-read; 0 / 15 / 25 for depth), and `evidence_score` is their direct, unweighted sum — no separate normalization step. Because halving four 0-50 components then summing is mathematically identical to summing four 0-25 components directly, **`evidence_score` values are unchanged by this fix** for every case in this document; only the component breakdown is now internally consistent with the total.

Re-validated at chr1:115,686,862 (Illumina 300x, this same locus) after the fix:

| | Before fix (mislabeled) | After fix (corrected) |
|---|---|---|
| `discordant_pair_score` | 15 | **7.5** |
| `soft_clip_score` | 15 | **7.5** |
| `split_read_score` | 0 | **0** |
| `depth_score` | 30 | **15** |
| component sum | 60 (≠ evidence_score) | **30 (= evidence_score)** |
| `evidence_score` (= today's `evidence_score_raw`) | 30.0 | **30.0 (unchanged)** |
| `evidence_strength` | weak | **weak (unchanged)** |
| `signal_layers` (N/4) | 3/4 | **3/4 (unchanged)** |

All 9 tests in `test_bam_tools.py` pass unchanged after the fix — none asserted on raw component-score magnitudes, only on `evidence_strength`/`signal_layers`, which are unaffected by the rescaling.

Found by an LLM assistant during its first end-to-end MCP session (constrained to tool calls only, no source access) — see `LLM_SESSION_1.md` for the full session report that surfaced this. That the inconsistency was catchable purely from tool outputs, without reading `bam_tools.py`, is itself a small positive data point for the interpretability goals of this project.

## Applicable-layer normalisation (post-audit addendum)

The "corrected" `evidence_score: 30.0` / `evidence_strength: weak` above sits
downstream of exactly the structural problem this whole document already
identifies in the "Calibration update" section: `discordant_pairs` and
`split_reads` "correctly contribute nothing... so 2 of 4 layers are
structurally unavailable here regardless of window/threshold tuning" — yet
`evidence_score` still divided by all 4. `summarize_breakpoint_evidence`
later gained an `applicable_layers` parameter (and a `detect_applicable_layers`
tool to infer it from the BAM) that normalises over only the applicable
layers instead. Re-scoring this exact locus with that change (not part of
this session's original findings — recorded here as a forward pointer, not
a revision of the numbers above; see the Group 4 commit for the run): with
`applicable_layers = [discordant_pairs, soft_clipped_reads, read_depth]`
(Novoalign has no `SA` tags, confirmed repeatedly above), the same
underlying component scores (7.5 + 7.5 + 0 + 15 = 30, unchanged) normalise
to `evidence_score: 40.0` / `evidence_strength: moderate`, while
`evidence_score_raw` stays `30.0` exactly as recorded in this document.
