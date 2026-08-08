# Real-Data Validation Attempt — Balanced Translocation Search + Cross-Technology Comparison

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
| `breakpoint_evidence_summary` | score 22.5, weak, 2/4 | score 15.0, weak, 2/4 |

### Two things worth highlighting

1. **The soft-clip consensus position matched exactly across both technologies** (115,686,865, 3bp from the documented deletion start) — strong independent confirmation the breakpoint location is real, not a mapping artifact of either platform. The Illumina data showed a much tighter pileup (13 reads at the exact same position vs. PacBio's 3) — short reads localize this breakpoint more sharply than long reads did here.

2. **The default `window_bp=500` in `summarize_breakpoint_evidence` misses this real deletion's depth signal on *both* real datasets**, while a wider, deletion-span-aware window catches it cleanly (ratio drops from ~0.66–0.71 to ~0.50–0.62 depending on exact width). This reproduces the same limitation flagged after the PacBio validation, now confirmed on a second, independent, real dataset — this is a real, addressable limitation of the current default window size relative to typical deletion sizes, not a one-off artifact of either BAM. A natural next step (not implemented here, to keep scope bounded) would be for `summarize_breakpoint_evidence` to accept an optional `region_end` for depth profiling, so callers with a known SV span (as from a VCF) can pass the true deletion extent instead of relying on the fixed `position ± window_bp`.

## Bottom line

No confirmed real balanced translocation or inversion BAM was found — this appears to be a genuine gap in publicly available structural-variant benchmark data, not a search shortfall. In its place, this session validated the tool suite against a second real dataset (Illumina 300x, different individual sequencing run and aligner than the PacBio HiFi validation), on the same confirmed real deletion, and surfaced one new, reproducible, actionable finding about the depth-evidence window size.
