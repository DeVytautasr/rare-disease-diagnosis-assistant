---
NOTE: This file records the demo session from the original Stage 1 run.
Two caveats listed below were fixed in subsequent commits:
- Caveat #1 (interpretation_template missing) → fixed in 385c77f
- Caveat #4 (no retry logic in get_gene_at_locus) → fixed in 385c77f
REAL_DATA_VALIDATION.md contains the current, fully up-to-date validation record.
---

# End-to-End Demo — Full Pipeline on a Synthetic Balanced Translocation

**Date:** 2026-08-09
**Case ID:** `demo_synthetic_translocation`
**BAM:** `stage1_igv_assistant/data/bam/synthetic_translocation_demo.sorted.bam` (gitignored — regenerate with `create_translocation_bam()` from `test_bam_tools.py`)

## Data note: synthetic, not real

Before running this demo, we searched for a publicly downloadable BAM/CRAM from a **confirmed germline balanced translocation carrier**, across:
- The three sources named for this search (HGSVC2 CRAM streaming, SRA via `esearch`/`efetch`, GIAB HG002 inversion search)
- Two further targeted web searches (general SRA/ENA balanced-translocation queries; the literature on `t(11;22)`, `t(8;22)`, `t(17;22)` — the best-characterized recurrent constitutional translocations in humans)

None produced a public raw-alignment file. Papers characterizing these translocations either predate NGS (GenBank junction sequences only) or explicitly restrict raw sequencing data to reasonable-request/institutional access (e.g. the pig translocation study, Topigs Norsvin). Unlike GIAB's curated deletion benchmark (validated against real HG002 PacBio HiFi data in the previous session), there does not appear to be an equivalent public benchmark for balanced translocations.

**This demo therefore uses the synthetic translocation BAM** built by `create_translocation_bam()` (chr1↔chr8, with proper `SA` tags, validated in TEST 6 and TEST 8) to prove the full 8-tool pipeline runs correctly end-to-end. It demonstrates pipeline correctness, not a real-world positive control.

---

## Step-by-step tool output

### BP1 — chr1:1,050,000

| Tool | Result |
|---|---|
| `bam_stats_at_locus` | 28 reads in ±500bp window, mean depth 2.73x, mean MAPQ 55.5, 0% low-MAPQ |
| `discordant_pairs` | 15/28 (54%) discordant, all mates on **chr8** |
| `soft_clipped_reads` | 5/13 (38%) clipped, consensus position 1,049,950 |
| `split_reads` | 8/13 (62%) split, all partners on **chr8** (e.g. chr8:47,000,000, chr8:47,000,020, ...) |
| `read_depth_profile` | min 0, max 15, mean 5.4 reads/window → `depth_ratio_min_to_mean: 0.0`, `likely_deletion: True` |
| `gene_at_locus` | **AGRN** (protein_coding), chr1:1,020,069–1,056,119 |
| `breakpoint_evidence_summary` | **evidence_score: 100.0 → "strong"**, `signal_layers: "4/4"` |

### BP2 — chr8:47,000,000

| Tool | Result |
|---|---|
| `discordant_pairs` | 15 discordant, all mates back on **chr1** |
| `soft_clipped_reads` / `split_reads` | 0 (this reciprocal side of the fixture only carries discordant-pair evidence, not soft-clips/split-reads — see caveats below) |
| `read_depth_profile` | also flags `depth_score: 50` |
| `gene_at_locus` | intergenic (`gene_count: 0`) |
| `breakpoint_evidence_summary` | evidence_score: 50.0 → "moderate", `signal_layers: "2/4"` |

### Reciprocal breakpoint check

```
verdict: "RECIPROCAL CONFIRMED — both breakpoints show concordant inter-chromosomal signal"
is_balanced: true
```

Both sides show discordant pairs pointing at each other (15 chr1→chr8, 15 chr8→chr1) — the expected signature of a true balanced translocation, as opposed to a one-sided artifact.

## Case summary (`BamCase.summary()`)

```
Case: demo_synthetic_translocation
BAM: .../synthetic_translocation_demo.sorted.bam
Technology: short_read_illumina
Suspected SV: balanced_translocation
Cytogenetics: t(1;8)(synthetic demo locus)
Breakpoints investigated: 2
  BP1_chr1: chr1:1050000 — strong (4/4)
  BP2_chr8: chr8:47000000 — moderate (2/4)
```

Full structured record saved to `stage1_igv_assistant/results/demo_case.json`.

---

## Caveats found while building this demo

1. **(FIXED in 385c77f)** ~~`breakpoint_evidence_summary`'s docstring (in `server.py`) claims it returns an `interpretation_template` field — it doesn't.~~ The actual keys at the time were: `chromosome, depth_profile, depth_score, discordant_pair_score, discordant_pairs, evidence_score, evidence_strength, label, locus_stats, position, signal_layers, soft_clip_score, soft_clips, split_read_score, split_reads, supporting_observations`. The docstring should be corrected in a future pass — left as-is here since Step 1 asked for the file verbatim. *(As of `385c77f`, `interpretation_template` is a real field and the docstring lists it accurately — see `REAL_DATA_VALIDATION.md` for current state.)*

2. **(FIXED 2026-08-11)** ~~The depth-profile layer flagged `likely_deletion: True` at a translocation locus, contributing a full 50-point `depth_score`.~~ This isn't because the depth tool is wrong about deletions — it's a known artifact of this synthetic fixture (documented in `test_bam_tools.py`'s TEST 6 comment): the fixture deliberately isolates its signal reads with read-free gaps to keep other fractions clean, and one of those gaps has 0 reads, which trips the ratio threshold. The `server.py` system instructions explicitly warn the LLM that "flat depth is expected [for balanced translocations]; do not interpret it as negative evidence" — but here the depth layer produced a *positive* (deletion-like) signal instead of a flat one, for a translocation. *(`get_read_depth_profile` gained a `focus_position` parameter and `summarize_breakpoint_evidence`'s `depth_score` now requires `dip_is_at_focus` — the region's true minimum, not just a nearby one, has to actually be near the queried position — before awarding points. This fixture's read-free gap sits well outside the default 1000bp tolerance of BP1's position, so it's now correctly excluded: `depth_score=0` for BP1, surfaced instead as an off-position observation. As a direct consequence, BP1's headline `evidence_score` in the table above is no longer 100.0/"strong" under the current code — re-running this exact case now gives `evidence_score: 50.0`/"moderate", `signal_layers: "2/4"` (soft_clip_score also drops to 0 for the same class of reason — see FIX 3 / `results/LLM_SESSION_3_BLIND.md`'s "Post-fix re-verification" section — this fixture's soft-clipped reads are staggered 3bp apart and never formed a genuine same-position pileup). The table above is left as originally recorded, not rewritten.)*

3. **Gene lookups (`gene_at_locus`) are independent of the synthetic BAM's read data.** AGRN is genuinely at chr1:1,050,000 on GRCh38 — but the synthetic reads there have nothing to do with AGRN biology; they were placed purely to exercise the discordant/clip/split-read tools. This lookup would carry real meaning in production, on real coordinates.

4. **(FIXED in 385c77f)** Ensembl's REST API was intermittently slow/unresponsive (multiple 15s timeouts across this session, both for this demo and earlier gene lookups) — genuine third-party flakiness, not a bug in `get_gene_at_locus`. ~~Production use of this tool should implement retry logic, which the current implementation does not have.~~ *(As of `385c77f`, `get_gene_at_locus` retries with exponential backoff on HTTP 429 and network/timeout errors, up to `max_retries` attempts — see `REAL_DATA_VALIDATION.md` for current state.)*
