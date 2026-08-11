# evidence_panel Validation — HCC1143 chr21 Locus

**Date:** 2026-08-11
**Tool:** `evidence_panel` (added 2026-08-11, commits `deb67c9`/`dec56c5`) — generates one screenshot per evidence layer instead of a single combined image, each with the IGV settings that isolate that layer visually.
**BAM:** `stage1_igv_assistant/data/bam/HCC1143.normal.21.19M-20M.bam`
**Locus:** chr21:19,089,694 — same locus already characterized in `RESULTS_HCC1143.md` (labeled "translocation" in the source workshop's tutorial script, but no real signal found there).

This session had no working example output for `evidence_panel` since it was added the same day. Purpose here: confirm the tool runs end-to-end against a local BAM (no network dependency) and that its screenshots visually agree with the numeric findings already recorded in `RESULTS_HCC1143.md`.

## Setup

`applicable_layers` on this BAM returns:

```json
{"applicable_layers": ["discordant_pairs", "soft_clipped_reads", "read_depth"],
 "split_reads": "not applicable — no SA tags observed in 1000 sampled reads"}
```

Consistent with the whole-BAM check in `RESULTS_HCC1143.md` §4 (0 supplementary-alignment reads across all 572,731 reads on chr21 — this is a 2018-era alignment with no chimeric/SA-tag output). `bam_stats_at_locus` on the full documented region (chr21:19,089,694–19,095,362) reproduced the same figures on file: 4,097 reads, mean depth 62.33x, mean MAPQ 31.09.

## Note on this session's environment

The first `evidence_panel` attempt failed for an infrastructure reason unrelated to the tool itself: the MCP server process (launched via `~/.claude.json`, pointed directly at the `rda` conda env's Python) had no `PATH` entry for that env's `java` binary, so IGV's batch script exited immediately with `java: not found`. Fixed by adding `PATH` (rda `bin/` first) to that server's `env` in `~/.claude.json` and reconnecting the MCP session. Not a bug in `bam_tools.py`/`server.py` — worth knowing if `evidence_panel`/`igv_screenshot` fail silently again in a fresh MCP session.

## Results

| Layer | Window used | Screenshot | Agrees with RESULTS_HCC1143.md? |
|---|---|---|---|
| `discordant_pairs` | chr21:19,088,194–19,091,194 (±1500bp) | `screenshots/hcc1143_chr21_discordant_pairs.png` | Yes — `UNEXPECTED_PAIR` coloring shows scattered single-colored pairs, no dominant partner-chromosome color cluster, matching the "5 discordant reads, one each to chr20/chr5/chr19/chr12/chr18, no partner clustering" finding. |
| `soft_clipped_reads` | chr21:19,089,544–19,089,844 (±150bp, dominant side auto-detected) | `screenshots/hcc1143_chr21_soft_clipped_reads.png` | Yes — visible staircase clip patterns at several distinct positions rather than one shared boundary, matching "22/287 soft-clipped, but only 1 read at the consensus position — no real pileup." `clip_side_determination` reported 7 right-clipped vs. 6 left-clipped reads in this window (near-even split; sorted by `RIGHT_CLIP`). |
| `read_depth` | chr21:19,086,694–19,092,694 (±3000bp, fixed scale 0–230) | `screenshots/hcc1143_chr21_read_depth.png` | Yes — coverage track is essentially flat across the window with no dip, consistent with this locus never being flagged as a depth-supported deletion anywhere in prior validation. |
| `split_reads` | — | skipped | Correctly skipped — excluded via `applicable_layers` (no SA tags on this BAM). |

## Bottom line

`evidence_panel` runs end-to-end on a local, real, network-free BAM and its per-layer screenshots are visually consistent with the numeric evidence already validated in `RESULTS_HCC1143.md` for this same locus — i.e., this remains a **weak/no-signal** locus, and the panel correctly shows *why* (scattered discordant partners, scattered clip positions, flat depth) rather than just asserting it. This is not a new positive-signal validation case (that would need a locus with real, strong signal — the synthetic translocation BAM used in `DEMO_END_TO_END.md` would be a better candidate for a follow-up demo of what the panel looks like on true signal); it validates that the tool's plumbing (window selection, layer skipping, dominant-clip-side detection) works correctly on non-trivial real data.
