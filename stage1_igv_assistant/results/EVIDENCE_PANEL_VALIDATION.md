# evidence_panel Validation

**Tool:** `evidence_panel` (added 2026-08-11, commits `deb67c9`/`dec56c5`) — generates one screenshot per evidence layer instead of a single combined image, each with the IGV settings that isolate that layer visually.

This session had no working example output for `evidence_panel` since it was added the same day. Two runs below cover both ends of the signal spectrum: a positive control (synthetic translocation, signal present on all 4 layers) and a negative control (real HCC1143 locus, no real signal on any layer) — confirming the tool's window selection, layer skipping, and dominant-clip-side detection all work correctly whether or not there's anything to find.

## Note on this session's environment

The first `evidence_panel` attempt failed for an infrastructure reason unrelated to the tool itself: the MCP server process (launched via `~/.claude.json`, pointed directly at the `rda` conda env's Python) had no `PATH` entry for that env's `java` binary, so IGV's batch script exited immediately with `java: not found`. Fixed by adding `PATH` (rda `bin/` first) to that server's `env` in `~/.claude.json` and reconnecting the MCP session. Not a bug in `bam_tools.py`/`server.py` — see `stage1_igv_assistant/README.md`'s Known limitations for the general note; worth knowing if `evidence_panel`/`igv_screenshot` fail silently again in a fresh MCP session.

---

## 1. Positive control — synthetic translocation (chr1↔chr8)

**Date:** 2026-08-11
**BAM:** `stage1_igv_assistant/data/bam/synthetic_translocation_demo.sorted.bam` (gitignored; regenerate with `create_translocation_bam()` from `test_bam_tools.py` — same fixture the automated TEST 8 real-panel check uses, see `test_bam_tools.py` line ~998).
**Locus:** chr1:1,050,000, mate breakpoint chr8:47,000,000 — 15 discordant pairs each direction, 8 SA-tagged split reads, 5 soft-clipped reads, per the fixture's docstring.

`applicable_layers` on this BAM returns all 4 layers applicable (243 reads sampled, at least one SA tag and one paired read observed) — as expected, since this fixture was built specifically to exercise every layer, unlike the real HCC1143 BAM below.

| Layer | Window used | Screenshot | What it shows |
|---|---|---|---|
| `discordant_pairs` | chr1:1,048,500–1,051,500 (±1500bp) | `screenshots/synthetic_translocation_discordant_pairs.png` | `UNEXPECTED_PAIR` coloring shows two distinct, tight clusters rather than the scattered single-read colors in the HCC1143 negative control below: a blue/orange cluster at chr1:1,049,950–1,050,014 (the fixture's overlapping soft-clip and local-mate split-read groups — not discordant, but flagged as not-proper-pair so they're colored as anomalous too) and a separate red cluster at chr1:1,050,250–1,050,320, ~236bp to the right — the actual 15-read discordant-pair group whose mates are on chr8. (The 15 *reciprocal* discordant pairs are recorded at chr8:47,000,000+ in the fixture, so they correctly do not appear in this chr1-restricted view at all.) |
| `soft_clipped_reads` | chr1:1,049,850–1,050,150 (±150bp) | `screenshots/synthetic_translocation_soft_clipped_reads.png` | A clean staircase pileup converging on one clip boundary at chr1:1,050,000 — `clip_side_determination`: "Sorted by LEFT_CLIP: 5 left-clipped vs 0 right-clipped reads observed in this window," a real dominant side (unlike HCC1143's near-even, no-pileup split). |
| `split_reads` | chr1:1,048,500–1,051,500 (±1500bp) | `screenshots/synthetic_translocation_split_reads.png` | This is the layer only this tool can show visually. Reads are grouped by `TAG SA` into individual labeled rows, each showing its exact SA partner locus text — all 8 rows read `chr8,47000000` through `chr8,47000140` (20bp apart, one per read, per the fixture's placement pattern), each ending `,+,30M70S,60,0;`, i.e. every read independently points to a partner position tightly clustered around chr8:47,000,000, visible directly in the screenshot. IGV logged non-fatal `WARNING [SupplementaryAlignment] ... couldn't be sorted unambiguously` for 3 of the 8 reads (`split_read_5/6/7`) — a display-grouping quirk when multiple supplementary alignments share a position, not a data or scoring problem; all 8 reads are still visible and correctly labeled in the image. |
| `read_depth` | chr1:1,047,000–1,053,000 (±3000bp, fixed scale 0–18) | `screenshots/synthetic_translocation_read_depth.png` | Two adjacent narrow coverage spikes (one orange/blue, one red) rather than a dip — this is the known fixture artifact already documented in `DEMO_END_TO_END.md` (the read-free flanking gaps that isolate other signal read groups trip the deletion-ratio check and register as depth peaks here, not evidence of an actual deletion at a translocation locus). |

## 2. Negative control — HCC1143 chr21 real-data locus

**Date:** 2026-08-11
**BAM:** `stage1_igv_assistant/data/bam/HCC1143.normal.21.19M-20M.bam`
**Locus:** chr21:19,089,694 — same locus already characterized in `RESULTS_HCC1143.md` (labeled "translocation" in the source workshop's tutorial script, but no real signal found there).

`applicable_layers` on this BAM returns only 3 layers (`discordant_pairs`, `soft_clipped_reads`, `read_depth`) — `split_reads` is correctly excluded, consistent with the whole-BAM check in `RESULTS_HCC1143.md` §4 (0 supplementary-alignment reads across all 572,731 reads on chr21; this is a 2018-era alignment with no chimeric/SA-tag output). `bam_stats_at_locus` on the full documented region (chr21:19,089,694–19,095,362) reproduced the figures already on file: 4,097 reads, mean depth 62.33x, mean MAPQ 31.09.

| Layer | Window used | Screenshot | Agrees with RESULTS_HCC1143.md? |
|---|---|---|---|
| `discordant_pairs` | chr21:19,088,194–19,091,194 (±1500bp) | `screenshots/hcc1143_chr21_discordant_pairs.png` | Yes — `UNEXPECTED_PAIR` coloring shows scattered single-colored pairs, no dominant partner-chromosome color cluster, matching the "5 discordant reads, one each to chr20/chr5/chr19/chr12/chr18, no partner clustering" finding. |
| `soft_clipped_reads` | chr21:19,089,544–19,089,844 (±150bp, dominant side auto-detected) | `screenshots/hcc1143_chr21_soft_clipped_reads.png` | Yes — visible staircase clip patterns at several distinct positions rather than one shared boundary, matching "22/287 soft-clipped, but only 1 read at the consensus position — no real pileup." `clip_side_determination` reported 7 right-clipped vs. 6 left-clipped reads in this window (near-even split; sorted by `RIGHT_CLIP`). |
| `read_depth` | chr21:19,086,694–19,092,694 (±3000bp, fixed scale 0–230) | `screenshots/hcc1143_chr21_read_depth.png` | Yes — coverage track is essentially flat across the window with no dip, consistent with this locus never being flagged as a depth-supported deletion anywhere in prior validation. |
| `split_reads` | — | skipped | Correctly skipped — excluded via `applicable_layers` (no SA tags on this BAM). |

## Bottom line

`evidence_panel` runs end-to-end on both a synthetic positive-control BAM and a real, network-free negative-control BAM, and in both cases the per-layer screenshots are visually consistent with evidence already validated elsewhere (the fixture's own docstring guarantees for the synthetic case; `RESULTS_HCC1143.md`'s numeric findings for the real case). The positive control additionally confirms the one thing only this tool's split-read layer can show: individual SA-tag partner-locus labels rendered directly in the image, each independently pointing to a partner position tightly clustered around chr8:47,000,000. The negative control shows the tool correctly reporting *why* a locus is weak (scattered discordant partners, scattered clip positions, flat depth) rather than just asserting it, and correctly skipping a structurally-inapplicable layer. Together they validate the tool's plumbing (per-layer window selection, layer skipping, dominant-clip-side detection, SA-tag grouping) on both signal-present and signal-absent data.
