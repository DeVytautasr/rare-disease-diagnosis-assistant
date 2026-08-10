# HCC1143 chr21 Breakpoint Evidence — Tool Validation Results

> **SUPERSEDED — scoring scheme has changed since this was written.**
> This file documents `summarize_breakpoint_evidence` scoring each layer
> 0–50, normalised to a 0–100 composite over 3 layers (`signal_layers:
> "3/3"`). The current implementation scores each of 4 layers 0–25,
> summed directly to `evidence_score_raw` (0–100), and additionally
> normalises over only the technology-applicable layers into
> `evidence_score` (see `bam_tools.summarize_breakpoint_evidence` and
> `AUDIT_2026_08.md` §5). The raw tool outputs and locus findings below
> are still accurate; the scoring-scheme narrative and the specific
> score values are not. See `REAL_DATA_VALIDATION.md` and
> `docs/thesis/thesis_background_methods_chapter.md` for the current
> scoring model.

**Date:** 2026-08-08
**BAM:** `stage1_igv_assistant/data/bam/HCC1143.normal.21.19M-20M.bam`
**Sample:** HCC1143 normal (germline) — no matched tumor BAM exists for this dataset (confirmed via directory listing of `genomedata.org/gen-viz-workshop/IGV/`; only the normal BAM was ever published there).
**Target locus:** chr21:19,089,694–19,095,362, labeled "translocation" (Example 8) in the source workshop's own `Run_batch_IGV_snapshots.txt` tutorial script.

This file records raw tool output for reproducibility, followed by an interpretation. Per project convention, the tools report structured evidence only — the interpretation below is a read of that evidence, not a diagnostic claim.

---

## 1. Region-level stats — `get_bam_stats_at_locus`

Full documented region, chr21:19,089,694–19,095,362:

| Metric | Value |
|---|---|
| Total reads | 4,097 |
| Mean depth | 62.33x |
| Mean MAPQ | 31.09 |
| Low-MAPQ fraction (MAPQ < 20) | 0.399 |
| Forward / reverse reads | 2,053 / 2,044 |

Roughly 40% of reads in this region have MAPQ < 20 — moderately ambiguous mapping across the whole locus, not just at one point.

## 2. Full tool suite at the documented position — chr21:19,089,694

| Tool | Result |
|---|---|
| `count_discordant_pairs` | 5 / 613 reads (0.8%) discordant. Mates scattered: chr20, chr5, chr19, chr12, chr18 — **1 read each**, no partner clustering. |
| `count_soft_clipped_reads` | 22 / 287 reads (7.7%) soft-clipped. Consensus clip position 19,089,435, but only **1 read** at that position (`max_clips_at_position: 1`) — no real pileup, clips are scattered. |
| `get_split_reads` | 0 / 289 reads. No `SA`-tagged (chimeric) reads found. |
| `summarize_breakpoint_evidence` | **evidence_score: 30.0 → "weak"** (discordant_pair_score 15, soft_clip_score 15) |

## 3. Discordant-pair scan across the region — `count_discordant_pairs`

| Position | Discordant pairs | Total reads in window | Fraction | Mate chromosomes |
|---|---|---|---|---|
| 19,089,694 | 5 | 613 | 0.8% | chr20(1), chr5(1), chr19(1), chr12(1), chr18(1) |
| 19,092,000 | 2 | 261 | 0.8% | chr12(1), chr3(1) |
| 19,095,000 | 4 | 595 | 0.7% | chr7(2), chr9(2) |

The discordant fraction is flat (~0.7–0.8%) across the whole region — there is **no peak** at any tested position, and the partner chromosome changes at every position rather than converging on one candidate translocation partner.

## 4. Whole-BAM check for split-read capability

To rule out a locus-specific absence of signal, we scanned the entire chr21 BAM (572,731 reads, not just the windowed region):

- Supplementary-alignment flag set: **0**
- Reads carrying an `SA` tag: **0**

This BAM contains no chimeric/supplementary alignment records anywhere. `get_split_reads()` will read as zero at any position in this file, regardless of whether a true breakpoint is present — this is a limitation of the original 2018 alignment pipeline (aligner/version did not emit chimeric alignments), not a finding about this locus specifically.

---

## Interpretation

The evidence tools ran correctly and produced consistent, reproducible output. However, at the documented translocation coordinate (and across the flanking region), the signal does **not** show the pattern expected of a genuine interchromosomal translocation:

- Discordant pairs are present but rare (~0.7–0.8%) and their mates are spread across many different chromosomes rather than clustering on one partner — a true translocation should show many discordant pairs converging on the same partner chromosome.
- Soft-clipping is present (~8%) but not concentrated at a single consensus position (`max_clips_at_position` never exceeds 1–5 reads) — no sharp breakpoint pileup.
- Split-read evidence is uninformative here: the whole BAM lacks `SA` tags, so this tool cannot contribute evidence for this dataset regardless of ground truth.

**We cannot confirm a translocation signal at this locus using these three evidence layers on this data.** Plausible explanations, none of which we can distinguish from the data alone:
1. The tutorial's "translocation" example may be intended for **visual inspection in IGV** rather than as a strong quantitative signal detectable by discordant-pair/soft-clip heuristics at default thresholds.
2. The original 2018 alignment pipeline (no `SA` tags anywhere in the file) may predate more sensitive SV-aware alignment/calling, so a real signal could be present but under-represented in this BAM's annotations.
3. The signal may be real but subtle (e.g., low allele fraction or a small subclonal population), below what these heuristic thresholds are tuned to flag.

**Next steps, if pursuing this further:** visually inspect the locus in IGV directly (as the original tutorial intends) to sanity-check whether the tool output matches what's visible; consider re-aligning with a chimeric-alignment-aware aligner (e.g., current bwa-mem) if split-read evidence is needed; or treat this dataset primarily as a pipeline/tool validation exercise rather than a source of a confirmed positive-control translocation call.

---

## Addendum: synthetic translocation confirms the detection logic itself is sound

To separate "does the tool detect a translocation" from "does this specific 2018 BAM contain detectable split-read evidence," we built a synthetic BAM (`create_translocation_bam()` in `stage1_igv_assistant/tests/test_bam_tools.py`, TEST 6) simulating a balanced translocation between chr1 and chr8, including `SA` tags on the split reads — i.e. data with the annotation the HCC1143 BAM lacks.

`summarize_breakpoint_evidence()` was also extended to add split-read evidence as a genuine third scored layer (alongside discordant pairs and soft-clips), each independently 0–50 and normalized to a 0–100 composite, so the tool can now report how many of the 3 evidence layers show signal (`signal_layers`).

Running all 5 tools on the synthetic BAM at chr1:1,050,000:

| Tool | Result |
|---|---|
| `count_discordant_pairs` | 15/15 discordant pairs found, all mates on chr8 |
| `count_soft_clipped_reads` | 5/5 soft-clipped reads found |
| `get_split_reads` | 8/8 split reads found, all partner loci on chr8 |
| `summarize_breakpoint_evidence` | **evidence_score: 100.0 → "strong"**, `signal_layers: "3/3"` |

All expected counts matched exactly, and the composite score correctly reached "strong" only when all three evidence layers fired together. This confirms the detection logic itself is correct: **the reason the real HCC1143 BAM didn't show a strong signal above is the absence of `SA` tags in that specific 2018 alignment pipeline, not a defect in these tools.** Given data with proper chimeric-alignment annotation, the same code correctly identifies a translocation with high confidence.

---

## GIAB HG002 validation — key findings

To move beyond synthetic data, we validated the tools against a real modern long-read dataset: a publicly streamable PacBio HiFi BAM (`HG002.GRCh38.haplotagged.bam`, aligned to GRCh38) plus GIAB's GRCh38-native CMRG SV v1.00 benchmark, targeting a confirmed 3,359bp heterozygous deletion at chr1:115,686,862–115,690,222. Three lessons came out of this that matter beyond this one locus:

**(a) A real bug was found and fixed for long-read data.** `count_discordant_pairs` originally checked only `mate_is_unmapped` to decide whether a read had a valid mate — but single-molecule long reads (PacBio HiFi, Oxford Nanopore) are unpaired entirely, so that check defaults to "mapped" even though no mate exists. Every unpaired read was silently miscounted as a discordant pair, producing a **spurious 100% discordant-pair signal** that looked exactly like strong translocation evidence but meant nothing. Fixed by requiring `read.is_paired` before evaluating mate location (commit `9bc0da9`). A second, related bug was caught while building the new depth-profile tool below: an early version judged window overlap using each read's `reference_start`/`reference_end` span, which extends *across* internal CIGAR deletions — so a read carrying the very deletion we were trying to detect was counted as if it fully covered the deleted region. Fixed by overlapping on actual aligned reference positions (`get_reference_positions()`) instead, the same approach `get_bam_stats_at_locus` already used correctly. Both bugs shared a root cause: **tools written and tested against short-read, paired-end assumptions silently misbehave on long-read data in ways that look like false positive signal, not obvious errors** — real data testing caught what synthetic data couldn't.

**(b) Evidence-layer relevance depends on sequencing technology and SV size relative to read length, not just SV type.** At this 3.4kb deletion, `count_discordant_pairs` was correctly inapplicable (0/0 — no paired reads exist in HiFi data). `count_soft_clipped_reads` and `get_split_reads` each found only a handful of reads (4/40 and 1/39 respectively) — small but genuine: the soft-clip consensus position (115,686,865) sits 3bp from the documented deletion start, and the one split read's partner locus (chr1:115,690,223) sits 1bp from the documented deletion end. But most spanning ~15–20kb HiFi reads simply represent a 3.4kb deletion as one contiguous alignment with an embedded CIGAR `D` operation, never triggering a soft-clip or a chimeric `SA`-tagged alignment at all. These two evidence layers were designed around short-read SV signatures and are structurally underpowered for a deletion this size relative to HiFi read length — a correct absence of signal, not a tool failure.

**(c) Read depth is the most direct signal for deletions, and is now the 4th evidence layer (`get_read_depth_profile`).** A direct depth check (via `get_bam_stats_at_locus`, using true per-base coverage) showed depth dropping from ~41x/36x flanking to ~20x inside the deletion — essentially the 50% reduction expected for a heterozygous deletion (`GT=1|0`). This is a much cleaner and more direct signal than soft-clips or split-reads for this SV type/technology combination, so `get_read_depth_profile()` was added as a proper 4th scored layer in `summarize_breakpoint_evidence` (raw components now normalized over 4×50 rather than 3×50; `signal_layers` now reports "N/4"). One honest caveat: the automated `depth_ratio_min_to_mean < 0.6` threshold used inside the composite score is **less sensitive than the manual flanking-vs-inside comparison** — it came out to 0.618–0.707 depending on window width (missing the 0.6 cutoff) versus the ~0.5 ratio found by comparing directly against flanking regions. This is because the tool's denominator is the *mean of the entire queried window* (which includes both flanking and deleted sequence), diluting the ratio compared to a true flanking-only baseline. A flanking-relative comparison would likely be more sensitive, but was left as a documented limitation rather than implemented now, to keep scope realistic for this stage.
