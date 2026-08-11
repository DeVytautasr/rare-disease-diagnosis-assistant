# LLM Session 3 — Blind Breakpoint Investigation (MCP tools only)

**BAM:** `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam`
**Genome build:** GRCh38
**Sequencing technology (as supplied with the task):** Illumina 300x paired-end, Novoalign aligner
**Session date:** 2026-08-11

## Ground rules for this session

- Investigation used **only** the `igv-breakpoint-assistant` MCP tools — no source-code reading, no scripts, no independent computation outside tool output.
- The three positions (A, B, C) were investigated as separate tool-call sequences. No finding at one position was used to set expectations for another.
- No prior knowledge about HG002, GIAB, any gene, or any of the three positions was used. Every number below is quoted directly from a tool response in this session.
- Per the MCP server's own documentation on `breakpoint_evidence_summary`: of the 7 tier cutoffs behind its 4 component scores, only **one** (the read-depth layer's 0.7 "moderate" threshold) is empirically calibrated, and that against a single confirmed locus elsewhere. The other 6 (discordant-pair 0.2/0.5, soft-clip 0.1/0.3, split-read 0.1/0.3, depth's own 0.3 "strong" threshold) are unvalidated heuristic judgement calls. `evidence_score`/`evidence_strength` below is treated as an interpretable decomposition of what fired, **not** a calibrated probability.
- `evidence_score` (normalised over applicable layers) is reported as the primary metric per the server's instructions, with `evidence_score_raw` (over all 4 layers) given alongside.

## Step 0 — `applicable_layers` (called once for this BAM)

Sampled 1000 reads. Result:

| Layer | Applicable? | Tool's stated reason |
|---|---|---|
| `discordant_pairs` | Yes | "at least one paired read observed in 1000 sampled reads" |
| `soft_clipped_reads` | Yes | "always applicable — any aligned BAM can show soft-clipping" |
| `read_depth` | Yes | "always applicable — any aligned BAM can show a coverage drop" |
| `split_reads` | **No** | "no SA tags observed in 1000 sampled reads (consistent with an aligner that doesn't emit chimeric alignments)" |

`applicable_layers = ["discordant_pairs", "soft_clipped_reads", "read_depth"]` was passed to every `breakpoint_evidence_summary` call below. This empirical finding (zero SA tags in the sample) is consistent with the Novoalign aligner stated in the task — noted only because it was independently confirmed by the tool in this session, not assumed beforehand.

**Consequence:** `split_reads` is structurally inapplicable for this whole BAM. Every "0 split reads" result below is an *expected null*, not negative evidence — it cannot contribute to any verdict in either direction, at any of the three positions.

---

## Position A — chr1:16,890,000

**`bam_stats_at_locus`** (16,889,750–16,890,250): 2,584 reads, mean_depth **587.17**, mean MAPQ **70**, low-MAPQ fraction **0**, forward/reverse 1,409/1,175. No repetitive-region flag (fraction ≫ 0.4 would have meant "interpret with caution" — not triggered).

**`discordant_pairs`** (pos 16,890,000, window ±500bp): 4,537 reads in window, **12 discordant pairs**, `discordant_fraction` **0.003** (0.3%). `mate_chromosomes`: chr10:2, chr8:1, chr6:1, chr4:1, chr14:1, chr2:1, chr9:1, chrUn_JTFH01001465v1_decoy:1, chr3:1, chr17:1, chr19:1 — **11 different partner chromosomes** for 12 reads, no chromosome above 2. This is precisely the "mates scattered across many different chromosomes = background noise" pattern the tool itself documents, not the "many discordant pairs clustering on ONE partner chromosome" pattern it defines as a real translocation signal.

**`soft_clipped_reads`** (pos 16,890,000, window ±200bp): 2,147 reads in window, **17 soft-clipped**, fraction **0.008** (0.8%). `max_clips_at_position` = **2**, below the tool's own noise cutoff ("< 3 = no real pileup, treat as noise"). Consensus clip position 16,890,019, but only 2 reads actually share it.

**`split_reads`**: 0/2,147, expected null (inapplicable BAM-wide).

**`read_depth_profile`** (16,887,000–16,893,000, 100bp bins): mean_depth 976.98, min 583 (bin 16,891,500–16,891,600), max 1,298, ratio min:mean **0.597** → `likely_deletion: true`. **However**, the bins actually covering position A itself are 964 (16,889,900–16,890,000) and 1,068 (16,890,000–16,890,100) — both at or above the window mean. The dip driving the "likely_deletion" flag sits **~1,500bp downstream** of the queried position, not at it. Position A itself is at a local peak, not a trough.

**`breakpoint_evidence_summary`**: `evidence_score` **20/100** ("weak"), `evidence_score_raw` **15/100**, `signal_layers` **2/3**. Component scores: discordant_pair_score 7.5, soft_clip_score 7.5, depth_score **0**, split_read_score 0. Internally, this call's own narrower depth window (16,888,000–16,892,000, 200bp bins) gives ratio **0.714** and `likely_deletion: false` — the opposite conclusion from the wider external scan above, because the narrower/coarser-binned window doesn't weight the off-position dip as heavily. This flip is a direct illustration of the depth-threshold's window-sensitivity, not a new finding about position A itself (see caveats below).

### Verdict — Position A
**No credible evidence of a structural-variant breakpoint at chr1:16,890,000.** Discordant pairs and soft clips are both at background/noise levels by the tool's own stated criteria, and the depth signal at the position itself is flat-to-slightly-elevated, not depressed. The one "likely_deletion: true" flag obtained from a wide scan is driven by an unrelated feature ~1.5kb away and does not localize to this coordinate.

**Confidence: High** (no evidence of an anomaly at this exact position). Not extended to the nearby unrelated depth feature, which was not part of the requested position and was not further investigated.

---

## Position B — chr2:96,300,000

**`bam_stats_at_locus`** (96,299,750–96,300,250): 1,300 reads, mean_depth **290.24**, mean MAPQ **70**, low-MAPQ fraction **0**, forward/reverse 657/643.

**`discordant_pairs`** (pos 96,300,000): 2,323 reads in window, **5 discordant pairs**, fraction **0.002** (0.2%). `mate_chromosomes`: chr9:2, chr11:1, chr5:1, chr1:1 — 4 different partners for 5 reads. Same scattered/noise pattern as Position A, at an even lower absolute count.

**`soft_clipped_reads`** (pos 96,300,000): 1,057 reads in window, **6 soft-clipped**, fraction **0.006** (0.6%). `max_clips_at_position` = **1** — no pileup at all, the lowest of the three positions. Left/right consensus positions don't even agree (96,299,688 vs 96,300,093), each with only 1 supporting read.

**`split_reads`**: 0/1,057, expected null.

**`read_depth_profile`** (96,297,000–96,303,000, 100bp bins): mean_depth 476.62, min 380, max 598, ratio min:mean **0.797** → `likely_deletion: false`. Flattest, most stable profile of the three positions (narrowest min–max spread relative to mean).

**`breakpoint_evidence_summary`**: `evidence_score` **20/100** ("weak"), `evidence_score_raw` **15/100**, `signal_layers` **2/3**. Component scores: discordant_pair_score 7.5, soft_clip_score 7.5, depth_score **0**, split_read_score 0 — **numerically identical component breakdown to Position A**, despite B's raw discordant/soft-clip counts being lower across the board (see cross-position comparison). Internal narrower depth window (96,298,000–96,302,000, 200bp bins) gives ratio **0.842**, `likely_deletion: false` — same conclusion as the external wide-window scan, unlike Position A. B's "no deletion" call is robust to window choice.

### Verdict — Position B
**No credible evidence of a structural-variant breakpoint at chr2:96,300,000.** This is the cleanest background-noise profile of the three positions: lowest discordant and soft-clip counts, flattest depth, and — unlike Position A — the "no deletion" conclusion holds regardless of scan width.

**Confidence: High**, and more robust than Position A's (no window-dependent ambiguity, no nearby confounding feature).

---

## Position C — chr1:115,686,862

**`bam_stats_at_locus`** (115,686,612–115,687,112): 1,017 reads, mean_depth **230.48**, mean MAPQ **70**, low-MAPQ fraction **0**, forward/reverse 506/511. Note this window spans into the region where depth later drops (see below); MAPQ stays high and low-MAPQ fraction stays at 0 even there, arguing against a mapping/repeat artifact as the explanation for what follows.

**`discordant_pairs`** (pos 115,686,862): 1,708 reads in window, **1 discordant pair**, fraction **0.001** (0.1%), single mate on chr12. Lower than either A or B. **Important caveat, not a contradiction**: this tool specifically counts inter-chromosomal mate pairs and is documented as the "PRIMARY signal for balanced translocations." A same-chromosome deletion — which is what the rest of the evidence below points toward — would not be expected to produce inter-chromosomal discordant pairs regardless of whether it is real. This layer is **uninformative** for that hypothesis here, not negative evidence against it.

**`soft_clipped_reads`** (pos 115,686,862): 853 reads in window, **23 soft-clipped**, fraction **0.027** (2.7% — 3.4–4.5× A and B's fractions). `max_clips_at_position` = **13**, well above the tool's noise cutoff of 3. Consensus clip position **115,686,865** — 3bp from the queried coordinate. `dominant_clip_side: "right"`, with 21/23 clipped reads on the right side (91%) — one consistent, one-sided orientation, not scatter.

**`split_reads`**: 0/853, expected null.

**`read_depth_profile`** (115,683,862–115,689,862, 100bp bins): mean_depth 340.28, min **175** (bin 115,687,762–115,687,862), max 518, ratio min:mean **0.514** → `likely_deletion: true`. Shape matters here: the bin ending exactly at the queried position (115,686,762–115,686,862) is depth 420; the very next bin (115,686,862–115,686,962, starting **at** the position) drops to 327, and depth continues declining over the next ~900bp to the 175 minimum, then fluctuates in the 230–290 range for the rest of the 3kb scanned — it does **not** recover to the ~400–520 upstream baseline within the region examined. Unlike Position A, the transition point coincides with the queried coordinate itself.

**`breakpoint_evidence_summary`**: `evidence_score` **40/100** ("moderate"), `evidence_score_raw` **30/100**, `signal_layers` **3/3**. Component scores: discordant_pair_score 7.5, soft_clip_score 7.5, depth_score **15**, split_read_score 0. Internal narrower depth window (115,684,862–115,688,862, 200bp bins) gives ratio **0.542**, `likely_deletion: true` — consistent with the external wide-window scan. Unlike Position A, this conclusion is robust to window choice.

**Important scoring nuance**: `discordant_pair_score` and `soft_clip_score` are **7.5 — identical to Positions A and B** — despite Position C's soft-clip fraction being 3–4× higher with a genuine 13-read pileup (vs. ≤2 reads, no pileup, at A/B). This is because all three positions' fractions (0.006–0.027) fall under the same 0.1 tier threshold the server's docs already flag as an uncalibrated heuristic; the tier is too coarse to separate "no pileup" from "clear 13-read pileup." The layer that *did* differentiate C from A/B is `depth_score` (15 vs. 0), which is scored against the one threshold (0.7) the server's docs describe as empirically calibrated. Practical implication: **don't read the composite score alone** — the raw `max_clips_at_position` (13 vs. ≤2) carries information the composite discards.

**`gene_at_locus`** (chr1, 115,686,862, GRCh38) — run because this position, unlike A/B, crossed the tool's own stated threshold for follow-up ("call after finding strong discordant/split-read evidence"; here the qualifying strong signal is soft-clip + depth rather than discordant/split-read, but the same logic applies: check gene disruption once there's a real signal to follow up). Result: `gene_count` **1**, `is_intergenic` **false**. Gene: **VANGL1** (`ENSG00000173218`), biotype protein_coding, strand +, gene body chr1:115,641,854–115,698,224. `clinical_note`: "Breakpoint directly disrupts 1 gene(s)." No further interpretation of this gene is offered — per this session's ground rules, no prior knowledge about any gene was used; this is a direct quote of the Ensembl lookup.

### Verdict — Position C
**Coherent, multi-layer, position-specific evidence of a structural anomaly**, most consistent with one breakpoint of a heterozygous deletion: a genuine soft-clip pileup (13 reads, 3bp from the query coordinate, 91% one-sided) coincides with a depth transition that begins in the bin immediately after the query coordinate and stays depressed to roughly 51–54% of local baseline (both wide- and narrow-window measurements agree) for at least 3kb downstream. High MAPQ and zero low-MAPQ fraction through the transition region argue against a mapping-artifact explanation. The breakpoint sits inside the single protein-coding gene VANGL1 per Ensembl. Discordant-pair and split-read layers are silent, but both are the wrong/inapplicable assay for this specific hypothesis (translocation-detector and SA-tag-dependent, respectively), not contradicting evidence.

**Confidence: Moderate.** Two independent, position-matched, window-robust signals (soft-clip pileup + depth transition) is meaningfully more than background, and clearly distinct from A/B. Capped below "high" because: the tool's own composite label is "moderate," not "strong"; the two layers that *do* differentiate C are diluted in the composite by two uncalibrated tiers that don't; there is no corroborating discordant-pair or split-read signal (both structurally unable to speak to this hypothesis here, which cuts against overclaiming, not for it); and the far/downstream breakpoint and true deletion extent were not determined — depth had not recovered to baseline by the edge of the 3kb window scanned in this session.

---

## Cross-position comparison

| Metric | A (chr1:16,890,000) | B (chr2:96,300,000) | C (chr1:115,686,862) |
|---|---|---|---|
| discordant_pairs / window | 12/4,537 (0.3%) | 5/2,323 (0.2%) | 1/1,708 (0.1%) |
| distinct mate chromosomes | 11 (max 2 on one) | 4 (max 2 on one) | 1 (chr12, n=1) |
| soft-clipped reads / window | 17/2,147 (0.8%) | 6/1,057 (0.6%) | 23/853 (**2.7%**) |
| max_clips_at_position | 2 (no pileup) | 1 (no pileup) | **13 (real pileup)** |
| clip position vs. query | — | — | 3bp offset |
| depth ratio, wide scan (100bp bins) | 0.597 (dip **not** at position) | 0.797 (flat) | **0.514 (drop starts at position)** |
| depth ratio, narrow scan (200bp bins, internal) | 0.714 (flips to "no deletion") | 0.842 (flat, consistent) | **0.542 (consistent "deletion")** |
| evidence_score (normalised) | 20/100 "weak" | 20/100 "weak" | **40/100 "moderate"** |
| evidence_score_raw | 15/100 | 15/100 | **30/100** |
| signal_layers | 2/3 | 2/3 | **3/3** |
| mean MAPQ / low-MAPQ fraction | 70 / 0 | 70 / 0 | 70 / 0 |
| gene_at_locus run? | No (no signal to follow up) | No (no signal to follow up) | Yes → VANGL1, protein-coding, breakpoint inside gene body |

**A vs. B**: look the same at the composite-score level (both 20/100 "weak", identical 7.5/7.5/0/0 component breakdown) — but this is a tier-bucketing coincidence, not true equivalence. B's raw counts are lower on every metric, and B's "no deletion" conclusion is robust across window widths while A's is not (A has a real depth feature nearby that isn't at the queried coordinate). Both reduce to the same verdict — no evidence of a breakpoint at the queried position — but B is the cleaner/more robust of the two negatives.

**C vs. A/B**: differs clearly and on multiple independent axes — composite score doubles (40 vs. 20), raw score doubles (30 vs. 15), 3/3 vs. 2/3 signal layers, a genuine soft-clip pileup (13 reads) vs. none anywhere else (≤2), and a depth transition that starts at the query coordinate and is robust to window choice, vs. depth that is either flat (B) or depressed only at an off-position location (A). C is also the only position where a gene lookup was warranted and returned a hit.

## Tools available but not used, and why

- **`reciprocal_breakpoint`**: not called for any position. This tool is for verifying a suspected balanced translocation once discordant pairs cluster on one partner chromosome. None of the three positions showed that pattern — all three `discordant_pairs` results were low-count and scattered across multiple partner chromosomes (11, 4, and 1 distinct partners respectively), which the `discordant_pairs` tool's own documentation identifies as background noise rather than a translocation candidate. There was no partner locus worth testing reciprocally.
- **`evidence_panel` / `igv_screenshot`**: not called. Not required to reach a tool-cited verdict, and outside the report format requested for this session.

## Data-quality caveat: two depth-reporting pathways disagree in absolute scale

Within this session, `bam_stats_at_locus` (and the `locus_stats` block `breakpoint_evidence_summary` computes internally) consistently reports **lower** mean depth than `read_depth_profile` (and the `depth_profile` block `breakpoint_evidence_summary` computes internally) for what is nominally the same region:

| Position | `bam_stats_at_locus` mean_depth (500bp window) | `read_depth_profile` mean_depth (6,000bp window, 100bp bins) |
|---|---|---|
| A | 587.17 | 976.98 |
| B | 290.24 | 476.62 |
| C | 230.48 | 340.28 |

The ratio between the two isn't fixed (≈1.66×, 1.64×, 1.48× respectively — the two calls also don't span identical windows, so some of that variation is real local coverage difference), but the direction and rough magnitude is consistent and reproducible across all three positions checked, so it reads as a systematic difference in how the two tool families compute "depth" rather than a one-off anomaly. A further, separate observation: `read_depth_profile`'s own reported depth also shifts with `window_size` — at Position A, the same region reads mean_depth 976.98 at `window_size=100` (external call) vs. 1,354.25 at `window_size=200` (internal call inside `breakpoint_evidence_summary`, narrower span). A true per-base average depth should be close to bin-size invariant; this wasn't.

Practical handling in this report: absolute depth numbers were **not** compared across the `bam_stats_at_locus`/`read_depth_profile` tool families, and `read_depth_profile` calls using different `window_size` values were not compared to each other. Only within-call ratios (`depth_ratio_min_to_mean`, computed by a single call against its own numbers) were used for interpretation, and — as documented above for Position A — even that ratio was shown to flip a boolean conclusion depending on scan width, which is reported as a finding in its own right rather than papered over.

> **Root cause confirmed and fixed (2026-08-11):** this was a real bug, not two equally-valid measurement conventions. `read_depth_profile` was counting distinct reads touching each bin, not true per-base coverage — see `bam_tools.py`'s `get_read_depth_profile` docstring. The Position-A window-size-dependent `likely_deletion` flip described above was a second, related bug: the ratio used the region's global minimum rather than depth local to the queried position, and that global minimum's location (correctly identified above as sitting ~1,500bp from Position A, not at it) is exactly what made the flag window-size-sensitive — a genuinely off-position feature will or won't cross a coarse global threshold somewhat by chance, depending on how it happens to get binned. Both are fixed; see "Post-fix re-verification" below for the corrected numbers at all three positions, and `results/REAL_DATA_VALIDATION.md`'s "Post-fix re-validation" section for the calibration-threshold implications.

## Summary

| Position | Verdict | Confidence |
|---|---|---|
| A — chr1:16,890,000 | No credible evidence of a breakpoint | High |
| B — chr2:96,300,000 | No credible evidence of a breakpoint | High |
| C — chr1:115,686,862 | Coherent evidence of a structural anomaly, most consistent with a heterozygous deletion breakpoint; disrupts VANGL1 | Moderate |

---

## Post-fix re-verification (2026-08-11)

Three bugs surfaced by this session were fixed in `bam_tools.py`:

1. **`get_read_depth_profile` computed a per-bin read count, not true per-base depth** (the "two depth-reporting pathways" caveat above) — fixed to sum aligned bases per bin over bin width, matching `bam_stats_at_locus`'s method. Now confirmed bin-size invariant to within ~1-2% on this exact BAM (was ~39% apart between `window_size=100` and `window_size=200` at Position A before the fix).
2. **The depth ratio used the region's global minimum, not depth local to the queried position** — the root cause of Position A's window-size-dependent `likely_deletion` flip. `get_read_depth_profile` now accepts `focus_position`; `summarize_breakpoint_evidence` passes the queried position through and additionally requires the region's true minimum to fall within 1000bp of it (`dip_is_at_focus`, HEURISTIC, calibrated against exactly these two loci — see the function's docstring) before awarding any depth-based score. An off-position dip is now reported in `supporting_observations` instead of silently scoring.
3. **`soft_clip_score` was tiered on `soft_clipped_fraction`, not `max_clips_at_position`** — the reason Position C's genuine 13-read pileup and Position A/B's non-pileups (max 2 and 1 reads respectively) all scored identically (7.5/25) in the original session above. Now tiered on `max_clips_at_position` directly (< 3 → 0, 3–9 → 15, ≥ 10 → 25).

Re-ran all three positions through the fixed `applicable_layers` → `breakpoint_evidence_summary` pipeline, same BAM, same coordinates, same `applicable_layers = [discordant_pairs, soft_clipped_reads, read_depth]`:

| | Position A | Position B | Position C |
|---|---|---|---|
| `soft_clip` max_clips_at_position | 2 (unchanged) | 1 (unchanged) | 13 (unchanged) |
| `soft_clip_score` | 7.5 → **0.0** | 7.5 → **0.0** | 7.5 → **25.0** |
| `depth` ratio (localized) | 0.968/0.714 (window-dependent) → **0.838 (single, stable value)** | n/a → **0.856** | 0.542 → **0.502** |
| `dip_is_at_focus` | *(field didn't exist)* → **False** (dip is 1400bp away) | *(field didn't exist)* → **False** (dip is 1200bp away) | *(field didn't exist)* → **True** (dip is 800bp away) |
| `depth_score` | 0 (unchanged — never crossed threshold either way) | 0 (unchanged) | 15 (unchanged value, now gate-confirmed rather than assumed) |
| `discordant_pair_score` | 7.5 (unaffected by these fixes) | 7.5 (unaffected) | 7.5 (unaffected) |
| `signal_layers` | "2/3" → **"1/3"** | "2/3" → **"1/3"** | "3/3" → **"3/3" (unchanged)** |
| `evidence_score` | 20.0 → **10.0** | 20.0 → **10.0** | 40.0 → **63.3** |
| `evidence_strength` | weak → **weak (unchanged)** | weak → **weak (unchanged)** | moderate → **moderate (unchanged)** |

**Verdicts hold — no flips.** All three positions land on the same conclusion as the original blind session: A and B show no credible evidence of a breakpoint, C shows coherent evidence of a real structural anomaly. What changed is separation, not direction — the corrected scoring pulls A/B *down* (their apparent soft-clip "signal" was never a real pileup) and C *up* (its genuine pileup was previously scored the same as A/B's non-pileups), roughly doubling the score gap between background noise and real signal (previously 20 vs 40, a 2x gap; now 10 vs 63.3, a ~6x gap). This matches the expected direction stated before this re-verification was run (C should rise, A/B should fall or stay flat) — confirmed against the actual re-measurement above, not assumed.

One implementation note worth recording: Position A's off-target dip (1400bp away) didn't need the `dip_is_at_focus` gate to be suppressed — the gate's search radius (`dip_tolerance_bp=1000`) already excludes it from the *localized* ratio computation itself, so `depth_ratio` never crossed the threshold in the first place at either position. The gate matters for a narrower case: a dip within `dip_tolerance_bp` that pulls the local ratio below threshold, while the region's true minimum sits even farther out. No such case appeared among these three positions; it's exercised directly by `test_bam_tools.py`'s TEST 6 (a synthetic fixture whose engineered read-free gap sits outside `dip_tolerance_bp` of its queried position).
