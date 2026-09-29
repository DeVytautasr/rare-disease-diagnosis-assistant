# Figure map — the thesis and the conference abstract against the records

Written 2026-09-25, after the implants were rebuilt, the evidence chain run on
them and the ceiling experiment rerun. For every numeric claim: where it
appears, the value it gives, the value the records give now, the record and key
the new value comes from, and a status.

## What was read

| Document | Found | Modified | Words | Read how |
|---|---|---|---|---|
| Conference abstract, `Rimas_tezes_DI_medicinoje_LT.docx` | `C:\Users\rvyta\Desktop` | 2026-09-24 12:06:05 (+03:00), 15,473 bytes | 356 | read only, text extracted from `word/document.xml`; nothing written under `/mnt/c` |
| Thesis, `Rimas_MSc_Thesis.docx` | **not found** — not on the Desktop, in Downloads or Documents, anywhere under `C:\Users` (searched to depth 7 for every `.docx`), or anywhere on `C:` outside the Windows, Program Files, ProgramData and Recycle Bin folders (searched by name) | — | — | — |

Because the thesis file was not available, its claims could not be extracted.
This map therefore covers: **(A)** every numeric claim in the abstract; **(B)**
the Phase 5, 6 and 10 figures that the task statement lists as surviving only as
text in the thesis and the abstract (their location in the thesis could not be
checked); **(C)** the patient figures, from the committed patient records. The
repository's own chapter draft (`docs/thesis/thesis_background_methods_chapter.md`,
2026-08-31, 15,101 words) is an earlier, different document: it predates the
implant and ceiling experiments, and its figures (literature, the audit, the GIAB,
HCC1143 and patient-data validations) were not re-measured in this session, so it
is not mapped here. Section A can be extended to the thesis as soon as the file is
available.

Paths below are relative to `stage1_igv_assistant/results/` unless they start
with `stage1_igv_assistant/` or `scripts/`. `AN` is
`synthetic_control_2026-09/analysis_2026-09-25/`; `P10` is
`stage1_igv_assistant/benchmark/runs/phase10_rerun_2026-09-25/`.

**Status:** *reproduced exactly* — derived again from the same data, or checked in
the code, with the same value. *Replaced by the new experiment* — measured again
on the rebuilt implants or in the rerun; the new value replaces the old one, even
where the two coincide. *Record lost, not regenerated* — the original record was
lost in the 2026-09-23 reinstall and this session did not measure it again.
*No longer supported* — the new measurement does not support the claim as stated.

## A. The conference abstract — every numeric claim

| # | Claim | Where | Old value | New value | Record, key | Status |
|---|---|---|---|---|---|---|
| A1 | the assistant is made of eleven tools reachable over the Model Context Protocol | Metodai, sentence 1 ("vienuolika įrankių") | 11 | 11 evidence tools (and 4 candidate-set tools on a second server) | `stage1_igv_assistant/server.py`; `ui.assert_tool_contract()` returns (11, 4); `tests/test_server.py` passes | reproduced exactly |
| A2 | four of them measure breakpoint evidence independently (discordant pairs, soft clips, split reads, depth) | Metodai, sentence 2 ("Keturi iš jų") | 4 | 4 | `stage1_igv_assistant/tools/bam_tools.py`, `EVIDENCE_LAYER_NAMES` | reproduced exactly |
| A3 | sensitivity measured by implanting twelve balanced translocations into real sequencing data | Metodai, last sentence ("dvylika") | 12 | 12, rebuilt 2026-09-25 | `synthetic_control_2026-09/implants_ground_truth.json`, `implants` | replaced by the new experiment |
| A4 | two clinical whole-genome samples | Rezultatai, sentence 1 ("Dviejuose") | 2 | 2 (SAMPLE_A, SAMPLE_B) | `patient_rerun_2026-09.json`, `samples` | reproduced exactly |
| A5 | about 9,000 interchromosomal junctions | Rezultatai, sentence 1 ("apie 9 000") | about 9,000 | 9,172 and 9,655 BND junctions after deduplication | `patient_rerun_2026-09.json`, `funnel.SAMPLE_A.bnd_after_dedup`, `funnel.SAMPLE_B.bnd_after_dedup` | reproduced exactly |
| A6 | reduced to 17 and 19 candidates | Rezultatai, sentence 1 ("iki 17 ir 19") | 17, 19 | 17, 19 | `patient_rerun_2026-09.json`, `funnel.SAMPLE_A.after_recurrence_500bp`, `funnel.SAMPLE_B.after_recurrence_500bp` | reproduced exactly |
| A7 | 14 of 24 implanted junctions detected | Rezultatai, sentence 2 ("Iš 24 … aptiktos 14") | 14 of 24 | **16 of 24** | `AN/figures.json`, `sensitivity.overall.junctions` | replaced by the new experiment |
| A8 | clean sequence 8 of 8 | Rezultatai, sentence 2 | 8 of 8 | 8 of 8 | `AN/figures.json`, `sensitivity.by_class.chosen_rule.clean_unique.junctions` | replaced by the new experiment (same value) |
| A9 | next to repeats 6 of 8 | Rezultatai, sentence 2 | 6 of 8 | **8 of 8** — the difference is IMP06, missed then and detected now; realigned without the `.alt` it is missed again, so ALT-aware alignment accounts for it | `AN/figures.json`, `sensitivity.by_class.chosen_rule.repeat_adjacent.junctions`; `AN/noalt/compare.json`, `implants.IMP06` | replaced by the new experiment |
| A10 | poorly mappable regions 0 of 8 | Rezultatai, sentence 2 | 0 of 8 | 0 of 8 | `AN/figures.json`, `sensitivity.by_class.chosen_rule.low_mappability.junctions` | replaced by the new experiment (same value) |
| A11 | every loss happened at discovery, none at filtering | Rezultatai, sentence 2 | all at discovery | all 8 at discovery: delly wrote no record for them; every PE × SR cell keeps all 16 detected | `AN/figures.json`, `loss_summary`, `grid` | replaced by the new experiment (same finding) |
| A12 | two fifths of the reads crossing a junction were invisible to the split-read signal … | Rezultatai, sentence 3 ("Dvi penktosios") | two fifths (121 of 304) | 149 of 390 (38 %) carry no SA tag | `AN/adjudication/summary.json`, `b_short_overhang.without_any_sa_tag`, `.truly_spanning_reads` | replaced by the new experiment (still about two fifths) |
| A13 | … and were recovered by the soft-clip signal | Rezultatai, sentence 3 ("atkurtos minkšto iškirpimo signalu") | 114 of 121 | **77 of 149**; of the 72 missed, first reason: 19 outside both breakpoint windows, 19 MAPQ under 20, 31 an overhang under 10 bp, 3 a clip under 10 bp | `AN/adjudication/summary.json`, `b_short_overhang.without_sa_caught_by_soft_clip_layer`, `.not_caught_first_reason` | **no longer supported** as stated: about half are recovered, not nearly all |
| A14 | the scoring function cannot rate a heterozygous balanced translocation as strong | Rezultatai, sentence 4 | — (the arithmetic) | attainable ceiling 57.5 at all 32 detected breakends, strong band 70, reachable nowhere | `AN/figures.json`, `ceiling_summary.attainable_here_range`, `.any_strong_band_reachable` | reproduced exactly (the arithmetic, re-measured on the rebuilt implants) |
| A15 | four models of different capability — two local, two cloud | Rezultatai, sentence 5 | 4 (2 + 2) | 2 local models run; the 2 API models could not run (stored key rejected, HTTP 401) | `P10/claude-sonnet-5/status.json`, `P10/claude-opus-5/status.json` | record lost, not regenerated (API half) |
| A16 | none explained the ceiling: 1 of 20 attempts | Rezultatai, sentence 5 ("1 iš 20") | 1 of 20 | local models without the fields: 0 of 10 | `P10/tally.json`, `scheduled_runs["local models WITHOUT"].C2` | replaced by the new experiment (local half); the 20-run figure's record is lost |
| A17 | with it added to one tool's return, 19 of 20 attempts explained the ceiling | Rezultatai, sentence 7 ("19 iš 20") | 19 of 20 | local models with the fields: **2 of 10** (the fields reached the model in 7 of the 10 runs) | `P10/tally.json`, `scheduled_runs["local models WITH"].C2`, `.ceiling_seen` | **no longer supported** for the local models; API half not rerun |
| A18 | the local 4-billion-parameter model explained it 5 of 5 | Rezultatai, sentence 7 ("5 iš 5") | 5 of 5 | **2 of 5** (qwen3.5:4b, full argument both times) | `P10/tally.json`, `scheduled_runs["qwen3.5:4b WITH"].C3`; quotes in `P10/scores.json` | **no longer supported** |

## B. The Phase 5, 6 and 10 figures (thesis and abstract, per the task statement)

| # | Claim | Old value | New value | Record, key | Status |
|---|---|---|---|---|---|
| B1 | translocations detected | 7 of 12 | 8 of 12, every one with both junctions | `AN/figures.json`, `sensitivity.overall.translocations_any_junction`, `.translocations_both_junctions` | replaced by the new experiment |
| B2 | junctions detected; by class | 14 of 24; 8/8, 6/8, 0/8 | 16 of 24; 8/8, 8/8, 0/8 | as A7–A10 | replaced by the new experiment |
| B3 | class robustness under the upper-decile rule | — | IMP01 moves to repeat-adjacent: clean 6/6, repeat-adjacent 10/10, low-mappability 0/8 — denominators change, rates (100 %, 100 %, 0 %) do not | `AN/figures.json`, `sensitivity.by_class.strict_rule` | new |
| B4 | where losses happen | all at discovery | all 8 at discovery | `AN/figures.json`, `loss_summary` | replaced by the new experiment (same finding) |
| B5 | background funnel | 894 → 513 → 155 → 27 | 894 → 513 → 155 → 27 (→ 27 → 27 after the primary and mask steps) | `AN/figures.json`, `background_funnel` | reproduced exactly (same background, same delly) |
| B6 | each detected implant added exactly 2 survivors | 2 | 2 in all 8 detected; the other 27 are the background's 27, identical by candidate_id, in all 12 BAMs | `AN/figures.json`, `precision` | replaced by the new experiment (same finding) |
| B7 | localisation | 12 of 14 at 0 bp, 2 of 14 at 2 bp | 6 of 16 at 0 bp, 10 of 16 at 1 bp (larger offset of the two ends) | `AN/figures.json`, `localisation_histogram_max_abs_bp`, `localisation` | replaced by the new experiment |
| B8 | reused coordinates | IMP01 detected; IMP06, IMP10 missed | IMP01 detected; **IMP06 detected**; IMP10 missed. IMP06 without the `.alt`: missed, junction pairs with both mates at MAPQ ≥ 20 37 → 0; IMP01 and IMP10 unchanged | `AN/figures.json`, `reused_coordinates`; `AN/noalt/compare.json` | replaced by the new experiment; IMP06 attributed to ALT-aware alignment |
| B9 | split-read tool vs delly, share of truly spanning reads | 60.5 % vs 54.4 % | 39.5 % (154) vs 35.1 % (137) of 390; the tool's raw count is closer to the truth in 10 of 12 implants (2 ties) | `AN/adjudication/summary.json`, `a_tool_vs_delly_vs_truth` | replaced by the new experiment |
| B10 | tool false-positive rate by class; r with low_mapq_fraction | 0 % / 3.9 % / 56.9 %, one locus; r = +0.142 | 11.1 % / 20.4 % / 87.7 %, r = 0.668 (41 counted reads span the junction but carry only a mapQ-0 SA to a wrong chromosome; counted as true: 0 % / 5.8 % / 60.0 %, r = 0.298); low-mappability FPs concentrated at IMP11 chr20:31,100,000 (40 of 41) | `AN/adjudication/summary.json`, `a_tool_vs_delly_vs_truth.by_class`, `.fp_rate_pearson_r_vs_low_mapq_fraction`, `.alternative_reading_wrong_partner_spanning_reads_counted_true` | replaced by the new experiment |
| B11 | short overhang: without SA; caught by soft clips | 121 of 304; 114 | 149 of 390; 77; no SA-less read under 30 bp lifted to -T 30 by chance matches | as A12–A13 | replaced by the new experiment |
| B12 | split_reads min_mapq 0 → 20 | removes 38 false positives, costs 35 true; zeroes the layer in low-mappability regions | removes 85 of 86 false positives, costs 30 genuine (154 → 124); zeroes all 8 low-mappability breakpoints | `AN/adjudication/summary.json`, `c_min_mapq_20` | replaced by the new experiment. The old figure is also in a comment in `stage1_igv_assistant/server.py`, a protected file not edited here |
| B13 | caller vs tool at the background's survivors | within ~1.5 reads where low_mapq_fraction ~0.004; up to 100-fold where ~0.457 | −5 to +9 reads where both breakends are under 0.01; up to 100-fold (1,700 vs 17, low_mapq_fraction 0.905 at that breakend); r = 0.689 | `AN/figures.json`, `caller_vs_tool_background`, `caller_vs_tool_pearson_r_absdiff_vs_low_mapq` | the 100-fold reproduced; the "~1.5 reads" summary's definition is lost, so that part is replaced by the new rows |
| B14 | PE × SR grid | 14/24 in every cell; SR ≥ 1 cut 156 → 28; SR ≥ 2 = SR ≥ 1 | 16/24 in every cell; 156 → 28 at PE ≥ 1 and 2 (155 → 27 at PE ≥ 3, 118 → 25 at PE ≥ 5); SR ≥ 2 = SR ≥ 1 | `AN/figures.json`, `grid` | replaced by the new experiment (156 → 28 same) |
| B15 | delly at -q 0 on the missed implants | 2 of 10 junctions called, 1 survived; raw BND ~7×; survivors 27 → 28 | 1 of 8 called (IMP11 J20, LowQual), 0 survived; raw BND ×7.15 (-q 0 -r 5) and ×6.54 (-q 0 -r 0); survivors 27 → 28 | `AN/ladder/analysis.json`, `summary` | replaced by the new experiment |
| B16 | scoring ceiling | 40.0–55.0 all moderate; discordant fraction ≤ 0.175; ceiling 57.5 at IMP01 chr20; IMP02, IMP04 at 55 via spurious depth | 32.5–55.0 (30 moderate, 2 weak); ≤ 0.164; 57.5 at every breakend; depth scored at 7 breakends of IMP04, IMP05, IMP07 | `AN/figures.json`, `ceiling_summary`, `ceiling` | replaced by the new experiment |
| B17 | rescue at a caller-missed breakpoint | IMP06 chr20:3,900,000: 44 discordant pairs, 20 soft clips, 14 split reads | IMP06 is no longer missed. Missed now: IMP09–IMP12; e.g. IMP09 chr20:33,700,000 10 / 7 / 9, 22.5 weak; IMP11 chr20:31,100,000 11 / 4 / 41 (1 from the junction), QUALITY-LIMITED; provenance `caller_supplied` on every call | `AN/figures.json`, `rescue`; `AN/rescue/` | **no longer supported** (the example is detected now); replaced |
| B18 | ceiling experiment, WITHOUT: assert / complete | 1 of 20 / 0 of 20 | local models: 0 of 10 / 0 of 10 | `P10/tally.json` | replaced by the new experiment (local half) |
| B19 | ceiling experiment, WITH: assert / complete | 19 of 20 / 17 of 20 | local models: 2 of 10 / 2 of 10 | `P10/tally.json` | **no longer supported** for the local models |
| B20 | qwen3.5:4b completed the argument | 5 of 5 | 2 of 5 | `P10/tally.json`, `scheduled_runs["qwen3.5:4b WITH"]` | **no longer supported** |
| B21 | qwen2.5:7b misused the ceiling | 2 of 5 | asserted it in 0 of 5 although the fields reached it in 4; one run names both limits yet implies strong is reachable, one contradicts itself | `P10/scores.json` | replaced by the new experiment |
| B22 | API spend | $0.96 for 10 API runs | no API run could be made; $0 spent; projected $5.77 for the planned 20 at current list prices | `P10/projection.json` | record lost, not regenerated |

## C. Patient figures

All from `patient_rerun_2026-09.json`, re-derived on 2026-09-25 from the same two
BAMs, verified byte for byte against the Phase 0 record, with delly v2.6.0.
Counts only; no coordinates or identifiers.

| # | Claim | Old value | New value | Key | Status |
|---|---|---|---|---|---|
| C1 | BAM sizes | 38,959,428,903 / 41,617,797,998 bytes | same | `task1_verification.SAMPLE_*.bytes` | reproduced exactly |
| C2 | "primary mapped reads" (idxstats, chr1-22, X, Y) | 631,618,015 / 677,604,873 | same | `task1_verification.SAMPLE_*.index_mapped_chr1_22_X_Y` | reproduced exactly (the only one of 1,024 definitions that gives both) |
| C3 | @PG records | 34 / 20 | same | `task1_verification.SAMPLE_*.pg_records_samtools_head` | reproduced exactly |
| C4 | delly records | 30,980 / 32,451 | same | `funnel.SAMPLE_*.total_records` | reproduced exactly |
| C5 | BND records | 9,187 / 9,676 | same | `funnel.SAMPLE_*.bnd_records` | reproduced exactly |
| C6 | BND after deduplication | 9,172 / 9,655 | same | `funnel.SAMPLE_*.bnd_after_dedup` | reproduced exactly |
| C7 | PASS BND | 896 / 923 | same | `funnel.SAMPLE_*.pass_bnd` | reproduced exactly |
| C8 | candidates before recurrence | 57 / 63 | same | `funnel.SAMPLE_*.before_recurrence` | reproduced exactly |
| C9 | candidates after recurrence (500 bp) | 17 / 19 | same | `funnel.SAMPLE_*.after_recurrence_500bp` | reproduced exactly |

The patient funnel filters svtype = BND first; the implant funnel in B5 does not
(it reproduces the first run's background funnel, which had no svtype step).

## D. Phase 12 — 2026-09-26

New figures, and every thesis figure marked † (supported by no committed record
when the thesis was revised at `62eb4d7`), each with its definition in the record
named. **Status:** *supported* — the new record gives the thesis value under the
stated definition; *replaced* — the new record gives a different value, which
replaces it; *still unsupported* — no record gives it and none could be made.
Nothing here maps the thesis text itself. `PP` is `results/patient_properties_2026-09-26.json`,
`SF` `results/supplementary_flags_2026-09-26.json`, `PB` `results/prefix_bridge_2026-09-26.json`,
`TC` `results/test_census_2026-09-26.json`, `P9` `stage1_igv_assistant/benchmark/runs/phase9_rerun_2026-09-26/`.

| # | Figure | Thesis value | Record value | Record, key | Status |
|---|---|---|---|---|---|
| D1 | bwa version † | 0.7.17-r1188 | no bwa `@PG` record carries `VN` (16 and 18 records: ID, PN, CL only); 0.7.17-r1188 occurs in the CL of every one | `PP`, `SAMPLE_*.bwa_version`, `.bwa_version_fields` | still unsupported as a declared version; supported only as the command line's executable |
| D2 | read length | — | 150 in every primary record examined | `PP`, `SAMPLE_*.read_length` | new |
| D3 | nominal depth † | 30.7 / 32.9-fold | 30.68 / 32.91 (index-mapped chr1-22, X, Y × 150 / 3,088,269,832) | `PP`, `SAMPLE_*.nominal_depth` | supported |
| D4 | duplicate fraction † | 12.96 / 12.97 % | 12.963 / 12.968 % of primary records; 12.73 / 12.76 % of all records | `PP`, `SAMPLE_*.duplicate_fraction` | supported, as primary duplicates over primary records |
| D5 | median insert size † | 317 / 311 bp | 314 / 312 bp (48,887 / 49,045 first-in-pair records, first 50 per seeded position) | `PP`, `SAMPLE_*.median_insert_size` | replaced |
| D6 | reference concordance † | 99.90 % in both | 951 of 952 (99.895 %) and 954 of 954 (100 %) positions meeting MAPQ 20, BQ 20, depth 10; seed 20260926 | `PP`, `SAMPLE_*.concordance` | replaced |
| D7 | supplementary records † | none (validation document; thesis: no read carries the flag) against flagstat 3,694,802 / 3,726,497 | 3,694,802 / 3,726,497, all on alt, decoy and HLA contigs, none on chr1-22, X, Y; 0 on chr21 by direct count and by every reconstruction of the old scan | `SF`, `SAMPLE_*.count`, `.oldscan_chr21` | flagstat supported; "no read carries the supplementary flag" replaced: true of the primary chromosomes, not of the files |
| D8 | `@PG` programs | — | bwa (per lane: 16 / 18), bamcat, bamsormadup; 16 records without PN in SAMPLE_A | `SF`, `SAMPLE_*.pg` | new |
| D9 | background merge groups DEL+DUP † | 33 | 33 groups exactly DEL+DUP; 34 hold at least one of each | `PB`, `public.background_prefix_dedup` | supported, under "exactly DEL+DUP" |
| D10 | background merge groups, both INV orientations † | 13 | 12 (one more group joins both orientations of a BND pair) | `PB`, `public.background_prefix_dedup.inv_both_groups` | replaced |
| D11 | pre-fix patient BND after dedup † | 8,716 / 9,144 | 8,716 / 9,144 | `PB`, `patient["prefix_dedup_prefix_compare …"]` | supported |
| D12 | pre-fix patient PASS † | 838 / 862 | 838 / 862 | as D11 | supported |
| D13 | pre-fix final survivors † | 15 / 16 | 15 / 16, whichever way the comparison is called | as D11 | supported |
| D14 | pre-fix two-direction difference † | 1,738 | 1,738 (SAMPLE_A: 15,640 scanned, 13,902 as target) only with the dedup already fixed; 53 and 89 before both fixes | `PB`, `patient[*].two_direction_diff` | supported only for the comparison defect alone |
| D15 | survivors unchanged by the comparison fix | 17 / 19 | 17 / 19 if each sample is compared as the scanning side; 17 / 30 with one compare(A, B) call, as the committed funnel makes it | `PB`, `patient["fixed_dedup_prefix_compare …"]`, `final_survivors`, `final_survivors_two_calls` | still unsupported: depends on how the old funnel called the comparison, which is not recorded |
| D16 | test suites † | seventeen | 24 | `TC`, `*.summary.suites` | replaced |
| D17 | assertions: no network / annotation service / remote BAM † | 306 / 308 / 309 | 660 / 662 / 685 held, 0 failed, under the recorded definition (the old count's definition is not recorded); conditions imposed from outside, no switches exist | `TC`, `*.summary` | replaced |
| D18 | IGV-dependent checks | — | not run in any condition: IGV is not installed on this machine | `TC`, `*.suites["test_bam_tools.py"].not_run` | new |
| D19 | frontier runs † | sonnet 30, opus 15 | 30 and 15 (reconstructed cases) | `P9/measures.json`, `per_model` | supported |
| D20 | malformed arguments † | 0 in 313 calls | 0 in 341 calls | `P9/measures.json` | replaced |
| D21 | refusals † | 0 | 0 | `P9/measures.json` | supported |
| D22 | frontier spend † | $3.86 | $3.5875 | `P9/measures.json`, `total_cost_usd` | replaced |
| D23 | tool calls written as prose † | none (no record) | 0 | `P9/measures.json` | supported |
| D24 | coordinate drift † | none (no record) | 0 runs | `P9/measures.json` | supported |
| D25 | invented findings † | none (no record) | 34 numbers unmatched to any tool return (the registered count); each traced to a tool description, the system prompt or a rounded returned position, so 0 after reading; no unrun check reported, no image described | `P9/reading_scores.json`, `per_model` | supported after reading, with the registered count stated beside it |
| D26 | false premise rejected (frontier) | every run | 8 of 8 (sonnet 5, opus 3) | `P9/reading_scores.json`, `case_a` | supported |
| D27 | Phase 10 API spend † | $0.96 for 10 runs (lost); projection $5.77 for 20 | $2.049 for 20 runs | `stage1_igv_assistant/benchmark/runs/phase10_rerun_2026-09-25/spend_2026-09-26.json` | replaced |
| D28 | ceiling, WITHOUT the fields: state / complete (four models × five runs) † | 1 of 20 / 0 of 20 | 0 of 20 / 0 of 20, scored blind (local runs 1-5, a replacement where one exists, and the API runs) | `P10/blind/blind_scores.json`, `cells.original_runs_1_5`, `cells.api_runs_1_5` | replaced |
| D29 | ceiling, WITH the fields: state / complete † | 19 of 20 / 17 of 20 | 12 of 20 / 12 of 20: `claude-sonnet-5` 5/5, `claude-opus-5` 5/5, `qwen3.5:4b` 2/5, `qwen2.5:7b` 0/5 | as D28 | replaced |
| D30 | `qwen3.5:4b` completes the argument † | 5 of 5 | 2 of 5 (runs 1-5); 4 of 15 (runs 1-15) | `P10/blind/blind_scores.json`, `cells` | replaced |
| D31 | `qwen2.5:7b` states the ceiling and misuses it † | 4 of 5, misused in 2 | states it 0 of 5 (runs 1-5), 3 of 15 (runs 1-15); completes it 0 of 15 | as D30 | replaced |
| D32 | registered test, 15-run cells, WITH vs WITHOUT (Fisher, two-sided) | — | `qwen3.5:4b`: C2 5/15 vs 0/15, p = 0.042; C3 4/15 vs 0/15, p = 0.100. `qwen2.5:7b`: C2 3/15 vs 0/15, p = 0.224; C3 0 vs 0, p = 1.0 | `P10/blind/blind_scores.json`, `fisher_exact_two_sided_15_run_cells` | new |
| D33 | blind against the unblinded scores of 2026-09-25 (23 local runs) | — | C1 22/23, C1b 22/23, C2 23/23, C3 23/23; both disagreements are `qwen2.5:7b` WITH run 2 | `P10/blind/blind_scores.json`, `agreement_with_unblinded_2026_09_25`, `disagreements` | new |
| D34 | `claude-sonnet-5` calls the evidence strong | — | with the fields 5 of 5 (C1b 0/5); without them 0 of 5 (C1b 5/5). Not a registered test | `P10/blind/blind_scores.json`, `cells.api_runs_1_5` | new |

`P10` is `stage1_igv_assistant/benchmark/runs/phase10_rerun_2026-09-25/`. The
blind scores come from a fresh verifier subagent that read only the packet
(prompt and final answer under random IDs); the key's salted hash was committed
before scoring (`blind/commitment.json`, `bcd7669`) and the key after it.

## E. The thesis and the abstract — 2026-09-27

Mapped from `docs/thesis/Rimas_MSc_Thesis.docx` (working draft of 26 September,
modified 2026-09-27 01:00; 24,776 words, 12 † marks) and the abstract
`Rimas_tezes_DI_medicinoje_LT.docx` (Desktop, same time; 520 words), read only:
text extracted from `word/document.xml`. Scope: Methods through both summaries of
the thesis (literature review and references skipped) and the whole abstract.
Every numeric claim is listed; a figure repeated in several places is listed once,
with its places. **Status:** *supported* — a committed record gives the value;
*differs* — a committed record gives another value; *unsupported* — no committed
record gives it. **†** marks a figure the thesis marks as resting on lost records.
"Derived" means computed here from the named record's own fields, the derivation
stated. Short names: `PP`, `SF`, `PB`, `TC`, `P9`, `P10`, `AN` as in sections A–D;
`PR` = `results/patient_rerun_2026-09.json`; `P8` = `results/phase8_final_record.json`;
`RPD` = `results/REAL_PATIENT_DATA_VALIDATION.md`; `SC` = `results/synthetic_control_2026-09/`;
`P9F` = `stage1_igv_assistant/benchmark/runs/phase9_case_f_2026-09-27/`.

### E.1 Flags — what the thesis or abstract should change

| # | Where | Claim | Record | Flag |
|---|---|---|---|---|
| E1 | Results, calibration | initial window "returning ratios of 0.66–0.71" | `GIAB_PUBLIC_DATA_VALIDATION.md` line 70: 0.618–0.707 (PacBio), 0.662 (Illumina) | **differs**: the range is 0.62–0.71 |
| E2 | Results, second round; Discussion; Conclusion 4 | the 7B model confirmed the false premise "in two of five runs" | `P8` has no false-premise judgement (only answered-by-case) | **unsupported, unmarked** |
| E3 | Results, second round | the 4B model "rejected the false premise in four of five" | as E2 (`P8`: case a answered 4 of 5) | **unsupported, unmarked** |
| E4 | Results, second round; Discussion; Conclusion 4 | frontier sweep: "no malformed arguments in 313 tool calls across forty-five runs, refused none, and rejected the false premise in every run" | the sweep's record was lost (`results/README.md`, 2026-09-26 note); the rerun gives 341 calls (`P9`) | **unsupported, unmarked** (the thesis's own provenance section lists this sweep's records as lost, but these figures carry no †) |
| E5 | Results, second round; Discussion | "Each sweep cost under four United States dollars" | first sweep $3.86 in lost record; rerun $3.5875 (`P9/measures.json`) | the first half **unsupported, unmarked** |
| E6 | Results, Deployability | "the same locus scored 43.3 instead of 47.5" | 47.5 and 43.3 are the first implant set's values (`docs/NAUDOJIMAS.md` note of 2026-09-25: "to rinkinio įrašai prarasti"); no committed record | **unsupported, unmarked** |
| E7 | Methods, gate | "the first construction's result of 82 per cent" | quoted in `SC/README.md` (50/61) from the lost construction | **unsupported, unmarked** (a first-construction figure without †) |
| E8 | Results, Deployability | bundle "of 1.59 megabytes" | `docs/DIEGIMAS.md` and `results/README.md`: 1.6 MB; no record of 1.59 | **unsupported** (rounding of an unrecorded size) |
| E9 | Results, The Instrument | at the missed implants manual entry returns "at best ten discordant pairs and a weak score of 22.5, otherwise a withheld score" | `AN/figures.json` `rescue`: IMP09 chr20 10 pairs, 22.5 weak; IMP11 chr20 **11** pairs, withheld | wording: the most pairs returned is 11 (at a withheld call); "at best ten" holds only for scored calls |
| E10 | Table 2 | first-construction column | cells not marked †; the table note says the records were lost | the † convention is carried by the note, not the cells |
| E11 | Discussion, methodological | the concealed gate firing was "a false positive, as it proved, and amended before anything was published" | no committed record (the push history shows no identifier; the amend itself is not recorded) | **unsupported** wording |
| E13 | Results, audit and reproducibility | the suite "now stands at nineteen tests"; "all nineteen tests passed" from a clean clone | git history of `tests/test_bam_tools.py`: 9 numbered TEST blocks at the audit commit `4b7f8de`, 18 at `2044852` (2026-08-11); "nineteen" occurs only in the earlier chapter draft (`9dce3bc`), not in a record | **unsupported, unmarked** |
| E12 | Results, GIAB | "deletion of 3,357 bases" | `stage1_igv_assistant/data/HG002_GRCh38_CMRG_SV.vcf.gz`, chr1:115686862: len(REF) − len(ALT) = 3,357 | **supported** — but `RESULTS_HCC1143.md`, `TUTORIAL.md` say 3,359; the documents, not the thesis, are off |

No † figure is supported by a committed record: every † (memory model, first
construction's 14 of 24 and 114 of 121 and r = 0.142, eight of eighteen prose
calls, the withheld-as-absence run, the first ceiling comparison's 1 and 19 of 20,
R² 0.9998) was checked and none has a record. Wording the records contradict: none
beyond E9 and E11; the "all", "every", "none" statements checked below hold.

### E.2 Methods

| # | Claim | Value | Record : key | Status |
|---|---|---|---|---|
| E20 | background slice aligned with BWA 0.7.15 mem -Y, ALT-aware | — | `SC/run_record_2026-09-25_continued.json` (implant alignment note); `AN/noalt/align.json` | supported |
| E21 | clinical BAM sizes | 36.3, 38.8 GiB; 39 and 42 GB | `PR` `task1_verification.SAMPLE_*.bytes` (38,959,428,903; 41,617,797,998) | supported (derived) |
| E22 | bwa 0.7.17-r1188 in every command line, no version field | — | `PP` `SAMPLE_*.bwa_version_fields` | supported |
| E23 | merged and duplicate-marked (bamcat, bamsormadup) | — | `SF` `SAMPLE_*.pg` | supported |
| E24 | every aligner record -M, not -Y | — | `PR` `bwa_flags_in_all_command_lines` | supported |
| E25 | 3,171 alternate contigs marked | 3,171 | `PR` `ah_contigs`, `ah_set_equals_hs38DH_alt` | supported |
| E26 | pa tag on 2.66 and 2.26 % of primary reads, chr20+21 | 2.66, 2.26 | `PR` `pa_tagged_chr20_chr21` / `primary_reads_chr20_chr21` (648,506/24,364,354; 593,172/26,296,677) | supported (derived) |
| E27 | reads 150 bases, ~400,000 examined per sample | 150; 412,826 / 446,318 | `PP` `SAMPLE_*.read_length` | supported |
| E28 | nominal depth, definition | 30.7 / 32.9 | `PP` `nominal_depth` | supported |
| E29 | duplicates of primary records | 12.96 / 12.97 % | `PP` `duplicate_fraction.primary_records` | supported |
| E30 | depth after duplicates "nearer 27 and 29-fold" | 26.7, 28.6 | derived from E28 × (1 − E29) | supported |
| E31 | median insert size | 314 / 312 | `PP` `median_insert_size` | supported |
| E32 | supplementary records: totals, ~2.9 M decoy, 0.76 / 0.82 M alt, ~3,000 HLA, none primary | 3,694,802 / 3,726,497 | `SF` `SAMPLE_*.count.flag_0x800.by_class` | supported |
| E33 | 3,366 contigs match; 1,000 positions, depth 10, MAPQ and BQ 20, seed | — | `PR` `task1_verification.*.contig_md5`; `PP` `concordance`, `floors`, `seed` | supported |
| E34 | eleven tools, four layers, seven others; startup assertion | 11, 4, 7 | A1, A2; `ui.assert_tool_contract()` | supported |
| E35 | 0–25 per layer, bands 70 and 40 | — | `stage1_igv_assistant/score_tiers.py` (derived from `bam_tools.py`); `P10/blind/packet.json` ceiling fields | supported |
| E36 | DELLY v2.6.0, exclude template, one thread | — | `PR` `delly.*.command` (`OMP_NUM_THREADS=1`, `-x`) | supported |
| E37 | six filters, seventh by comparison, 500 bp | — | `PR` `funnel_filters` | supported |
| E38 | twelve implants chr20–chr21, 24 junctions, three contexts | 12, 24 | `SC/implants_ground_truth.json` | supported |
| E39 | reinstall on 23 September 2026; generator never committed | — | `results/README.md` (reinstall note); git history has no generator before the rebuild | supported |
| E40 | three first-construction coordinates reused | 3 | `AN/figures.json` `reused_coordinates` | supported |
| E41 | ART HiSeq 2500, five profiles, 0.001905 vs 0.002104 | — | `SC/art_profiles.json`, `SC/background_measurements.json` | supported |
| E42 | fragment length mean 442.6, sd 104.6 | 442.64, 104.58 | `SC/background_measurements.json` `fragment_length` | supported |
| E43 | eight-value base qualities | 8 | `SC/background_measurements.json` `quality_bins.distinct` | supported |
| E44 | one alignment call, BWA 0.7.15 mem -Y, whole-genome index with .alt | — | `SC/run_record_2026-09-25_continued.json` | supported |
| E45 | 978 scanned candidate breakends; stricter rule (dinucleotide run 8) | 978 | `SC/selection.json`; `AN/adjudication/summary.json` `strict_classes` | supported |
| E46 | gate: 70 %, 18 of 33 = 54.5 %, fifteen explained by score 30 / seed 19, 61.1 %, P = 0.27, −7.3 / +1.6 % | — | `SC/IMP01_gate_revised_2026-09-25.json`; `SC/README.md` | supported |
| E47 | wrong-partner criterion fails across twelve | 41 reads | `AN/adjudication/summary.json` `spanning_reads_with_a_wrong_partner_sa` | supported |
| E48 | 5 runs per local model and per sonnet; opus 15 | — | `P8` `models.*.runs` (30 = 6 × 5); `P9/registration.json` | supported |
| E49 | blind rescoring: 23 earlier local runs; agreement all on C2 and C3, 22 of 23 on the others | — | `P10/blind/blind_scores.json` `agreement_with_unblinded_2026_09_25` | supported |
| E50 | 8 GB accelerator; fallbacks disabled | — | `P8` `hardware`; `P9/registration.json` `config` | supported |
| E51 | versions: Python 3.11 / 3.14, pysam 0.24.0, FastMCP 3.4.6, samtools 1.21 + htslib 1.24, ART 2.5.8, IGV 2.17.4, Ollama 0.32.9 | — | `requirements.txt`; `results/README.md` (3.14.4); `SC/art_profiles.json`; `scripts/install_igv.sh`; `P8` `hardware.ollama` | supported |

### E.3 Results

| # | Claim | Value | Record : key | Status |
|---|---|---|---|---|
| E60 | fixture chr1↔chr8: 15 / 5 / 8; 100/100 STRONG, 4/4 | — | `tests/test_bam_tools.py` fixture; `DEMO_END_TO_END.md` | supported |
| E61 | HCC1143 normal (germline) at a tutorial translocation locus; 0.7–0.8 %; zero SA in 572,731 reads | — | `RESULTS_HCC1143.md` | supported |
| E62 | GIAB deletion, VCF position chr1:115,686,862 | 3,357 bases | see E12 | supported |
| E63 | calibration: ±500 → ±2,000 bp, 0.6 → 0.7; 0.609 / 0.542 | — | `GIAB_PUBLIC_DATA_VALIDATION.md` lines 70–97 | supported (see E1 for 0.66–0.71) |
| E64 | 16 thresholds (13 scoring, 3 text), 2 empirical, 14 judgement; eleven in the scoring documentation | — | `bam_tools.py` threshold inventory (lines 28–40); the summary tool's description ("only ONE … The other 10") | supported |
| E65 | 42 control loci; localisation 1 kb, 1,400 bases, 800–900 bases | — | `RPD` findings 4–5; `bam_tools.py` localisation note | supported |
| E66 | first session: 100-base bin, 115,686,865, 13-read pileup; 15 + 15 + 0 + 30 = 60 vs 30 | — | `LLM_SESSION_1.md` | supported |
| E67 | second session: nine tools; 410–520 to 175, ratio 0.547; true depth 100–320; one mate on chr12 | — | `LLM_SESSION_2_WITH_VISUAL.md` lines 82–92 | supported |
| E68 | audit: five critical, twelve substantive, four housekeeping; seven of nine tools; two wrappers | 5 / 12 / 4; 7/9; 2 | `AUDIT_2026_08.md` SUMMARY; commit `4b7f8de` message (split_reads and the summary's label default) | supported |
| E69 | nineteen tests after the audit; nineteen passed from a clean clone | 19 | see E13 | unsupported |
| E70 | bin width moved depth 39 %; invariant within 1–2 % after the fix; pileups 13 / 2 / 1; dip 1,500 bases away; separation twofold → sixfold | — | `LLM_SESSION_3_BLIND.md` lines 156, 174 | supported |
| E71 | patient validation: 261 alt, 2,512 unplaced and decoy, 525 HLA; ~31-fold; 70 calls, 16 probe groups; ≤ 0.096 s; annotation up to 3.6 s; 68 MB; 100-kb window; ~470,000 entries | — | `RPD` lines 40–60, 445–472 (3.584 s), 418–421 | supported |
| E72 | ~2.4 million decoy supplementary records per sample with a primary on the main chromosomes | 2,358,228 / 2,362,639 | `SF` `SAMPLE_*.count.sa_first_class.decoy.primary` | supported |
| E73 | 28 cells, 16 strong; up to 452 reads; 71 %; −29 / −60 % and −66 / −87 %; 10–12 % strand; 16 of 44, 102 entries; MHC moderate → weak, mostly from the MAPQ filter | — | `RPD` lines 120–141, 249–262, 555–557, 238, 569–586 | supported |
| E74 | false-positive rate: 120 positions, 26 %, 0.739 / 0.753, 56–58 %, 3× less scatter, 0.523–0.868 crossing between 200 and 500 | — | `RPD` lines 172–187, 313 | supported |
| E75 | support minimum: 33 → 71 / 60 %, 60 → 17 / 21 %, 7.5 → 0; 43 / 50 %; 17 / 26 % | — | `RPD` finding 4 and re-measurement | supported |
| E76 | fourteen suites added, seventeen total; census 24 suites, 660 / 662 / 685, none failing | — | git history (17 test files at `96b5e25`, 2026-08-31); `TC` `*.summary` | supported |
| E77 | twelve fully resolved, one with caveat, one part-resolved | 12 / 1 / 1 | `RPD` status lines | supported |
| E78 | concordance 951 / 952 and 954 / 954 | — | `PP` `concordance` | supported |
| E79 | first delly run: 46 min, 2 h 4 min, 1,184 / 1,187 MiB | — | `PR` `phase3_targets` | supported (as the Phase 3 record quoted there) |
| E80 | 30,980 / 32,451; 9,172 / 9,655; 896 / 923; rerun reproduced every count; rerun memory 1,158 / 1,187; runs concurrent | — | `PR` `funnel`, `all_counts_match`, `delly.*.peak_rss_mib`, `delly` note | supported |
| E81 | memory model: 0.9998, 13.19 GiB, 467 / 444 Mb, 1,109 / 461 MiB, 1,121 MiB, 5.6 / 5.9 % † | — | none (lost) | unsupported, marked † — correct |
| E82 | bridge: 33 DEL+DUP groups; 13 = 12 inversion + 1 breakend; 840 → 894, 24 → 27; 8,716 / 9,144 → 9,172 / 9,655; 838 / 862 → 896 / 923; 15 / 16 → 17 / 19; rerun reproduces | — | `PB` `public`, `patient` (D9–D13) | supported |
| E83 | 54 of 56; 1,738 and 1,795; ~70 % removed at the last step; 17 / 19 two calls, 17 / 30 one call | — | `PB` `two_direction_diff`, `final_survivors*`; `PR` funnel (57 → 17, 63 → 19) | supported |
| E84 | 9,172 → 17, 99.8 % | — | `PR` funnel | supported |
| E85 | first implant depth −7.3 / +1.6 %, −4.6 / +3.7 % without duplicates | — | `SC/README.md` line 70 | supported |
| E86 | 16 of 24; 8 of 12 each with both; 8/8, 8/8, 0/8; strict rule 6 and 10, rates 100 / 100 / 0 | — | `AN/figures.json` `sensitivity` (A7–A10, B1–B3) | supported |
| E87 | first construction 14 of 24 †, 6 of 8 | — | none (lost) | unsupported, marked † (see E10) |
| E88 | IMP06 realigned without .alt: 37 → 0 pairs; the two control implants change in no read and no call | — | `AN/noalt/compare.json` `implants.*.reads.changed`, `detection_identical` | supported |
| E89 | no filter removed a true positive in 12 grid cells | — | `AN/figures.json` `grid` | supported |
| E90 | missed implants: 114 spanning reads, 106 MAPQ 0, 49 with the primary in a breakpoint window, 5 MAPQ ≥ 20; detected classes 276, all in a window, none MAPQ 0 | — | derived from `AN/adjudication/IMP*.json.gz` `reads` (mapq, breakpoint_window) | supported (derived) |
| E91 | MAPQ floor: none at default or 1; at 0 one called, LowQual, not surviving; raw BND 13 → 85 / 93, 6.5–7.2×; survivors 27 → 28 | — | `AN/ladder/analysis.json` `summary` | supported |
| E92 | localisation 6 of 16 at 0 bp, 10 within 1 bp | — | `AN/figures.json` `localisation_histogram_max_abs_bp` | supported |
| E93 | 27 identical survivors in all 12; each detected implant adds two; 156 → 28 and 155 → 27; SR ≥ 2 = SR ≥ 1 | — | `AN/figures.json` `precision`, `grid` | supported |
| E94 | Table 2, rebuild column | — | B1–B16 | supported |
| E95 | enumeration reproduced every tool count at all 24 breakpoints | — | `AN/adjudication/summary.json` `enumeration_reproduces_every_tool_count` | supported |
| E96 | 390 spanning, 149 without SA; seed 19 / score 30; shorter side 1–30; none alignable under 30 | — | `AN/adjudication/summary.json` `b_short_overhang` (`reads[*].shorter_side`, `flagged_alignable_despite_side_rule` empty) | supported |
| E97 | soft-clip layer 77 of 149; misses 31 / 19 / 19 / 3; about a fifth of all spanning reads | 77 / 390 = 19.7 % | `b_short_overhang.not_caught_first_reason` | supported |
| E98 | first construction 114 of 121 † | — | none (lost) | unsupported, marked † — correct |
| E99 | tool 154 (39.5 %), caller 137 (35.1 %); closer in 10 of 12, 2 ties | — | `a_tool_vs_delly_vs_truth` | supported |
| E100 | 41 wrong-partner reads; FP 11.1 / 20.4 / 87.7 %, r 0.668; 0 / 5.8 / 60.0 %, r 0.298; one breakpoint 40 of 57 | — | `a_tool_vs_delly_vs_truth.by_class`, `per_breakpoint` (IMP11 chr20:31,100,000: 40 of 57) | supported |
| E101 | first construction r = 0.142 † | — | none (lost) | unsupported, marked † — correct |
| E102 | MAPQ 20: removes 85 of 86 FPs, costs 30, 154 → 124, silences 8 breakpoints | — | `c_min_mapq_20` | supported |
| E103 | background: −5 to +9 below 0.01; 1,700 vs 17 at 0.905; r = 0.689 | — | `AN/figures.json` `caller_vs_tool_background` | supported |
| E104 | ceiling: 32.5–55.0, 30 moderate and 2 weak; 24–46 pairs; 48 of 1,094 mates elsewhere; 6–24 clips; 0.164; 57.5 at all 32 | — | `AN/figures.json` `ceiling`, `ceiling_summary`; the 1,094 / 48 derived from `AN/chain/IMP*.json.gz` (discordant_pairs call at each breakend, `mate_chromosomes`) | supported (derived) |
| E105 | 0.25 at a 200-base window, 65 | — | `P10/proof_prepared_payload_2026-09-26.json` `loci[1]` (chr21:14,100,000) | supported |
| E106 | depth at 7 breakends of 3 implants, 15 points; 5 at 55.0 | — | `AN/figures.json` `ceiling` (depth_score, evidence_score) | supported |
| E107 | first round: 3 models × 3 cases × 3 runs; 3 of 9 at the turn cap; 6–10 calls; 5 of 6; 38 characters; 3 of 3; one further site and a soft-clip one | — | `BENCHMARK_LOCAL_MODELS.md` lines 74, 140–162, 294; `BENCHMARK_CLAUDE_BASELINE.md` line 15; commit `ae4a02c` message | supported |
| E108 | eight of eighteen prose runs †; the withheld-as-absence run † | — | none | unsupported, marked † — correct |
| E109 | 4B: 5 of 30 unanswered, 26.3 %; 9B: 2,048 vs 2,690, 11 of 30 answered | — | `P8` `models`, `vram_measurements` | supported |
| E110 | rerun: 0 malformed in 341; no refusal, prose call or drift; 8 of 8; 34 numbers all traced | — | `P9/measures.json`, `reading_scores.json` | supported |
| E111 | local malformed 13–26 % | 13.4 / 17.1 / 26.3 | `P8` `models.*.schema_invalid_pct` | supported |
| E112 | reasoning off 16.4 → 35.8 %; low 9.3 % with one fewer answered of 18; trim 16.4 → 31.8 %; 5,114, 2,690, 866 tokens | — | `P8` `thinking.on_task_18_runs_each`, `schema_token_cost` | supported |
| E113 | fifteen tool descriptions; 2,600 characters; 71 checks; one replaced for the check, two for server errors | — | `P10/proof_prepared_payload_2026-09-26.json`; `P10/*/meta_replacement_*.json` | supported |
| E114 | first comparison 1 of 20 / 19 of 20 † | — | none | unsupported, marked † — correct |
| E115 | Table 3, all cells; Fisher p = 0.042, 0.100, 0.224 | — | `P10/blind/blind_scores.json` `cells`, `fisher_exact_two_sided_15_run_cells` | supported |
| E116 | 0 of 20 without, 12 of 20 with; frontier every run; local 7 of 10 called the summary; qwen2.5:7b 3 of 15, fields reached it in 11; one run implied strong reachable, one strong-then-moderate | — | `P10/blind/blind_scores.json` `cells`; `P10/scores.json` notes (qwen2.5:7b WITH runs 3 and 5) | supported |
| E117 | sonnet strong in 5 of 5 with, moderate in 5 of 5 without, quoted sentence; opus strong in 2 of 5 in each condition | — | `P10/blind/blind_scores.json` `items` (C1b), `P10/blind/packet.json` | supported |
| E118 | eleven and four tools at startup | — | A1 | supported |
| E119 | three packages, a fourth optional; third installation attempt | — | `requirements.txt`, `requirements-api.txt`; `results/README.md` ("three iterations") | supported |
| E120 | demo locus 40.0, moderate, three of four layers | — | `AN/figures.json` `rescue.IMP01.ends.chr20` (the demo locus); `docs/NAUDOJIMAS.md` | supported |
| E121 | shortlists 17 and 19, reproduced | — | `PR` | supported |

### E.4 Discussion, conclusions, summaries

| # | Claim | Record | Status |
|---|---|---|---|
| E130 | 16 of 24, eight losses at discovery; PE and SR minima remove 95 % of the PASS background (513 → 27) | `AN/figures.json` `background_funnel` | supported |
| E131 | 106 of 114 MAPQ 0, fewer than half at the locus; one of eight called, removed; sevenfold | E90, E91 | supported |
| E132 | soft-clip recovers about half; r 0.668 / 0.298; construction check at one implant | E96–E100 | supported |
| E133 | every detected breakend moderate or weak | E104 | supported |
| E134 | 5 of 15 and 3 of 15 with the fields; frontier every run; 19 of 20 † | E115, E114 | supported / † correct |
| E135 | "confirmed a false clinical premise in two of five", "313 tool calls", "under four US dollars each time" | E2, E4, E5 | **unsupported, unmarked** |
| E136 | "failed to answer in five of thirty" | `P8` | supported |
| E137 | seven apparatus defects; the post-reinstall instances | the defects' own records (`RPD`, `BENCHMARK_*`, `AN/noalt/align_first_check_defective.json`, `SF`); the gate firing: see E11 | supported except E11 |
| E138 | limitations: 70 % recurrence; 14 of 16; 26 %; 11 vs 16; 5 frontier and 15 local runs per condition | E64, E74, E83, E115 | supported |
| E139 | Conclusions 1–5 | E86, E90, E96, E97, E111, E115, E119, E121; Conclusion 4's "two runs of five" and "313": E2, E4 | supported except E2, E4 |
| E140 | Summary and Lithuanian summary | as E121, E86, E90, E96, E115 | supported |

### E.5 Abstract

| # | Claim | Record | Status |
|---|---|---|---|
| E150 | eleven tools, four measuring; twelve implants; four models, two cloud and two local; 71 checks; blind scoring | A1–A3, E113, E49 | supported |
| E151 | ~9,000 → 17 and 19; reproduced exactly | `PR` | supported |
| E152 | 16 of 24; 8/8, 8/8, 0/8; all 8 lost at discovery; 106 of 114 MAPQ 0 | E86, E90 | supported |
| E153 | 149 of 390 without SA; soft clip about half | E96, E97 | supported |
| E154 | 57.5 at all 32 breakends; strong from 70 | E104, E35 | supported |
| E155 | 0 of 20 → 12 of 20; both cloud models all 10; local 2 of 10; 5 of 15 (0 of 15, p = 0.042); 3 of 15 (p = 0.22); local models often did not call the tool | E115, E116 | supported |
| E156 | one cloud model called the evidence strong in all 5 with the limit, 0 of 5 without | E117 | supported |
| E157 | three local models malformed 13–26 % | E111 | supported |
| E158 | "capable models used the exposed quantity always" | E115 (10 of 10 stated it) | supported |
| E159 | "sensitivity limited by candidate discovery in ambiguously mapped regions, not by filter thresholds" | E89, E90 | supported |

The abstract has no † and needs none: every figure in it has a record.

## F. Phase 14 — the decoy-partner fix, 2026-09-27

`DF` is `results/decoy_fix_2026-09-27/`. Before = `bam_tools.py` at `7bb4379`,
after = the fix (`a6c887f`), registered first (`DF/registration.json`, `8ed46e1`).
Definitions beside every figure in `DF/remeasure.json` and `DF/suites.json`.

| # | Figure | Before | After | Record : key | Status |
|---|---|---|---|---|---|
| F1 | chr21:10,770,078 (NA12878 background): score, band | 72.5, strong | 55.0, moderate | `DF/remeasure.json` `public.positions[…]` | new; as registered |
| F2 | same: split reads (MAPQ 20); split component | 269; 25.0 | 3; 7.5 (266 decoy-only, `chrUn_KN707891v1_decoy`) | as F1 | new |
| F3 | same: split-read sentence | "269 split reads … predominantly to chrUn_KN707891v1_decoy (266/269)" | "3 split reads … predominantly to chr21 (3/3)" | as F1 | new |
| F4 | background survivors, 54 breakends: positions changed; split reads at MAPQ 20 / 0 | — | 1; 643 → 377 / 2,994 → 2,617 | `DF/remeasure.json` `public.aggregate.background` | new |
| F5 | 32 detected implant breakends: changed; split reads MAPQ 20 / 0 | — | 0; 250 → 250 / 350 → 346 | `…aggregate.detected` | new; nothing changes, as registered |
| F6 | 8 missed breakpoints: changed; split reads MAPQ 20 / 0 | — | 0; 0 → 0 / 65 → 55 | `…aggregate.missed` | new |
| F7 | Phase 10 locus IMP01 chr20:200,000 | — | every earlier field identical; model-visible payload identical (4,111 characters), 9 ceiling keys present | `DF/remeasure.json` `public.phase10_locus` | new; the ceiling experiment stands |
| F8 | patients, 100 Phase 13 positions each: split reads MAPQ 20 | 20 / 7 | 18 / 7; no score or band changes; one SAMPLE_A split sentence changes | `DF/remeasure.json` `SAMPLE_*` | new; totals only |
| F9 | regression suite on the pre-fix code | — | exit 1: every positive-control check fails | `DF/test_on_prefix_code_7bb4379.log.txt` | new |
| F10 | all suites after the fix | 24 suites | 25 suites, all passed; 734 assertions held, 0 failed, none NOT RUN (IGV available) | `DF/suites.json` `summary` | new |
| F11 | Phase 13 figure: background strong because of decoy partners | 72.5 | superseded by F1 | `decoy_partners_2026-09-27.json` | historical (the pre-fix measurement) |

## G. Phase 15 — 2026-09-27

`SC2` is `results/same_chrom_partners_2026-09-27.json`; `GC` is `results/guard_check_2026-09-27.json`.

| # | Figure | Value | Record : key | Status |
|---|---|---|---|---|
| G1 | guard base | a6c887f (was 9ae73bc); clean before, scratch edit to bam_tools.py flagged (exit 1), clean after revert | `GC` | new |
| G2 | same-chromosome-only split reads, 32 detected / 8 missed implant breakends (MAPQ 20; MAPQ 0) | 0 of 250 / 0 of 0 (0 of 346 / 0 of 55) | `SC2` `public.aggregate` | new |
| G3 | same, 54 background breakends | 365 of 377 (512 of 2,617 at MAPQ 0) | `SC2` `public.aggregate.background` | new |
| G4 | background breakends whose summary would change without same-chromosome partners (counterfactual) | 36 scores, 9 bands, 40 split components, 47 sentences; every one an intra-chromosomal call (DEL 34 of 40, DUP 9 of 10, INV 4 of 4), where a same-chromosome partner is the event's own split evidence | `SC2` `public.background_by_svtype`, `positions[*].svtype` | new; not a defect |
| G5 | patients, 100 positions each: same-chromosome-only split reads at MAPQ 20 | 13 of 18 (SAMPLE_A), 0 of 7 (SAMPLE_B); counterfactual changes 1 score and 1 split component, no band | `SC2` `SAMPLE_*` | new; totals only |

**The thesis, 2026-09-27.** The file in `docs/thesis/` was not replaced: it is
byte-identical to the version mapped in section E (74,617 bytes, modified
2026-09-27 01:00, sha256 2fdf7c27aa284139…), as is the Desktop copy; no newer
version was found. It has no Results subsection on decoy partners. Every flag of
E.1 still applies to it (E1–E11, E13). Two sentences are now overtaken by records
made after it: Results, "how the split-read layer weighs a decoy partner has not
yet been examined", and Limitations, "how the split-read layer weighs partners on
decoy contigs … has not been examined" — measured in
`decoy_partners_2026-09-27.json` and fixed in `a6c887f`
(`decoy_fix_2026-09-27/remeasure.json`).

## H. The final thesis draft (27 September), 2026-09-27

`docs/thesis/Rimas_MSc_Thesis.docx`, sha256 eeccb653…30ed7d (76,231 bytes; the
Desktop file `Rimas_MSc_Thesis_2026-09-27.docx`, copied; `thesis_install_2026-09-27.json`),
read only. 25,637 words, 20 † marks. Methods through both summaries, including the
Results subsection "Decoy Contigs as Split-Read Partners" and its same-chromosome
paragraph, against `decoy_partners_2026-09-27.json`, `decoy_fix_2026-09-27/` and
`same_chrom_partners_2026-09-27.json`. Every flag of section E.1 is resolved in this
draft (E1 0.62–0.71; E2–E7 and the Table 2 column now †; E8 1.6 MB — the bundle is
1,590,676 bytes; E9 "at most eleven"; E11 removed; E13 "nineteen tests" removed),
and the decoy and same-chromosome figures agree with their records.

Flagged:

| # | Where | Claim | Record | Flag |
|---|---|---|---|---|
| H1 | Results, ceiling section; the decoy subsection; Discussion | the decoy-driven background junction was "the only position rated strong anywhere in the synthetic data" / "the one strong rating anywhere in the synthetic data" | `decoy_fix_2026-09-27/remeasure.json` `public.positions`: before the fix three positions were strong — chr20:2,822,501 and chr20:17,122,806 (deletion breakends of the background, 72.5) and chr21:10,770,078 (72.5); after it the two deletion breakends remain strong | **differs**. The claim was first made in the Phase 13 report to the user, not in any record; the records never supported it. The two remaining strong ratings are deletion breakends whose split partners lie on their own chromosome — genuine evidence for a deletion (`same_chrom_partners_2026-09-27.json`) |
| H2 | Results, Deployability | "scored 43.3 instead of 47.5†" | 43.3 is, like 47.5, a figure of the first demonstration bundle (`docs/NAUDOJIMAS.md`, 2026-09-25 note); no record | minor: the † follows 47.5 only; 43.3 has no record either |

No † figure is supported by a committed record.

*Re-check, 2026-09-27 13:20 draft* (sha256 b519469a…a413115a, replacing eeccb653…30ed7d; the thesis file is not committed): H2 resolved (43.3† and 47.5†). H1 resolved in the count — the ceiling section, decoy subsection and Discussion now state three strong positions before the fix, two of them deletion breakends — but one clause still differs: the decoy subsection calls chr21:10,770,078 "the only one that was not a deletion breakend" (and the ceiling section contrasts it with "two deletion breakends"), whereas it too is a breakend of a deletion call (`same_chrom_partners_2026-09-27.json` `public.positions[*].svtype` = DEL, from the evidence chain's get_candidate); what set it apart is that its split partners lay on a decoy rather than on its own chromosome.

*Re-check, 2026-09-27, draft sha256 be6c9d72…b99709245* (replacing b519469a…a413115a; the thesis file is not committed): H1 resolved. The ceiling section, the decoy subsection and the Discussion now describe three strong positions, all deletion-call breakends (`same_chrom_partners_2026-09-27.json` `public.positions[*].svtype` = DEL), two with every split partner on the breakend's own chromosome (`mapq20.same_chrom_only` 11 of 11 and 17 of 17; `decoy_partners_2026-09-27.json` `synthetic.positions[*].mapq20.with_decoy_partner` 0) and one, chr21:10,770,078, with 266 of 269 split reads partnered only on a decoy (`mapq20.decoy_only`); every clause agrees with these fields. No flag remains open in section H.

## I. The revised conference abstract, 2026-09-28

`Rimas_tezes_DI_medicinoje_LT_2026-09-27.docx` (Desktop, read only with python-docx
1.2.0 from a scratch directory; sha256 12b8b3a9…d7a63ba1; 867 words), not yet sent.
The English copy (`Rimas_abstract_EN_2026-09-27.docx`, sha256 e77df205…467dd49a)
carries the same figures (only the thousands separator differs) and renders
"keliolikos" as "under twenty". `AN` = `results/synthetic_control_2026-09/analysis_2026-09-25/`,
`P10B` = `benchmark/runs/phase10_rerun_2026-09-25/blind/blind_scores.json`.

| # | Claim | Value | Record : key | Status |
|---|---|---|---|---|
| I1 | tools called through the Model Context Protocol with arguments, returning structured results; genome data reached only through tools | — | `stage1_igv_assistant/server.py`, `candidate_server.py` (FastMCP); `ui._chat_exec` (datasets by label, returns scrubbed) | supported |
| I2 | four measuring tools: discordant pairs, soft clips, split reads, depth | 4 | `ui_check_2026-09-27.json` `stdout` ("11 evidence tools + 4 bridge tools"); A1–A2 | supported |
| I3 | each layer 0–25; sum over applicable layers normalised to 0–100; strong from 70, moderate from 40 | — | `P10/proof_prepared_payload_2026-09-26.json` `loci[0].with_result.score_bands` (strong ≥ 70, moderate ≥ 40); `score_tiers.py` | supported |
| I4 | seven more tools: locus quality, layer applicability, score, gene, reciprocal check, IGV images; four load, filter and compare candidate sets | 7; 4 | as I2; the four are load_candidate_set, list_candidates, get_candidate, compare_candidate_sets | supported (the fourth bridge tool, get_candidate, retrieves one candidate — not named, not contradicted) |
| I5 | input: indexed BAM (positions, MAPQ, mate and supplementary locations); DELLY VCF/BCF with breakend coordinates and supporting read counts | — | `tools/bam_tools.py` (read fields used); `tools/vcf_tools.py` `Junction` (pos1, pos2, pe, sr) | supported |
| I6 | sensitivity measured in public NA12878 whole-genome data [5] | — | `na12878_provenance_2026-09-27.json`: the background is the chr20–chr21 slice (`source.from_header_last_pg.cl`, regions chr20 chr21) of `1000G_2504_high_coverage/data/ERR3239334/NA12878.final.cram`; all 12 bwa records match the NYGC settings (`nygc_high_coverage_match.all_match`) | supported; the background is a two-chromosome slice of the whole-genome alignment |
| I7 | twelve heterozygous balanced t(20;21) at known coordinates; spanning reads simulated with ART | 12 | `results/synthetic_control_2026-09/implants_ground_truth.json`; `…/art_profiles.json` (ART 2.5.8) | supported |
| I8 | models: Claude Sonnet 5, Claude Opus 5, Qwen2.5 7B, Qwen3.5 4B | — | `P10/*/meta.json`, run files `model` | supported |
| I9 | local models: 7 and 4 billion parameters, 4-bit quantisation, 8 GB graphics card, Ollama | 7, 4; 4-bit; 8 GB | `local_models_2026-09-28.json`: `models.*.ollama_show.Model.parameters` 7.6B and 4.7B, `quantization` Q4_K_M for both, IDs 845dbda0ea48 / 2a654d98e6fb (the IDs in the Phase 10 runs), `gpu` RTX 5060 Laptop 8,151 MiB, Ollama 0.32.9 | **differs slightly**: 7 and 4 billion are the model names; Ollama reports 7.6 and 4.7 billion parameters. 4-bit and 8 GB supported |
| I10 | the one-variable difference confirmed by 71 automatic checks; answers scored blind | 71 | `P10/proof_prepared_payload_2026-09-26.json` (PROVEN 71 of 71); `P10B` | supported |
| I11 | the interface: load a candidate set, see each filter's effect, open a candidate or enter a coordinate and see the four layers with the tool returns behind them; a language model in the same window with the same tools; the ceiling shown to the user | — | `ui_check_2026-09-27.json` (exit 0 at HEAD 599e742; every MINIMAL and FULL condition PASS); `ui.py` `/api/load`, `/api/funnel` (per-step and standalone filter counts), `assess()` via `/api/assess` and the call log `/api/calls`, `/api/chat` through the same recorder, ceiling rendered at `ui.py:1063–1074` | supported |
| I12 | 16 of 24 junctions; 8/8 unique, 8/8 repeat-adjacent, 0/8 low-mappability | — | `AN/figures.json` `sensitivity` | supported |
| I13 | all 8 lost at discovery; 106 of 114 spanning reads MAPQ 0 | — | `AN/figures.json` `loss_summary`, `grid`; E90 (derived from `AN/adjudication/IMP09–12.json.gz`) | supported |
| I14 | 149 of 390 spanning reads without SA; about half caught by soft clips | 77 of 149 | `AN/adjudication/summary.json` `b_short_overhang` | supported |
| I15 | ceiling 57.5 at all 32 detected breakends; strong from 70 | — | `AN/figures.json` `ceiling_summary.attainable_here_range` [57.5, 57.5] | supported |
| I16 | "cannot be rated strong"; the ceiling follows from the event: depth adds 0, and the discordant fraction (at most 0.164) does not reach the top tier (0.5) | 0.164; 0.5 | `AN/figures.json` `ceiling_summary.max_observed_discordant_fraction` 0.164; the tool's `attainable_basis` in `P10/proof_prepared_payload_2026-09-26.json` `loci[0].with_result` (top band 0.5; depth "contributing nothing, because a balanced rearrangement gains and loses no DNA") | supported under that stated premise. By arithmetic a spurious depth contribution could still lift a balanced event past 70: depth added 15 at 7 breakends (`AN/figures.json` `ceiling[*].depth_score`), and 7.5 + 25 + 25 + 15 = 72.5; the highest observed was 55.0 |
| I17 | without the ceiling in any tool return 0 of 20 runs stated it; with it 12 of 20 — cloud 10 of 10, local 2 of 10 | — | `P10B` `cells.original_runs_1_5`, `cells.api_runs_1_5` (C2) | supported |
| I18 | pre-registered extension, 15 runs: Qwen3.5 4B 5 (0 without; p = 0.042), Qwen2.5 7B 3 (p = 0.22) | — | `P10B` `cells.runs_1_15`, `fisher_exact_two_sided_15_run_cells` (0.042146, 0.224138) | supported |
| I19 | with the ceiling, Claude Sonnet 5 called the evidence strong in 5 of 5 runs although the tool rated it moderate; 0 of 5 without | — | `P10B` `cells.api_runs_1_5` (C1b 0 of 5 with, 5 of 5 without); the summary returned 40.0 "moderate" (window 500) in all five runs (`claude-sonnet-5/WITH__run*.json` `recorder_calls_unablated`) | supported |
| I20 | local models relayed the ceiling rarely | 2 of 10; 5 and 3 of 15 | I17, I18 | supported |
| I21 | about 9,000 inter-chromosomal junctions reduced to "keliolika" (11–19) candidates per genome | 9,172 / 9,655 → 17 / 19 | `patient_rerun_2026-09.json` `funnel` | supported |
| I22 | balanced translocations will score at most moderate; some low-mappability breakpoints will go undetected; sensitivity limited by discovery, not by filter thresholds | — | I15–I16 (with I16's qualification); `AN/figures.json` `grid` (no threshold removed a true positive in 12 cells), `loss_summary` | supported (a forecast from the records) |

Flags: I9 (parameter counts 7.6 and 4.7 billion, not 7 and 4); I16 (a qualification,
not a contradiction: "cannot be rated strong" holds when the depth layer behaves as
it should, which the Conclusions state). No other claim differs from a record.

*Re-check, 2026-09-28, LT sha256 5ad80c07…c5b24788 and EN 120aa565…62c6847d7709* (replacing 12b8b3a9…d7a63ba1 and e77df205…467dd49a; read only): a paragraph diff of the 27 and 28 September files shows changes only in the Methods, Results and Conclusions paragraphs, and only these, identically in both languages — I4 the fourth bridge tool named ("pateikia pasirinktą kandidatą" / "retrieve a single candidate"; get_candidate); I6 the background stated as chromosomes 20 and 21 of the public NA12878 data (`na12878_provenance_2026-09-27.json` `source.from_header_last_pg.cl`); I9 "7,6 ir 4,7 mlrd." (`local_models_2026-09-28.json` `models.*.ollama_show.Model.parameters`); I16 "gali gauti tik dėl klaidingo signalo (didžiausias stebėtas įvertis — 55,0)" (`AN/figures.json` `ceiling_summary.evidence_score_range` [32.5, 55.0], `attainable_here` 57.5 with depth at 0); I22 the forecast qualified "nesant klaidingo signalo"; and the removal of "sulygintų skaitinių" and "ir modelių nekeičiant", which leaves every claim true. All five changed clauses agree with their records; no flag remains open in section I.

## J. The plain-language abstract and the meeting sheet, 2026-09-29

Three Desktop files, read only, not yet sent; each sha256 matches the value given
with the task (`SH` `files`): `Rimas_abstract_EN_2026-09-29.docx` (ca3b3303…2387ebcf;
1,053 words), `Rimas_tezes_DI_medicinoje_LT_2026-09-29.docx` (d2471702…3ab00732; 894
words) and the meeting sheet `Rimas_tezes_paaiskinimai_2026-09-29.docx`
(c6bb20a6…51d66dc9; 1,676 words, 3 tables, 2 images). The two abstracts carry the same
figures (compared number by number over the Methods, Results and Conclusions) and the
same references, so one row covers both. The abstract was rewritten in plain language
and nearly every clause is reworded: each row gives the new wording and cites the
section I row whose record it rests on. `SH` = `results/meeting_sheet_checks_2026-09-29.json`,
`DR` = `results/demo_dry_run_2026-09-29.json`, `CH` = `AN/chain/IMP01.json.gz` (call ids),
`P10` = `benchmark/runs/phase10_rerun_2026-09-25/`, `PF` = `P10/proof_prepared_payload_2026-09-26.json`;
`AN` and `P10B` as in section I. The two scripts (`scripts/abstract_records.py sheet`,
`scripts/demo_dry_run.py`) were committed before they ran.

### J.1 The abstracts (EN wording; the LT says the same)

| # | Claim (29 September) | Value | Record : key | Status |
|---|---|---|---|---|
| J1 | tools are "small programs that the model runs on a patient's aligned sequencing reads (a BAM file) …, such as 'count the read pairs at chromosome 20, position 200,000', and that return exact measurements" | — | the MCP tool `discordant_pairs` (`count_discordant_pairs`, `bam_tools.py:784`: pairs whose mate maps to a different chromosome, ±500 bases); at IMP01 chr20:200,000 it returns 44 (`DR` `steps.4_tool_return`; `CH` call 323) | supported; the example tool counts only the pairs whose mate lies on another chromosome |
| J2 | "The model sees the data only through these tools … so every number it reports can be checked against a tool output" | — | I1; `chat.verify_numbers` (every number in the final answer checked against the turn's tool returns) | supported |
| J3 | four tools measure "independent signs": pairs whose ends map to different chromosomes, reads cut off at the breakpoint, reads split between the two chromosomes, depth | 4 | I2; `bam_tools.py:1808` (each layer "scored independently on a 0-25 scale") | supported for the scores. The reads can overlap: a split read's primary alignment carries a soft clip where it splits |
| J4 | each sign 0–25 points; total scaled to 0–100; 70 or more strong, 40 or more moderate | — | I3 | supported |
| J5 | other tools "check data quality, name the gene at the breakpoint, check the partner breakpoint and draw IGV images" | — | I4 | supported (the applicability and score tools and the four candidate-set tools are not named) |
| J6 | "The assistant does not scan the genome itself: … DELLY scans the BAM file and lists candidate breakpoints, which the assistant filters and assesses, as it does any position the user enters" | — | `vcf_tools.load_candidate_set`, `list_candidates` (the filters); `ui.assess` runs the same four layers for a candidate and a typed position (`DR` steps 3 and 5); the limits panel's first item, "It checks positions. It does not search for them." (`ui.build_limits`) | supported |
| J7 | "twelve balanced translocations between chromosomes 20 and 21 were inserted at known positions into public whole-genome data of NA12878 …, with reads across the new junctions simulated by ART" | 12 | I6, I7 | supported with I6's note: the background is the chr20–chr21 slice of the whole-genome alignment (the Limitations say "two chromosomes of one public genome") |
| J8 | four models; the open ones "run through Ollama on an ordinary computer with an 8 GB graphics card, so the data never leave it" | 8 GB | I8, I9 (`local_models_2026-09-28.json` `gpu` 8,151 MiB) | supported; the parameter counts are no longer given |
| J9 | "Each model answered the same question about one breakpoint five times with each of two versions of the tool output, identical except for one piece of information (confirmed by 71 automated checks)" | 5 × 2; 71 | `P10B` `cells.original_runs_1_5`, `cells.api_runs_1_5` (`runs` 5 per cell); I10 | supported; the one piece is the ceiling, carried in 10 fields (`PF` `ceiling_keys`) |
| J10 | "an AI evaluator that did not know which version it was reading graded the answers against fixed criteria" | — | `P10B` `scorer` ("verifier subagent, fresh context, blind: read the packet only"), `method` (the packet holds the prompt and the final answer only: no model, condition, file name, tool call or tool return), `criteria` C1, C1b, C2, C3, fixed before the extension (`P10/registration_local_extension_2026-09-26.json` `analysis_fixed_in_advance`) | supported; the evaluator did not know the model either |
| J11 | the interface: the user "loads the candidate list, sees how many candidates each filter removes, opens a candidate or types a position, and sees the four measurements, each linked to the tool output behind it; a language model can be asked questions in the same window" | — | I11; walked on 29 September (J.3) | supported |
| J12 | 16 of 24 junctions found: all 8 in unique sequence, all 8 next to repeats, none of the 8 in low-mappability regions "where reads cannot be placed unambiguously" | — | I12 | supported |
| J13 | "missed by the candidate search, not removed by the filters: 106 of the 114 reads crossing them could not be assigned to a single location" | 106 / 114 | I13 (MAPQ 0: an equally good alignment elsewhere) | supported |
| J14 | "at all 32 breakpoint ends found, the highest score the tool could give was 57.5, below the 70 needed for 'strong' (only a false copy-number signal from the depth measurement could lift it higher; the highest score observed was 55.0)" | 57.5; 55.0 | I15, I16 and the 28 September re-check | supported at the default ±500 window; the ceiling moves with the window (65.0 at ±200 at chr21:14,100,000, `PF` `loci[1]`; 65.0 at ±50 at chr20:200,000 in the dry-run chat, `DR` `observations.chat_vs_panel`), still below 70 |
| J15 | without the ceiling 0 of 20 answers said strong was impossible; with it 12 of 20, "all 10 from the cloud models and 2 of 10 from the open models" | — | I17 | supported |
| J16 | "A balanced translocation neither adds nor removes DNA, so the depth measurement correctly scores 0" | 0 | `PF` `loci[0].with_result.attainable_basis`; `DR` step 3 (depth component 0 at both IMP01 ends) | supported as what the layer should do; in the synthetic data it added 15 at 7 of the 32 detected breakends (I16), and 15 of the 22.5 at the hand-entered IMP09 end (`DR` `steps.5_hand_entry`) |
| J17 | "because only one of the two copies of each chromosome is rearranged, about half of the reads at the breakpoint come from the normal copy and look normal, so the share of abnormal read pairs (at most 0.164) stays far below the 0.5 needed for the top score" | 0.164; 0.5 | construction: half the pairs spanning each breakpoint were removed and replaced by junction fragments (`implants_ground_truth.json` `implants[0].compensation.removal`: 34 of 68 at chr20:200,000, 24 of 49 at chr21:14,100,000); `AN/figures.json` `ceiling_summary.max_observed_discordant_fraction` 0.164; top tier 0.5 (`score_tiers`) | **qualification**: the numbers hold at ±500, but the share also depends on the window. At chr21:14,100,000 it is 0.164 at ±500 (`CH` call 303) and 0.25 at ±200 (`PF` `loci[1]`). The normal copy explains why the share cannot pass about 0.5, and the window explains why it sits at 0.16–0.25; "far below 0.5" holds at both windows. The tool's own `attainable_basis` gives the same one-factor account |
| J18 | Novelty: an assistant that "ties every number it reports to a tool output and shows its own limits to the user" | — | evidence panel: each layer count and the combined score carries a call chip (`ui.py` `layerBlock`, `renderEvidence`; `DR` step 4: 44 and 40.0 equal their returns); chat panel: `chat.verify_numbers`, `ui.py` `markUnsupported`; limits: as J42 | supported with J41's qualification: in the chat panel a number is matched by value to any return of the turn, and one with no match is shown, marked |
| J19 | Novelty: "an experiment showing that a language model can explain a limit of the analysis only when a tool reports it" | 0/20 → 12/20 | I17, I18 | supported for this experiment (one question at one locus; 5 or 15 runs per cell) |
| J20 | Limitations: simulated translocations on two chromosomes of one public genome; low-mappability breakpoints missed by the candidate search | — | I6, I7, I12, I13 | supported |
| J21 | Limitations: "10 of the 11 scoring thresholds are judgement calls not yet calibrated on confirmed cases" | 10 of 11 | `bam_tools.py:1883–1900` (`breakpoint_evidence_summary`, THRESHOLD PROVENANCE: "of the 11 tier-cutoff values … only ONE … is empirically calibrated, and against a single confirmed real locus … The other 10 … are HEURISTIC"); `server.py:318–319` | supported. The interface's limits panel counts the whole inventory instead: "14 of the 16 cut-offs are judgement calls" (E64; `DR` `steps.1_load.limits_panel_headings`); a viewer of the demo sees 14 of 16 |
| J22 | Limitations: "the open models used the tools unreliably" | — | `phase8_final_record.json` `models.*.schema_invalid_pct` (qwen2.5:7b 13.4, qwen3.5:4b 26.3; 30 runs each); `P10B` `cells.*.ceiling_fields_reached_model` (WITH runs whose ceiling reached the model: 4 and 3 of 5, 11 and 8 of 15); 1 of 7 calls rejected in the dry run (`DR` `steps.7_chat`) | supported |
| J23 | Limitations: "one cloud model, given the ceiling, called the evidence strong in all 5 answers although the tool rated it moderate" | 5 of 5 | I19 (Claude Sonnet 5: `P10B` `cells.api_runs_1_5` C1b 0 of 5 WITH) | **incomplete**: the other cloud model, Claude Opus 5, called the evidence strong in 2 of 5 answers with the ceiling and 2 of 5 without (`P10B` `items`, claude-opus-5 WITH__run1, WITH__run5, WITHOUT__run3, WITHOUT__run5: C1b false; e.g. "strong in substance but labelled moderate"). I19 missed this too |
| J24 | "in these genomes the candidate search and filters already reduced about 9,000 possible joins between chromosomes to fewer than twenty per genome" | 9,172 / 9,655 → 17 / 19 | I21 | supported ("about 9,000" for 9,172 and 9,655) |
| — | "Interpretable AI systems should therefore report the facts needed to check their conclusions, not only the conclusions" | — | — | a recommendation, not a factual claim |

### J.2 The meeting sheet

| # | Claim | Value | Record : key | Status |
|---|---|---|---|---|
| J25 | DELLY reads the whole BAM "once per sample, about 1–2 hours" | 1–2 h | `patient_rerun_2026-09.json` `delly.SAMPLE_*.wall_clock` 1:01:48 and 1:44:54 | supported |
| J26 | filter table: DELLY junctions between chromosomes 9,172 / 9,655; PASS 896 / 923; ≥ 3 discordant pairs 664 / 642; ≥ 1 split read 57 / 63; none found (±500) "in another, unrelated genome" 17 / 19 | as stated | `patient_rerun_2026-09.json` `funnel.SAMPLE_*.per_step` (svtype, filter_pass, min_pe, min_sr, after_recurrence_500bp) | every count matches. **First row**: 9,172 / 9,655 are the inter-chromosomal junctions after the bridge merged records within 500 bases; DELLY wrote 9,187 / 9,676 BND records (`bnd_records`). Two steps that removed none (both ends on primary contigs; exclude template) are not listed. **Last row**: the other genome is the other patient's (`survivors_recurrent_in_other_sample` 40 / 44); "unrelated" has no record |
| J27 | "Only aggregate counts"; the candidates (coordinates, genes) not yet looked at | — | `patient_rerun_2026-09.json` `blinding` | supported |
| J28 | a known position can be typed in, but "the evidence tools look only at a ±200–500 base window around it, so it must be accurate to within a few hundred bases" | ±200–500 | `ui.assess`: discordant pairs and the summary ±500, soft clips and split reads ±200 (tool defaults), read depth ±2,000 (`read_depth_profile` start and end; `summarize_breakpoint_evidence` `depth_window_bp=2000`), a dip counting only within 1,000 bases of the position (`dip_tolerance_bp`) | **differs** for the depth layer (±2,000); the conclusion (a few hundred bases) holds for the three read layers |
| J29 | the candidate list cannot at present be filtered by region | — | `vcf_tools.list_candidates(set_id, svtype, filter_pass, min_pe, min_sr, primary_only, mask_path, limit, offset)`: no region parameter; the page offers none | supported |
| J30 | if DELLY misses a breakpoint, as 8 of 24 synthetic junctions in low-mappability regions, the assistant does not see it unless the exact position is typed in | 8 of 24 | I12, I13; typed in on 29 September, IMP09 chr20:33,700,000 returned 22.5 weak (`DR` `steps.5_hand_entry`, as `AN/figures.json` `rescue`); the other seven ends are withheld (QUALITY-LIMITED, `rescue`) | supported |
| J31 | the IGV tool makes one image per evidence layer; two of the four at IMP01 chr20:200,000 (public data); "a vertical line marks the position" | 4 | `imp01_panels_2026-09-27.json` (four PNGs, one per layer); `SH` `sheet_images`: the first image is rows 0–499 and the second rows 0–469 of the committed `discordant_pairs.png` and `soft_clipped_reads.png`, pixel-exact at full width (`image_search_controls.pass`; no display crop in Word); the line at the 200,000 tick is visible in the committed panel | supported |
| J32 | discordant-pair caption: coloured reads are pairs whose mate maps to another chromosome, the colour naming it; they cluster on both sides of the breakpoint because the translocation is reciprocal, with reads from one derivative on one side and from the other on the other | — | `SH` `captions.discordant_pairs_panel`: in chr20:198,500–201,500, 25 reads end left of the breakpoint, all from der(20), all forward, and 19 start right of it, all from der(21), all reverse; every mate is on chr21, MAPQ ≥ 20. The 6 other inter-chromosomal reads (mates elsewhere) are also there before implanting (negative control) | supported; the colours are IGV's (`color_by` UNEXPECTED_PAIR) |
| J33 | soft-clip caption: the clipped, coloured parts begin at exactly the same place, the breakpoint | — | `SH` `captions.soft_clipped_reads_panel`: 13 of 14 primary non-duplicate clipped alignments in chr20:199,850–200,150 clip after base 200,000; the tool (±200, MAPQ ≥ 20, clips ≥ 10 bases, supplementary alignments counted): 19 of 20 at 200,000 (`CH` call 324 `consensus_clip_position`, `max_clips_at_position`) | supported |
| J34 | section 3: the ablation (one question, two tool-output versions differing only in the ceiling); the evaluator, a "Claude Code verifier subagent" without earlier context, blind to the version, agreeing with the earlier unblinded scoring in 22–23 of 23; the 15-run extension (5 and 3; p = 0.042 and 0.22) | 22–23 of 23 | J9, J10; `P10B` `agreement_with_unblinded_2026_09_25` (C1 and C1b 22 of 23, C2 and C3 23 of 23: the 23 local runs scored unblinded on 25 September); I18 | supported |
| J35 | section 3, bands: "70 and more strong, 40 and more moderate, less weak" | — | I3 | **differs slightly**: a score of 0 is "none", not weak (the glossary row has it right: "> 0 weak") |
| J36 | section 3: "57.5 = discordant pairs 7.5 + clipped reads 25 + split reads 25 + depth 0"; above 70 only with a false depth signal (7.5 + 25 + 25 + 15 = 72.5); in the synthetic data depth wrongly added 15 at 7 breakends, but no score passed 55.0 | 57.5; 72.5; 7; 55.0 | `PF` `loci[0].with_result` (`attainable_here` 57.5, `attainable_basis`: discordant held at its band, 7.5 in `DR` step 3's components); `AN/figures.json` `ceiling_summary` (7 breakends with a depth contribution; `evidence_score_range` [32.5, 55.0]) | supported at the default window (J14) |
| J37 | section 3: about half the reads covering the breakpoint come from the normal copy, so the share of unusual pairs stays small ("at most 0.164") and "does not reach even the 0.2 tier" | 0.164; 0.2 | as J17 | supported at ±500; at ±200 the same breakend reaches 0.25, above the 0.2 tier (`PF` `loci[1]`) |
| J38 | section 3: DELLY found about 9,000 junctions between chromosomes in each genome "(mostly artefacts)"; 17 and 19 after the filters | — | I21 | "mostly artefacts": **no record**, and it reads as a judgement on the patient calls, which stay blinded. What is recorded is that only 896 and 923 pass DELLY's own quality filter (`per_step.filter_pass`) |
| J39 | section 3: three forecasts (few candidates; no more than moderate unless depth errs; some low-mappability breakpoints undetected, lost at discovery) | — | I22 | supported (forecasts from the records) |
| J40 | Novelty: read-level evidence at breakpoints "including balanced translocations, which most existing interpretation tools do not examine" | — | — | **no record** (a literature claim) |
| J41 | Novelty: "every number the assistant reports is linked to the tool return behind it; the interface does not show a number without such a return" | — | evidence panel as J18; chat panel: the model's prose is shown in full ("model — prose, not verified except where marked"), a number with no matching return is shown, marked "no tool call returned this number" (`ui.py` `markUnsupported`), and listed in the verification pass; matching is by value against every return of the turn (`chat.verify_numbers`: within 0.011, roundings, percentages), not a link to one call | **differs** for the chat panel: such numbers are shown, marked, not withheld |
| J42 | Novelty: the system shows its own limits: the highest attainable score, the withheld score (QUALITY-LIMITED) and the limits of the candidate search | — | `ui.py` `ceilBlock` ("Highest score reachable at this position", `DR` step 3), the "withheld — quality limited" badges (`layerBlock`, `renderEvidence`), the limits panel ("It checks positions. It does not search for them.") | supported |
| J43 | Novelty: 0 of 20 → 12 of 20; Limitations: sensitivity from 12 simulated translocations on two chromosomes of one public genome, no confirmed real cases yet; DELLY misses low-mappability breakpoints (0 of 8); 10 of 11 thresholds set by the author's judgement; blind scoring by an AI agent, not a person | — | I17; I7; I12 (and `AN/ladder/analysis.json` `summary.q0_r0`: with DELLY's mapping-quality floors at 0 it called 1 of the 8, and none survived the filters); J21; `P10B` `scorer` | supported |
| J44 | Limitations: a balanced translocation cannot score "strong": "the scoring function does not distinguish a layer that does not apply to this event from a layer whose data are missing" | — | `applicable_layers` decides from the BAM, not the event: depth is "always applicable — any aligned BAM can show a coverage drop" (`CH` call 20); a layer the data cannot show is dropped from the normalisation; `max_with_flat_depth` 75 (`DR` step 3) | supported in substance; **wording**: the function does set aside a layer whose data are missing (no paired reads, no SA tags). What it cannot tell apart is a layer the event cannot move (depth, for a balanced event) from a layer that looked and found nothing |
| J45 | Limitations: few runs (5 or 15 per condition); local models called tools with wrong arguments in 13–26 % of calls; "one cloud model re-rated the tool's score" | 5, 15; 13–26 % | `P10B` `cells`; `phase8_final_record.json` `models.*.schema_invalid_pct` 13.4 and 26.3 (Phase 8; qwen3.5:9b's 17.1 lies between); J23 | first two supported; the third **incomplete** as J23 (the other cloud model did so in 4 of its 10 answers) |
| J46 | glossary: reads of 150 bases; pair ends a few hundred bases apart; a discordant pair has its ends on different chromosomes; a clip marks the breakpoint; a split read carries an SA tag; MAPQ 0 means several equally good places | 150 | `implants_ground_truth.json` `implants[*].simulation` (read length 150, fragment mean 443, sd 105); `count_discordant_pairs`; J33; `get_split_reads` (SA) | supported (150 is recorded for the simulated reads; the background's read length is in no record) |
| J47 | glossary: "12 translocations = 24 junctions = 48 breakends (32 detected)" | 12, 24, 48, 32 | I12; `AN/figures.json` `ceiling_summary.breakends` 32 | supported |
| J48 | glossary: bands "≥ 70 strong, ≥ 40 moderate, > 0 weak" | — | I3 | supported |
| J49 | glossary: the ceiling is the highest score the tool can give at a particular position; "for a heterozygous balanced translocation, 57.5" | 57.5 | `AN/figures.json` `ceiling_summary.attainable_here_range` [57.5, 57.5] (the 32 detected breakends, default window); `PF` `loci[1].with_result.attainable_here` 65.0 (±200) | **differs** as a general statement: 57.5 is the default-window value at every detected breakend of the synthetic control; the ceiling is computed per position and window |
| J50 | demo: `python -m stage1_igv_assistant.ui`, http://127.0.0.1:8765; the interface finds only public data and shows no patient data | — | `DR` (J.3) | supported; autodiscovery reads only ~/public_data, but a config file can register other data (`ui.discover_public` reads `CFG.registered` first), and none exists on this machine |
| J51 | demo: load a synthetic implant's candidate set (IMP01) and show how many candidates each filter removes | — | `DR` steps 1–2 | supported; the page's default type is "any", so the chain starts from all 896 IMP01 junctions, not from BND as in the patient table |
| J52 | demo: the IMP01 candidate (chr20:200,000 and chr21:14,100,001), four layers at both ends, the chr20 end 40 "moderate", ceiling 57.5 and its explanation | 40; 57.5 | `DR` `steps.3_candidate` (CAND_52179e966345; both ends 40.0 moderate, components 7.5, 25, 7.5, 0; `attainable_here` 57.5 with `attainable_basis`) | supported |
| J53 | demo: "clicking any number opens the tool return behind it" | — | `ui.py:940` `chip`: a "tool call #N" button beside each layer count and the combined score (`showCall`); one chip for the whole filter table; none in the ceiling block or the limits panel; "Where every number came from" lists every call | **differs slightly**: one clicks the chip beside a number, and not every number has one |
| J54 | demo: a missed implant's position typed in; the interface shows what the reads show and marks the position as entered by the user | — | `DR` `steps.5_hand_entry` ("This position was entered by hand."; `position_provenance` caller_supplied) | supported |
| J55 | demo: IGV images, four layers | 4 | `DR` `steps.6_igv` | **differs on 29 September**: three of four produced; the read-depth panel timed out (J.3) |
| J56 | demo: a question to a local model or to Claude, with its tool calls | — | `DR` `steps.7_chat` (qwen2.5:7b); Claude not exercised (no API call in this run) | supported for the local model |
| — | section 7 (questions for the meeting); the authors, affiliations and e-mail line | — | — | not mapped (no repository claim; e-mail addresses stay out of the repository) |

### J.3 Demo dry run, 2026-09-29 (`DR`)

The interface was started as the sheet describes (`python -m stage1_igv_assistant.ui`,
repository root, no flags). There was no config file and no SV_* variable, and
ANTHROPIC_BASE_URL pointed at a closed local port (positive control passed). It
answered in 1.7 s and was stopped by its PID (SIGTERM; no child process left; port
8765 free after). The run took 404 s, 378 s of it in the IGV step. The seven steps:

1. Listed: datasets IMP01–IMP12 and NA12878.chr20_chr21 (13); candidate sets IMP01–IMP12
   and background (13). Loaded: IMP01 (918 records, 896 junctions after the merge).
2. Filters at the page's defaults: 896 → PASS 515 → PE ≥ 3 157 → SR ≥ 1 29 → primary 29
   → exclude template 29.
3. CAND_52179e966345 (chr21:14,100,001 ↔ chr20:200,000; PASS, PE 25, SR 6): both ends
   40.0 moderate (7.5, 25, 7.5, 0), ceiling 57.5 with its basis, and the page's sentence
   that the position cannot reach strong.
4. The discordant count 44 (call 13) and the combined score 40.0 (call 17) equal their
   tool returns.
5. IMP09 chr20:33,700,000 typed in: 22.5 weak (7.5, 0, 0, 15), marked "This position was
   entered by hand."; 10 discordant pairs and 7 clipped reads, as the committed rescue
   record. 15 of its 22.5 points are depth, where 36.4 % of reads have MAPQ < 20 and no
   copy-number change was built in.
6. IGV: the discordant-pair, soft-clip and split-read panels were produced (the
   discordant-pair PNG is byte-identical to the committed one; the other two differ in
   size from theirs). The read-depth panel timed out at 180 s. IGV's log shows each
   panel's IGV loading the hg38 genome from igv.org and the RefSeq track from UCSC over
   the internet. The fourth stopped while still loading RefSeq, before the BAM
   (`observations.igv_log_during_run`).
7. qwen2.5:7b through Ollama, one question, 22.5 s, 7 tool calls (1 rejected: `window_bp`
   passed to `read_depth_profile`). The verification pass found 7 numbers, none
   unsupported. The model chose a 50-base window (47.5 moderate, ceiling 65.0) where the
   panel uses 500 (40.0, ceiling 57.5). Its answer calls the evidence "strong" while
   quoting the tool's "moderate", says the position "could theoretically reach" strong,
   and does not give the ceiling (`observations.final_text_mentions`: strong 4,
   moderate 1, 65.0 0).

No patient dataset was listed or loaded. The listed labels equal those derived from
~/public_data; no label contains an identifier (positive control caught); the 31 files
named in the call log are all listed labels; and no open file lay outside the allowed
locations in 8 samples (`no_patient_data`).

### J.4 Flags

1. J28: "±200–500 bases": the depth layer uses ±2,000 (`ui.assess` `read_depth_profile`; `depth_window_bp=2000`).
2. J26, first row: 9,172 / 9,655 are counted after the bridge's 500-base merge; DELLY's BND records are 9,187 / 9,676 (`bnd_records`).
3. J26, last row: "unrelated" has no record; the other genome is the other patient's.
4. J38: "(mostly artefacts)" has no record and judges blinded patient calls; the recorded fact is 896 / 923 PASS.
5. J35: the bands row omits "none" (a score of 0).
6. J49: the glossary gives 57.5 as the ceiling of any heterozygous balanced translocation; it is the default-window value here (65.0 at ±200).
7. J41: "the interface does not show a number without such a return": the chat panel shows such numbers, marked.
8. J40: "most existing interpretation tools do not examine" balanced translocations: no record.
9. J23, J45: "one cloud model": Claude Opus 5 also called the evidence strong, in 4 of its 10 answers (2 with the ceiling, 2 without); I19 did not note it.
10. J17, J37: the account of the 0.164 share leaves out the window (0.25 at ±200 at the same breakend); "does not reach even the 0.2 tier" holds at ±500 only.
11. J44: "a layer whose data are missing": the function does set such a layer aside; the case it cannot tell apart is a layer that found nothing.
12. J53: "clicking any number": one clicks the chip beside a number, and not every number has one.
13. J55: in the dry run three of four IGV panels were made; the read-depth panel timed out while IGV was still downloading RefSeq. Every panel downloads the genome and RefSeq over the internet.
14. J.3 step 7: the local model's numbers matched its tool returns, but its words did not ("strong", "could theoretically reach"); the number check cannot catch this.

Notes, not flags: J3 ("independent" holds for the scores, not the reads); J16 (depth
"correctly scores 0" is its intent: it added 15 at 7 of the 32 detected breakends); J21
(the limits panel says 14 of 16 cut-offs, the abstract 10 of 11 thresholds; both are
recorded).

### J.5 Re-check of the revised files (v2), 2026-09-29

`Rimas_abstract_EN_2026-09-29_v2.docx` (9def8b94…c659328f), `Rimas_tezes_DI_medicinoje_LT_2026-09-29_v2.docx`
(0c583660…57a06c74) and `Rimas_tezes_paaiskinimai_2026-09-29_v2.docx` (343642d8…04b5dee6):
Desktop, read only, hashes as given with the task. Two paragraph diffs against the
v1 files agree. The first compares body paragraphs, table rows, images, headers and
footers; the second compares every paragraph in `word/document.xml`, table cells
included. In each file only `word/document.xml` and the file metadata
(`docProps/core.xml`) differ; there are no text boxes, tracked changes, comments or
fields, and the sheet's two images are byte-identical to v1. In the EN and LT only the
Methods, Results and Conclusions paragraphs changed, with the same edits and identical
figures in both languages. In the sheet the 13 changed regions are the intended
changes, including the new question in section 7 (J40: compare the novelty wording with
the literature, naming four published tools; a question, not a claim). The only other
edits are the version labels ("v2") in the header and in the heading of the section 3
table. Nothing else changed.

| J.4 | Revised wording (v2) | Record : key | Status |
|---|---|---|---|
| 1 (J28) | "three read layers look only at a ±200–500 base window … (the depth layer ±2,000 bases)" | `ui.assess` (`read_depth_profile` ±2,000); `depth_window_bp=2000` | agrees |
| 2 (J26) | two rows: DELLY's BND records 9,187 / 9,676; "the same, with duplicate records merged (±500 bases)" 9,172 / 9,655 | `patient_rerun_2026-09.json` `funnel.SAMPLE_*.bnd_records`, `bnd_after_dedup`; `funnel_filters.dedup_tolerance_bp` 500 | agrees |
| 3 (J26) | "found (±500 bases) in the other patient's genome" | `survivors_recurrent_in_other_sample` 40 / 44; `recurrence_tolerance_bp` 500 | agrees |
| 4 (J38) | "only 896 and 923 passed its own quality filter" | `per_step.filter_pass` | agrees |
| 5 (J35, J48) | "above 0 weak, 0 no evidence", in the section 3 table and the glossary | `bam_tools.py:2243–2250` (> 0 weak, otherwise "none") | agrees |
| 6 (J49; J36 row) | the ceiling is the highest score "at a given position with a given window … at the default ±500 window 57.5 at all 32 ends; at a narrower window 65.0; below 70 in both cases" | `AN/figures.json` `ceiling_summary.attainable_here_range` [57.5, 57.5]; 65.0 at ±200 at chr21:14,100,000 (`PF` `loci[1]`) and at ±50 at chr20:200,000 (`DR` `observations.chat_vs_panel`) | 57.5 agrees. **Partly open**: 65.0 is recorded at those two IMP01 ends only, not at the other 30 ends or at every narrower window. No Phase 10 summary call used a window below 500 |
| 7 (J41; J18) | sheet: "in the chat window every number of the model is checked against the tool returns; a number no tool returned is shown, but marked"; abstracts: "checks every number it shows against a tool output" | `chat.py:495`, `chat.py:767` (`verify_numbers(final_text, …)`); `ui.py` `sendChat`, `markUnsupported` | the flagged point agrees (shown, marked). **Partly open**: the check covers the model's final answer only; text from its intermediate turns is shown unmarked ("model — prose, not verified except where marked"), and the limits panel's figures come from committed records, not tool outputs. "every number in the model's answer" would be exact |
| 8 (J40) | the literature clause removed from Novelty; a question in section 7 | — | agrees |
| 9 (J23, J45) | "both cloud models at times called the evidence strong although the tool rated it moderate (Claude Sonnet 5 in all 5 answers given the ceiling, Claude Opus 5 in 4 of 10)", in the EN, the LT and the sheet | `P10B` `cells.api_runs_1_5` (Sonnet WITH C1b 0 of 5); `items` (Opus WITH__run1, WITH__run5, WITHOUT__run3, WITHOUT__run5: C1b false); the summary returned 40.0 "moderate" at the default window in all nine runs (`recorder_calls_unablated`) | agrees |
| 10 (J17, J37) | abstracts: "(at most 0.164 at the default window)"; sheet: "at the default ±500 window at most 0.164 (not reaching even the 0.2 tier), at ±200 at the same place 0.25" | `AN/figures.json` `ceiling[*].observed_discordant_fraction` (the maximum, 0.164, is at IMP01 chr21:14,100,000 / 14,100,001); 0.25 at ±200 there (`PF` `loci[1]`) | agrees |
| 11 (J44) | "the scoring function drops a layer its data cannot show (no pairs, no SA tags) but does not tell a layer this event cannot change (depth, for a balanced translocation) from a layer that looked and found nothing" | `CH` call 20 (`applicable_layers` decided from the BAM; depth "always applicable") | agrees |
| 12 (J53) | "the 'tool call #N' button beside a layer's number or the combined score …; not every number has one; every call is listed under 'Where every number came from'" | `ui.py:1057` (layer chips "tool call #N"); `ui.py:1125` (the combined score's chip reads "breakpoint_evidence_summary #N"); `showCalls` | **minor, open**: the combined score's button is labelled "breakpoint_evidence_summary #N"; the rest agrees |
| 13 (J55) | "show the pre-made images (stage1_igv_assistant/screenshots/imp01_2026-09-27/); live, each image downloads the genome and the RefSeq track from the internet; in the dry run the depth image did not finish within 180 s" | the four committed PNGs (`imp01_panels_2026-09-27.json`); `DR` `observations.igv_log_during_run`, `steps.6_igv.panel_errors.read_depth` | agrees |
| 14 (J.3 step 7) | a new limitation: "the number check checks only numbers, not words: in the dry run qwen2.5:7b gave all 7 numbers as the tools returned them but called the evidence 'strong' although the tool said 'moderate'"; the demo line: window 50, 47.5, ceiling 65.0 | `DR` `steps.7_chat.numbers_matched` (7, none unsupported), `final_text`, `observations.chat_vs_panel` | agrees |
| note J3 | "measure, and score separately" | `bam_tools.py:1808` | agrees |
| note J16 | depth "should score 0"; "such a signal added 15 points at 7 ends, while the highest score observed was 55.0" | `AN/figures.json` `ceiling[*].depth_score` (15 at 7 of the 32 ends, 0 at 25); `ceiling_summary.evidence_score_range` [32.5, 55.0] | agrees |
| note J21 | abstracts: "14 of the tool's 16 thresholds are judgement calls, not calibrated values"; the sheet adds "(the interface's limits panel says so too); of the 11 tier cut-offs only one is calibrated, against one confirmed locus" | `bam_tools.py:32–39` (16 thresholds, 2 empirically derived); `ui.build_limits`; `bam_tools.py:1885–1889` | agrees |

The other intended changes agree with their records:
- J1: the tool example is removed.
- J10: the evaluator knew neither the version nor the model (`P10B` `method`).
- J14: "at its default window".
- J50–J52 and J54: the demo lines quote `DR`. Datasets IMP01–IMP12 and the NA12878 background; no config file; 896 → 515 → 157 → 29; 40.0 = 7.5 + 25 + 7.5 + 0 with ceiling 57.5; 44 and 40.0 equal their returns; IMP09 22.5 weak, 15 points of it depth; 404 s.

**New, open (found while checking J14 and J49).** The ceiling of 57.5, and "only a false
copy-number signal could lift it higher", hold when all four measurements are counted,
as the interface and the evidence chain count them. But the summary tool also accepts a
caller's own list of applicable layers. In one Phase 10 run the model declared two:
`P10/qwen3.5-4b/WITH__run15.json`, `recorder_calls_unablated`, the call with
applicable_layers [discordant_pairs, soft_clipped_reads]. It returned 65.0 "moderate" at
the default window, normalised over those two layers ((7.5 + 25) / 50). The same return
reported `attainable_here` 57.5 and "the top band is UNREACHABLE here". By the same
arithmetic, soft clips and split reads counted alone could reach 100. The recorded
figures stand, since they count four layers. The ceiling field does not follow a
caller's layer list; that is a matter for the protected tool, which is unchanged.

Still open: the new item (a qualifier such as "counting all four measurements" would
cover it); 6 (65.0 recorded at IMP01's two ends only); 7 (the check covers the final
answer only); 12 (the combined score's button label). Flags 1–5, 8–11, 13 and 14 and
all three notes are resolved.

### J.6 Re-check of the v3 files, 2026-09-29

`Rimas_abstract_EN_2026-09-29_v3.docx` (5ad2fc7e…b6b9ee57), `Rimas_tezes_DI_medicinoje_LT_2026-09-29_v3.docx`
(fa4f3255…1f7530fe) and `Rimas_tezes_paaiskinimai_2026-09-29_v3.docx` (4e7098eb…09a3c8d1):
Desktop, read only, hashes as given with the task. Two paragraph diffs against v2 agree
(body blocks; every paragraph of `word/document.xml`, table cells included). In each
file only `word/document.xml` and `docProps/core.xml` differ; there are no text boxes,
tracked changes, comments or fields, and the sheet's images are byte-identical. The EN
changed in its aim sentence (shortened; no figure), the Results and the Conclusions; the
LT in the Results and the Conclusions; the figures are identical in both languages. The
listed removals ("each linked to the tool output behind it", "(two per translocation)",
"already", and their LT forms) carried no figure, and the first survives in the Novelty.
The sheet changed in 8 regions: the items below, the new limitation and the version
labels ("v3"). Nothing else changed.

| J.5 open item | v3 wording | Record : key | Status |
|---|---|---|---|
| new: the ceiling and the layer list | abstracts: "the highest score the tool could give at its default window, counting all four measurements, was 57.5"; sheet: the section 3 row ("at the default ±500 window and counting all four layers"), the glossary ("with a given window, counting all four layers") and the J44 sentence ("at the default window, counting all four layers") | `AN/figures.json` `ceiling_summary`, counted over four layers (`applicable` lists all four at all 12 implants, `AN/chain/IMP*.json.gz`) | agrees |
| the new limitation | "the ceiling field ignores which layers the model asked to count: when qwen3.5:4b in one run asked for discordant pairs and soft clips only, the tool returned 65.0 at the default window although the ceiling read 57.5. The interface always counts all four layers, so its numbers are correct; the tool's bug will need fixing" | `P10/qwen3.5-4b/WITH__run15.json` `recorder_calls_unablated` (applicable_layers [discordant_pairs, soft_clipped_reads], window 500, evidence_score 65.0, attainable_here 57.5); `ui.py:238` (the evidence panel passes the layers `applicable_layers` finds in the BAM); `chat.resolve_args`, `ui._chat_exec` (a model's own layer list reaches the tool unchanged); `benchmark/phase10_rerun.py:39–40, 198` (the Phase 10 runs, run 15 among them, called their tools through `ui._chat_exec`) | the run-15 facts agree. **Open**: "the interface always counts all four layers, so its numbers are correct", and the section 3 row's "(the interface always counts all four)", hold for the evidence panel only. In the chat panel the model sets the tool's arguments, and run 15 went through that same executor |
| 6 | glossary: "at IMP01 with ±200 and ±50 base windows 65.0; below 70 in every case"; section 3 row: "at IMP01 at narrower windows (±200 and ±50 bases) the ceiling was 65.0, still below 70" | `PF` `loci[1]` (±200 at chr21:14,100,000: 65.0); `DR` `observations.chat_vs_panel` (±50 at chr20:200,000: 65.0) | agrees |
| 7 | abstracts: "shows each measurement with the tool output behind it, checks every number in the model's answer against the tool outputs"; sheet: "every number in the model's final answer is checked (intermediate model text is shown unchecked)" | `ui.py:1057`, `ui.py:1125` (a chip beside each measurement and the combined score); `chat.py:495`, `chat.py:767` (`verify_numbers(final_text, …)`); `ui.py` `sendChat` (only the final text is marked) | agrees |
| 12 | "the button beside a layer's number ('tool call #N') or the combined score ('breakpoint_evidence_summary #N')" | `ui.py:1057`, `ui.py:1125` | agrees |

Still open: the new limitation's last sentence and the section 3 row's parenthesis,
which are true of the evidence panel, not of the chat panel. Wording such as "the
evidence panel always counts all four layers, so its figures are right; in the chat
window the model sets the tool's arguments, so the bug can show there" would match the
records. Flags 6, 7 and 12 are resolved, and so is the J.5 ceiling item as worded in
the abstracts, the glossary and the J44 sentence. The tool's ceiling field itself is
unchanged; fixing it is a code change outside this re-check.

*Follow-up, 2026-09-29, meeting sheet v4* (`Rimas_tezes_paaiskinimai_2026-09-29_v4.docx`, 161f7040…eecd9a5d, read only; the abstracts stay at v3): two paragraph diffs against v3 agree (body blocks, and every paragraph of `word/document.xml`). Only `document.xml` and `docProps/core.xml` differ; there are no hidden elements and the images are byte-identical. What changed is the header's version label ("lapo versija v4") and the two clauses J.6 left open, nothing else. The limitation now reads "Įrodymų skydelis visada skaičiuoja visus keturis sluoksnius, todėl jo skaičiai teisingi; pokalbio lange įrankio argumentus parenka modelis, todėl ten ši klaida gali pasirodyti", and the section 3 row's bracket "(įrodymų skydelis visada skaičiuoja visus; pokalbio lange sluoksnius gali parinkti modelis)". Both agree with the records J.6 cited: `ui.py:238` (the evidence panel's layers come from `applicable_layers`, all four at all 12 implants), and `chat.resolve_args`, `ui._chat_exec` and `benchmark/phase10_rerun.py:39–40, 198` for the chat path. No flag remains open in section J; the tool's ceiling field is unchanged, and fixing it is a code change outside the map.

*Fix, 2026-09-30 (commit 25d12bb):* the ceiling-field defect of J.5 is fixed. The ceiling is now computed over the layers the score counts (`results/ceiling_layers_fix_2026-09-30/`: registration, before, after, suites). At IMP01 chr20:200000 the default window now gives 65.0 with discordant pairs and soft clips, 100 with soft clips and split reads, 43.3 with discordant pairs, soft clips and depth, and 57.5 with all four. Before the fix every list reported 57.5. At the 94 public breakends every ceiling field is byte-identical to before, and the prepared-payload proof's model-visible payload is unchanged. The thesis sentences that call the defect uncorrected (K3(b), K7) are now out of date.

## K. The thesis draft of 30 September, 2026-09-30

`Rimas_MSc_Thesis_2026-09-30.docx` (32d5ff13…9abbd121) was copied from the Desktop over
`docs/thesis/Rimas_MSc_Thesis.docx`, replacing be6c9d72…b99709245 (copy verified by sha256;
not committed). Two paragraph diffs against the 27 September draft agree: one compares body
blocks, the other every paragraph of `word/document.xml`. Only `document.xml` and
`docProps/core.xml` differ, in 22 regions, and every one is listed below together with the
date line. There are no tracked changes, comments or text boxes; the two field codes were
already there. The English summary body has 298 words (≤ 300; the title and byline add 24).
`PR` = `results/na12878_provenance_2026-09-27.json`; `R15` = `P10/qwen3.5-4b/WITH__run15.json`.

| # | Claim (30 September) | Record : key | Status |
|---|---|---|---|
| K1 | the background is "the 1000 Genomes Project high-coverage alignment of NA12878, run ERR3239334 [48], produced with BWA-MEM 0.7.15 (mem -Y) [49]" | `PR` `source` (`1000G_2504_high_coverage/data/ERR3239334/NA12878.final.cram`), `nygc_high_coverage_match` (all 12 bwa records 0.7.15-r1140, `mem`, `-Y`, GRCh38_full_analysis_set_plus_decoy_hla) | supported. The unchanged clause "alternate-haplotype-aware mode" is not in `PR` |
| K2 | "the model reaches genomic data only through the tools, so that every figure it reports can be checked against a tool return": Task 1, Methods, Conclusion 1; both summaries | I1 (`ui._chat_exec`, datasets by label); `chat.verify_numbers` on the final answer (`chat.py:495`, `767`) | supported. Both summaries carry only the first half ("…only through the tools") |
| K3 | the ceiling paragraph: the figures assume all four layers counted, as the evidence panel always does; in one extension run a local model counted discordant pairs and soft clips and received 65.0 at the default window with a reported ceiling of 57.5; with soft clips and split reads only, 100 would be attainable; "the defect lies in the protected scoring code and is recorded here rather than corrected" | J.5; `R15` `recorder_calls_unablated` (applicable_layers [discordant_pairs, soft_clipped_reads], window 500, 65.0, `attainable_here` 57.5); run 15 is in the registered extension (`registration_local_extension_2026-09-26.json`); `ui.py:238` (the evidence panel's layers come from `applicable_layers`); 100 is arithmetic, not observed | supported, with two flags: (a) the ceiling fields are added by `server.py` `_with_ceiling` (protected), but derived in `score_tiers.ceiling_from_observed` (not protected), and the tool's score itself is right, so "the protected scoring code" is imprecise; (b) "rather than corrected" goes stale once Phase 20's fix is committed |
| K4 | The Instrument: no build step; the interface and the four evidence layers need no network access; for every image the viewer loads the genome and a gene track from the internet, and one of four panels timed out on 29 September while the track was loading; the gene lookup queries an external service | the page loads nothing from off the machine (no external URL in `ui.py` `PAGE`); the layers read local BAMs with pysam; `DR` `observations.igv_log_during_run`, `steps.6_igv.panel_errors`; `bam_tools.py:2484–2486` (Ensembl REST) | supported. The cloud-model chat option also needs the network, which the text does not say |
| K5 | every tool call passes through one recorder; each layer count and the combined score has a button to its call, and a list gives every call; the ceiling is not a tool return but is derived from the scoring tiers read at startup and the observed fractions | `ui.py` `RECORDER.call` (assess, `_chat_exec`); `ui.py:1057`, `1125`, `showCalls`; `ui.ceiling_for` (tiers from `score_tiers` at import; observed values from the four layer calls) | supported. Two notes: one observed value is a count (`max_clips_at_position`), not a fraction; and the summary tool also returns the same ceiling (`server.py` `_with_ceiling`) |
| K6 | chat: every number in the final answer checked, unmatched numbers shown marked, intermediate text unchecked; in the dry run the model's seven numbers matched while it called moderate evidence strong | J.5 row 7; `DR` `steps.7_chat` | supported |
| K7 | Limitations, two new sentences: the reported ceiling ignores the caller's layers, and the panel always counts every layer; the chat check cannot detect a wrong verbal rating. Recommendation: compute the ceiling over the layers actually counted | as K3, K6 | supported ("always counts every layer": the layers the BAM shows, all four here); the first sentence goes stale after Phase 20 |
| K8 | references [48]–[62], each checked against its source | Crossref: [48] Cell 2022;185(18):3426–3440.e19, [51], [52], [53] (Crossref issued 2011, print 2012;28(4)), [59]: titles, first authors and pages match. arXiv: [49] 1303.3997, [60] 2412.15115, [61] 2407.21783 match. [54] and [55]: the PDFs' first pages read "System Card: Claude Sonnet 5, June 30, 2026" and "System Card: Claude Opus 5, July 24, 2026". [50], [56], [57], [58]: the URLs resolve (HTTP 200). [62]: the title matches secondary listings (Qwen Team, February 2026); the qwen.ai page renders by script and could not be read | supported |
| K9 | citations: each at the first use of its dataset, tool or model | [48] [49] at the background's description; [50] Methods (architecture); [51] [52] Methods (discovery); [53] simulation; [54] [55] the adversarial cases; [56] local serving; [57]–[59] the software list; [60] [61] the frontier sweep; [62] qwen3.5:4b and 9b; [43][44] the case object | all at first use in Methods. Earlier mentions in the abbreviations and the literature review remain uncited: MCP (MARRVEL-MCP), DELLY (a list of callers), VCF, Phenopackets |
| K10 | [43] [44] (existing entries, newly cited) | NCBI: PMID 35705716 is Jacobsen JOB et al., "The GA4GH Phenopacket schema defines a computable representation of clinical data", Nat Biotechnol 2022; "GA4GH Phenopackets: A Practical Introduction" is by **Ladewig MS** et al., Adv Genet 2023, doi:10.1002/ggn2.202200016. PMC11564936 is "A corpus of GA4GH phenopackets: Case-level phenotyping for genomic diagnostics and discovery", HGG Adv 2025 | **differs**: [43] pairs the Ladewig title (with initial "E") with the Jacobsen link; [44] shortens the title and gives no authors, journal or year |

Flags: K10 ([43] title and link point to different papers; [44] incomplete); K3(a) (the
defect lies in `server.py`'s wrapper and `score_tiers.py`, not the scoring itself); K3(b) and
K7 (to be updated after the Phase 20 fix); K2 (the summaries carry half the restated claim);
K4 (the cloud chat needs the network); K9 (earlier uncited mentions in the review). Nothing
else changed.

## L. The thesis after the ceiling fix (30 September, v2) and the meeting sheet v5, 2026-09-30

`Rimas_MSc_Thesis_2026-09-30_v2.docx` (1df2f62f…d0946fc2) was copied over `docs/thesis/`
(verified by sha256, not committed), replacing 32d5ff13…9abbd121. The sheet
`Rimas_tezes_paaiskinimai_2026-09-30_v5.docx` (a632f2ec…2983168b) replaces v4. Both were
read only from the Desktop. `FX` = `results/ceiling_layers_fix_2026-09-30/`.

| # | Change | Record : key | Status |
|---|---|---|---|
| L1 | K10: [43] → Ladewig MS, Jacobsen JOB, Wagner AH, et al., Advanced Genetics 2023;4(1):2200016, doi:10.1002/ggn2.202200016; [44] → Danis D, Bamshad MJ, Bridges Y, et al., HGG Advances 2025;6(1):100371, doi:10.1016/j.xhgg.2024.100371 (now [21] and [22]) | Crossref for both DOIs: title, first three authors, volume, issue and article number match (Ladewig issued online 2022, in the 2023 issue) | resolved |
| L2 | K9: citations at the first mentions in the literature review: VCF [13] (Exomiser paragraph), MCP [20] and Phenopackets [21][22] (MARRVEL-MCP paragraph), DELLY [39] (list of callers) | the new numbers map back to the VCF, MCP, Phenopackets and DELLY entries | resolved |
| L3 | renumbering by order of first citation | 60 entries matched between the drafts by identical text; they give exactly the stated mapping (old→new), and the two rewritten entries are 43→21 and 44→22 as stated. All 165 earlier body citations point at the same work after mapping; the only additions are L2's four. First citations run 1, 2, …, 62 in order; no entry is uncited; there are no citations after the list | agrees |
| L4 | K2: both summaries carry the whole claim ("…so that every figure it reports can be checked against a tool return" / "…todėl kiekvieną jo pateiktą skaičių galima patikrinti pagal įrankio atsakymą"); EN rewording: "four signals", "A separate caller performs discovery", "while no tool returned the number it required; with the number added, …, often never calling the tool that returned it" | K2's records | resolved. The EN summary is 300 words (≤ 300), the LT 268 |
| L5 | K4: "the image tools, the gene lookup and the chat panel's cloud models" need the network | K4's records; the cloud chat goes to the Anthropic API (`chat.run_turn_api`) | resolved |
| L6 | K3 and the fix: until 30 September the reported ceiling did not follow a restricted list; the score was right, and the ceiling fields were derived in a separate module and attached by the wrapper; 72 summary calls, 5 restricted, all from qwen3.5:4b, 1 above its ceiling; a registered fix; at IMP01 chr20:200000 57.5 / 65.0 / 100 / 43.3; the test failed first; the 94 breakends identical; the payload unchanged; 26 suites, 759 assertions. The interface's ceiling is "over the layers the score counts" | `FX/before.json` `phase10.total` (72, 5, 1; all restricted calls from qwen3.5:4b); `FX/registration.json` (commit 95e1927, before the fix 25d12bb); `FX/after.json` `imp01_agreement`, `identity.all_identical`, `proof`; `FX/test_before_fix.log`, `test_head_copies_prefix.log`; `FX/suites.json` `summary`; `ui.py` `ceiling_for(obs, counted_layers(sm))` | agrees |
| L7 | K7: the limitation sentence on the ceiling and the recommendation clause removed | — | resolved (both were stale after the fix) |
| L8 | sheet v5: the limitation now reads "found and fixed": until 30 September the ceiling ignored the model's layer list; one case (65.0 against 57.5) in 72 summary calls; now computed over the score's layers (IMP01: 57.5 / 65.0 / 100), with all four layers every ceiling field unchanged. Section 3 bracket: "from 09-30 the ceiling is computed over the same layers as the score". Glossary: the ceiling is per position, window and counted layers | as L6 (`FX/before.json`, `after.json`); `test_ceiling_layers.py` negative control (160 grid cases) and `after.json` `identity` (94 breakends) for "unchanged" | agrees |
| L9 | nothing else changed | thesis: two paragraph diffs against the morning draft, with citation numbers masked (body and summaries), show 9 changed paragraphs, all under L2 and L4–L7; apart from those, only citation numbers and the reference list's order changed. Only `document.xml` and `docProps/core.xml` differ; no tracked changes, comments or text boxes. Sheet: only the version label ("v5") and the three L8 items changed | agrees |

No flag remains open from section K, and none is raised by these files.
