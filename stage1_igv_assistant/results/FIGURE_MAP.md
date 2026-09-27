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
