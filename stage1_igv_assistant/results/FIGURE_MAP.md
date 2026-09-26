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
