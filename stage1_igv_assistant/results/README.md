# Results index

Prose only. Run logs live under `benchmark/runs/` (see its README for the
stage layout); the demo case object is at `data/demo_case.json`.

## Reading order

Start at 1 and 2 — they carry correction notices that supersede parts of the
older documents, so reading them first prevents taking a retracted finding at
face value.

### Current

| # | Document | Lines | What it is |
|---|---|---:|---|
| 1 | `BENCHMARK_LOCAL_MODELS.md` | 411 | qwen2.5:7b vs llama3.1:8b across three cases. Opens with the correction notice covering both retracted findings, and states why `correct_verdict` is a screening aid rather than a verdict. |
| 2 | `BENCHMARK_CLAUDE_BASELINE.md` | 362 | claude-sonnet-5 on the same cases and code. Carries the second correction notice, the cost accounting, and the FIX C architectural result. |
| 3 | `REAL_PATIENT_DATA_VALIDATION.md` | 813 | Validation of the Stage 1 tools against the two real patient BAMs: 14 findings and 3 crashes, each with a status column recording what was done about it, plus a post-fix re-measurement addendum and two later findings (15, 16) from an LLM session, with the control grid re-run at two window widths. Not to be confused with `GIAB_PUBLIC_DATA_VALIDATION.md`, which covers public data and reaches different conclusions about threshold calibration.  *Correction, 2026-09-26: its statement that `bwa mem -M` produces zero `0x800` reads holds for the chromosomes it scanned (chr21 counts zero again, directly and under every reconstruction of the scan) but not for the files: 3,694,802 and 3,726,497 records carry `0x800`, all of them on alt, decoy and HLA contigs (`supplementary_flags_2026-09-26.json`).* |
| 4 | `LLM_SESSION_5_PATIENT_DATA_qwen.md` | 503 | First LLM session on the post-fix tools, and the first on patient data. qwen2.5:7b, MCP tools only, 15 runs over 5 positions. Records that the honest observation sentences work as intended, and finds two new defects — both since fixed, one of them older than the fix it was blamed on. |
| 5 | `EVIDENCE_PANEL_VALIDATION.md` | 71 | **Short, but the authority for the visual-interpretation trap cited in three other documents** — the discordant-pairs panel looks visually busy while the underlying count is 1 of 1,708, because IGV colours anomalous insert size and inter-chromosomal mates alike. Its line count badly understates its role; most of its content is in long table rows. Read it before interpreting any panel image. |
| 6 | `GIAB_PUBLIC_DATA_VALIDATION.md` | 205 | Tool behaviour against the GIAB HG002 BAM at the confirmed CMRG deletion. The ground-truth reference the benchmark cases are built on. |
| 7 | `LLM_SESSION_3_BLIND.md` | 176 | The instruction-blind methodology the benchmark's blind arm follows. |
| 8 | `LLM_SESSION_4_VISUAL_claude-sonnet-5.md` | 152 | Visual-tool session. Documents the stale-file false-success bug and the image-description finding. |
| 9 | `LLM_SESSION_4_VISUAL_qwen2.5-7b.md` | 139 | Visual-tool session. Opens with the retraction of the "fabricated predominantly claim" accusation. |

### History — retained, superseded, not current

| Document | Lines | Why kept |
|---|---:|---|
| `AUDIT_2026_08.md` | 442 | The systematic audit that produced the scoring-threshold work. Largest document here; its findings are folded into the current tool behaviour, so it reads as the record of how the pipeline got here rather than as current guidance. |
| `LLM_SESSION_1.md` | 85 | First LLM session. Predates the depth-ratio, localisation and soft-clip fixes. |
| `LLM_SESSION_2_WITH_VISUAL.md` | 200 | First session with screenshots. Same caveat; its `session2_chr1_deletion*.png` images predate the `evidence_panel` tool entirely. |
| `RESULTS_HCC1143.md` | 146 | HCC1143 negative-control validation, on an older scoring pipeline. |
| `DEMO_END_TO_END.md` | 85 | End-to-end demo on the synthetic translocation fixture. References `data/demo_case.json`. |

**Why the history set is not deleted:** these predate the depth-ratio,
localisation and soft-clip scoring fixes, so their numbers describe an older
system. That makes them unusable as current results and useful as the record
of what changed — `BENCHMARK_CLAUDE_BASELINE.md` explains at length why the
Claude baseline was regenerated rather than scored from these transcripts.

## Figure map — 2026-09-25

`FIGURE_MAP.md` sets every numeric claim of the conference abstract, the Phase 5,
6 and 10 figures that survive only as text, and the patient figures against the
records that now back them, with a status for each (reproduced exactly /
replaced by the new experiment / record lost, not regenerated / no longer
supported). The thesis file itself was not available on this machine and is not
yet mapped.

*2026-09-26:* section D adds the Phase 12 figures and the status of each thesis figure
marked †. The thesis file is now in `docs/thesis/` (untracked, read only); it is still
not mapped, as instructed.

## Not written up here — phases 5 to 11

Everything indexed above predates a run of later work. That work produced
measurements but not prose documents, so this index would mislead anyone who
read it as complete. What exists, and where:

| Work | What was measured | Where the record is |
|---|---|---|
| Candidate-set bridge | Deduplication keyed on position alone collapsed reciprocal junctions; keying on orientation took the public set from 840 to 894 junctions and survivors 24 → 27 (rerun on 2026-09-25: the same 894 and 27, `analysis_2026-09-25/figures.json`). `compare_candidate_sets` was asymmetric by construction — 54 of 56 reported "unmatched" junctions were mis-reported. *Addition, 2026-09-26: the pre-fix code, never committed, was reconstructed in a scratch module and reproduces 840/24 and 894/2 against 840/56. On the patient rerun BCFs it gives the pre-fix BND after dedup 8,716/9,144, PASS 838/862 and final survivors 15/16. Background merge groups: 33 are exactly DEL+DUP (34 hold at least one of each); 12, not 13, hold INV junctions of both orientations. The 1,738 two-direction difference appears only with the dedup already fixed; before both fixes it is 53 and 89.* | `tests/test_vcf_tools.py`, and the commit message for the bridge; the reconstruction: `prefix_bridge_2026-09-26.json` |
| Controlled positive test | 12 heterozygous balanced translocations implanted into real NA12878 reads at known positions. Rebuilt 2026-09-25 as a new experiment: **16 of 24 junctions** recovered (8 of 12 translocations): 8/8 clean, 8/8 repeat-adjacent, 0/8 low-mappability; every loss at discovery (no delly record), none at any filter; each detected implant adds exactly its 2 junctions to the background's 27 survivors. *Correction, 2026-09-25: the first run found 14 of 24 (6/8 repeat-adjacent); its records were lost. The difference is IMP06 — missed then, detected now — and realigning it without the `.alt` file loses it again, so ALT-aware alignment accounts for it* | `synthetic_control_2026-09/analysis_2026-09-25/` (figures.json, adjudication/, ladder/, noalt/), with the implants and delly in `synthetic_control_2026-09/`. The first run's `~/public_data/sim/implants.json` was lost in the reinstall |
| Local front end | Filter chain, four layers at both breakends, hand-entered coordinates, two-sample comparison, call log | `stage1_igv_assistant/ui.py` |
| Patient delly rerun, 2026-09-25 | Both patient BAMs re-verified against the Phase 0 record after the reinstall and rerun through delly v2.6.0, one sample per process: every Phase 3 funnel count reproduced exactly (total, BND, BND after dedup, PASS BND, before and after recurrence). Counts only | `results/patient_rerun_2026-09.json` |
| Patient BAM properties, 2026-09-26 | The bwa `@PG` records carry no `VN` field (the string 0.7.17-r1188 is in their command lines); read length 150; nominal depth 30.68 / 32.91-fold; duplicates 12.963 % / 12.968 % of primary records (12.73 % / 12.76 % of all records); median insert size 314 / 312 bp; reference concordance 951 of 952 and 954 of 954 positions meeting the floors. Definitions beside every figure; totals only | `results/patient_properties_2026-09-26.json` |
| Supplementary flags, 2026-09-26 | Every `0x800` record in both BAMs (3,694,802 / 3,726,497, equal to flagstat) lies on an alt, decoy or HLA contig; none on chr1-22, X or Y; all carry SA, none is also `0x100`. The `@PG` chain is bwa per lane, then bamcat and bamsormadup | `results/supplementary_flags_2026-09-26.json` |
| Test census, 2026-09-26 | All 24 suites under three network conditions (no suite has a switch for them, so they were imposed from outside): 660 / 662 / 685 assertions held without network / with Ensembl only / with the network; none failed; `test_bam_tools` incomplete in all three because IGV is not installed on this machine | `results/test_census_2026-09-26.json` |
| Local model comparison | 3 models × 6 cases × 5 runs. Malformed tool arguments: `qwen3.5:4b` 51/194, `qwen2.5:7b` 19/142, `qwen3.5:9b` 29/170. Completion 25/30, 30/30, 11/30 | `results/phase8_final_record.json` — a restored copy of the original; see below |
| Claude API control | `claude-sonnet-5` 30 runs, `claude-opus-5` 15 runs, same six cases. 45/45 answered, 0 refusals, 0 malformed arguments in 313 calls. Spend $3.86 *Correction, 2026-09-26: that record is lost. The sweep was rerun on reconstructed cases (the original prompts were never committed) with the same models and run counts: 45/45 answered, 0 refusals, 0 malformed arguments in 341 calls, no tool call written as prose, no coordinate drift, the false premise rejected 5/5 (sonnet) and 3/3 (opus); 34 numbers matched no tool return, every one traced to a tool description, the system prompt or a rounded returned position; spend $3.59* | `~/public_data/sim/phase9_final_record.json` (lost); the rerun: `benchmark/runs/phase9_rerun_2026-09-26/` (every trace, `measures.json`, `reading_scores.json`) |
| The ceiling result | Before: 1 of 20 runs across four models stated that a balanced translocation cannot reach "strong". After adding nine derived fields to one tool return, with no model change: 19 of 20, and 17 of 20 gave the full argument *Correction, 2026-09-25: both counts come from records lost in the reinstall. Rerun with the fields as the only difference (proven), local models only (API key rejected): without the fields 0 of 10 runs state the ceiling, with them 2 of 10 (qwen3.5:4b 2 of 5, qwen2.5:7b 0 of 5); the tool's attainable ceiling is 57.5 against 70 at all 32 detected breakends of the rebuild* *Correction, 2026-09-26: scored blind, with the API arm run: over the same twenty runs 0 of 20 state the ceiling without the fields and 12 of 20 with them, all 12 with the full argument (sonnet 5/5, opus 5/5, qwen3.5:4b 2/5, qwen2.5:7b 0/5); at 15 runs per cell qwen3.5:4b 5/15 and 4/15 (p = 0.042, 0.100), qwen2.5:7b 3/15 and 0/15 (p = 0.224, 1.0); with the fields sonnet also called the evidence strong in 5 of 5 runs* | `~/public_data/sim/phase10_final_record.json`, `tests/test_ceiling_echo.py`; the rerun: `benchmark/runs/phase10_rerun_2026-09-25/` (proof, every trace, `scores.json`, `tally.json`); the blind scores: `blind/blind_scores.json` |
| Installability | Installed into a clean environment on Python 3.14.4 with no inherited packages; three iterations. A 1.6 MB demo bundle regenerates from public data via `make_demo_bundle.py` | `install.sh`, `docs/DIEGIMAS.md` |

**Two things to be aware of before citing any of this.** The JSON records under
`~/public_data/sim/` were outside the repository and not committed, and the
2026-09-23 reinstall destroyed them: that directory no longer exists, and the
implant records and the Phase 9 and 10 records are gone. Only the Phase 8
record survived, as a copy, and is now committed (below). And the measurements
above have not been through the same read-and-correct pass that produced the
correction notices on documents 1 and 2 — the caution below applies to them
with more force, not less.

### `phase8_final_record.json` — restored copy

A restored copy of the original Phase 8 final record (its `generated` field:
2026-09-09T16:45:31), put back into the repository on 2026-09-25 because it is
the only Phase 5–11 primary record that survived the reinstall. It was not
regenerated and has not been edited; sha256
`8dd84a65f61922c53006e590b0d4f86bf012f68284c6ea976a61c69e4144cf61`. It
covers the local model comparison:

- the hardware it ran on (an 8 GB laptop GPU, 16 CPU cores) and the chosen
  configuration — context 32,768, `max_iters` 40, untrimmed tool descriptions,
  `think='low'` where the model supports it — with the reason for each;
- how much of each model stays on the GPU at that context, and the token cost of
  the tool schemas, including a trim that was measured and reverted;
- the iteration-cap and thinking-level experiments;
- per-model results for `qwen3.5:4b`, `qwen2.5:7b` and `qwen3.5:9b` over the
  six cases (answered runs, empty finals, tool calls, schema-invalid calls,
  speed, transport errors) — the figures in the table above;
- the Phase 8 and 8b figures it supersedes, each with the reason, and four open
  questions.

## One caution that applies to every document here

Three of the five scoring criteria are regex heuristics over free text.
`correct_verdict` has required three separate corrections, each found by a
human reading full reports and none by the metric itself. Where a document
reports an adversarial verdict, that verdict has been read individually;
where you extend this work, do the same. The score flags a run for
inspection — it does not settle it.
