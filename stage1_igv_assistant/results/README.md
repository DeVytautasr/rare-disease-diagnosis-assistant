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
| 3 | `EVIDENCE_PANEL_VALIDATION.md` | 71 | **Short, but the authority for the visual-interpretation trap cited in three other documents** — the discordant-pairs panel looks visually busy while the underlying count is 1 of 1,708, because IGV colours anomalous insert size and inter-chromosomal mates alike. Its line count badly understates its role; most of its content is in long table rows. Read it before interpreting any panel image. |
| 4 | `REAL_DATA_VALIDATION.md` | 197 | Tool behaviour against the GIAB HG002 BAM at the confirmed CMRG deletion. The ground-truth reference the benchmark cases are built on. |
| 5 | `LLM_SESSION_3_BLIND.md` | 176 | The instruction-blind methodology the benchmark's blind arm follows. |
| 6 | `LLM_SESSION_4_VISUAL_claude-sonnet-5.md` | 152 | Visual-tool session. Documents the stale-file false-success bug and the image-description finding. |
| 7 | `LLM_SESSION_4_VISUAL_qwen2.5-7b.md` | 139 | Visual-tool session. Opens with the retraction of the "fabricated predominantly claim" accusation. |

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

## One caution that applies to every document here

Three of the five scoring criteria are regex heuristics over free text.
`correct_verdict` has required three separate corrections, each found by a
human reading full reports and none by the metric itself. Where a document
reports an adversarial verdict, that verdict has been read individually;
where you extend this work, do the same. The score flags a run for
inspection — it does not settle it.
