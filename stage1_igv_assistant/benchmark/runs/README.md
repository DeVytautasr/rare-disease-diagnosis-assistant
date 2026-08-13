# Benchmark run logs

Every run of the benchmark harnesses, grouped by which tool and scorer fixes
were in place when it executed. Each log also carries that information in its
own `stage` field, so a file remains self-describing if it is ever moved or
read in isolation — the directory name is a convenience, not the record.

## Stages

| Directory | Runs | `stage` |
|---|---:|---|
| `v1_prefix/` | 27 | *pre-fix: predominantly bug and 30-char regex active* |
| `v2_post_fixA/` | 18 | *post fix A: predominantly corrected* |
| `v3_post_fixACD/` | 12 | *post fixes A, C, D: paths server-assigned, sentence-scoped scoring* |
| `blind_arm/` | 3 | *instruction-blind arm, self-reported trace* |
| `visual_sessions/` | 9 | per-file — these span four stages; see below |
| `unsorted/` | — | default output of a fresh harness run, before a stage is assigned |

The fixes referenced:

- **Fix A** — `summarize_breakpoint_evidence` asserted a "predominantly"
  partner chromosome that did not exist (it took the first-inserted dict key
  and never checked for dominance), so a single discordant read produced
  *"mates mapping predominantly to chr12"*.
- **Fix C** — `igv_screenshot` and `evidence_panel` stopped accepting an
  output path and stopped returning one; the server assigns the location and
  returns an opaque handle.
- **Fix D** — `correct_verdict`'s negation detector moved from a fixed
  30-character window to sentence-scoped, nearest-cue classification.

(Fix B — harness-side redaction of paths from tool results — was superseded
by Fix C and has no directory of its own. It survives only in
`visual_sessions/`, where two logs were recorded under it.)

## What each stage supports

- **`v1_prefix/`** — the originally published numbers, and the evidence that
  the tool bug reached committed results. Both of the corrections issued in
  `results/BENCHMARK_LOCAL_MODELS.md` are demonstrated from these files.
- **`v2_post_fixA/`** — isolates the effect of the tool fix. Because the
  input data is identical to `v1_prefix/` for the adversarial case
  (`mate_chromosomes = {'chr12': 1}` in all six runs), the pair forms a
  controlled comparison in which only the tool's sentence changed.
- **`v3_post_fixACD/`** — the current numbers, and the basis for the FIX C
  verification (grepping a full message history for any path-shaped string).
- **`blind_arm/`** — the Claude Code arm. Its tool-call trace is
  self-reported by the agent rather than captured from the wire, so it is
  read as a plausibility-checked transcript, not a ground-truth log.
- **`visual_sessions/`** — the one-off visual-tool sessions. These span four
  stages (original, advisory-note, harness-redaction, server-assigned paths),
  which is why they are grouped by kind rather than by stage; each file's
  `stage` field records its own. Reading them in filename order shows the
  image-path constraint tightening and the model behaviour changing with it.

## Why `v1_prefix/` is kept, given it contains results now known to be wrong

Because "known to be wrong" is a finding, and deleting the evidence would
leave the correction unverifiable.

Two published claims turned out to be measurement artifacts: qwen2.5:7b was
accused of fabricating a "predominantly chr4" claim it had in fact quoted
verbatim from its own tool, and claude-sonnet-5 was scored 2/3 on the
adversarial case by a regex too narrow to see a rejection phrased across 38
characters. Both corrections are checkable only against the runs that
produced them — the exact tool output the model received, and the exact
report text the scorer misread.

Keeping the superseded logs also keeps the controlled comparison intact: the
strongest statement available about the tool fix is that identical input data
produced different model behaviour once the sentence changed, and that
requires both halves. Discarding `v1_prefix/` would reduce the correction to
an assertion.

Superseded numbers are marked as superseded — in the stage field, in this
file, and in the correction notices at the top of both benchmark documents.
They are not presented as current anywhere.
