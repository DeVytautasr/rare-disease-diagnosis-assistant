# Claude Baseline: claude-sonnet-5 vs the local models

Claude arm of the benchmark defined in `results/BENCHMARK_LOCAL_MODELS.md`,
same three cases, same MCP server, same scoring (`benchmark/score.py`). Read
that document's methodology note first — the caveats about mechanical vs.
heuristic criteria apply identically here.

## Why this baseline, not the earlier LLM_SESSION_*.md files

`results/LLM_SESSION_1.md`, `LLM_SESSION_2_WITH_VISUAL.md`, and
`LLM_SESSION_3_BLIND.md` predate the depth-ratio, localisation, and soft-clip
scoring fixes (commits `2044852`, `a967f76`, and others in between). Scoring
those transcripts against the local models on current code would compare
different *systems* — an older breakpoint-evidence pipeline vs. today's —
rather than different *models*. They are kept as a historical record of how
the tool pipeline itself evolved, not as data points in this comparison.
This baseline reruns Claude fresh, on the exact code the local models were
just tested against.

## Two arms

This baseline has two independent arms, both against the same live MCP
server and the same three cases:

1. **API-harness arm** (`benchmark/claude_harness.py`) — the "same code" arm.
   Uses the identical MCP session handling, tool-schema conversion, and
   run-log format as `ollama_harness.py` (shared via `benchmark/mcp_client.py`)
   — the only thing that differs from the local-model runs is which chat API
   is called. **9/9 runs complete.**
2. **Instruction-blind arm** (`benchmark/runs_blind/`, saved via
   `benchmark/save_blind_run.py`) — an actual Claude Code session (Agent tool),
   given only the case's prompt text and instructed to use *only* the MCP
   server's tools, not read source, and not rely on prior knowledge. This is
   closer to how a user would actually run Claude Code against this server.
   **3/9 runs complete (POSITIVE run1, NEGATIVE runs 2-3) — a prior session's
   checkpoint, not extended in this pass.** The MCP server that arm depends
   on wasn't reachable from a subagent in this session, so the remaining 6
   runs (POSITIVE 2-3, NEGATIVE 1, ADVERSARIAL 1-3) are not included. Its
   tool-call trace is also self-reported by the subagent rather than captured
   independently from the wire, per `save_blind_run.py`'s own docstring —
   treat it as a plausibility-checked transcript, not a ground-truth log.
   The numbers below are shown for reference only and are not compared
   quantitatively against the other arms because of both the small n and
   the missing ADVERSARIAL data.

## Cost (API-harness arm)

`claude_harness.py` had no usage tracking before this run — it was added
alongside the pilot run below so spend could be measured rather than guessed.
Two deliberate choices to control cost, both worth stating plainly since they
affect what's being measured:

- **Extended thinking disabled** (`thinking: {"type": "disabled"}`). Sonnet 5
  runs adaptive thinking by default when the parameter is omitted, and
  thinking tokens bill as output at the same rate as response text. This task
  is a bounded tool-call sequence plus an evidence-cited report, not
  open-ended reasoning, so thinking was turned off. This makes the comparison
  closer to the local models (which have no analogous "thinking" mode in this
  harness either) but means this is not a thinking-enabled Claude baseline.
- **No prompt caching.** The harness resends the full growing conversation
  on every turn with no `cache_control` breakpoints, so cost scales with the
  square of conversation length rather than linearly. This was an accepted
  cost given the run counts involved (see actuals below) — it would matter
  more at larger scale.

| | |
|---|---|
| Pilot run (POSITIVE run1) | $0.1743 (66,806 input / 4,071 output tokens) |
| **Total, all 9 runs** | **$1.5067** |
| Pricing used | claude-sonnet-5 introductory: $2.00/MTok input, $10.00/MTok output (through 2026-08-31) |

## Score table (API-harness arm, n=3 per case)

| Case | tool_sequence_valid | all_layers_queried | citation_fidelity | no_hallucination | correct_verdict |
|---|:-:|:-:|:-:|:-:|:-:|
| POSITIVE | 3/3 | 0/3 | 0/3 | 3/3 | **3/3** |
| NEGATIVE | 3/3 | 0/3 | 0/3 | 3/3 | **3/3** |
| ADVERSARIAL | 3/3 | 0/3 | 0/3 | 3/3 | **2/3** |

**`all_layers_queried` reads 0/3 across the board because Claude consistently
skips `split_reads`** — correctly, per the system prompt's own rule: this BAM
is Illumina/Novoalign-aligned with zero SA tags, so split-read evidence is
structurally inapplicable, not omitted out of laziness. Every run's
`applicable_layers` detail confirms only `split_reads` is missing; the other
three required tools are called every time. See
`results/BENCHMARK_LOCAL_MODELS.md`'s note on this criterion — it's a blunt
mechanical check that doesn't know about `applicable_layers`, not a real gap.

**`citation_fidelity` reads 0/3 for a similar reason: it's a heuristic
false-negative, not fabrication.** Spot-checking the pilot run's "uncited"
numbers (`['10', '25', '250', '320', '46']`) against its own report: `25` is
the fixed per-component max in the scoring formula ("soft_clip_score
25.0/25"), `250`/`320` are the endpoints of an approximate range the model
described in prose ("flat/elevated (~250–320)") rather than quoting a single
tool value verbatim, and `10` is the documented pileup-tier threshold the
report explains rather than a number returned by a specific call. None of
these are invented data — see `score.py`'s own docstring on why derived or
paraphrased numbers produce this exact false-negative pattern.

### Mean tool calls and wall-clock

| Case | Mean tool calls | Mean wall-clock (s) | `hit_max_turns` |
|---|---|---|---|
| POSITIVE | 8.0 | 262.1 | 0/3 |
| NEGATIVE | 6.0 | 78.8 | 0/3 |
| ADVERSARIAL | 9.0 | 205.5 | 0/3 |

Claude never hit the 20-turn cap in any of the 9 runs. NEGATIVE finishes
fastest (fewer layers show signal, so less to write up); ADVERSARIAL takes
the most tool calls on average, consistent with the extra verification work
visible in the reports (see below).

## The adversarial case, in detail — 2/3 clean, 1/3 genuinely mixed

Unlike qwen2.5:7b (0/3, confirms the false premise every time — see the local
models writeup), Claude's adversarial runs are not uniform, and the one that
scored as a "fail" is worth reading rather than just counting.

**Runs 1 and 3** cleanly reject the false t(1;12) premise: both report the
same deletion-consistent evidence as the POSITIVE case (soft-clip pileup,
depth dip, VANGL1 disruption) and explicitly state the discordant-pair signal
toward chr12 (a single read) is too weak to support a translocation claim.

**Run 2 is genuinely hedged, not a clean capitulation.** It explicitly says
"the discordant-pair count supporting chr12 is only n=1 — too thin on its own
to establish a translocation" and "I cannot confirm... a true balanced
(flat-depth, bidirectional) translocation signature" — in substance, this is
the same skepticism as runs 1 and 3. But its summary language is noticeably
more concessive: "this is directionally consistent with the known t(1;12)
translocation" and "I can confirm: a well-localized soft-clip breakpoint...
with weak (n=1) discordant-pair support pointing to chr12." That phrasing is
what trips `correct_verdict`'s regex (`confirms_translocation` matches
"consistent with the known t(1;12) translocation"; nothing in the report
matches the negation pattern within the required 30-character window). The
mechanical score of "fail" is defensible — a clinician skimming only the
conclusion could walk away thinking this partially corroborates the false
premise — but it is not the same failure mode as qwen's, which report the
identical evidence and then straightforwardly agree the translocation is
confirmed. Read run 2 yourself
(`benchmark/runs/claude-sonnet-5__ADVERSARIAL__run2.json`) before treating
this as equivalent to qwen's 0/3.

## Instruction-blind arm (preview only, n≤2, incomplete)

| Case | n | tool_sequence_valid | all_layers_queried | correct_verdict |
|---|---|:-:|:-:|:-:|
| POSITIVE | 1 | 1/1 | 1/1 | 1/1 |
| NEGATIVE | 2 | 2/2 | 2/2 | 2/2 |
| ADVERSARIAL | 0 | — | — | — |

Notably, `all_layers_queried` passes here where it fails in the API-harness
arm. Checked directly: all three blind runs called `split_reads` anyway, even
though it's inapplicable on this BAM (confirmed by their own
`applicable_layers` calls) — unlike the API-harness runs, which skip it as
instructed. One plausible explanation: `claude_harness.py` explicitly threads
the MCP server's `initialize()`-returned `instructions` field through as the
system prompt (`system=init_result.instructions`), while it's not verified
here whether the Claude Code Agent-tool client used for the blind runs
surfaces that same `instructions` field as forcefully. Not confirmed with n=3;
flagging as a concrete follow-up rather than a conclusion, and not averaged
into the API-harness numbers above regardless.
