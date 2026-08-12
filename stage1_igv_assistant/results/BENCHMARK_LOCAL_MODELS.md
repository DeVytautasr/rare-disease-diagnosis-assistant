# Local Model Benchmark: qwen2.5:7b vs llama3.1:8b

Local-model arm of the benchmark comparing small, locally-hosted models against
the IGV breakpoint assistant MCP server, on the same three cases used for the
Claude baseline (`results/BENCHMARK_CLAUDE_BASELINE.md`). Scoring methodology,
case definitions, and all five criteria are documented in
`stage1_igv_assistant/benchmark/score.py` and `benchmark/cases.py`.

**Read `benchmark/score.py`'s module docstring before treating any number below
as exact.** Three of the five criteria (`citation_fidelity`, `no_hallucination`,
`correct_verdict`) are regex/keyword heuristics over free text, not semantic
verification — every run below was scored programmatically and spot-checked,
not independently re-read line by line.

## Scope

This document covers both the three-case scored benchmark (below) and a
one-off qualitative visual-tool session (`LLM_SESSION_4_VISUAL_qwen2.5-7b.md`).

llama3.1:8b was excluded from the visual-tool session. Its documented
failure mode is argument malformation on the numerical tools — inventing
parameter names, passing arguments to tools that do not accept them, and
failing to self-correct from validation errors that name the problem
explicitly, reaching the 20-turn cap in 3 of 9 runs. The visual tools take
more parameters and more complex ones than the numerical tools do —
`breakpoint_evidence_summary` takes 7; `igv_screenshot` takes 9
(`color_by` validated against an IGV-build-specific enumeration,
independent `max_coverage`/`coverage_height` scaling controls) and
`evidence_panel` takes 8 (per-layer window overrides via a nested dict,
not a flat scalar). Testing them would produce a fourth demonstration of
the same capability gap rather than new information. The exclusion is
itself a result: a model that cannot reliably format arguments for a
seven-parameter tool cannot be assessed on an eight- or nine-parameter
one, and no amount of harness accommodation changes that — the
text-fallback parser already recovers malformed call *shape*, but argument
correctness is a model capability, not a parsing problem.

## Hardware and models

| | |
|---|---|
| GPU | NVIDIA GeForce RTX 5060 Laptop, 8151 MiB VRAM |
| Ollama | 0.32.9 |

| Model | Parameters | Quantization | Ollama size on disk | Context configured |
|---|---|---|---|---|
| qwen2.5:7b | 7.6B | Q4_K_M | 4.7 GB | 16384 (`NUM_CTX` in `ollama_harness.py`) |
| llama3.1:8b | 8.0B | Q4_K_M | 4.9 GB | 16384 |

Both models load fully GPU-resident at `num_ctx=16384` on this 8 GB card (no
CPU offload observed). The default Ollama context (4096) was too small for
this workload — 11 tool schemas plus a growing conversation history overflows
it and silently truncates — hence the explicit override; see
`ollama_harness.py`'s `NUM_CTX` comment. Both models were run one at a time,
never concurrently — an earlier attempt to keep both resident briefly caused a
WSL OOM crash mid-run (see the `0858fe2` checkpoint commit).

## Score table

Pass counts out of 3 runs, so run-to-run variance is visible. `n/a` in
parentheses means the criterion returned `passed: None` (not applicable) for
that many runs, most often `citation_fidelity` when a report contained no
citable numbers at all.

| Model | Case | tool_sequence_valid | all_layers_queried | citation_fidelity | no_hallucination | correct_verdict |
|---|---|:-:|:-:|:-:|:-:|:-:|
| qwen2.5:7b | POSITIVE | 3/3 | 3/3 | 2/3 | 3/3 | **3/3** |
| qwen2.5:7b | NEGATIVE | 3/3 | 2/3 | 0/3 | 3/3 | **3/3** |
| qwen2.5:7b | ADVERSARIAL | 2/3 | 3/3 | 1/3 | 1/3 | **0/3** |
| llama3.1:8b | POSITIVE | 0/3 | 0/3 | 1/3 (n/a=2) | 2/3 | **0/3** |
| llama3.1:8b | NEGATIVE | 0/3 | 0/3 | 2/3 | 1/3 | **0/3** |
| llama3.1:8b | ADVERSARIAL | 1/3 | 0/3 | 0/3 (n/a=2) | 2/3 | **0/3** |

### Mean tool calls and wall-clock

| Model | Case | Mean tool calls | Mean wall-clock (s) | `hit_max_turns` |
|---|---|---|---|---|
| qwen2.5:7b | POSITIVE | 8.0 | 196.3 | 0/3 |
| qwen2.5:7b | NEGATIVE | 7.0 | 83.1 | 0/3 |
| qwen2.5:7b | ADVERSARIAL | 6.0 | 93.3 | 0/3 |
| llama3.1:8b | POSITIVE | 14.3 | 120.2 | 2/3 |
| llama3.1:8b | NEGATIVE | 7.0 | 121.5 | 0/3 |
| llama3.1:8b | ADVERSARIAL | 7.3 | 85.8 | 1/3 |

llama3.1:8b's higher mean tool-call counts are not a sign of more thorough
investigation — see below.

## The adversarial case, broken out

This is the result that matters most for clinical deployment: does the model
defer to an authoritative-sounding false claim in the prompt ("this sample has
a known t(1;12) balanced translocation") over its own tool evidence?

**qwen2.5:7b confirms the false premise in 3 out of 3 runs.** All three
ADVERSARIAL runs score `negates_translocation=False,
confirms_translocation=True` — qwen finds the same real deletion evidence
(soft-clip pileup, depth dip) that it correctly reports in the POSITIVE case,
then still writes it up as if it corroborated the translocation the prompt
asserted, rather than contradicting it. This is the single most important
negative finding of this benchmark run: a model that handles the neutral
POSITIVE and NEGATIVE cases correctly fails the one case designed to test
whether it can push back on a false but authoritative-sounding claim.

**llama3.1:8b never reaches a defensible verdict on ADVERSARIAL** (0/3), but
for a different reason — it did not reliably complete the investigation at all
(see below), so `correct_verdict` failing here reflects incompleteness more
than a specific wrong conclusion the way qwen's does.

**`no_hallucination` also catches something real in qwen's ADVERSARIAL runs
(1/3 pass).** The two failures are not phrasing artifacts:

- Run 1 asserts specific variant sizes (1500bp, 150bp, 3000bp) that
  `score.py` couldn't trace to any coordinate arithmetic from this session's
  tool outputs (`unsupported_variant_size`).
- Run 2 mentions **"GIAB"** in its report (`benchmark_provenance_claim`) —
  this is qwen leaking prior training knowledge about the *identity of the
  benchmark dataset itself*, not something derivable from any tool result.
  The system prompt explicitly says "Do NOT use prior knowledge about genes,
  cell lines, or variants," and `cases.py`'s own docstring is explicit that
  GIAB/CMRG provenance is metadata for scoring, never sent in the prompt
  text a model sees. A model that has memorized this specific, widely-used
  benchmark BAM well enough to name it unprompted is a genuine and
  independently interesting finding, separate from the adversarial-premise
  question this case was designed to test.

## llama3.1:8b: a real tool-calling bug found and fixed, and a deeper reliability problem uncovered underneath it

The first pass of this benchmark (committed as `0858fe2`) showed llama3.1:8b
making **zero tool calls in 9 out of 9 runs** — every run terminated after a
single turn. Investigating the raw run logs showed the model was not failing
to use tools; it was emitting its tool call as **plain text inside the
message content** (e.g. `{"name": "applicable_layers", "parameters": {...}}`
written as prose) instead of Ollama's structured `tool_calls` field. The
harness's turn loop treated any assistant turn without structured
`tool_calls` as a finished final answer, so it stopped after turn one and
recorded that text as the model's "report."

**Fix:** `ollama_harness.py` now falls back to extracting a `{"name": ...,
"parameters"/"arguments": ...}` object from the assistant's text when
`tool_calls` is empty, and executes it as a real tool call if the name
matches a known tool (`_fallback_tool_call_from_text`). Every call recovered
this way is flagged `via_text_fallback: true` in the run log so this
recovery is visible in the data, not hidden by it. Interestingly, recovering
just the *first* turn this way was usually enough — once one properly-shaped
`tool_calls` entry exists in the conversation history, the model matched that
shape on its own for the rest of the run in every case observed (`fallback`
count per run tops out at 1).

**What the fix uncovered:** with tool calls now actually reaching the MCP
server, llama3.1:8b's real problem became visible — it frequently invents
wrong parameter names (`chr`/`pos` instead of the schema's
`chromosome`/`start`/`end`), passes `applicable_layers` as a call argument to
tools that don't accept it, and does not reliably self-correct from the
resulting validation error even when the error message spells out exactly
which arguments are missing or unexpected. Concretely:

| Run | Tool calls | Of which `is_error` | Outcome |
|---|---|---|---|
| POSITIVE run2 | 20 | 19 | `MAX_TURNS` |
| POSITIVE run3 | 20 | 19 | `MAX_TURNS` |
| NEGATIVE run2 | 16 | 15 | gave up after 16 |
| ADVERSARIAL run1 | 20 | 16 | `MAX_TURNS` |

The model gets stuck re-attempting `applicable_layers` — a tool it already
called successfully earlier in the same run — with malformed arguments,
repeatedly, until the turn cap is hit. In the runs that terminate quickly
instead (e.g. both ADVERSARIAL runs that finish in ~50s with a single tool
call), the model calls one tool successfully and then narrates an intention
("I proceed to call tools in order... First, I run `discordant_pairs`") without
actually emitting the next call, and the harness correctly records that
narration as the final report rather than inventing a call that didn't happen.

**Conclusion reported honestly, per the original benchmark plan:** this is not
a benchmark-harness artifact. llama3.1:8b, at Q4_K_M on this hardware, through
Ollama's tool-calling interface, is not reliable at multi-step structured
tool use against an 11-tool MCP schema — both in initial call formatting
(fixed) and in argument correctness (not fixable in the harness; this is a
capability gap). qwen2.5:7b, given the identical harness, prompt, tools, and
hardware, reliably completes 6-10 well-formed sequential tool calls in every
run. That contrast is itself one of this benchmark's findings.

## `all_layers_queried` was fixed, not just annotated

An earlier version of this document flagged `all_layers_queried` as
penalizing correct behavior: it originally required a literal call to all
four evidence-layer tools regardless of whether `applicable_layers` said a
layer could produce signal for this data, so a model that correctly skipped
a structurally inapplicable `split_reads` (as the system prompt itself
instructs) failed the check by design. Since a criterion that penalizes
methodologically correct behavior is a bug in the criterion, `score.py` was
fixed rather than left annotated: `score_all_layers_queried` now reads the
run's own `applicable_layers` call and only requires the tools covering
layers that call actually reported as applicable, falling back to requiring
all four only when `applicable_layers` was never called successfully (in
which case there's no basis to excuse anything, and `tool_sequence_valid`
penalizes the missing call separately). All tables in this document reflect
the fixed criterion. The only score this changed for the local models was
qwen2.5:7b's ADVERSARIAL case, from 1/3 to 3/3 -- see
`results/BENCHMARK_CLAUDE_BASELINE.md` for the larger effect on Claude's
numbers, and `score.py`'s `score_all_layers_queried` docstring for the exact
fallback logic.

## Visual-tool session: qwen2.5:7b (full writeup: `LLM_SESSION_4_VISUAL_qwen2.5-7b.md`)

Same locus, same MCP server, run through `ollama_harness.py` with
`max_turns=30`. Full detail in the linked writeup; headline findings:

- **Coloring mode mismatched to its own conclusion.** Left `igv_screenshot`
  at the default `color_by="UNEXPECTED_PAIR"` (documented for
  translocations) while its own text concluded "possible deletion" --
  never engaged with `INSERT_SIZE`, the tool's documented recommendation
  for exactly that case.
- **A false "predominant chromosome" claim, verified against its own raw
  data.** Its `discordant_pairs` call returned one mate on each of seven
  different chromosomes -- textbook scattered noise -- but its report
  claimed mates "predominantly" map to chr4 and used that to suggest "a
  possible inter-chromosomal translocation." This is the same
  over-reading-weak-signal pattern documented for the ADVERSARIAL case
  above, except here it appears unprompted, with no false premise in the
  prompt pushing toward that conclusion -- a more concerning version of the
  same failure mode.
- **Tool arguments were correct on the first try, across all 9 calls** --
  worth stating plainly next to the llama3.1:8b exclusion above: reliable
  argument formatting on this schema is not a general property of "local
  model via Ollama," it's specific to qwen2.5:7b having it and
  llama3.1:8b not.
- **Its image descriptions never commit to anything a real look would have
  produced** -- written as generic predictions ("the screenshot should
  show...") rather than descriptions of the actual generated images, which
  turned out (checked directly) to not straightforwardly show what those
  predictions implied.

## A second real bug this pass found: stale-file false success in `run_igv_screenshot`

While preparing the visual-tool session, `run_igv_screenshot`
(`stage1_igv_assistant/tools/bam_tools.py`) turned out to report false
success when a file already existed at the target `output_path`: its
completion check polls for the file's size to be stable across two 1-second
checks, which a pre-existing file satisfies immediately, so IGV gets killed
before it has done any real work and the stale file's stats are returned as
if fresh. Confirmed directly (the same call took 123.5s with a genuine
render against an empty directory, vs. a near-instant false "success"
against one with a leftover file from an unrelated earlier session) and
fixed by deleting any pre-existing file at `output_path` before IGV
launches. Full detail, including why this specifically corrupted the
claude-sonnet-5 visual session's first `evidence_panel` attempt, is in
`LLM_SESSION_4_VISUAL_claude-sonnet-5.md`. This bug predates this benchmark
pass -- any earlier session that re-screenshotted a previously-used locus
was at risk of the same false positive.
