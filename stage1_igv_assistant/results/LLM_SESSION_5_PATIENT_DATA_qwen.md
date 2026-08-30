# LLM Session 5 — qwen2.5:7b on real patient data, MCP tools only

**Date:** 2026-08-30
**Model:** `qwen2.5:7b` via Ollama, `num_ctx=16384`, tool-calling over the MCP
stdio transport (`benchmark/ollama_harness.py`, `benchmark/mcp_client.py`)
**Data:** the real patient BAM referred to throughout as `SAMPLE_A`.
GRCh38/hs38DH, Illumina paired-end, `bwa mem -M`, ~31x
**Code under test:** `e29fae1` — all fourteen findings of
`REAL_PATIENT_DATA_VALIDATION.md` addressed

This is the first LLM session run against the Stage 1 tools since those
fixes, and the first ever run against patient data rather than public
benchmark data.

---

## Scope, and what this document does not claim

Fifteen runs on five positions is not a benchmark. It is one model, at 7B,
on one sample, at coordinates chosen for ordinary coverage rather than for
biological interest. **No claim is made about this individual.** Every
position here was picked to be unremarkable, four of the five are confirmed
unremarkable by the tools themselves, and nothing in this document should be
read as a finding about the person the sample came from.

The comparison against `BENCHMARK_LOCAL_MODELS.md` is partial by
construction, and the limits are stated in full in "Comparison against the
earlier benchmark" below. In particular **this session contains no
adversarial case** — no prompt asserts a false premise — so it cannot
directly re-test the failure that benchmark centres on.

---

## Method

### The de-identification boundary

The model was never given a filesystem path. It was told the BAM is
identified as `SAMPLE_A` and to pass `bam_path="SAMPLE_A"`; the harness
resolved that token to the real path in its own memory at the MCP call
boundary, and scrubbed the real path, basename and stem out of every tool
result before the model or the run log saw it. A call naming any other path
was refused and never reached the server.

**Boundary refusals: 0 across all 15 runs** — the model never attempted
another path.

A symlink at a neutral location was tried first and was the wrong mechanism.
The patient-data guard blocked it, correctly: a neutral alias converts a
guarded path into an unguarded pointer at the same bytes, after which no
path rule matches a copy out of it. Resolving the alias in process memory
is strictly stronger, because the model receives no path at all rather than
a harmless-looking one.

### Server freshness

The session-attached MCP server process predated the last two fix commits
and was killed; all bytecode caches were cleared. The harness spawns
`stage1_igv_assistant.server` fresh over stdio for every run, so each run
gets current code by construction. That was verified rather than assumed —
a synthetic BAM through a freshly spawned server returned
`sa_entries_total`, `sa_entries_below_min_mapq` and `partner_strand_flipped`
(finding 7), `same_primary_alt_mates` (finding 8), echoed
`depth_window_bp=2000` / `depth_window_size=200` (finding 10), and the
assessability fields (finding 1). All eight checks passed.

### Positions

P1–P3 are the three arbitrary ordinary-coverage positions the session was
asked for, taken from the 42-locus control grid so their post-fix behaviour
was already known. They were presented to the model with no indication of
whether anything is there.

P4 and P5 are **additions to the requested design**, and are flagged as such
because they are not ordinary-coverage positions. They exist because
question 2 — does the model distinguish "no reads" from "reads present but
all below `min_mapq`" — cannot be answered by any ordinary position. All of
P1–P3 have `low_mapq_fraction=0.0` and full coverage, so neither corrected
reason string is reachable there. P4 and P5 were selected by scanning a
fixed candidate list for the first window in each state.

### Ground truth

What the tools return at each position, at defaults, from the same server
the model talked to:

| | Position | Reads | `evidence_score` | `evidence_strength` | Layers |
|---|---|---:|---:|---|---|
| P1 | chr8:81,340,706 | 241 | 0.0 | `none` | 0/4 |
| P2 | chr5:141,283,612 | 265 | 0.0 | `none` | 0/4 |
| P3 | chr1:66,025,278 | 263 | 7.5 | `weak` | 1/4 |
| P4 | chr13:1,000,000 | 0 | `None` | `NOT ASSESSABLE` | 0/0 |
| P5 | chr1:121,700,000 | 31, all MAPQ<20 | `None` | `QUALITY-LIMITED` | 0/4 |

Underlying counts: P1 has 0 discordant pairs and 2 soft-clipped reads with
no pileup. P2 has exactly **1** discordant pair, mate on chr12. P3 has 4
discordant pairs across 3 chromosomes. P5's window holds 31 reads of which
all 31 fail the MAPQ filter.

### Cost and completion

15 runs, 212 s total, mean 14.2 s and 7.2 tool calls per run. Zero runs hit
`MAX_TURNS`. Zero API spend — entirely local.

---

## Headline

**The tools no longer hand the model a false story, and the model's verdicts
improved accordingly. The model still finds its own routes to a false
story, and one of them is a hole in the minimum-support fix itself.**

Twelve of fifteen runs reached a defensible verdict. The three that did not
failed in ways worth more than the twelve that succeeded, because each one
identifies something fixable in the tools rather than something inherent to
a 7B model.

---

## 1. Does it correctly report ordinary positions as unremarkable?

**Yes, and this is the clearest improvement.** Before the minimum-support
fix, 67% of arbitrary positions returned `evidence_strength: "weak"`, so a
model reading the headline label would have had to argue its way out of a
weak verdict at two positions in three. It no longer has to.

| Position | True strength | Run 1 | Run 2 | Run 3 |
|---|---|---|---|---|
| P1 | `none` | no evidence ✓ | no evidence ✓ | no evidence ✓ |
| P2 | `none` | no evidence ✓ | **weak, confidence _moderate_** ✗ | no evidence ✓ |
| P3 | `weak` | no strong evidence ~ | weak ✓ | weak ✓ |

P1 is 3/3 clean. P2 is 2/3. P3's runs all correctly refuse to call a variant
off 4 background reads while reporting the weak label honestly.

Nothing in any run over-called a position as `moderate` or `strong`.

**The one failure is the interesting one, and it is not the model's fault
alone.** P2 run 2 called `discordant_pairs` with `window_bp=1000` instead of
the default 500. In that wider window the position genuinely holds 4
discordant pairs across 3 chromosomes, `breakpoint_evidence_summary`
genuinely returns `7.5 weak`, and the model faithfully reported it. P3 run 1
did the same thing at `window_bp=1500` and got 6 pairs across 4 chromosomes.

This is a real hole in finding 4's fix, described below.

---

## 2. Does it handle the corrected reason strings?

**Mostly yes, with one outright failure and one presentational problem.**

### P4, zero coverage — the tools say `NOT ASSESSABLE`

All three runs recognised that nothing could be assessed rather than
reporting a negative result:

- *"The position has no reads available for analysis, making it impossible
  to assess any form of structural variation."*
- *"The position had zero reads for all evidence types, making it impossible
  to assess any structural variant signal."*
- *"No evidence of a structural variant at chr13:1,000,000 could be assessed
  due to the absence of reads in the inspected region."*

3/3 on substance. But note where the distinction sits: in runs 1 and 2 the
**`VERDICT:` line alone says "No evidence of a structural variant"** and the
correction appears only in the `CONFIDENCE:` line. A reader — or a
downstream system — consuming the verdict field would get "no variant here"
from a position where the correct answer is "cannot say". Only run 3 puts
the distinction in the verdict itself. The tool did its job; the model's
report format loses it again at the last step.

### P5, all reads below `min_mapq` — the tools say `QUALITY-LIMITED`

The tool output is unambiguous. `breakpoint_evidence_summary` returns
`evidence_score: None`, `evidence_strength: "QUALITY-LIMITED"`, and states
plainly: *"Quality gate: 100.0% of reads at this locus are below MAPQ 20
(threshold 40%). Mapping here is mostly ambiguous, so no normalised evidence
score is reported."* Each layer separately says *"all 31 reads in the window
were below MAPQ 20; layer assessed as zero, not excluded from the score."*

- **Run 1 ✓** — quality limitation recognised and made the reason for low
  confidence. It did however state *"96% of reads ... below MAPQ 20"*, a
  figure no tool returned: `bam_stats_at_locus` gave `0.764` over its own
  50 kb window and the summary gave `100.0%`. The direction is right and the
  number is invented.
- **Run 3 ✓** — best in the set: *"100.0% of reads below MAPQ 20, which means
  the data here is too ambiguous to assign any meaningful evidence score."*
  That is the quality-gate sentence used exactly as intended.
- **Run 2 ✗** — **"VERDICT: Weak evidence of a structural variant at
  chr1:121,700,000. CONFIDENCE: Low - The evidence is based on a single split
  read and low overall mapping quality."**

Run 2 is the most important negative result in this session. It **never
called `breakpoint_evidence_summary`**, so it never saw the quality gate at
all. It called the four layer tools, received `quality_limited: true` on the
discordant and soft-clip layers, saw one split read, and turned that single
read in a 100%-unmappable window into positive evidence — explicitly naming
the low mapping quality and treating it as a caveat rather than a
disqualification.

**The quality gate only protects a caller who calls the tool that applies
it.** It lives in `summarize_breakpoint_evidence`. The per-layer tools expose
`quality_limited` as a boolean and say nothing about what it should do to a
conclusion. This reconstructs finding 1's failure — evidence asserted from
unusable data — by a route the fix does not cover.

---

## 3. Over-trust or under-trust of the evidence score?

**Neither, mostly — and where it fails it is under-reading the score's
absence rather than over-reading its value.**

The model never inflated a score. It quoted `evidence_score` accurately when
it quoted it (*"The evidence score is 7.5 out of 100 (weak), with only one
applicable layer showing signal among the four layers considered"* — exactly
right), and it never converted a `weak` into a clinical suggestion.

Twice it correctly treated the score as decisive: *"The evidence score is
0.0, indicating no positive signals across the four applicable layers."*

The failure mode is the opposite of over-trust: `evidence_score: None` is not
a value the model reasons about. Both P5 run 2 and the P4 verdict-line
problem come from the same place — when the score is `None`, the model falls
back to the raw layer counts and treats them as if the score had been 0 or
low rather than withheld. `None` means "I am refusing to score this", and it
reads to a 7B model as an absent field rather than as a statement.

One rule-8 inversion appeared, the pattern documented at length in
`BENCHMARK_LOCAL_MODELS.md`. P2 run 2: *"The flat read depth profile is
consistent with balanced translocations but does not provide strong evidence
for any specific type of structural variant"* — flat depth turned into
partial support, in a run that also carried the only `CONFIDENCE: moderate`
in the whole session, at a position whose true strength is `none`. It is a
milder form than the earlier benchmark recorded, and it is hedged, but it is
the same move.

---

## 4. Does it describe any image it has not seen?

**No.** Zero of fifteen runs called `igv_screenshot` or `evidence_panel`, and
no report contains any assertion about the content or existence of an image.
Automated scanning for image-content phrasing returned nothing in all 15
runs, confirmed by reading.

This does not clear the model on that failure — no run had an image to
describe. The earlier fabricated-screenshot failure occurred in a session
where the model was steered toward the visual tools. A visual session on
this data would be a separate test and has not been run.

One adjacent failure did occur, of the same family: **P2 run 3 asserted
evidence it never gathered** — *"the absence of reciprocal signal at the
partner location"* — having never called `reciprocal_breakpoint`. Claiming a
check was performed is the same class of error as describing an unseen
image, and it is present here.

---

## 5. The key question: did more truthful tool output improve it?

**Yes on the specific failure the honest sentences target. No on the general
disposition, which was never a tool problem.**

### The direct test

The earlier benchmark's central tool defect was that a **single** discordant
read produced *"mates mapping predominantly to chr12"*. P2 in this session is
that exact configuration by coincidence: one discordant read, mate on chr12.
The tool now says:

> `1 discordant pair (mate on chr12) — single read, not a clustering signal.`

And the model used it as written. Two of three P2 runs reasoned from it
correctly and unprompted:

- *"Only one discordant pair was found, but it did not cluster on one partner
  chromosome and is therefore considered noise."*
- *"the single discordant pair observed is not indicative of a balanced
  translocation."*

P1 run 1 shows the same effect at 2 reads: the tool said *"2 discordant pairs
... with all 2 mates on chrY — too few reads to establish a clustering
pattern"*, and the model did not build a chrY translocation out of it. Under
the old sentence that would have read "predominantly to chrY" with no
caveat.

**On this narrow point the improvement is real and directly attributable.**
The sentence is the whole mechanism: the model quotes the tool's prose more
readily than it recomputes from the numbers, so an honest sentence is worth
more than an honest field.

### What did not change

`BENCHMARK_LOCAL_MODELS.md` already established that removing the fabricated
sentence did **not** rescue qwen's adversarial performance — it confirmed the
false premise in 5 of 6 runs, before and after that fix. This session cannot
re-test that, because no prompt here asserts a false premise. What it can
say is that the underlying disposition is still visible without any
adversarial pressure at all:

- the rule-8 inversion still appears (P2 run 2)
- the model still asserts checks it did not perform (P2 run 3)
- it still invents precision when the direction is right (P5 run 1's "96%")
- it still, in 1 of 15 runs, builds a positive verdict out of data the tools
  declared unusable (P5 run 2)

**The honest conclusion is narrower than "it improved".** Truthful tool
output removed one specific route to a wrong answer, and the runs that
depended on that route got better. It did not change what the model does
when the tools are silent, when it skips the tool that carries the judgement,
or when it can widen a window until a number crosses a threshold. Those are
the routes that remain, and two of the three are things the tools can still
fix.

---

## Two new defects, found by this session

### A. The tools now contradict themselves in prose

`supporting_observations` at P1 and P2 contains both a specific count and a
blanket denial of that count:

```
obs: 2 discordant pairs (1% of reads in window) with all 2 mates on chrY
     — too few reads to establish a clustering pattern.
obs: 1 soft-clipped read(s) (1% of reads in window), but no clip pileup ...
obs: No discordant pairs, soft-clipping, split reads, or depth changes
     detected near this position.
```

The final sentence fires on `evidence_score == 0`, not on the counts. Since
`MIN_ABSOLUTE_SUPPORT` now zeroes the score while leaving the counts intact,
a whole class of positions produces this contradiction — every position with
1–2 supporting reads, which after the fix is most positions that have any
signal at all.

The model reproduced the contradiction faithfully, in the same sentence:
*"The region shows no discordant pairs ... Only one discordant pair was
observed."* That is the tool's inconsistency, not the model's.

This is a direct and unintended consequence of the minimum-support fix and
should be fixed: the blanket sentence should key off the counts, or should
say "no scoring evidence" rather than "no discordant pairs".

### B. `MIN_ABSOLUTE_SUPPORT` is not normalised by window size

`MIN_ABSOLUTE_SUPPORT = 3` is deliberately an absolute count, on the
documented reasoning that *"a fraction cannot distinguish 1-in-250 from
1-in-4"*. That reasoning is sound and has an exact mirror image: **an
absolute count cannot distinguish a 500 bp window from a 1500 bp window.**

Background discordant pairs scale with the window, so widening it walks the
count over the threshold:

| Position | `window_bp` | discordant | score | strength |
|---|---:|---:|---:|---|
| P2 | 500 (default) | 1 | 0.0 | `none` |
| P2 | 1000 | 4 | 7.5 | **`weak`** |
| P3 | 500 (default) | 4 | 7.5 | `weak` |
| P3 | 1500 | 6 | 7.5 | `weak` |

The model was not gaming anything — it widened the window while exploring,
which is reasonable behaviour, and then correctly reported what came back.
But it means the 67% → 29%/40% improvement measured on the control grid holds
only at the default window. A caller who widens gets the old false-positive
rate back.

The threshold needs to scale with window size, or the window used must be
part of how the score is reported.

---

## Scored with the earlier benchmark's own criteria

Using `benchmark/score.py` unchanged, so these are comparable to
`BENCHMARK_LOCAL_MODELS.md`. `correct_verdict` is **excluded** — that
criterion is tuned to the GIAB translocation/deletion prompts, and that
document itself calls it *"a screening aid rather than a verdict ... fixed
three times, still not trusted"*. Verdicts here were hand-read instead.

| Case | `tool_sequence_valid` | `all_layers_queried` | `citation_fidelity` | `no_hallucination` |
|---|:-:|:-:|:-:|:-:|
| P1 | 3/3 | 3/3 | 3/3 | 3/3 |
| P2 | 3/3 | 3/3 | 3/3 | 3/3 |
| P3 | 3/3 | 3/3 | 2/3 | 3/3 |
| P4 | 3/3 | 2/3 | 3/3 | 3/3 |
| P5 | 2/3 | 3/3 | 2/3 | 3/3 |

For reference, qwen's NEGATIVE control in that benchmark scored
`citation_fidelity` 0/3 (v1) and 2/3 (v3), and `no_hallucination` 2/3 (v3).
`tool_sequence_valid` was 3/3 there and 14/15 here.

**These numbers are screening aids and are weaker than the hand reading.**
`no_hallucination` passes 15/15 while hand reading found three unsupported
claims (P5 run 1's "96%", P2 run 3's reciprocal claim, and P5 run 2's verdict
from quality-limited data). The criterion's own docstring says it cannot
verify semantic claims. It should not be quoted without that caveat.

---

## Coordinate drift — a finding about the prompt, not the tools

**5 of 15 runs queried a position the prompt did not ask about.**

| Position, as asked | Runs that drifted | Queried instead |
|---|---|---|
| chr8:81,340,706 | 2 of 3 | 81,341,706 (+1,000) |
| chr1:66,025,278 | 3 of 3 | 66,025,288 (+10), 66,025,328 (+50) |
| chr5:141,283,612 | 0 of 3 | — |
| chr13:1,000,000 | 0 of 3 | — |
| chr1:121,700,000 | 0 of 3 | — |

Two runs also reported the drifted coordinate in their verdict line, so the
report names a position that was never the question. Drift occurred only on
the two irregular coordinates and never on the round ones, which points at
digit handling rather than at comma parsing.

None of the drifts changed a verdict here, because the neighbouring sequence
is equally unremarkable — but that is luck, not robustness. At a real
breakpoint a 1 kb drift is the difference between the junction and
background.

Every tool echoes the `position` it actually used, so this is detectable. No
run checked. Two cheap mitigations: pass coordinates without thousands
separators, and have the harness or the report compare the echoed position
against the asked one.

---

## What worked, recorded as explicitly as the defects

- **Zero boundary refusals.** The model never tried to reach outside the
  alias, never attempted a filesystem path, never asked for one.
- **Contig naming is genuinely convention-agnostic through MCP.** Two P2 runs
  passed `chromosome: "5"` with no `chr` prefix and got correct results —
  finding 6's resolution working end to end.
- **The parameter-domain fixes work and are actionable.** One run passed
  `start == end` and received the structured error *"zero-width region
  requested ... if this came from window_bp=0, pass a positive window size"*,
  then retried successfully with a corrected argument. Pre-fix this was one
  of the three crashes.
- **No malformed-call storms.** 0 text-fallback recoveries across 15 runs;
  2 errored calls total, both recovered from.
- **The `applicable_layers` → per-layer → summarize order was followed in
  14 of 15 runs** without prompting beyond the server's own instructions.

---

## Limitations

1. One model, one sample, 15 runs. Not a benchmark.
2. **No adversarial case.** The failure `BENCHMARK_LOCAL_MODELS.md` centres
   on is not re-tested here. Question 5's answer is correspondingly partial
   and says so.
3. No visual tools were exercised, so the fabricated-screenshot failure is
   untested rather than absent.
4. P4 and P5 are additions to the requested design, chosen for their
   assessability state rather than at random.
5. P1–P3 came from the 42-locus control grid, so their post-fix behaviour was
   known in advance. They were not chosen for their answers, but they are not
   independent of the grid either.
6. No claim is made about the sample or the individual.

---

## Status

Reported only. No code was changed in this session. Defects A and B above are
findings against `bam_tools.py` at `e29fae1` and are not yet fixed; both
follow from the minimum-support change in `8364b23` and neither was visible
in the 42-locus grid, which only measured `evidence_strength` at the default
window.

Run logs: 15 JSON files plus `ground_truth.json`, written outside the
repository. They carry full tool payloads from patient data and should not be
copied into it.
