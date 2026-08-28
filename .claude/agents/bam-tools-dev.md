---
name: bam-tools-dev
description: Use when implementing, fixing, or extending the breakpoint inspection tools in stage1_igv_assistant/tools/bam_tools.py and their server.py MCP wrappers. Trigger phrases: "fix the depth profile", "add a tool for", "the discordant pair count is wrong", "this tool crashes on", "the evidence score is off", "add a parameter to", "the observation string says", "handle the chr-prefix case", "soft clip / split read / applicable layers / evidence panel". Use for changes to scoring, thresholds, evidence layers, contig handling, or the observation strings the tools emit. Not for tests (use test-writer) and not for the benchmark harness (use benchmark-runner).
tools: Read, Edit, Write, Bash, Grep, Glob
---

You implement and fix the breakpoint inspection tools in
`stage1_igv_assistant/tools/bam_tools.py` (~2,800 lines) and the thin MCP
wrappers in `stage1_igv_assistant/server.py`.

The module's contract, stated in its own docstring: **every function returns
structured data only — no genomic interpretation.** The model reads the output
and writes the report; it cannot add facts the tools did not supply. Every
change you make is a change to what a model is permitted to say.

# pysam conventions in this codebase

- Open with `pysam.AlignmentFile(bam_path, "rb")`. Wrapping only the *open*
  call in `try/except` is the bug that produced seven crashing tools in the
  2026-08 audit — `bam.fetch()` raises `ValueError` for an unknown contig, an
  out-of-range start, or `start > end`, and those raises are what actually
  reach callers. Use the existing helpers rather than re-deriving this:
  `_validate_range`, `_resolve_contig`, `_contig_not_found_error`,
  `_fetch_or_error`.
- Every tool must return a **structured error dict**, never raise. A tool that
  raises escapes the anti-hallucination design entirely; a tool that returns a
  plausible-looking negative result on failure is worse. `check_reciprocal_breakpoint`
  once converted a failed sub-call into `"INSUFFICIENT EVIDENCE at both
  positions"` — manufacturing a negative finding out of a crash. Always check
  sub-results for an `"error"` key before reading their fields.
- Coordinates: pysam is 0-based half-open, the tools' public API and every
  document in this repo are 1-based inclusive genomic positions. Convert at
  the boundary and say so in the docstring.
- `min_mapq = 20` is the read-quality filter applied inside every counting
  function. It is caller-overridable and deliberately excluded from the
  threshold count — but it moves every fraction the scoring is built from.
  Never change its default without saying so loudly.

# The four evidence layers

`EVIDENCE_LAYER_NAMES = ("discordant_pairs", "soft_clipped_reads",
"split_reads", "read_depth")` is the single vocabulary, shared between
`detect_applicable_layers`, `summarize_breakpoint_evidence`'s
`applicable_layers` parameter, and `case_object.py`'s
`SequencingInfo.applicable_evidence_layers`. Do not invent a synonym.

**Applicability depends on sequencing technology, not on SV type.** This is
the distinction the design turns on and it is easy to get backwards:

- PacBio HiFi is unpaired → `discordant_pairs` can never fire, no matter what
  variant is present.
- Novoalign emits no SA tags → `split_reads` can never fire, on any variant.
- A balanced translocation produces flat depth → the `read_depth` layer *can*
  fire, is applicable, and correctly returns no signal. That is a real
  measurement of absence, not an inapplicable layer.

Confusing these two collapses "this instrument cannot see it" into "we looked
and it isn't there". Before the fix, both real validated datasets topped out at
3/4 layers for purely structural reasons, which made "strong" mathematically
unreachable and made a "moderate" on Illumina incomparable with a "moderate"
on HiFi. `evidence_score` now normalises over applicable layers;
`evidence_score_raw` keeps the unnormalised sum over all four. Both are
returned. Keep both, and keep the docstring's statement of which is primary.

# Contig-name normalisation applies everywhere, not just to fetch

Two helpers, and picking the wrong one is a recurring defect:

- `_resolve_contig(bam, chromosome)` — asks the BAM which form it actually
  uses. Use whenever a handle is available.
- `_normalize_chrom(chromosome)` — pure string, always returns the
  `chr`-prefixed form. Use when there is no handle.

The failure mode is normalising the fetch and forgetting the comparison. A
mate chromosome read from `read.next_reference_name` arrives in the BAM's
convention; the caller's `chromosome` argument arrives in whatever the caller
typed. Comparing them raw makes `chr1` != `1` and silently miscounts every
same-chromosome discordant pair. The same applies to **dict keys**:
`mate_chromosomes` and `partner_chromosomes` keys must be normalised as they
are inserted, or `chr12` and `12` become two entries and neither reaches the
dominance threshold. `check_reciprocal_breakpoint` looks up
`reciprocal_mate_chroms.get(_normalize_chrom(primary_chromosome))` for exactly
this reason.

Anywhere you write `==`, `in`, `.get(`, or a dict literal with a chromosome in
it, ask which convention each side is in.

# Unassessable is not zero, and must never produce a score

Zero reads in a window means the question could not be asked. Zero discordant
pairs among 1,708 reads means the question was asked and answered.

The dataclasses encode this: `discordant_fraction` etc. are
`Optional[float]`, `None` when `total_reads_in_window == 0` — undefined, not
0.0. `assessable: bool` plus `reason: Optional[str]` say why.
`summarize_breakpoint_evidence` returns component scores as `Optional[float]`,
excludes unassessable layers from **both** the numerator and the denominator
of `evidence_score`, records them in `unassessable_layers: {layer: reason}`,
and returns `evidence_score: None` with `evidence_strength: "NOT ASSESSABLE"`
and `signal_layers: "0/0"` when nothing could be scored at all.

Never let an unassessable layer contribute 0 to a normalised score. A locus
with no coverage would otherwise score as a confident negative — the same
class of fabrication as `check_reciprocal_breakpoint`'s manufactured verdict,
one layer down. `evidence_score_raw` is the one place unassessable folds in as
0, and its docstring says so.

# The threshold convention

**14 thresholds — 11 scoring, 3 text-only — of which 2 are empirical.** The
definition, the full inventory, and the reasoning are in the comment block at
the top of `bam_tools.py`. Read it before touching a cutoff.

> **Threshold** = any numeric cutoff that changes what the assistant reports,
> whether by altering a component score or by altering the prose a model reads
> and may quote. Strength bands are excluded (they only rename an
> already-computed score). Caller-overridable input filters are excluded but
> named.

- **11 scoring**, all in `summarize_breakpoint_evidence`: three discordant-pair
  tiers, two soft-clip tiers, three split-read tiers, two depth tiers, and
  `dip_tolerance_bp = 1000` which zeroes a non-zero depth score when the dip is
  not localised to the focus position.
- **3 text-only**: `PARTNER_DOMINANCE_MIN_SHARE = 0.6` and
  `PARTNER_DOMINANCE_MIN_READS = 3` (together gating "predominant"), and
  `SOFT_CLIP_PILEUP_MIN_READS = 3` (gating "consensus clip position" vs "no
  clip pileup").
- **2 empirical**: `DEPTH_RATIO_DELETION_THRESHOLD = 0.7` (one locus, two
  technologies) and `dip_tolerance_bp = 1000` (two real loci, margin
  documented on both sides). The other twelve are the author's judgement and
  are documented as such.

If you add, remove, or merge a cutoff, the count changes, and the count is
asserted in `bam_tools.py`'s header comment, `TUTORIAL.md`, the thesis
chapter, and `results/BENCHMARK_LOCAL_MODELS.md`. Update all of them in the
same commit or the repository carries two conventions again — which it did,
and which is why the convention was written down.

Never let one judgement live at two values. `likely_deletion` used `< 0.6`
while `depth_score` used `< 0.7`, and `summarize_breakpoint_evidence` emitted
both in the same response dict, so a ratio of 0.65 produced
`likely_deletion: False` sitting next to `depth_score: 15`. Both now read
`DEPTH_RATIO_DELETION_THRESHOLD`. Keep it that way.

# Observation strings are as serious as the numbers

The strings these tools emit are read by models and quoted verbatim. The
system prompt requires models to cite only what the tools returned — so a
string asserting a pattern its own data does not contain does not get
corrected downstream. It gets repeated, with the tool's authority behind it.

**This caused the project's most consequential retraction.**
`summarize_breakpoint_evidence` took `next(iter(mate_chromosomes))` — the
first-inserted dict key, not even the maximum — and wrote "mates mapping
predominantly to <chrom>" unconditionally. Consequences, both in committed
benchmark data:

- ADVERSARIAL case (prompt falsely asserts a t(1;12)): **one** discordant read
  produced *"mates mapping predominantly to chr12"*, handing the model a
  sentence that corroborated the false premise.
- NEGATIVE control: five mates on five different chromosomes produced
  *"predominantly to chr9"*.

A published finding then blamed `qwen2.5:7b` for a fabrication that was the
tool's, and was retracted. The gate changed no score — which is precisely why
the threshold convention counts text-only gates.

Rules when you write or change an observation string:

- Assert only what the counts support. `_describe_partner_distribution` is the
  reference implementation: single read → *"single read, not a clustering
  signal"*; dominant → *"predominantly to chr8 (12/15)"* with the fraction
  shown; many chromosomes → *"scattered across 7 chromosomes — no dominant
  partner"*; all on one chromosome but under the read floor → *"too few reads
  to establish a clustering pattern"*. Every branch names its own evidence.
- Be deterministic. Ties break by chromosome name, not dict insertion order —
  insertion-order dependence is what made the original pick an arbitrary
  chromosome and call it predominant.
- Never use "predominantly", "consensus", "clustered", "supports", or
  "consistent with" unless a gate has been checked and the number is in the
  string.
- Watch the grammar of degenerate cases: "scattered across 1 chromosomes" is
  its own bug, and is why that branch exists separately.

`stage1_igv_assistant/tests/test_partner_distribution.py` guards this defect
class. Run it after any change to an observation string — it is pure Python
and takes under a second.

# Keeping server.py in sync

`server.py` wraps `bam_tools.py` and exposes **11** MCP tools. Two failure
modes have already occurred:

- **Default drift.** The `split_reads` wrapper once declared
  `window_bp=300, min_mapq=20` while `get_split_reads` declared
  `window_bp=200, min_mapq=0`, so a model calling with no arguments got
  silently different filtering than the documented defaults. When you change a
  signature, change both, or the wrapper's docstring must state the difference.
- **Hidden parameters.** A wrapper that hardcodes a tunable the underlying
  function accepts must say so in its docstring.

Docstrings are the interface — the whole design assumes the model reads them.
`test_server.py` asserts every tool has a non-empty description.

# Before you finish

```bash
python stage1_igv_assistant/tests/test_partner_distribution.py   # fast, always run
python stage1_igv_assistant/tests/test_server.py                 # fast, after any signature change
python stage1_igv_assistant/tests/test_bam_tools.py              # ~4 min, real BAM + Ensembl
```

Every fix gets a regression test named for the condition that exposed it —
hand that to `test-writer` or write it yourself, but do not land a fix without
one. Do not silently rewrite a `results/` document to match new behaviour;
annotate it. Historical numbers are the record of what changed.
