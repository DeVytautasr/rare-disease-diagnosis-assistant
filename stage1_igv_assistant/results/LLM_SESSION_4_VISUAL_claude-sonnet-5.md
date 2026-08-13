# LLM Session 4 (Visual): claude-sonnet-5

One-off qualitative session testing visual-tool choice and, most importantly,
whether the model's written description of a generated image matches what
the image actually shows. Run via `benchmark/claude_harness.py` (API
harness, not Claude Code) so context is identical to the qwen2.5:7b session
in this same pair — see `LLM_SESSION_4_VISUAL_qwen2.5-7b.md`. Full run log:
`LLM_SESSION_4_VISUAL_claude-sonnet-5.json`. Case definition:
`benchmark/cases.py`'s `"VISUAL"` entry. `max_turns=30` (raised from the
standard 20 for this session only, via a new optional parameter on
`run_once` — the three-case benchmark's default is unchanged).

**Session stats:** 8 tool calls, 0 errors, 0 malformed arguments, 129.3s
wall-clock, **$0.1897** (thinking disabled, per the existing cost-control
rationale in `BENCHMARK_CLAUDE_BASELINE.md`).

## A real bug this session surfaced: stale-file false success in `run_igv_screenshot`

Before any assessment of Claude's behavior, a tool bug has to be disclosed,
because it initially made this session's own images untrustworthy.

Claude's `evidence_panel` call used `output_dir="./chr1_115686862_evidence"`
— a directory that already existed with PNGs from the original benchmark's
POSITIVE-case pilot run, over an hour earlier.

> **That directory no longer exists, and that call can no longer be made.**
> `output_dir` was removed from the tool signature entirely (FIX C): output
> paths are now assigned by the server under its own session directory, and
> a caller-supplied path is not accepted. The collision described here — a
> model choosing an output directory that another run had already written
> to — is therefore no longer reachable through the MCP tools. The path is
> retained in this paragraph because it is what the run actually did; see
> `benchmark/runs/README.md` for the stage layout and
> `BENCHMARK_CLAUDE_BASELINE.md`'s FIX C section for why the parameter was
> removed rather than the behaviour merely discouraged.
 `run_igv_screenshot`'s
completion check polls `output_path` for a stable (unchanged) file size
across two 1-second checks, then kills IGV. If a file already exists at that
path when the poll starts, its size is trivially "stable" from the very
first check — IGV gets SIGTERM'd before it has done more than start the
JVM, and the **stale** file's size is reported back as a fresh success.
Confirmed directly, not inferred: the same call re-run against a genuinely
empty directory took **123.5s** with `shutdown_method: clean_exit`; against
the pre-existing directory it returned `shutdown_method:
terminated_after_snapshot` in a small fraction of that time. The images
Claude's own tool call "generated" in this session were therefore leftover
files from an unrelated earlier run, coincidentally at the same locus —
not proof the call was fabricated, since content happened to still be
locus-correct, but not evidence the call worked, either.

**Fixed** in `stage1_igv_assistant/tools/bam_tools.py`'s
`run_igv_screenshot`: any pre-existing file at `output_path` is now deleted
before IGV launches, so the file's mere existence during polling is
unambiguous proof of a fresh write. Verified the fix twice: once against an
empty directory (123.5s, `clean_exit`), once by re-running against the now
non-empty directory the fix had just populated (40.0s, `clean_exit`, new
mtime, slightly different file size) — confirming it doesn't get fooled by
its own prior output either.

**The images assessed below are the fixed, freshly-regenerated ones**, using
Claude's exact call arguments (same `chromosome`, `position`, `output_dir`,
`applicable_layers`, `start`, `end` as its own tool call — reproducible from
`batch_script` in the run log). Claude's original report text is unedited —
it was written against what the model actually received back (which, per
the harness's own construction, was never the image pixels either way; see
below) — but the report is graded against the correct images, not the
accidentally-reused stale ones.

## Visualization choice

Claude used `evidence_panel` only, explicitly declining `igv_screenshot`:

> "Given the numbers point to a deletion signature ... not a translocation, I
> used the per-layer evidence panel rather than a single UNEXPECTED_PAIR
> screenshot ... a translocation-oriented single view would emphasize
> inter-chromosomal pairing, which the numbers show is not the dominant
> signal here."

This is a stronger choice than picking a `color_by` value: `evidence_panel`
sidesteps the single-coloring-mode decision entirely, since each layer gets
its own tool-managed settings (documented in `evidence_panel`'s own
docstring — a single shared window can't serve every layer well). It also
means this session doesn't directly test whether Claude knows
`igv_screenshot`'s documented recommendation (`INSERT_SIZE` for deletions,
not the `UNEXPECTED_PAIR` default) — see the qwen session for that test,
where the model did use `igv_screenshot` and left the mismatched default.

## Assessment

**1. Panel/coloring choice appropriate to a suspected deletion?** Yes —
see above. `split_reads` was correctly skipped (excluded from the panel as
inapplicable, matching `applicable_layers`' own determination for this
aligner).

**2. Tool arguments correct on the first try?** Yes, across all 8 calls —
0 errors, 0 malformed arguments (see full log).

**3. Did it notice the discordant-pairs "looks busy but n=1/1708" trap
(`EVIDENCE_PANEL_VALIDATION.md`)?** Partially, and the partial failure is
specific and worth stating precisely. Claude's *numeric* verdict never falls
for it — it correctly treats discordant_fraction=0.001 (1 of 1,708 reads) as
"background-level noise, not a clustered translocation signal" throughout.
But its *description of the image itself* asserts the opposite of what the
image shows (see #4) — so the trap was avoided at the number-interpretation
level and walked into at the image-description level, by a model that, per
#4, never saw the image to check either way.

**4. MOST IMPORTANT — does the written description match the actual image,
or is it plausible narration of an unviewed image?** Verified by generating
the corrected images (above) and viewing them directly.

For `discordant_pairs.png`, Claude wrote:

> "the image should show essentially uniform, non-anomalous pairing,
> corroborating that this is not a translocation."

**The actual image does not show this.** Viewed directly: a substantial
cluster of reads — roughly a third to half of those visible — render in
IGV's anomalous-pair red, forming a clear wedge shape in the left-center of
the frame. This is not subtle. The almost-certain explanation is IGV's
`UNEXPECTED_PAIR` coloring also flagging same-chromosome pairs with
anomalous **insert size** — exactly what a deletion between paired mates
produces — not just cross-chromosome mates. That's consistent with the
numbers (a real deletion is present; true cross-chromosome discordance is
noise-level), but it means the image is visually busy for a reason Claude's
report doesn't mention, and directly contradicts the "essentially uniform"
claim. A clinician who read this report and then glanced at the actual
image would notice the mismatch immediately.

For `soft_clipped_reads.png`, Claude's description ("a vertical wall of
aligned read starts/ends at ~115,686,865") **does** match the image — there
is a clear, tight vertical marking at that column. For `read_depth.png`,
Claude's description ("a visible step-down ... beginning right around the
center") is **directionally consistent but oversold** — the image shows a
gentle, undulating decline rather than a sharp step; "step-down" implies
more visual abruptness than is actually there.

**Why this happened, confirmed by reading the code, not inferred:**
`claude_harness.py` constructs every tool result as
`{"type": "tool_result", ..., "content": str(result["payload"])}` — a plain
string. `mcp_client.py`'s `call_mcp_tool` only ever extracts
`structuredContent` (a JSON dict) or parses text content; the underlying
`evidence_panel`/`igv_screenshot` MCP tools return a plain dict with a
`screenshot_path` string, never an MCP image content block. **The model
never receives the image's pixels through this pipeline — only a file path
as text.** Every visual "description" in this session (and, by the same
construction, in every other run across this whole benchmark) is
necessarily synthesized from the numeric tool data already in context, not
from observation. In this instance, that synthesis produced one accurate
description (soft-clip), one oversold-but-directionally-right description
(depth), and one description that flatly contradicts the real image
(discordant pairs) — because nothing in the pipeline could have corrected
it against the actual pixels. This is a harness-architecture gap, not a
claude-sonnet-5-specific failure: `ollama_harness.py`'s tool-result
construction has the same shape, and qwen2.5:7b isn't a vision model in the
first place (see the qwen session).

**Not attempted in this pass:** actually wiring a vision-capable image
content block into `claude_harness.py`'s tool-result construction (feasible
in principle — Claude models support image input, and the file is right
there on disk) so a follow-up test could show whether a *genuinely* shown
image changes what gets written. That's a real, scoped follow-up, not done
here to keep this session to what was asked: run it, view the images,
compare.
