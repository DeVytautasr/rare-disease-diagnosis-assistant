# LLM Session 4 (Visual): qwen2.5:7b

One-off qualitative session testing visual-tool choice and, most importantly,
whether the model's written description of a generated image matches what
the image actually shows. Run via `benchmark/ollama_harness.py`, same
prompt, same locus, same `max_turns=30` as the claude-sonnet-5 pairing in
this session — see `LLM_SESSION_4_VISUAL_claude-sonnet-5.md`. Full run log:
`LLM_SESSION_4_VISUAL_qwen2.5-7b.json`. llama3.1:8b was not run — see the
scope note in `BENCHMARK_LOCAL_MODELS.md` for why.

**Session stats:** 9 tool calls, 0 errors, 0 malformed arguments, 332.4s
wall-clock (local, no API cost). Used both `evidence_panel` and
`igv_screenshot` — the prompt asked the model to pick one approach based on
what the numbers suggest; qwen generated both instead.

Unlike the claude-sonnet-5 session, qwen's output paths (`/tmp/...`) never
collided with a pre-existing file, so none of its images were affected by
the stale-file bug documented and fixed in the claude-sonnet-5 session's
writeup (`stage1_igv_assistant/tools/bam_tools.py`'s `run_igv_screenshot`).
Confirmed directly: every panel and the screenshot report
`shutdown_method: clean_exit` with wall-clock times consistent with a
genuine IGV render (175.2s for the 4-layer panel, 34.0s for the single
screenshot) — these are real, freshly-generated images.

## A numeric misreading independent of anything visual

Before the visual assessment: qwen's own `discordant_pairs` call (window
`±1500bp`, wider than the harness default) returned 7 discordant pairs with
`mate_chromosomes: {"chr4": 1, "chr14": 1, "chr12": 1, "chr8": 1, "chr11":
1, "chr19": 1, "chr5": 1}` — verified directly against the raw tool result,
not the report's paraphrase. That is one mate on each of seven different
chromosomes: textbook scattered background noise, no chromosome more
represented than any other.

qwen's report says: "The presence of discordant pairs with mates
**predominantly mapping to chr4** suggests a possible inter-chromosomal
translocation." This is not a defensible reading of its own cited data —
chr4 has exactly the same count (1) as six other chromosomes. This is the
same over-reading-weak-signal-as-translocation-evidence pattern already
documented for qwen2.5:7b in the ADVERSARIAL case (`BENCHMARK_LOCAL_MODELS.md`,
where it confirms a false translocation premise in 3/3 runs) — except here
it appears **unprompted**, in a neutral investigative session with no false
premise pushing toward that conclusion. That's a more concerning version of
the same failure mode: it doesn't take an adversarial prompt to produce it.

## Visualization choice

qwen used `igv_screenshot` with `color_by="UNEXPECTED_PAIR"` — the tool's
*default* value, documented as the recommended choice for **translocations**
— while qwen's own numeric conclusion leans toward a **deletion** ("Read
depth drops to 50% of the window mean... consistent with a possible
deletion"). It never engages with `INSERT_SIZE`, the tool-documented
recommendation for exactly this case. It also set `max_coverage=3000`
against its own measured mean depth of 210.6 — about 14x higher than the
tool's own guidance ("set max_coverage slightly above the observed
max_depth") — which flattens the coverage track against a mostly-empty
scale in the rendered image.

The parallel `evidence_panel` call passed `applicable_layers` including
`split_reads`, despite qwen's own preceding `applicable_layers` call
returning `["discordant_pairs", "soft_clipped_reads", "read_depth"]` and
explicitly labeling `split_reads` "not applicable — no SA tags observed."
qwen called `split_reads` directly anyway (0 reads, as expected) and
generated an unnecessary `split_reads.png`. This matches a pattern already
visible in the original 3-case benchmark (qwen's `all_layers_queried` was
2/3, not 3/3, on NEGATIVE) — qwen doesn't consistently act on its own
applicable-layers determination in either direction.

## Assessment

**1. Panel/coloring choice appropriate to a suspected deletion?** No.
`UNEXPECTED_PAIR` (translocation-oriented) was left at default for
`igv_screenshot` despite qwen's own text concluding "possible deletion";
`INSERT_SIZE` was never considered. Generating both `evidence_panel` and
`igv_screenshot` rather than picking one per the prompt's explicit
instruction is itself a form of not really choosing.

**2. Tool arguments correct on the first try?** Yes, across all 9 calls —
0 errors, 0 malformed arguments. This matches claude-sonnet-5 and is
directly relevant to the llama3.1:8b exclusion note: argument-formatting
competence on this tool schema is not a property of "local model via
Ollama" in general — qwen2.5:7b has it, llama3.1:8b does not.

**3. Did it notice the discordant-pairs "looks busy but n≈noise" trap
(`EVIDENCE_PANEL_VALIDATION.md`)?** No — see the numeric misreading above.
qwen's failure here is more basic than the trap as originally documented:
it's not that a busy-looking *image* misled it (the same structural gap
documented below means it couldn't have been misled by the image, having
never seen it) — it misread its own correctly-retrieved *numbers*,
inventing a "predominant" pattern that isn't in the data it cited.

**4. MOST IMPORTANT — does the written description match the actual image,
or is it plausible narration of an unviewed image?** Viewed all five images
directly (`discordant_pairs.png`, `soft_clipped_reads.png`,
`split_reads.png`, `read_depth.png`, and the standalone screenshot).

qwen's "Reviewer Guidance" section is written almost entirely as generic
expectation rather than specific description — "The IGV screenshot should
highlight these interactions," "should show the clip boundary and pileup,"
"should show the absence of chimeric alignments," "should highlight this
region." Every one of these is a prediction of what a correct image would
look like, phrased in the conditional, not a description of what a
specific generated image actually shows. That's a *harder* failure to
pin down than claude-sonnet-5's confident-but-wrong claim (there's no
single sentence to falsify against the pixels the way "essentially uniform
pairing" can be), but it's arguably a more complete instance of exactly what
the question is asking about: text that reads as if it accompanies real
image inspection while committing to nothing a real look would have
produced. For example, the `split_reads.png` image actually shows a busy
frame with red-highlighted segments and purple insertion markers on several
reads — not empty, not obviously "absent chimeric alignments" to a glancing
reader — but qwen's "should show the absence of chimeric alignments" reads
as confirmed rather than predicted, without ever describing what's actually
in the frame.

**Same root cause as the claude-sonnet-5 session, confirmed the same way:**
`ollama_harness.py` constructs tool-result content as
`json.dumps(result["payload"], default=str)` — plain text, never image
data — so qwen never received image pixels through this pipeline either.
Additionally, and unlike claude-sonnet-5, qwen2.5:7b is not a vision-capable
model at all (`ollama show qwen2.5:7b` lists `completion` and `tools` under
Capabilities, no `vision`) — so even a harness fix that fed the PNG back as
an image content block would have nothing to show it to. For qwen2.5:7b
specifically, "did it examine the image" is not just unanswered by this
harness, it's structurally unanswerable by this model.
