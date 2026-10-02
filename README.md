# Rare Disease Diagnosis Assistant

Interpretable phenotype-driven rare disease diagnosis assistant for an MSc thesis in Systems Biology.

## Aim

This project aims to prioritize candidate genes and inherited disease diagnoses by integrating:
- SNVs/indels,
- structural variants,
- transcriptome evidence near breakpoints,
- and HPO phenotypes.

## Modules

- Variant/gene + phenotype prioritization
- Structural variant / chromothripsis interpretation with transcriptome integration

## Setup

Two paths, depending on what you want.

**To run the breakpoint tool** (this is the supported path — no reference
genome, no samtools, no delly, no GPU):

```
bash install.sh                                     # venv + 3 Python packages
.venv/bin/python -m stage1_igv_assistant.ui --check # says what is available
.venv/bin/python -m stage1_igv_assistant.ui         # http://127.0.0.1:8765
```

`install.sh` needs Python 3.10+ and internet at install time only. Verified on
a clean environment against Python 3.14.4 with `pysam==0.24.0`,
`fastmcp==3.4.6`, `requests==2.34.2`. If `python3 -m venv` fails because
Debian ships `ensurepip` separately and you have no root, the script builds the
environment without pip and bootstraps it over the network instead.

Paths are configurable — copy `sv-assistant.conf.example` to
`sv-assistant.conf`, or set `SV_DATA_DIR` / `SV_EXCLUDE_TEMPLATE` / `IGV_PATH`,
or pass `--dataset LABEL=PATH` / `--candidates LABEL=PATH`. The startup banner
reports which of three tiers is available (MINIMAL, FULL with IGV, COMPLETE
with a local model) and what is missing.

Lithuanian install and usage documentation: `docs/DIEGIMAS.md`,
`docs/NAUDOJIMAS.md`.

**To reproduce the development environment** (adds htslib, IGV, the
benchmark tooling):

```
conda env create -f environment.yml -n rda
conda activate rda
bash scripts/install_igv.sh   # IGV 2.17.4, needed for the screenshot tools
```
`environment.yml` is version-pinned and solves on linux-64 and macOS
(Intel and Apple Silicon); see the comment at the top of the file for what
was pruned/adjusted to make that true, and why. It does not support
Windows (some bioconda tools aren't published there). To reproduce the
linux-64 development environment exactly, build hashes included, use
`environment-linux64-exact.yml` instead.

## Status

**Stage 1 (`stage1_igv_assistant/`) — SV/breakpoint inspection module.**
Two MCP servers exposing 15 tools: 11 evidence tools over BAM files
(`server.py` — discordant pairs, soft clips, split reads, depth, quality,
layer applicability, gene lookup, reciprocal check, integrated summary, two
IGV image tools) and 4 candidate-set tools (`candidate_server.py` — load a
caller's VCF/BCF, filter it down a reported chain, open one junction, compare
two sets). Validated on synthetic translocation data, HCC1143 (real
short-read, 2018 pipeline), and GIAB HG002 — cross-technology, on both real
PacBio HiFi (long-read) and real Illumina 300x (short-read) alignments of the
same confirmed deletion. See `stage1_igv_assistant/README.md` for tool
details, `stage1_igv_assistant/results/` for validation write-ups, and
`TUTORIAL.md` for a guided walkthrough (written for external reviewers).

**Local front end (`stage1_igv_assistant/ui.py`, page `ui_page.html`).** A
dependency-free browser interface on 127.0.0.1 in three steps. *Candidates*: the
filter chain, each threshold labelled with its provenance, the count left after
every step, and an optional comparison that drops junctions also found in another
sample. *Evidence*: the four layers at both ends of a junction, the combined score
against the score reachable at that position, hand-entered positions (marked as
such), and IGV images on request. *Ask the assistant*: below. The combined score
and every measurement link to the tool call that returned them (the "reachable
here" mark is worked out by the page from the scoring tiers and the measured
values; it is not a tool return), and a call log lists every call of the session;
"Limits" states what the tool cannot tell you. Each sample is marked *Test data*
(autodiscovered under the data directory, or declared in the config file's
`[test_data]`) or *Private data* (registered explicitly), and the page always opens
on test data. *Note, 2026-10-02:* the page was rebuilt in Phase 24; the previous
page is served at `/classic`.

**Assistant (`stage1_igv_assistant/chat.py`).** A model calls the same tools
through the same recorder. It is a second consumer, never a second path to the
data: it receives dataset labels rather than file paths, and any number in its
final answer that no tool returned is marked on screen (its intermediate text is
shown unchecked). The page offers the models the local Ollama serves and, when an
Anthropic API key is present, `claude-sonnet-5` and `claude-opus-5` (the two
evaluated). A cloud model receives private data only on the person's explicit
confirmation for that question: without it the model is offered the test-data
labels alone, any call naming a private label, a set loaded from a private file,
or a file by its path is refused before it runs, and a question that names a
private sample or contains a position read from private data is not sent at all
(`ui.private_refs_in_text`, `ui._private_hits`; tests:
`tests/test_interface_page.py`). The confirmation is cleared when the question is
sent; while it holds, the model is offered every private label, not only the ones
the question names.

**Controlled positive test.** Twelve heterozygous balanced translocations were
implanted into real NA12878 reads at positions known in advance. In the rebuild
of 2026-09-25, 16 of 24 junctions were recovered: 8 of 8 in clean unique
sequence, 8 of 8 next to repeats, 0 of 8 where the surrounding sequence maps
ambiguously. Every loss occurred at variant calling — delly reported nothing at
those junctions — and no filter setting tested discarded a single true junction.
Record: `stage1_igv_assistant/results/synthetic_control_2026-09/analysis_2026-09-25/`.

> *Correction, 2026-09-25.* This paragraph gave the first run's result: 14 of 24,
> with 6 of 8 next to repeats. That run's records were lost in the 2026-09-23
> reinstall and it could not be regenerated; the test was rebuilt as a new
> experiment (new random draws, an ALT-aware index, partly new breakpoints). The
> class that changed is repeat-adjacent. Its implant at the reused coordinates,
> IMP06, was missed then and is detected now; realigned without the `.alt` file,
> it loses both junctions again, so ALT-aware alignment rather than the random draw
> accounts for that difference (`analysis_2026-09-25/noalt/compare.json`).

**Model comparison.** Six adversarial cases, 5 runs each, same server, same
tools, same scoring. Local models on an 8 GB GPU (`qwen3.5:4b`, `qwen2.5:7b`,
`qwen3.5:9b`) produced malformed tool arguments in 13-26% of calls and did not
always finish; `claude-sonnet-5` and `claude-opus-5` through the same harness
produced 0 malformed arguments in 313 calls (*2026-09-26:* that record is lost; rerun on
reconstructed cases, 0 in 341 calls — `stage1_igv_assistant/benchmark/runs/phase9_rerun_2026-09-26/`;
`claude-opus-5` ran 15 runs, three or two per case). No local model is usable
unsupervised on that hardware. Findings are in
`results/BENCHMARK_LOCAL_MODELS.md` and `results/BENCHMARK_CLAUDE_BASELINE.md`,
both of which open with correction notices — two published findings turned out
to be measurement artifacts rather than model behaviour, and the documents say
so before they say anything else. `results/README.md` indexes all of it with a
reading order. Run logs from the earlier stages are grouped under
`benchmark/runs/` (see its README).

**The result worth knowing.** With all four measurements counted, at the default
window, a balanced translocation cannot score "strong" here unless the depth
layer errs: depth correctly contributes nothing when no DNA is gained or lost,
and the paired-read layer is capped because roughly half the reads at the
breakpoint come from the intact homolog. At all 32 breakends of the 16 detected
junctions the attainable ceiling is 57.5, below the 70 at which "strong" begins.
Scored blind, with four models and five runs per condition: while no tool
return carried the ceiling, 0 of 20 runs stated it; with nine fields added to
one tool's return and no change to any model, 12 of 20 did, all with the full
argument — `claude-sonnet-5` 5 of 5, `claude-opus-5` 5 of 5, `qwen3.5:4b` 2 of
5, `qwen2.5:7b` 0 of 5. In the pre-registered extension to 15 runs,
`qwen3.5:4b` stated it in 5 of 15 (0 of 15 without; Fisher's exact test,
two-sided, p = 0.042). With the fields, `claude-sonnet-5` also called the
evidence strong in 5 of 5 runs (in none without), and `claude-opus-5` in 2 of 5
(2 of 5 without). Exposing a quantity in a tool return is necessary for any
model to reason from it and sufficient only for capable ones. Record:
`stage1_igv_assistant/benchmark/runs/phase10_rerun_2026-09-25/blind/blind_scores.json`.

> *Correction, 2026-09-25.* The run counts above (1 of 20 before, 19 of 20 after)
> come from Phase 9 and 10 records that were lost in the 2026-09-23 reinstall and
> could not be regenerated. The arithmetic stands, re-measured on the rebuilt
> implants: at all 32 breakends of the 16 detected junctions the tool's own
> attainable ceiling is 57.5, below the 70 at which "strong" begins. The model
> comparison was rerun with the nine fields as the only difference between
> conditions — stripped in the harness, proven at four loci, everything else
> identical. (The original comparison set runs before the commit against runs after
> it, and that commit added two other fields to the same return; measured here, the
> nine fields alone also push the return across the 2,600-character budget at which
> the model loop truncates it, which changes what else the model sees — the rerun
> holds that fixed.) Only the two local models could run; the stored
> API key was rejected. Without the fields 0 of 10 runs stated the ceiling; with
> them 2 of 10 did, and gave the full argument — both from `qwen3.5:4b` (2 of 5;
> the lost record said 5 of 5). `qwen2.5:7b` stated it in none of its 5 runs,
> including the 4 in which the fields reached it. Record:
> `stage1_igv_assistant/benchmark/runs/phase10_rerun_2026-09-25/` (`tally.json`,
> `scores.json` with every deciding sentence quoted).

> *Correction, 2026-09-26.* The API arm has now run, and every run of the experiment was scored blind: a fresh scorer saw only the prompt and each final answer, under random IDs, and every sentence it quoted was checked against the traces. Over the same twenty runs (four models, five each): without the fields 0 of 20 state the ceiling; with them 12 of 20 state it, and all 12 give the full argument — `claude-sonnet-5` 5 of 5, `claude-opus-5` 5 of 5, `qwen3.5:4b` 2 of 5, `qwen2.5:7b` 0 of 5. Extended to 15 runs per cell as pre-registered: `qwen3.5:4b` states it in 5 of 15 and completes it in 4 of 15 (0 and 0 of 15 without; Fisher's exact test, two-sided, p = 0.042 and 0.100), `qwen2.5:7b` in 3 and 0 of 15 (p = 0.224 and 1.0). One further result, not a registered test: with the fields, `claude-sonnet-5` also called the evidence strong in all five runs (without them, in none) — it took the ceiling as a reason to re-rate the evidence above its band. Record: `stage1_igv_assistant/benchmark/runs/phase10_rerun_2026-09-25/blind/blind_scores.json`.

Stage 2 (variant/gene + phenotype prioritization) has not been started.
