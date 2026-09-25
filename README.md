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

**Local front end (`stage1_igv_assistant/ui.py`).** A dependency-free browser
interface on 127.0.0.1: the filter chain with every threshold's provenance,
the four evidence layers at both breakends of a candidate, hand-entered
coordinates, two-sample comparison, IGV panels, and a log of every tool call
behind every number on screen. Two panels state what the tool cannot do:
"What this tool cannot tell you" and "Where every number came from".

**Model chat panel (`stage1_igv_assistant/chat.py`).** An optional side panel
where a model calls the same tools through the same recorder. It is a second
consumer, never a second path to the data: it receives dataset labels rather
than file paths, and any number in its prose that no tool returned is marked
on screen. Backends: a local model via ollama, or the Anthropic API.

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
produced 0 malformed arguments in 313 calls. No local model is usable
unsupervised on that hardware. Findings are in
`results/BENCHMARK_LOCAL_MODELS.md` and `results/BENCHMARK_CLAUDE_BASELINE.md`,
both of which open with correction notices — two published findings turned out
to be measurement artifacts rather than model behaviour, and the documents say
so before they say anything else. `results/README.md` indexes all of it with a
reading order. Run logs from the earlier stages are grouped under
`benchmark/runs/` (see its README).

**The result worth knowing.** A balanced translocation can never score
"strong" here: depth correctly contributes nothing when no DNA is gained or
lost, and the paired-read layer is capped because roughly half the reads at
the breakpoint come from the intact homolog. Initially no model at any tier
could state this — 1 of 20 runs — because the score at which "strong" begins
appeared in no tool return and no tool description. After adding nine fields
and one sentence to one tool's return, with no change to any model, 19 of 20
runs state it. What a model can reach determines what it can say.

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

Stage 2 (variant/gene + phenotype prioritization) has not been started.
