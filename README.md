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

```
conda env create -f environment.yml -n rda
conda activate rda
bash scripts/install_igv.sh   # installs IGV 2.17.4, needed for screenshot tools
```
`environment.yml` is version-pinned and solves on linux-64 and macOS
(Intel and Apple Silicon); see the comment at the top of the file for what
was pruned/adjusted to make that true, and why. It does not support
Windows (some bioconda tools aren't published there). To reproduce the
linux-64 development environment exactly, build hashes included, use
`environment-linux64-exact.yml` instead.

## Status

**Stage 1 (`stage1_igv_assistant/`) — SV/breakpoint inspection module.**
An MCP server exposing 11 tools for interpretable breakpoint evidence
(discordant pairs, soft clips, split reads, depth, gene lookup, IGV
screenshots) over BAM files. Validated on synthetic translocation data,
HCC1143 (real short-read, 2018 pipeline), and GIAB HG002 — cross-technology,
on both real PacBio HiFi (long-read) and real Illumina 300x (short-read)
alignments of the same confirmed deletion. See
`stage1_igv_assistant/README.md` for tool details,
`stage1_igv_assistant/results/` for validation write-ups, and
`TUTORIAL.md` for a guided walkthrough (written for external reviewers).

**Model comparison benchmark (`stage1_igv_assistant/benchmark/`).**
Three models run against the same server, cases, and scoring:
`claude-sonnet-5` via an API harness, `qwen2.5:7b` and `llama3.1:8b`
hosted locally. Three cases — a confirmed deletion, a control locus, and
an adversarial variant whose prompt asserts a translocation the data does
not support. Run logs are grouped by fix stage under `benchmark/runs/`
(see its README); findings are in `results/BENCHMARK_LOCAL_MODELS.md` and
`results/BENCHMARK_CLAUDE_BASELINE.md`, both of which open with
correction notices — two published findings turned out to be measurement
artifacts rather than model behaviour, and the documents say so before
they say anything else. `results/README.md` indexes all of it with a
reading order.

Stage 2 (variant/gene + phenotype prioritization) has not been started.