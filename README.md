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
HCC1143 (real short-read), and GIAB HG002 (real long-read). See
`stage1_igv_assistant/README.md` for tool details and
`stage1_igv_assistant/results/` for validation write-ups.

Stage 2 (variant/gene + phenotype prioritization) has not been started.