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

## Status

**Stage 1 (`stage1_igv_assistant/`) — SV/breakpoint inspection module.**
An MCP server exposing 11 tools for interpretable breakpoint evidence
(discordant pairs, soft clips, split reads, depth, gene lookup, IGV
screenshots) over BAM files. Validated on synthetic translocation data,
HCC1143 (real short-read), and GIAB HG002 (real long-read). See
`stage1_igv_assistant/README.md` for tool details and
`stage1_igv_assistant/results/` for validation write-ups.

Stage 2 (variant/gene + phenotype prioritization) has not been started.