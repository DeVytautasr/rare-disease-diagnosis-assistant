# Project: Phenotype-driven Rare Disease Diagnosis Assistant

## High-level description

This repository contains the code for my MSc thesis in Systems Biology at Vilnius University.

The goal is to build an interpretable, phenotype-driven rare disease diagnosis assistant that takes:
- raw variant lists (SNVs/indels),
- structural variants (including complex rearrangements and chromothripsis-like events),
- transcriptome signals near breakpoints,
- and patient phenotypes encoded as HPO terms,

and produces a ranked, well-explained list of candidate genes and diagnoses.

This is a Master's thesis prototype, not a production system.

## Main modules

1. Variant/gene + phenotype module
- Combine SNV/indel evidence, gene annotations, inheritance hints, and HPO phenotypes.
- Prioritize candidate genes and inherited disease diagnoses.
- Keep scoring interpretable and decomposed into evidence components.

2. Structural variant / chromothripsis + transcriptome module
- Connect SV breakpoints to nearby or disrupted genes.
- Integrate transcriptome changes near breakpoints.
- Link structural events to phenotype and disease-gene knowledge.

## Project principles

- Interpretability first.
- Prefer simple and robust methods over complex black-box models.
- Keep everything within realistic MSc thesis scope.
- Use explicit evidence layers: variant-level, gene-level, phenotype-level, transcriptome-level, SV-level.
- Write clean, documented, reproducible code.

## Expected repository structure

- `src/` for implementation
- `docs/` for design notes and method descriptions
- `tests/` for unit tests
- `notebooks/` for exploratory work
- `data/` for small synthetic or benchmark examples only

## How Claude should help

- Help design data structures for a unified case object.
- Help implement interpretable scoring functions.
- Help write Python code for annotation, phenotype processing, SV interpretation, and integration.
- Warn when ideas drift beyond realistic thesis scope.
- Prefer step-by-step plans and clear code over overly clever solutions.

## Constraints

- Do not silently invent biological results or fake benchmark outcomes.
- Ask for clarification when input format or intended output is unclear.
- Prefer reusable simple pipelines over fragile complex architectures.