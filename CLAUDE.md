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

## Repository structure

Work is organized by stage, each in its own top-level folder (e.g.
`stage1_igv_assistant/` for the Stage 1 SV/breakpoint module). Within a
stage folder:

- `tools/` for implementation (`bam_tools.py` — the evidence tools;
  `vcf_tools.py` — candidate-set loading, filtering and comparison)
- `tests/` for unit tests — one file per defect class, each runnable alone
- `data/` for small synthetic or benchmark examples only (large sequencing
  files are gitignored, not committed)
- `results/` for session reports, validation write-ups, and audits
- `server.py` and `candidate_server.py` — the stage's two MCP entrypoints
  (11 evidence tools and 4 candidate-set tools)
- `ui.py` — the local browser front end, which runs both servers in-process
- `chat.py` — the optional model panel (ollama or the Anthropic API)
- `config.py` — path resolution: flag > env var > config file > default
- `score_tiers.py` — derives the scoring tiers and band boundaries from
  `bam_tools.py`'s own source, so displayed ceilings cannot go stale

`docs/` at the repo root holds cross-stage design notes, the thesis
chapter drafts (`docs/thesis/`), and the Lithuanian install and usage
documentation (`DIEGIMAS.md`, `NAUDOJIMAS.md`). At the repo root:
`install.sh`, `requirements.txt` and `sv-assistant.conf.example` are the
portable install path; `make_demo_bundle.py` regenerates the demo data.
There is no top-level `src/` or `notebooks/` — exploratory work and
implementation both live inside the relevant stage folder.

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