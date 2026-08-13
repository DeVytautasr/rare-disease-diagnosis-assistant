# Thesis Chapter: Background and Methodological Foundations

This directory contains the thesis background and methods chapter.

## Current version
- thesis_background_methods_chapter.md — full chapter, 47 references
- thesis_chapter_updated.docx — Word format for supervisor submission

The chapter's Stage 1 section now also covers the three-model comparison,
the two retracted findings, and the scoring-criteria reliability analysis —
see `stage1_igv_assistant/results/` for the underlying write-ups.

## Chapter sections
1. Introduction and Motivation
2. Related Work (includes Eilbeck 2017, AI-MARRVEL 2024, MARRVEL-MCP 2026)
3. Structural Variants and Chromothripsis
4. IGV and Breakpoint Inspection
5. Synthetic Cases and Benchmark Data
6. Prototype Design and Evaluation
7. Stage 1 Implementation and Validation
8. Scope and Limitations
9. Concluding Remarks
10. References (47 total)

## Stage 1 implementation summary (as of 9dce3bc — see stage1_igv_assistant/README.md for current counts)
- 11 tools in bam_tools.py, including the applicable-layers normalisation
  (`applicable_layers`) and the per-layer `evidence_panel` screenshot tool
  added 2026-08-11
- MCP server in server.py
- BamCase schema in case_object.py
- Test suite in stage1_igv_assistant/tests/, all passing (see
  stage1_igv_assistant/results/AUDIT_2026_08.md for a reproduced run)
- Validated: synthetic translocation (4/4 STRONG), HCC1143 2018 BAM
  (pipeline limitation documented), GIAB HG002 deletion on PacBio HiFi
  and Illumina 300x (depth threshold calibrated; normalised over
  applicable layers, e.g. 3/3 on Illumina 300x — evidence_score
  63.3/100, split_reads structurally excluded since novoalign emits no
  SA tags)
- Outstanding: real balanced translocation BAM not found. Real
  balanced translocation validation data (e.g. a known patient case)
  needed to complete this validation. `evidence_panel` (added
  2026-08-11) is now validated end-to-end across all three cases that
  matter — synthetic positive control, HCC1143 real negative control,
  and the real confirmed GIAB deletion — see
  `stage1_igv_assistant/results/EVIDENCE_PANEL_VALIDATION.md` §§1–3.
