# Stage 1 IGV Breakpoint Assistant

## Running the MCP server

From the repo root with rda environment active:
```
conda activate rda
python -m stage1_igv_assistant.server
```

## Running all tests
```
python stage1_igv_assistant/tests/test_bam_tools.py
```

## Tools available (8 total)
1. bam_stats_at_locus — quality check, always call first
2. discordant_pairs — inter-chromosomal translocation signal
3. soft_clipped_reads — breakpoint precision
4. split_reads — chimeric junction evidence
5. read_depth_profile — copy-number changes
6. breakpoint_evidence_summary — integrated 4-layer report
7. gene_at_locus — which gene is disrupted (Ensembl REST)
8. reciprocal_breakpoint — both sides of a balanced translocation

## Anti-hallucination design
The LLM receives only tool output. It cannot add genomic claims
from training data. Evidence must be stated with the tool that
produced it.

## Known limitations
- discordant_pairs: only valid for paired-end data (not PacBio HiFi)
- split_reads: requires modern aligner (BWA-MEM, minimap2). Zero SA
  tags in 2018-era BAMs means this tool cannot contribute.
- gene_at_locus: queries Ensembl REST API, requires internet,
  may be slow. Retry logic added (see bam_tools.py).
- Real balanced translocation BAM not yet found for validation.
  Demo used synthetic data. See DEMO_END_TO_END.md.
- run_igv_screenshot: IGV re-downloads genome annotation from igv.org 
  on every invocation, causing variable startup time and occasional 
  failures on repeated rapid calls. Observed non-deterministic 
  success/failure on the synthetic-BAM test case across repeated runs. 
  For batch use, consider pre-downloading the genome with IGV's 
  genome cache or calling with retries. Single interactive calls 
  are reliable.
- run_igv_screenshot requires a working X display. On WSL2 this is 
  provided by WSLg (DISPLAY=:0). Do not override DISPLAY — IGV's 
  AWT thread will crash before rendering.
