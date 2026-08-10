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

## Tools available (9 total)
1. bam_stats_at_locus — quality check, always call first
2. discordant_pairs — inter-chromosomal translocation signal
3. soft_clipped_reads — breakpoint precision
4. split_reads — chimeric junction evidence
5. read_depth_profile — copy-number changes
6. breakpoint_evidence_summary — integrated 4-layer report
7. gene_at_locus — which gene is disrupted (Ensembl REST)
8. reciprocal_breakpoint — both sides of a balanced translocation
9. igv_screenshot — headless IGV batch mode, generates PNG visual 
   evidence of the breakpoint region

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
- run_igv_screenshot: after the batch script's `snapshot` command 
  writes the PNG, IGV's JVM/AWT thread has been observed to hang on 
  the following `exit` command instead of terminating (not yet 
  root-caused; may be related to the genome re-download timing 
  above). The tool handles this itself: it polls for the output file 
  and, once its size is stable across two checks 1s apart, terminates 
  IGV directly (SIGTERM, then SIGKILL after 5s if needed) rather than 
  waiting for `exit` to work. Because igv.sh runs `java` as its last 
  command without `exec`, the shell stays alive as java's parent — the 
  tool signals the whole process group (`start_new_session=True` + 
  `killpg`), not just that wrapper shell, so the actual IGV/java GUI 
  process is reliably killed too. No manual window closing is needed. 
  The returned dict's `shutdown_method` field records what happened: 
  `"clean_exit"` (IGV exited on its own), `"terminated_after_snapshot"` 
  (the tool had to kill it after confirming the PNG was written), or 
  `"timeout"` (no output file appeared within `timeout_sec`, which 
  still acts as a hard ceiling and kills the process either way).
