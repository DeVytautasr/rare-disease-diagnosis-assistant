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
python stage1_igv_assistant/tests/test_server.py
```
`test_bam_tools.py`'s real-IGV assertions and `igv_screenshot`'s
functionality both require Java (present in the `rda` conda env, not the
base env) — run under `conda run -n rda python ...` or
`conda activate rda` first to exercise those paths; both test files skip
or fail cleanly rather than hanging if Java/IGV aren't available.

## Tools available (11 total)
1. applicable_layers — samples the BAM to determine which evidence layers
   apply (pairing, SA-tag support); call once per BAM, first
2. bam_stats_at_locus — quality check, call first per locus
3. discordant_pairs — inter-chromosomal translocation signal
4. soft_clipped_reads — breakpoint precision
5. split_reads — chimeric junction evidence
6. read_depth_profile — copy-number changes
7. breakpoint_evidence_summary — integrated evidence report, normalised
   over applicable_layers (evidence_score) with the unnormalised sum also
   available (evidence_score_raw)
8. gene_at_locus — which gene is disrupted (Ensembl REST)
9. reciprocal_breakpoint — both sides of a balanced translocation
10. igv_screenshot — headless IGV batch mode, generates a single PNG with
    the requested coloring/window
11. evidence_panel — one PNG per informative evidence layer (discordant
    pairs, split reads, read depth, soft clips), each with the IGV
    settings that actually isolate that layer visually; skips and
    explains layers detect_applicable_layers finds inapplicable

## Anti-hallucination design
The LLM receives only tool output. It cannot add genomic claims
from training data. Evidence must be stated with the tool that
produced it.

## Known limitations
- The MCP server must be registered with a `PATH` that includes the conda
  environment's `bin` directory. IGV requires `java`, which is only present
  in the `rda` environment. A server registered with an empty `env` will
  report screenshot failures with no clear cause — the tools work when
  called directly from an activated shell but not through MCP.
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
