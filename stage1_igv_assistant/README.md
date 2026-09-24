# Stage 1 IGV Breakpoint Assistant

## Installing IGV

`igv_screenshot` and `evidence_panel` need the IGV desktop app, which is
not a conda package. Install it with:
```
bash scripts/install_igv.sh
```
This downloads IGV 2.17.4 to `~/IGV_2.17.4` (the default path
`run_igv_screenshot` auto-detects) and checks that `java` and `DISPLAY`
are set up correctly, warning clearly if not. Safe to re-run. To point
the tools at a different IGV install instead, set `IGV_PATH=/path/to/igv.sh`
(checked before the hardcoded fallback locations).

## Running

The browser front end is the normal way in — it starts both servers in-process:
```
python -m stage1_igv_assistant.ui --check   # what is available, then exit
python -m stage1_igv_assistant.ui           # http://127.0.0.1:8765
```
`--check` prints the capability banner: which of MINIMAL / FULL (IGV) /
COMPLETE (local model) this install can deliver, where the exclude template
and data directory resolved from, and whether an Anthropic key is present.
It then prints every condition behind MINIMAL and FULL as PASS or FAIL: MINIMAL
is exercised on a small generated fixture, not assumed, and FULL checks the
IGV prerequisites (executable igv.sh, a new-enough Java, a display) without a
test render. The exit code is 0 only if every MINIMAL condition passed.

Either MCP server can also be run on its own, for an MCP client:
```
python -m stage1_igv_assistant.server            # the 11 evidence tools
python -m stage1_igv_assistant.candidate_server  # the 4 candidate-set tools
```

## Running the tests

There are 23 test files in `tests/`, each runnable on its own:
```
python stage1_igv_assistant/tests/test_bam_tools.py    # the largest; needs BAMs, and Java for the IGV assertions
python stage1_igv_assistant/tests/test_server.py       # MCP layer, tool contract
python stage1_igv_assistant/tests/test_vcf_tools.py    # candidate sets, dedup orientation, comparison symmetry
python stage1_igv_assistant/tests/test_ceiling_echo.py # the attainable-ceiling echo; exit 2 = INCOMPLETE (no IMP01.bam)
python stage1_igv_assistant/tests/test_api_leak.py     # nothing path-like reaches an external API
```
The remaining eighteen are focused pure-Python regression suites (no BAM, no
IGV, each under three seconds), one per defect class found during development —
`test_partner_distribution.py`, `test_quality_gate.py`,
`test_subthreshold_observations.py`, `test_minimum_support.py`,
`test_contig_naming.py` and others. Each is named for the condition that
exposed the defect it guards.

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
   available (evidence_score_raw). Also returns the band ladder and the
   attainable-score analysis: `score_bands`, `strong_band`, `attainable_here`,
   `strong_band_reachable_here`, `max_all_layers`, `max_with_flat_depth`,
   `attainable_basis`, `attainable_note`, `attainable_ceiling_derivable`. All
   are derived from this module's own source by `score_tiers.py`, so a revised
   scoring ladder moves them with it and an unparseable one is reported rather
   than guessed. They exist because the score at which "strong" begins was
   previously in no return and no description, which made the ceiling argument
   for a balanced event unstatable at every model tier.
8. gene_at_locus — which gene is disrupted (Ensembl REST)
9. reciprocal_breakpoint — both sides of a balanced translocation
10. igv_screenshot — headless IGV batch mode, generates a single PNG with
    the requested coloring/window
11. evidence_panel — one PNG per informative evidence layer (discordant
    pairs, split reads, read depth, soft clips), each with the IGV
    settings that actually isolate that layer visually; skips and
    explains layers detect_applicable_layers finds inapplicable

## Candidate-set tools (4 total, `candidate_server.py`)
1. load_candidate_set — register a caller's VCF/BCF for the session; reports
   record counts, the caller's breakend convention (`chr2_pos2` or `mateid`),
   and how many records merged during deduplication
2. list_candidates — filtered list with every applied threshold echoed back,
   each carrying its provenance (`tool-defined`, `reference-defined` or
   `author judgement`), the cumulative survivors after that step, and what
   that filter alone would remove from the unfiltered set
3. get_candidate — one junction in full, shaped so the evidence tools can
   consume its two breakends directly
4. compare_candidate_sets — junction-level recurrence between two sets
   (±500 bp), for two-sample artifact screening or two-caller concordance

Deduplication keys on orientation as well as position. It did not originally,
so the two halves of a reciprocal junction — same coordinates, opposite
orientation — collapsed into one; keying on orientation took the public test
set from 840 to 894 junctions and its survivors from 24 to 27.

## Anti-hallucination design
The LLM receives only tool output. It cannot add genomic claims
from training data. Evidence must be stated with the tool that
produced it.

This extends to images. `igv_screenshot` and `evidence_panel` take no
output path from the caller and return none: the server assigns the
location and hands back an opaque `image_ref` plus region, coloring mode,
dimensions, and success/failure. Since the tool call carries text only, a
model never receives pixel data — and with no path either, it has nothing
to describe. Benchmarking found models describing image contents they had
not been shown, and an advisory "you have not seen this image" instruction
stopped one model but not another; removing the parameter was the only
constraint that did not depend on the model choosing to comply. Handles
resolve to real files via `manifest.json` in the server's image session
directory (`IGV_IMAGE_SESSION_DIR`, else `screenshots/sessions/`).

## Known limitations
- The MCP server must be registered with a `PATH` that includes the conda
  environment's `bin` directory. IGV requires `java`, which is only present
  in the `rda` environment. A server registered with an empty `env` will
  report screenshot failures with no clear cause — the tools work when
  called directly from an activated shell but not through MCP. The same
  applies to `IGV_PATH` if you're using it to point at a non-default IGV
  install: it must be set in whatever `env` the MCP server is registered
  with, not just in your interactive shell.
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
  Any pre-existing file at the output path is deleted before IGV
  launches: a stale file's size is trivially "stable" from the first
  poll, which previously caused the tool to kill IGV before it had
  rendered anything and report the leftover file as a fresh success.
