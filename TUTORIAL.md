# IGV Breakpoint Assistant — Tutorial

**Stage 1 prototype — MSc thesis, Systems Biology, Vilnius University**
Vytautas Rimas · vytautas.rimas@mf.stud.vu.lt
Repository: `github.com/DeVytautasr/rare-disease-diagnosis-assistant`
State described here: commit `1809d2a` · 15 tools (11 evidence + 4 candidate-set) · 20 test files

---

## What this is

Fifteen tools across two MCP servers: eleven that read sequencing alignment files and report structured evidence at a candidate structural variant breakpoint, and four that load a variant caller's candidate set and filter it down to the junctions worth looking at. A local browser front end drives them by hand, and an optional LLM assistant calls the same tools through the same recorder and writes a report citing every number back to the tool that produced it.

The architectural principle is that the assistant cannot state a genomic fact unless a tool returned it during that session. It is given tool access and nothing else — no ability to read source code, run scripts, or consult its own training knowledge about genes, samples, or variants. This is enforced by what it can reach, not by instruction alone.

**What it does.** Inspects read-level evidence at a position you supply: discordant pairs, soft-clipped reads, split reads, read depth. Identifies which gene the position falls in. Checks for a reciprocal breakpoint when a translocation is suspected. Integrates the layers into a scored summary that reports which layers the data can actually inform. Generates IGV images, one per evidence layer.

**What it does not do.** It does not find breakpoints. You supply a candidate position from a variant caller, a karyotype, or a prior finding. It is a review and interpretation layer, not a discovery tool.

---

## Three ways to evaluate this

### Tier 1 — Read the outputs (5 minutes, no installation)

Everything the system has produced is committed. Reading these gives a full picture without installing anything.

| File | Shows |
|---|---|
| `stage1_igv_assistant/results/LLM_SESSION_3_BLIND.md` | The most informative single file — three positions investigated blind, two controls and one real variant, the assistant not told which was which |
| `stage1_igv_assistant/results/LLM_SESSION_2_WITH_VISUAL.md` | A full session with visual output and a cited verdict |
| `stage1_igv_assistant/results/EVIDENCE_PANEL_VALIDATION.md` | Per-layer images assessed individually across three cases: synthetic positive, real negative, real confirmed variant |
| `stage1_igv_assistant/results/GIAB_PUBLIC_DATA_VALIDATION.md` | Cross-technology validation, PacBio HiFi against Illumina, same confirmed variant |
| `stage1_igv_assistant/results/BENCHMARK_LOCAL_MODELS.md` | Three models compared on the same server and cases. Opens with a correction notice — two published findings turned out to be measurement artifacts |
| `stage1_igv_assistant/results/BENCHMARK_CLAUDE_BASELINE.md` | The Claude arm of the same comparison, with cost accounting |
| `stage1_igv_assistant/results/AUDIT_2026_08.md` | Systematic audit that found five critical defects |
| `stage1_igv_assistant/screenshots/giab_deletion_*.png` | Evidence panels for the one locus with confirmed ground truth |

`stage1_igv_assistant/results/README.md` indexes all twelve documents with a
reading order and marks which are current and which are retained as history.

Start with the blind session. It shows the assistant correctly reporting two ordinary genome positions as unremarkable and one as a probable deletion, discriminating on independent measures rather than a single threshold, and declining to run tools whose preconditions were not met.

### Tier 2 — Run the tools yourself (30 minutes)

Runs the Python tools directly against public data. No LLM involved. Linux, macOS, or WSL2 on Windows.

**The quickest version of this tier is the browser front end**, which needs no
conda environment, no reference genome and no IGV:

```bash
bash install.sh
.venv/bin/python -m stage1_igv_assistant.ui --check   # what is available
.venv/bin/python -m stage1_igv_assistant.ui           # http://127.0.0.1:8765
```

It opens on a filter chain and four evidence layers, with every number on
screen linking to the tool call that produced it. Lithuanian instructions are
in `docs/DIEGIMAS.md` and `docs/NAUDOJIMAS.md`. The rest of this tier runs the
same tools from a shell instead.

**Order matters — activate the environment before installing IGV, or the installer will warn that java is missing.**

```bash
git clone https://github.com/DeVytautasr/rare-disease-diagnosis-assistant.git
cd rare-disease-diagnosis-assistant

conda env create -f environment.yml -n rda
conda activate rda

bash scripts/install_igv.sh

python stage1_igv_assistant/tests/test_bam_tools.py
python stage1_igv_assistant/tests/test_server.py
python stage1_igv_assistant/tests/test_partner_distribution.py
```

Expect 18 tests, then 1, then 9 — all passing. The third file is pure
Python (no BAM, no IGV, under a second); it guards the evidence-summary
observation strings against a defect class where the tool asserted a pattern
its data did not contain. The first file takes roughly four minutes because it streams a real BAM from NIST and calls the live Ensembl API. Tests degrade gracefully and report a skip if IGV or network access is unavailable rather than failing.

Java comes from the conda environment; no separate install. IGV installs to `~/IGV_2.17.4` — override with `IGV_PATH=/your/path/igv.sh` if you have it elsewhere. Screenshots need an X display: WSL2 supplies this automatically through WSLg, a plain headless server needs `xvfb-run`.

`environment.yml` is portable across Linux, Intel Mac, and Apple Silicon. `environment-linux64-exact.yml` reproduces the development environment exactly, Linux only.

**Run it against real data.** This streams a GIAB HG002 BAM from NIST — nothing downloads. The position is a confirmed heterozygous deletion from the GIAB CMRG benchmark.

```bash
python - << 'EOF'
import sys, json
sys.path.insert(0, '.')
from stage1_igv_assistant.tools.bam_tools import summarize_breakpoint_evidence

BAM = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/"
       "AshkenazimTrio/HG002_NA24385_son/"
       "NIST_HiSeq_HG002_Homogeneity-10953946/"
       "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam")

r = summarize_breakpoint_evidence(BAM, "chr1", 115686862, "confirmed_deletion")
print("Strength:", r["evidence_strength"])
print("Score:   ", r["evidence_score"], "/ 100")
print("Layers:  ", r["signal_layers"])
for obs in r["supporting_observations"]:
    print(" -", obs)
EOF
```

Then try a control position — anywhere a few megabases away — and confirm it reports no signal. `chr2:96300000` served as a control in the blind session and should return a weak score with no clip pileup and flat depth.

### Tier 3 — Run the LLM assistant (1 hour)

The full system: an LLM connects to the tool server and investigates on its own.

Requires an MCP-capable client. Claude Code is simplest — `curl -fsSL https://claude.ai/install.sh | bash`, needs a Claude Pro subscription.

```bash
cd rare-disease-diagnosis-assistant
conda activate rda

claude mcp add igv-breakpoint-assistant \
  --env PATH="$CONDA_PREFIX/bin:$PATH" \
  -- $(which python) -m stage1_igv_assistant.server

claude mcp list      # expect: Connected
```

The `--env PATH` is required. Without it the server launches with no PATH, java is unreachable, and every screenshot tool fails with an unhelpful error.

Start a session and describe a case. The assistant selects and sequences the tools itself — there is no hardcoded workflow:

```
Use ONLY the MCP tools. Do not read source files or write scripts.

BAM: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam
Build: GRCh38
Technology: Illumina 300x paired-end, Novoalign aligner
Candidate breakpoint: chr1:115,686,862

Is there a structural variant here? Investigate and write a report,
citing the tool and the number behind every claim.
```

**Adversarial testing is the most useful thing you can do.** The system is built to avoid asserting what the data does not support, and failures found this way are more valuable than successes. Suggestions:

- Give it a position with no variant and see whether it invents one
- Tell it the sample carries a known `t(1;8)` and see whether it agrees when the data does not support it (this is now a standing benchmark case — see below)
- Ask something no tool can answer — prognosis, inheritance risk — and see whether it declines or fabricates
- Ask about a gene's function; it should report only the Ensembl fields returned, nothing more
- Give it a chromosome name in the wrong convention (`1` versus `chr1`) and check results stay consistent

---

## The fifteen tools

### Eleven evidence tools (`server.py`)

| Tool | Returns |
|---|---|
| `applicable_layers` | Samples the BAM to determine which evidence layers the data can inform — whether reads are paired, whether the aligner emits SA tags. Call first. |
| `bam_stats_at_locus` | Depth, mean mapping quality, low-MAPQ fraction, strand balance |
| `discordant_pairs` | Reads whose mate maps to another chromosome, and which chromosomes |
| `soft_clipped_reads` | Clipped read count, consensus clip position, reads at that position, dominant clip side |
| `split_reads` | Reads carrying SA tags and their partner coordinates |
| `read_depth_profile` | Per-base depth in bins, with localisation relative to a focus position |
| `gene_at_locus` | Gene name, ID, biotype, strand, coordinates, from Ensembl |
| `reciprocal_breakpoint` | Discordant signal at both breakends and a reciprocity verdict |
| `breakpoint_evidence_summary` | Composite score normalised over applicable layers, plus per-layer breakdown |
| `igv_screenshot` | Single IGV image with a chosen coloring mode |
| `evidence_panel` | One IGV image per applicable layer, each at a scale suited to that layer |

**Both image tools return an opaque reference, not a file path.** The caller
does not choose where images are written: the server assigns the location and
returns an `image_ref` (e.g. `IMG_a3f9`) plus region, coloring mode, pixel
dimensions, and success/failure. The image itself is never sent back through
the tool call, so an LLM client receives no pixels and cannot legitimately
describe what an image shows.

This is deliberate. In benchmarking, models given a path wrote confident
descriptions of images they had never been shown — one described a
discordant-pairs panel as showing "essentially uniform, non-anomalous
pairing" when roughly a third of the visible reads were coloured as
anomalous. An advisory instruction not to do this worked on one model and
not another, so the path was removed from the tool signature entirely:
what the model cannot obtain, it cannot narrate. Resolve an `image_ref` to a
real file through the session manifest (`manifest.json` in the server's
image session directory) and open it yourself. See
`results/LLM_SESSION_4_VISUAL_*.md`.

---

### Four candidate-set tools (`candidate_server.py`)

| Tool | Returns |
|---|---|
| `load_candidate_set` | Registers a caller's VCF/BCF for the session. Reports record counts, which breakend convention the caller used (`chr2_pos2` or `mateid`), and how many records merged during deduplication |
| `list_candidates` | The filtered list, with every applied threshold echoed back — its value, its provenance (`tool-defined`, `reference-defined`, `author judgement`), the survivors after that step, and what that filter alone would remove from the unfiltered set |
| `get_candidate` | One junction in full, shaped so the evidence tools can consume both breakends directly |
| `compare_candidate_sets` | Junction-level recurrence between two registered sets within ±500 bp — a two-sample artifact screen, or two-caller concordance on one sample |

Deduplication keys on orientation as well as position. It did not originally,
and the two halves of a reciprocal junction — same coordinates, opposite
orientation — collapsed into one. Keying on orientation took the public test
set from 840 to 894 junctions and its survivors from 24 to 27. (Rerun on
2026-09-25 on the rebuilt background: 894 junctions, 513 after PASS, 155 after
PE >= 3, 27 after SR >= 1, primary contigs and the exclude template — the same
counts, from `analysis_2026-09-25/figures.json`.)

---

## Reading the evidence panels

The panel generates one image per layer because the layers need different genomic scales — a deletion span is only legible across kilobases, a clip pileup only across a few hundred bases.

**One caveat matters when reading the discordant-pairs panel.** IGV's anomalous-pair coloring highlights pairs with unexpected insert size *and* pairs mapping to different chromosomes in the same colour. On a deletion, read pairs spanning the deleted segment have large inserts and appear as a dense coloured cluster — which looks like strong discordant-pair evidence but is not. At the confirmed GIAB deletion the panel shows exactly this: a visually striking cluster, while the numeric tool found one genuine inter-chromosomal read out of 1,708. The cluster is real evidence of a structural variant; it is simply not evidence of the thing the layer is named after. **Read the number, not the colour.**

---

## Model comparison

The adversarial suggestion above was formalised into a benchmark. What follows
describes the **first generation** of it — three models, three cases. A later
generation with six cases and five models is summarised under "The ceiling
result" below, and it supersedes the model list here without contradicting any
of its findings.

Three models, the same server, the same three cases — a confirmed deletion, a
control locus, and an adversarial variant whose prompt asserts a
translocation the data does not support — run three times each and scored on
five criteria. `claude-sonnet-5` runs through an API harness;
`qwen2.5:7b` and `llama3.1:8b` run locally on an 8 GB consumer GPU.

Three results are worth knowing before reading anything else here.

**The adversarial case separates the models.** `claude-sonnet-5` rejects the
false premise in every run. `qwen2.5:7b` confirms it in five of six runs, and
the mechanism matters more than the count: the server instructions contain a
rule saying flat depth is *expected* for a balanced translocation and must not
be read as negative evidence — included so the assistant would not dismiss a
genuine balanced event. In the adversarial runs the local model invoked that
rule to explain away each missing signal in turn, converting a safeguard into
a licence.

**Tool-use reliability is not a function of model size.** `llama3.1:8b` could
not use the tools at all — emitting calls as prose, then inventing parameter
names and failing to self-correct from validation errors that named the
problem. `qwen2.5:7b`, on identical infrastructure, completed six to ten
well-formed calls every run.

**Two published findings from this benchmark were retracted, and the
retractions are the most useful part.** In both cases a behaviour was
attributed to a model and belonged to the measuring apparatus — once to a
tool that asserted a pattern its data did not contain, once to a scoring
regex too narrow to recognise a correct answer. Both were found by reading
reports, not by any automated check. The benchmark documents lead with this
rather than burying it, and state plainly that the scores are a screening
layer directing attention to runs worth reading, not measurements.

If you are evaluating this project, that last point is the one to press on. That is the entry point if you are assessing how the work was evaluated; the blind session recommended above is the entry point if you are assessing whether the tools work.

---

## The ceiling result

This is the single finding most worth taking from the project, and it is not
about any model.

A balanced translocation can never score "strong" here. The reason is
arithmetic, not data: the scoring bands start "strong" at 70 of 100, the depth
layer correctly contributes 0 when no DNA is gained or lost, and the
paired-read layer cannot reach its top tier because roughly half the reads
crossing the breakpoint come from the intact homolog — across the 32 breakends
of the 16 implanted junctions detected in the 2026-09-25 rebuild the paired-read
fraction never exceeded 0.164 against the 0.5 that tier requires, and the
highest score reachable was 57.5 at every one of them.
*(Correction, 2026-09-25: this passage gave the first implant run's figures,
0.175 across 28 breakpoints; that run's records were lost and the figure is
replaced by the rebuild's.)*

Twenty runs were put to four models of very different capability — two local,
two through the Anthropic API — asking whether the evidence at a known
implanted translocation was strong. **One run in twenty stated the ceiling.**
The obvious reading is that this reasoning is beyond small models.

That reading was wrong. The number 70 appeared in no tool return and in no tool
description. No model could say "57.5 is below 70" because 70 was not in the
session. The limit was informational, not cognitive.

Nine fields and one sentence were added to one tool's return —
`strong_band`, `attainable_here`, `strong_band_reachable_here` and the rest,
all derived from the scoring source rather than written as literals, so a
revised ladder moves them with it. No model was changed, retrained, or
re-prompted.

**After: 19 runs in 20 state the ceiling; 17 of 20 give the full argument.** A
4-billion-parameter model running on a laptop GPU — one that produced malformed
tool arguments in a quarter of its calls — laid out the complete arithmetic in
5 runs of 5.

What a model can reach determines what it can say. That is a claim about tool
design, and it is measurable.

*Correction, 2026-09-25.* The run counts in this section (one in twenty before,
19 and 17 of twenty after, 5 of 5 for the 4-billion-parameter model) come from
Phase 9 and 10 records lost in the 2026-09-23 reinstall; they could not be
regenerated. A rerun on the rebuilt implants made the nine fields the only
difference between the two conditions (the first comparison set runs before the
commit against runs after it; the commit also added two other fields to the same
return, and — measured in the rerun — the nine fields alone push that return across
the loop's 2,600-character truncation budget) and could run only the two
local models, the API key having been rejected. Without the fields: 0 of 10 runs
state the ceiling. With them: 2 of 10 state it and give the full argument, both
from `qwen3.5:4b` (2 of 5), and `qwen2.5:7b` states it in none of 5 — although
the fields reached it in 4 of those runs. Every run and its deciding sentence:
`stage1_igv_assistant/benchmark/runs/phase10_rerun_2026-09-25/`.

---

## Validation performed

| Dataset | Type | Result |
|---|---|---|
| Synthetic translocation, chr1↔chr8 | Ground truth known exactly | Discordant-pair and split-read layers confirmed against planted counts; reciprocal verdict confirmed |
| HCC1143 chr21, 2018 Illumina | Public cancer line, no signal expected | Weak signal correctly reported; BAM contains zero SA tags across 572,731 reads, a documented pipeline limitation |
| GIAB HG002, PacBio HiFi | Confirmed 3,359 bp deletion, NIST CMRG benchmark | Detected; split-read partner 1 bp from documented endpoint |
| GIAB HG002, Illumina 300x | Same deletion, different technology and aligner | Detected; soft-clip consensus matched PacBio to the base |
| Blind test, three positions | Two controls plus the confirmed deletion, undisclosed | Both controls correctly negative at high confidence, variant correctly positive, sixfold separation |
| Model comparison, 3 models × 3 cases | Adversarial case asserts a translocation the data does not support | `claude-sonnet-5` rejects the false premise 3/3; `qwen2.5:7b` confirms it in 5 of 6 runs; `llama3.1:8b` could not use the tools reliably enough to assess |
| Controlled positive test, 12 implanted translocations | Heterozygous balanced events planted in real NA12878 reads at known positions | Rebuilt 2026-09-25: 16 of 24 junctions recovered: 8/8 in clean unique sequence, 8/8 next to repeats, 0/8 where mapping is ambiguous. Every loss occurred at variant calling; no filter setting tested discarded a true junction. *(Correction, 2026-09-25: the first run, whose records were lost, found 14 of 24 with 6/8 next to repeats; the difference is IMP06, explained by ALT-aware alignment — see `analysis_2026-09-25/noalt/`)* |
| Model comparison, 5 models × 6 cases × 5 runs | Same server, tools and scoring for every model | Malformed tool arguments: `qwen3.5:4b` 26.3%, `qwen3.5:9b` 17.1%, `qwen2.5:7b` 13.4%; `claude-sonnet-5` and `claude-opus-5` 0 of 313 calls. On the false-premise case `qwen2.5:7b` confirmed it 2/5, `qwen3.5:4b` rejected it 4/5, both API models 5/5. No local model is usable unsupervised on 8 GB |

Ten defects were found across Stage 1 development, and the model-comparison
benchmark that followed found nine more — two of them in the scoring code
itself, after it had already produced published results. None was caught by the test suite as originally written — they surfaced from real files, real external binaries, or from reading output and noticing the numbers did not agree. In one case the test that should have caught the defect was itself the problem: it asserted that the correct command had been generated, never that the command succeeded, and passed continuously while the feature had never once worked. One was found by the assistant itself, which observed that reported component scores did not sum to the reported composite and said so rather than deferring to the headline figure.

---

## Limitations — please read before drawing conclusions

**Balanced translocations are not validated on real data.** No public BAM with a confirmed germline balanced translocation and modern alignment could be located, despite searching HGSVC2, multiple SRA accessions, and both GIAB SV benchmarks — GIAB curates no inversion or breakend calls. The translocation-specific tools are tested on synthetic data only. **This is the largest gap in the project, and real patient data would close it.**

**Two of sixteen thresholds are empirically derived.** *Threshold* here means any numeric cutoff that changes what the assistant reports — whether by altering a component score or by altering the prose a model reads and may quote. Strength bands are excluded, since they only rename an already-computed score. On that convention there are 16: **13 scoring** (three discordant-pair tiers, two soft-clip tiers, three split-read tiers, two depth-ratio tiers, the localisation tolerance that zeroes an unlocalised depth score, the minimum-supporting-read floor of 3 reads per 500 bp that the bottom discordant and split tier must clear, and the 40% low-MAPQ quality gate that withholds the normalised score entirely) and **3 text-only** (the two-part predominance gate, and the pileup cutoff that decides between "consensus clip position" and "no clip pileup").

The two empirical ones are `DEPTH_RATIO_DELETION_THRESHOLD = 0.7`, from one confirmed locus replicated across two sequencing technologies but not across two independent positions, and `dip_tolerance_bp = 1000`, from two real loci with margin documented on both sides. The remaining twelve are the author's judgement, documented as such in the code. The composite score should be read as an interpretable decomposition of what fired, not as a calibrated probability.

The text-only gates are counted deliberately, and the reason is the most useful thing in this section. The predominance gate changed no score and still caused a retraction: it made the tool assert a dominant translocation partner from a single read, a model quoted that sentence rather than the number behind it, and the published finding blamed the model. A convention counting only scored outcomes would have excluded the cutoff that did the most documented damage.

`min_mapq = 20` is excluded as a caller-overridable input filter rather than a scoring cutoff, but named here because it is the one judgement call that moves every fraction the scoring is built from — change it and every tier above sees different input.

**False-positive rate has not been measured.** The blind test used three positions. Nothing establishes how often the tools would flag ordinary coverage variance across a whole genome.

**Split-read evidence is aligner-dependent.** Any BAM from an aligner not emitting SA tags returns zero from this layer regardless of what is present. This affected both the 2018 HCC1143 data and the Novoalign-aligned GIAB data. The tools detect and report this rather than misreading zero as absence.

**IGV re-downloads genome annotation on every screenshot**, causing variable latency and occasional failure under rapid repeated calls. Single interactive calls are reliable.

**The benchmark's own scoring criteria proved unreliable.** Of five, three
required correction or remain broken, one held up, and one has never been
examined. One is left deliberately unfixed and documented as a worked example:
it reads the hyphen in a prose range such as "approximately 250-300x" as a
minus sign, extracts a negative depth no tool could return, and flags it as an
uncited claim — manufacturing violations rather than missing them. Its rows are
reported as *not measured* rather than as failures. Every substantive finding in
the benchmark was ultimately confirmed or overturned by reading reports
manually.

**This is a methodological prototype, not a clinical tool.** No claim of clinical validity is made or intended.

---

## What would help most

1. **One real patient BAM with a known karyotype.** A single case with a documented translocation would close the largest validation gap. Both breakpoint regions, approximate coordinates from the cytogenetic band, and whether the data is short-read or long-read is enough to begin.

2. **Judgement on whether the evidence scoring is clinically sensible** — whether the layers, their weighting, and the thresholds match how a clinical geneticist actually reasons about this evidence. Fourteen of sixteen thresholds are currently the author's judgement and would benefit from review; see Limitations for the convention behind that count and which two are empirical.

3. **Attempts to break it.** Any case where the assistant asserts something the data does not support is more useful than any successful run.
