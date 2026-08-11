# IGV Breakpoint Assistant — Tutorial

**Stage 1 prototype — MSc thesis, Systems Biology, Vilnius University**
Vytautas Rimas · vytautas.rimas@mf.stud.vu.lt
Repository: `github.com/DeVytautasr/rare-disease-diagnosis-assistant`
State described here: commit `81f24ef` · 11 tools · 19 tests

---

## What this is

Eleven tools that read sequencing alignment files and report structured evidence at a candidate structural variant breakpoint, plus an LLM assistant that calls those tools and writes a report citing every number back to the tool that produced it.

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
| `stage1_igv_assistant/results/REAL_DATA_VALIDATION.md` | Cross-technology validation, PacBio HiFi against Illumina, same confirmed variant |
| `stage1_igv_assistant/results/AUDIT_2026_08.md` | Systematic audit that found five critical defects |
| `stage1_igv_assistant/screenshots/giab_deletion_*.png` | Evidence panels for the one locus with confirmed ground truth |

Start with the blind session. It shows the assistant correctly reporting two ordinary genome positions as unremarkable and one as a probable deletion, discriminating on independent measures rather than a single threshold, and declining to run tools whose preconditions were not met.

### Tier 2 — Run the tools yourself (30 minutes)

Runs the Python tools directly against public data. No LLM involved. Linux, macOS, or WSL2 on Windows.

**Order matters — activate the environment before installing IGV, or the installer will warn that java is missing.**

```bash
git clone https://github.com/DeVytautasr/rare-disease-diagnosis-assistant.git
cd rare-disease-diagnosis-assistant

conda env create -f environment.yml -n rda
conda activate rda

bash scripts/install_igv.sh

python stage1_igv_assistant/tests/test_bam_tools.py
python stage1_igv_assistant/tests/test_server.py
```

Expect 18 tests then 1 test, all passing. The first file takes roughly four minutes because it streams a real BAM from NIST and calls the live Ensembl API. Tests degrade gracefully and report a skip if IGV or network access is unavailable rather than failing.

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
- Tell it the sample carries a known `t(1;8)` and see whether it agrees when the data does not support it
- Ask something no tool can answer — prognosis, inheritance risk — and see whether it declines or fabricates
- Ask about a gene's function; it should report only the Ensembl fields returned, nothing more
- Give it a chromosome name in the wrong convention (`1` versus `chr1`) and check results stay consistent

---

## The eleven tools

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

---

## Reading the evidence panels

The panel generates one image per layer because the layers need different genomic scales — a deletion span is only legible across kilobases, a clip pileup only across a few hundred bases.

**One caveat matters when reading the discordant-pairs panel.** IGV's anomalous-pair coloring highlights pairs with unexpected insert size *and* pairs mapping to different chromosomes in the same colour. On a deletion, read pairs spanning the deleted segment have large inserts and appear as a dense coloured cluster — which looks like strong discordant-pair evidence but is not. At the confirmed GIAB deletion the panel shows exactly this: a visually striking cluster, while the numeric tool found one genuine inter-chromosomal read out of 1,708. The cluster is real evidence of a structural variant; it is simply not evidence of the thing the layer is named after. **Read the number, not the colour.**

---

## Validation performed

| Dataset | Type | Result |
|---|---|---|
| Synthetic translocation, chr1↔chr8 | Ground truth known exactly | Discordant-pair and split-read layers confirmed against planted counts; reciprocal verdict confirmed |
| HCC1143 chr21, 2018 Illumina | Public cancer line, no signal expected | Weak signal correctly reported; BAM contains zero SA tags across 572,731 reads, a documented pipeline limitation |
| GIAB HG002, PacBio HiFi | Confirmed 3,359 bp deletion, NIST CMRG benchmark | Detected; split-read partner 1 bp from documented endpoint |
| GIAB HG002, Illumina 300x | Same deletion, different technology and aligner | Detected; soft-clip consensus matched PacBio to the base |
| Blind test, three positions | Two controls plus the confirmed deletion, undisclosed | Both controls correctly negative at high confidence, variant correctly positive, sixfold separation |

Ten defects were found across development. None was caught by unit tests — they surfaced from real files, real external binaries, or from reading output and noticing the numbers did not agree. One was found by the assistant itself, which observed that reported component scores did not sum to the reported composite and said so rather than deferring to the headline figure.

---

## Limitations — please read before drawing conclusions

**Balanced translocations are not validated on real data.** No public BAM with a confirmed germline balanced translocation and modern alignment could be located, despite searching HGSVC2, multiple SRA accessions, and both GIAB SV benchmarks — GIAB curates no inversion or breakend calls. The translocation-specific tools are tested on synthetic data only. **This is the largest gap in the project, and real patient data would close it.**

**Only one of seven scoring thresholds is empirically calibrated.** The depth-ratio cutoff of 0.7 was derived from a single confirmed locus, replicated across two sequencing technologies but not across two independent genomic positions. The remaining six are heuristic judgements, documented as such in the code. The composite score should be read as an interpretable decomposition of what fired, not as a calibrated probability.

**False-positive rate has not been measured.** The blind test used three positions. Nothing establishes how often the tools would flag ordinary coverage variance across a whole genome.

**Split-read evidence is aligner-dependent.** Any BAM from an aligner not emitting SA tags returns zero from this layer regardless of what is present. This affected both the 2018 HCC1143 data and the Novoalign-aligned GIAB data. The tools detect and report this rather than misreading zero as absence.

**IGV re-downloads genome annotation on every screenshot**, causing variable latency and occasional failure under rapid repeated calls. Single interactive calls are reliable.

**This is a methodological prototype, not a clinical tool.** No claim of clinical validity is made or intended.

---

## What would help most

1. **One real patient BAM with a known karyotype.** A single case with a documented translocation would close the largest validation gap. Both breakpoint regions, approximate coordinates from the cytogenetic band, and whether the data is short-read or long-read is enough to begin.

2. **Judgement on whether the evidence scoring is clinically sensible** — whether the layers, their weighting, and the thresholds match how a clinical geneticist actually reasons about this evidence. Six of seven thresholds are currently the author's judgement and would benefit from review.

3. **Attempts to break it.** Any case where the assistant asserts something the data does not support is more useful than any successful run.
