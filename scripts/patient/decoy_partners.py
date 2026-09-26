#!/usr/bin/env python3
"""Phase 13 Task 2: how many counted split reads have a partner on a decoy contig,
and whether any score, band or partner statement depends on them. Measure only:
the evidence tools are not modified.

    decoy_partners.py synthetic          the 32 detected implant breakends, the 8
                                         breakpoints of the 4 missed implants, and the
                                         breakends of the background's 27 survivors
    decoy_partners.py patient LABEL      100 seeded random autosomal positions of one
                                         patient BAM -- aggregate totals only
    decoy_partners.py record             both results -> the committed record

HOW THE TOOL TREATS A DECOY PARTNER (read from bam_tools.py, not changed): in
get_split_reads a primary read that passes the read filters and carries an SA tag
counts as a split read if at least one of its SA entries survives the entry's own
mapQ filter (sa_mapq >= min_mapq). There is no filter on the entry's contig: a
decoy (name ending _decoy) is kept, added to partner_chromosomes under its own name,
and the read counts toward split_reads and split_read_fraction. No distinction is
made between inter- and intra-chromosomal entries. summarize_breakpoint_evidence
scores the split layer from split_read_fraction and split_reads (tiers 0.1 / 0.3,
min_support) and builds its partner sentence with _describe_partner_distribution
over partner_chromosomes, so a decoy can be the named partner.

METHOD. mirror_split_reads repeats get_split_reads' read loop (bam_tools.py lines
1165-1245) and records, for each counted split read, the set of partner contigs
it contributed. Control: its split_reads, total and partner_chromosomes must equal
the tool's on every position, or the position is reported as a mismatch and not
used. The counterfactual summary calls summarize_breakpoint_evidence unchanged,
with get_split_reads replaced IN THIS PROCESS ONLY by a wrapper that removes the
decoy partners from each read and drops reads left with none; bam_tools.py is not
edited. At the patient positions only aggregate totals are kept.

FIGURES AND THEIR DEFINITIONS (written into the record beside each figure)
  decoy_contig         a contig whose name ends in _decoy
  split_reads          the tool's own count at the summary's geometry: window
                       min(500, 200) = 200 bp, min_mapq 20, and also at min_mapq 0
  with_decoy_partner   counted split reads with at least one surviving SA entry on a
                       decoy contig
  decoy_only           counted split reads whose surviving SA entries are all on
                       decoy contigs; removing these is the counterfactual
  change               at a synthetic position, whether the evidence score, the
                       evidence strength, the split-read component or the split-read
                       observation sentence differs between the summary as it is and
                       the counterfactual (window 500, min_mapq 20, all four layers)
  patient_positions    100 positions per BAM drawn uniformly over chr1-22 (seed
                       20260927); a position whose split-read call fails is counted,
                       not retried. Kept: totals over the positions, never a position
"""
import gzip
import json
import os
import random
import sys
from collections import Counter

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
from common import DEID_DIR, RUN_DIR, die  # noqa: E402
from stage1_igv_assistant.tools import bam_tools as bt  # noqa: E402

HOME = os.path.expanduser("~")
BAMS = os.path.join(HOME, "public_data", "sim", "bams")
BG = os.path.join(HOME, "public_data", "NA12878.chr20_chr21.bam")
AN = os.path.join(REPO, "stage1_igv_assistant", "results", "synthetic_control_2026-09", "analysis_2026-09-25")
LAYERS = ["discordant_pairs", "soft_clipped_reads", "split_reads", "read_depth"]
SEED = 20260927
N_POS = 100
AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
RECORD = os.path.join(REPO, "stage1_igv_assistant", "results", "decoy_partners_2026-09-27.json")
OUT = {"synthetic": os.path.join(HOME, "public_data", "sim", "logs", "phase13_2026-09-27", "decoy_synthetic.json")}


def definitions():
    doc = __doc__.split("FIGURES AND THEIR DEFINITIONS (written into the record beside each figure)")[1]
    out, key = {}, None
    for line in doc.splitlines():
        if line.startswith("  ") and not line.startswith("    ") and line.strip():
            key, _, rest = line.strip().partition(" ")
            out[key] = rest.strip()
        elif key and line.strip():
            out[key] += " " + line.strip()
    return {k: " ".join(v.split()) for k, v in out.items()}


def is_decoy(name):
    return name.endswith("_decoy")


def mirror_split_reads(bam_path, chromosome, position, window_bp=200, min_mapq=20):
    """get_split_reads' read loop, returning per-read partner sets (see METHOD)."""
    start, end = max(0, position - window_bp), position + window_bp
    reads = []
    total = 0
    with pysam.AlignmentFile(bam_path) as bam:
        refs = set(bam.references)
        chrom = bt._resolve_contig(bam, chromosome)
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            total += 1
            if not read.has_tag("SA"):
                continue
            partners = set()
            for entry in read.get_tag("SA").rstrip(";").split(";"):
                f = entry.split(",")
                if len(f) < 2:
                    continue
                if len(f) >= 5:
                    try:
                        q = int(f[4])
                    except ValueError:
                        q = None
                    if q is not None and q < min_mapq:
                        continue
                partners.add(bt._canonical_chrom(f[0], refs))
            if partners:
                reads.append(partners)
    return total, reads


def tally(reads):
    return {"split_reads": len(reads),
            "with_decoy_partner": sum(1 for p in reads if any(is_decoy(c) for c in p)),
            "decoy_only": sum(1 for p in reads if p and all(is_decoy(c) for c in p))}


def control(bam_path, chrom, pos, window_bp, min_mapq):
    tool = bt.get_split_reads(bam_path, chrom, pos, window_bp=window_bp, min_mapq=min_mapq)
    total, reads = mirror_split_reads(bam_path, chrom, pos, window_bp, min_mapq)
    pc = Counter(c for p in reads for c in p)
    ok = ("error" not in tool and tool["split_reads"] == len(reads) and tool.get("total_reads_in_window") in (total, None)
          and dict(pc) == dict(tool["partner_chromosomes"]))
    return ok, reads, tool


ORIGINAL = bt.get_split_reads


def without_decoys(bam_path, chromosome, position, window_bp=200, min_mapq=20):
    r = ORIGINAL(bam_path, chromosome, position, window_bp=window_bp, min_mapq=min_mapq)
    if "error" in r:
        return r
    total, reads = mirror_split_reads(bam_path, chromosome, position, window_bp, min_mapq)
    kept = [{c for c in p if not is_decoy(c)} for p in reads]
    kept = [p for p in kept if p]
    r = dict(r)
    r["split_reads"] = len(kept)
    n = r.get("total_reads_in_window") or 0  # the tool's own denominator
    r["split_read_fraction"] = round(len(kept) / n, 3) if n else r.get("split_read_fraction")
    r["partner_chromosomes"] = dict(sorted(Counter(c for p in kept for c in p).items(), key=lambda x: -x[1]))
    return r


def summary_pair(bam_path, chrom, pos):
    a = bt.summarize_breakpoint_evidence(bam_path, chrom, pos, applicable_layers=LAYERS)
    bt.get_split_reads = without_decoys
    try:
        b = bt.summarize_breakpoint_evidence(bam_path, chrom, pos, applicable_layers=LAYERS)
    finally:
        bt.get_split_reads = ORIGINAL
    split_obs = lambda s: [o for o in s.get("supporting_observations", []) if "split read" in o]  # noqa: E731
    comp = lambda s: s.get("split_read_score")  # noqa: E731
    return {"score": [a.get("evidence_score"), b.get("evidence_score")],
            "strength": [a.get("evidence_strength"), b.get("evidence_strength")],
            "split_component": [comp(a), comp(b)],
            "changed": {"score": a.get("evidence_score") != b.get("evidence_score"),
                        "strength": a.get("evidence_strength") != b.get("evidence_strength"),
                        "split_component": comp(a) != comp(b),
                        "split_sentence": split_obs(a) != split_obs(b)},
            "split_sentence": [split_obs(a), split_obs(b)]}


def synthetic_positions():
    f = json.load(open(os.path.join(AN, "figures.json")))
    pos = []
    for r in f["ceiling"]:
        imp = r["junction"].split()[0]
        c, p = r["breakend"].split(":")
        pos.append(("detected", imp, os.path.join(BAMS, f"{imp}.bam"), c, int(p)))
    for imp, v in f["per_implant"].items():
        if not any(j.get("detected") for j in v["junctions"].values()):
            for c in ("chr20", "chr21"):
                pos.append(("missed", imp, os.path.join(BAMS, f"{imp}.bam"), c, int(v["breakpoints"][c])))
    d = json.load(gzip.open(os.path.join(AN, "chain", "background.json.gz")))
    seen = set()
    for e in d["evidence"]:
        for be in e["breakends"].values():
            k = (be["chromosome"], be["position"])
            if k not in seen:
                seen.add(k)
                pos.append(("background", "background", BG, be["chromosome"], be["position"]))
    return pos


def step_synthetic():
    rows, mism = [], []
    for group, label, bam, c, p in synthetic_positions():
        row = {"group": group, "label": label, "position": f"{c}:{p}"}
        for q in (20, 0):
            ok, reads, tool = control(bam, c, p, 200, q)
            if not ok:
                mism.append(f"{label} {c}:{p} mapq{q}")
            row[f"mapq{q}"] = {**tally(reads), "control_matches_tool": ok,
                               "decoy_partner_names": sorted({x for s in reads for x in s if is_decoy(x)})[:5]}
        row["summary"] = summary_pair(bam, c, p)
        rows.append(row)
        print(group, label, row["position"], row["mapq20"]["split_reads"], row["mapq20"]["with_decoy_partner"],
              row["mapq0"]["with_decoy_partner"], row["summary"]["changed"], flush=True)
    agg = {}
    for g in ("detected", "missed", "background"):
        rr = [r for r in rows if r["group"] == g]
        agg[g] = {"positions": len(rr),
                  **{f"mapq{q}_{k}": sum(r[f"mapq{q}"][k] for r in rr)
                     for q in (20, 0) for k in ("split_reads", "with_decoy_partner", "decoy_only")},
                  "positions_with_any_change": sum(1 for r in rr if any(r["summary"]["changed"].values())),
                  "changes": {k: sum(1 for r in rr if r["summary"]["changed"][k]) for k in
                              ("score", "strength", "split_component", "split_sentence")}}
    rec = {"aggregate": agg, "control_mismatches": mism, "positions": rows}
    os.makedirs(os.path.dirname(OUT["synthetic"]), exist_ok=True)
    json.dump(rec, open(OUT["synthetic"], "w"), indent=1)
    print(json.dumps(agg, indent=1), "\ncontrol mismatches:", mism)
    return 0


def step_patient(label):
    bam_path = os.path.join(DEID_DIR, f"{label}.bam")
    if not os.path.exists(bam_path):
        die(f"{label}: the de-identified link is missing")
    with pysam.AlignmentFile(bam_path) as bam:
        lengths = {c: bam.get_reference_length(c) for c in AUTOSOMES}
    rng = random.Random(f"{SEED}:{label}")
    total_len = sum(lengths.values())
    agg = Counter()
    for _ in range(N_POS):
        x = rng.randrange(total_len)
        for c in AUTOSOMES:
            if x < lengths[c]:
                break
            x -= lengths[c]
        agg["positions"] += 1
        for q in (20, 0):
            ok, reads, tool = control(bam_path, c, x, 200, q)
            if "error" in tool:
                agg[f"mapq{q}_failed_calls"] += 1
                continue
            agg[f"mapq{q}_control_mismatches"] += (not ok)
            t = tally(reads)
            for k, v in t.items():
                agg[f"mapq{q}_{k}"] += v
            agg[f"mapq{q}_positions_with_split_reads"] += t["split_reads"] > 0
            agg[f"mapq{q}_positions_with_decoy_partner"] += t["with_decoy_partner"] > 0
    rec = {"label": label, "seed": f"{SEED}:{label}", **dict(agg)}
    for q in (20, 0):
        n = rec.get(f"mapq{q}_split_reads", 0)
        rec[f"mapq{q}_fraction_with_decoy_partner"] = round(rec.get(f"mapq{q}_with_decoy_partner", 0) / n, 4) if n else None
    os.makedirs(RUN_DIR, exist_ok=True)
    json.dump(rec, open(os.path.join(RUN_DIR, f"decoy_{label}.json"), "w"), indent=1)
    print(json.dumps(rec, indent=1))
    return 0


def step_record():
    out = {"what": "Decoy partners in the split-read layer, measured, not fixed (Phase 13 Task 2, 2026-09-27)",
           "treatment_in_code": __doc__.split("HOW THE TOOL TREATS A DECOY PARTNER (read from bam_tools.py, not changed):")[1]
                                .split("METHOD.")[0].strip().replace("\n", " "),
           "method": __doc__.split("METHOD.")[1].split("FIGURES AND")[0].strip().replace("\n", " "),
           "definitions": definitions(), "code": "scripts/patient/decoy_partners.py",
           "blinding": "patient figures are totals over 100 seeded positions per BAM; no position is kept"}
    s = json.load(open(OUT["synthetic"]))
    out["synthetic"] = {"aggregate": s["aggregate"], "control_mismatches": s["control_mismatches"],
                        "positions": s["positions"]}
    for label in ("SAMPLE_A", "SAMPLE_B"):
        p = os.path.join(RUN_DIR, f"decoy_{label}.json")
        if not os.path.exists(p):
            die(f"{label}: not run")
        out[label] = json.load(open(p))
    if os.path.exists(RECORD):
        die("the record exists; a committed record is never overwritten")
    json.dump(out, open(RECORD, "w"), indent=1)
    print("written:", os.path.relpath(RECORD, REPO))
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if a == ["synthetic"]:
        sys.exit(step_synthetic())
    if len(a) == 2 and a[0] == "patient" and a[1] in ("SAMPLE_A", "SAMPLE_B"):
        sys.exit(step_patient(a[1]))
    if a == ["record"]:
        sys.exit(step_record())
    die("usage: decoy_partners.py synthetic | patient SAMPLE_A|SAMPLE_B | record")
