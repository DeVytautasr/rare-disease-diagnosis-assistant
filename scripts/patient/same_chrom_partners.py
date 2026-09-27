#!/usr/bin/env python3
"""Phase 15 Task 2: same-chromosome SA partners in the split-read layer -- measure only.

    same_chrom_partners.py public          the 54 background breakends and 40 implant breakends
    same_chrom_partners.py patient LABEL   the 100 seeded Phase 13 positions -- totals only
    same_chrom_partners.py record          -> results/same_chrom_partners_2026-09-27.json

HOW THE CODE TREATS THEM (bam_tools.py after the decoy fix a6c887f, unchanged here):
get_split_reads keeps every SA entry that passes its own mapQ filter and is not on a
decoy contig, whatever its contig: an entry on the queried chromosome itself (or on
one of that chromosome's own alt/random contigs) is a partner like any other. The
read counts toward split_reads and split_read_fraction and the contig is a key of
partner_chromosomes. summarize_breakpoint_evidence scores the split layer from that
fraction and names partners through _describe_partner_distribution. By contrast,
count_discordant_pairs sets aside mates on the queried chromosome's own alt
contigs (same_primary_alt_mates, bam_tools.py line 885) and counts only
inter-chromosomal mates.

METHOD. As Phase 13's decoy_partners.py, on the fixed code: a mirror of
get_split_reads' loop (with the decoy exclusion) records each counted read's partner
set, and must reproduce the tool's split_reads and partner_chromosomes exactly at
every position (the control). The counterfactual calls summarize_breakpoint_evidence
unchanged with get_split_reads replaced IN THIS PROCESS ONLY by a wrapper that drops
same-chromosome partners and the reads left without one. No code is changed.

FIGURES AND THEIR DEFINITIONS (written into the record beside each figure)
  same_chromosome   an SA partner whose primary contig (bam_tools._primary_contig of
                    its canonical name) equals the queried chromosome's
  same_chrom_only   counted split reads whose every partner is same-chromosome
  with_same_chrom   counted split reads with at least one same-chromosome partner
  split_reads       the tool's count at the summary's geometry (window 200,
                    min_mapq 20), and at min_mapq 0 (the MCP tool's default)
  changed           whether the evidence score, band, split component or split-read
                    sentence differs between the summary as it is and the
                    counterfactual (window 500, min_mapq 20, all four layers)
  patient           the 100 positions per BAM of Phase 13 (seed '20260927:<label>');
                    totals only, no position kept
"""
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
import decoy_partners as dp  # noqa: E402
from stage1_igv_assistant.tools import bam_tools as bt  # noqa: E402

LAYERS = dp.LAYERS
LOGDIR = os.path.expanduser("~/public_data/sim/logs/phase15_2026-09-27")
RECORD = os.path.join(REPO, "stage1_igv_assistant", "results", "same_chrom_partners_2026-09-27.json")
ORIGINAL = bt.get_split_reads


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


def mirror(bam_path, chromosome, position, window_bp=200, min_mapq=20):
    start, end = max(0, position - window_bp), position + window_bp
    reads = []
    with pysam.AlignmentFile(bam_path) as bam:
        refs = set(bam.references)
        chrom = bt._resolve_contig(bam, chromosome)
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < min_mapq:
                continue
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
                c = bt._canonical_chrom(f[0], refs)
                if c.endswith("_decoy"):
                    continue
                partners.add(c)
            if partners:
                reads.append(partners)
        qprim = bt._primary_contig(bt._canonical_chrom(chrom, refs))
    return reads, qprim


def same(c, qprim):
    return bt._primary_contig(c) == qprim


def tally(bam, c, p, q):
    tool = ORIGINAL(bam, c, p, window_bp=200, min_mapq=q)
    if "error" in tool:
        return None
    reads, qprim = mirror(bam, c, p, 200, q)
    ok = tool["split_reads"] == len(reads) and dict(Counter(x for s in reads for x in s)) == dict(tool["partner_chromosomes"])
    return {"split_reads": len(reads), "same_chrom_only": sum(1 for s in reads if all(same(x, qprim) for x in s)),
            "with_same_chrom": sum(1 for s in reads if any(same(x, qprim) for x in s)), "control_matches_tool": ok}


def without_same(bam_path, chromosome, position, window_bp=200, min_mapq=20):
    r = ORIGINAL(bam_path, chromosome, position, window_bp=window_bp, min_mapq=min_mapq)
    if "error" in r:
        return r
    reads, qprim = mirror(bam_path, chromosome, position, window_bp, min_mapq)
    kept = [{x for x in s if not same(x, qprim)} for s in reads]
    kept = [s for s in kept if s]
    r = dict(r)
    n = r.get("total_reads_in_window") or 0
    r["split_reads"] = len(kept)
    r["split_read_fraction"] = round(len(kept) / n, 3) if n else r.get("split_read_fraction")
    r["partner_chromosomes"] = dict(sorted(Counter(x for s in kept for x in s).items(), key=lambda kv: -kv[1]))
    return r


def summary_pair(bam, c, p):
    sent = lambda s: [o for o in s.get("supporting_observations", []) if "split read" in o]  # noqa: E731
    a = bt.summarize_breakpoint_evidence(bam, c, p, applicable_layers=LAYERS)
    bt.get_split_reads = without_same
    try:
        b = bt.summarize_breakpoint_evidence(bam, c, p, applicable_layers=LAYERS)
    finally:
        bt.get_split_reads = ORIGINAL
    va = {"score": a.get("evidence_score"), "band": a.get("evidence_strength"),
          "split_component": a.get("split_read_score"), "split_sentence": sent(a)}
    vb = {"score": b.get("evidence_score"), "band": b.get("evidence_strength"),
          "split_component": b.get("split_read_score"), "split_sentence": sent(b)}
    return {"as_is": va, "counterfactual": vb, "changed": {k: va[k] != vb[k] for k in va}}


def step_public():
    rows, mism = [], []
    for group, label, bam, c, p in dp.synthetic_positions():
        row = {"group": group, "label": label, "position": f"{c}:{p}"}
        for q in (20, 0):
            t = tally(bam, c, p, q)
            row[f"mapq{q}"] = t
            if t and not t["control_matches_tool"]:
                mism.append(f"{row['position']} mapq{q}")
        row["summary"] = summary_pair(bam, c, p)
        rows.append(row)
        print(group, label, row["position"], row["mapq20"], row["summary"]["changed"], flush=True)
    agg = {}
    for g in ("detected", "missed", "background"):
        rr = [r for r in rows if r["group"] == g]
        agg[g] = {"positions": len(rr),
                  **{f"mapq{q}_{k}": sum((r[f"mapq{q}"] or {}).get(k, 0) for r in rr)
                     for q in (20, 0) for k in ("split_reads", "same_chrom_only", "with_same_chrom")},
                  "changed": {k: sum(r["summary"]["changed"][k] for r in rr)
                              for k in ("score", "band", "split_component", "split_sentence")}}
    os.makedirs(LOGDIR, exist_ok=True)
    json.dump({"aggregate": agg, "control_mismatches": mism, "positions": rows},
              open(os.path.join(LOGDIR, "same_chrom_public.json"), "w"), indent=1)
    print(json.dumps(agg, indent=1), "\ncontrol mismatches:", mism)
    return 0


def step_patient(label):
    bam = os.path.join(DEID_DIR, f"{label}.bam")
    if not os.path.exists(bam):
        die(f"{label}: the de-identified link is missing")
    with pysam.AlignmentFile(bam) as f:
        lengths = {c: f.get_reference_length(c) for c in dp.AUTOSOMES}
    rng = random.Random(f"{dp.SEED}:{label}")
    total = sum(lengths.values())
    agg = Counter()
    for _ in range(dp.N_POS):
        x = rng.randrange(total)
        for c in dp.AUTOSOMES:
            if x < lengths[c]:
                break
            x -= lengths[c]
        agg["positions"] += 1
        for q in (20, 0):
            t = tally(bam, c, x, q)
            if t is None:
                agg[f"mapq{q}_failed_calls"] += 1
                continue
            agg[f"mapq{q}_control_mismatches"] += not t["control_matches_tool"]
            for k in ("split_reads", "same_chrom_only", "with_same_chrom"):
                agg[f"mapq{q}_{k}"] += t[k]
        s = summary_pair(bam, c, x)
        for k, v in s["changed"].items():
            agg[f"positions_changed_{k}"] += v
    rec = {"label": label, "seed": f"{dp.SEED}:{label}", **dict(agg)}
    os.makedirs(RUN_DIR, exist_ok=True)
    json.dump(rec, open(os.path.join(RUN_DIR, f"same_chrom_{label}.json"), "w"), indent=1)
    print(json.dumps(rec, indent=1))
    return 0


def step_record():
    pub = json.load(open(os.path.join(LOGDIR, "same_chrom_public.json")))
    out = {"what": "Same-chromosome SA partners in the split-read layer, measured, not fixed (Phase 15 Task 2)",
           "treatment_in_code": " ".join(__doc__.split("HOW THE CODE TREATS THEM")[1].split("METHOD.")[0].split()),
           "method": " ".join(__doc__.split("METHOD.")[1].split("FIGURES AND")[0].split()),
           "definitions": definitions(), "code": "scripts/patient/same_chrom_partners.py", "public": pub,
           "blinding": "patient figures are totals over the 100 seeded positions per BAM"}
    for label in ("SAMPLE_A", "SAMPLE_B"):
        p = os.path.join(RUN_DIR, f"same_chrom_{label}.json")
        if not os.path.exists(p):
            die(f"{label}: not run")
        out[label] = json.load(open(p))
    if os.path.exists(RECORD):
        die("the record exists")
    json.dump(out, open(RECORD, "w"), indent=1)
    print("written:", os.path.relpath(RECORD, REPO))
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if a == ["public"]:
        sys.exit(step_public())
    if len(a) == 2 and a[0] == "patient" and a[1] in ("SAMPLE_A", "SAMPLE_B"):
        sys.exit(step_patient(a[1]))
    if a == ["record"]:
        sys.exit(step_record())
    die("usage: same_chrom_partners.py public | patient SAMPLE_A|SAMPLE_B | record")
