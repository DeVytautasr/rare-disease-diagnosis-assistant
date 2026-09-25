#!/usr/bin/env python3
"""The PRE-REGISTERED gate, checked on IMP01 before anything downstream.

    gate_imp01.py [IMPxx]        (default IMP01)

  (a) at least 70% of truly spanning reads (per ART ground truth) carry an SA tag
      pointing at the partner chromosome                  Phase 6: 50/61 = 82%
  (b) zero SA tags pointing at a wrong partner            Phase 6: 0
  (c) local depth at both breakpoints within +/-10% of the pre-implant value
                                                          Phase 6: -0.8% and +5.3%

DEFINITIONS -- fixed here, committed before the gate is run:
  truly spanning read  a simulated read whose own ART template interval crosses its
                       junction, with at least one base on each side
  its record           its PRIMARY record in the implant BAM (not secondary, not
                       supplementary), found among the SIM read group's records
  partner chromosome   chr21 if that primary record is on chr20, chr20 if it is on
                       chr21; a primary record on any other contig, or unmapped,
                       fails (a)
  (a) numerator        truly spanning reads whose primary record has an SA tag with
                       at least one entry on the partner chromosome
  (b) wrong partner    SA entries, on those primary records, naming any contig other
                       than the partner chromosome. The same count over EVERY
                       simulated record is reported too, not gated.
  (c) local depth      the evidence tools' own mean_depth over +/-500 bp -- the
                       quality gate's window -- via bam_tools.get_bam_stats_at_locus:
                       pre = the untouched background, post = the implant BAM, at
                       chr20:A and chr21:B. Depth without duplicates is reported too.

Exit 0 all three hold; 1 at least one fails -- STOP; 2 could not evaluate.
"""
import json
import os
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from synth_common import BAMS, BG, SAMTOOLS, WORK, tilde, write_json  # noqa: E402
from stage1_igv_assistant.tools import bam_tools  # noqa: E402

MIN_SA_FRACTION = 0.70
MAX_WRONG_PARTNER = 0
MAX_DEPTH_CHANGE = 0.10
DEPTH_WINDOW = 500


def nondup_depth(bam_path, chrom, x):
    a, b = x - DEPTH_WINDOW, x + DEPTH_WINDOW
    depth = 0
    with pysam.AlignmentFile(bam_path) as f:
        for r in f.fetch(chrom, a, b):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                continue
            depth += sum(1 for p in r.get_reference_positions() if a <= p < b)
    return depth / (b - a)


def sa_entries(r):
    if not r.has_tag("SA"):
        return []
    return [e.split(",")[0] for e in r.get_tag("SA").split(";") if e]


def main(iid):
    prep = json.load(open(os.path.join(WORK, iid, "prepare.json")))
    A, B = prep["chr20"], prep["chr21"]
    bam = os.path.join(BAMS, f"{iid}.bam")
    sim_only = os.path.join(WORK, iid, "gate_sim_records.bam")
    subprocess.run([SAMTOOLS, "view", "-b", "-r", "SIM", "-o", sim_only, bam], check=True)
    primary, all_sa_other = {}, 0
    n_sim = 0
    with pysam.AlignmentFile(sim_only) as f:
        for r in f:
            n_sim += 1
            all_sa_other += sum(1 for c in sa_entries(r) if c not in ("chr20", "chr21"))
            if not (r.is_secondary or r.is_supplementary):
                primary[(r.query_name, 1 if r.is_read1 else 2)] = r
    spanning = [t for t in prep["truth_reads"]
                if t["spans_junction"] and t["bases_left"] >= 1 and t["bases_right"] >= 1]
    with_partner, wrong, missing, off_target = 0, 0, 0, 0
    per_read = []
    for t in spanning:
        r = primary.get((t["read"], t["mate"]))
        if r is None:
            missing += 1
            continue
        if r.is_unmapped or r.reference_name not in ("chr20", "chr21"):
            off_target += 1
            per_read.append({"read": t["read"], "mate": t["mate"], "primary": "unmapped" if r.is_unmapped
                             else r.reference_name, "sa": sa_entries(r), "partner_sa": False})
            continue
        partner = "chr21" if r.reference_name == "chr20" else "chr20"
        sa = sa_entries(r)
        hit = partner in sa
        with_partner += hit
        wrong += sum(1 for c in sa if c != partner)
        per_read.append({"read": t["read"], "mate": t["mate"], "primary": r.reference_name,
                         "bases_left": t["bases_left"], "bases_right": t["bases_right"],
                         "sa": sa, "partner_sa": hit})
    n = len(spanning)
    frac = with_partner / n if n else 0.0
    depth = {}
    for chrom, x in (("chr20", A), ("chr21", B)):
        pre = bam_tools.get_bam_stats_at_locus(BG, chrom, x - DEPTH_WINDOW, x + DEPTH_WINDOW)["mean_depth"]
        post = bam_tools.get_bam_stats_at_locus(bam, chrom, x - DEPTH_WINDOW, x + DEPTH_WINDOW)["mean_depth"]
        pre_nd, post_nd = nondup_depth(BG, chrom, x), nondup_depth(bam, chrom, x)
        depth[chrom] = {"position": x, "pre_mean_depth": pre, "post_mean_depth": post,
                        "change": round(post / pre - 1, 4) if pre else None,
                        "pre_nondup": round(pre_nd, 2), "post_nondup": round(post_nd, 2),
                        "change_nondup": round(post_nd / pre_nd - 1, 4) if pre_nd else None}
    a_ok = n > 0 and frac >= MIN_SA_FRACTION
    b_ok = wrong <= MAX_WRONG_PARTNER
    c_ok = all(d["change"] is not None and abs(d["change"]) <= MAX_DEPTH_CHANGE for d in depth.values())
    rec = {"implant": iid, "bam": tilde(bam), "breakpoints": {"chr20": A, "chr21": B},
           "a": {"truly_spanning_reads": n, "with_partner_sa": with_partner,
                 "fraction": round(frac, 4), "threshold": MIN_SA_FRACTION,
                 "primary_off_chr20_chr21": off_target, "not_found": missing, "holds": a_ok},
           "b": {"wrong_partner_sa_entries": wrong, "threshold": MAX_WRONG_PARTNER, "holds": b_ok,
                 "info_sa_entries_outside_chr20_chr21_over_all_simulated_records": all_sa_other,
                 "info_simulated_records": n_sim},
           "c": {"depth": depth, "threshold": MAX_DEPTH_CHANGE, "holds": c_ok},
           "verdict": "PASS" if a_ok and b_ok and c_ok else "FAIL", "per_read": per_read}
    write_json(os.path.join(WORK, iid, "gate.json"), rec)
    print(f"GATE on {iid} (chr20:{A:,} / chr21:{B:,})")
    print(f"  (a) {with_partner}/{n} truly spanning reads carry an SA tag to the partner = {frac:.1%}"
          f" (need >= {MIN_SA_FRACTION:.0%}; {off_target} primary off chr20/21, {missing} not found)"
          f"  -> {'HOLDS' if a_ok else 'FAILS'}")
    print(f"  (b) {wrong} SA entries point at a wrong partner (need {MAX_WRONG_PARTNER})"
          f"  -> {'HOLDS' if b_ok else 'FAILS'}   [info: {all_sa_other} SA entries outside chr20/21 "
          f"across all {n_sim} simulated records]")
    for chrom, d in depth.items():
        print(f"  (c) {chrom}:{d['position']:,} depth {d['pre_mean_depth']} -> {d['post_mean_depth']} "
              f"({d['change']:+.1%}); without duplicates {d['pre_nondup']} -> {d['post_nondup']} "
              f"({d['change_nondup']:+.1%})")
    print(f"  (c) both within +/-{MAX_DEPTH_CHANGE:.0%}  -> {'HOLDS' if c_ok else 'FAILS'}")
    print(f"GATE {rec['verdict']}")
    sys.exit(0 if rec["verdict"] == "PASS" else 1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "IMP01")
