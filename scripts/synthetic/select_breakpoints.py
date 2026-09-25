#!/usr/bin/env python3
"""Classify the scanned candidates and select the 12 implant breakpoints.

    select_breakpoints.py [--out FILE]

CLASS THRESHOLDS -- fixed on 2026-09-25 after seeing the distribution of 978
candidates (scan_breakpoints.py summary); a breakend's class is:

  low_mappability  low_mapq_fraction_200 > 0.4          the tools' own quality gate
  repeat_adjacent  low_mapq_fraction_200 <= 0.4 and at least one repeat signal:
                     max_homopolymer >= 10 bp            (the 90th percentile; a common
                                                          "long homopolymer" cut)
                     max_dinucleotide_run >= 10 bp       (5 units, a common STR minimum;
                                                          the 95th percentile)
                     entropy2 < 3.6148                   (the 5th percentile)
                     low_mapq_fraction_200 >= 0.05       (partial MAPQ loss below the gate)
                     mappability < 0.9                   (partial empirical ambiguity)
  clean_unique     low_mapq_fraction_200 == 0 and mappability == 1.0 and
                   max_homopolymer <= 8 and max_dinucleotide_run <= 8 and
                   entropy2 >= 3.6897                    (the 10th percentile)
  unclassified     anything else -- a buffer between the classes, never selected

An implant's two breakends must both be in its class.

The recorded Phase 6 coordinates were displayed together with that
distribution, so they were visible when the thresholds were chosen. Under a
strict upper-decile rule (dinucleotide run >= 8 bp, the 90th percentile) the
chr21 end of IMP01 (a run of exactly 8 bp) would be repeat_adjacent, not clean.

SELECTION -- 12 implants: IMP01-04 clean_unique, IMP05-08 repeat_adjacent,
IMP09-12 low_mappability. The recorded coordinates are reused where they
survive: IMP01 chr20:200,000 / chr21:14,100,000; IMP06 chr20:3,900,000 /
chr21:17,200,000; IMP10 chr20:25,800,000 / chr21:7,600,000; IMP05
chr21:13,300,000; IMP11 chr20:31,100,000. Each is checked against its recorded
class; one that no longer qualifies is REPORTED and replaced by a draw of the
right class, never forced. Every other breakend is drawn uniformly from its
class pool on its chromosome with a seeded RNG, at least 1 Mb from every
breakend already selected on that chromosome, in implant order, chr20 before
chr21. The seed of every draw is recorded.
"""
import argparse
import hashlib
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth_common import SIM, write_json  # noqa: E402

MASTER_SEED = 20260925           # the same master seed as make_implants.py
MIN_SPACING = 1_000_000
T = {"low_mapq_gate": 0.4, "homopolymer_repeat": 10, "dinucleotide_repeat": 10,
     "entropy_repeat_below": 3.6148, "lmf_partial": 0.05, "mappability_partial_below": 0.9,
     "homopolymer_clean_max": 8, "dinucleotide_clean_max": 8, "entropy_clean_min": 3.6897}
CLASS_OF = {**{f"IMP{i:02d}": "clean_unique" for i in range(1, 5)},
            **{f"IMP{i:02d}": "repeat_adjacent" for i in range(5, 9)},
            **{f"IMP{i:02d}": "low_mappability" for i in range(9, 13)}}
RECORDED = {"IMP01": {"chr20": 200_000, "chr21": 14_100_000},
            "IMP06": {"chr20": 3_900_000, "chr21": 17_200_000},
            "IMP10": {"chr20": 25_800_000, "chr21": 7_600_000},
            "IMP05": {"chr21": 13_300_000},
            "IMP11": {"chr20": 31_100_000}}


def seed_for(name):
    return int(hashlib.sha256(f"{MASTER_SEED}:{name}".encode()).hexdigest()[:8], 16) & 0x7FFFFFFF


def classify(r):
    l, m = r["low_mapq_fraction_200"], r["mappability"]
    if l > T["low_mapq_gate"]:
        return "low_mappability"
    if (r["max_homopolymer"] >= T["homopolymer_repeat"] or r["max_dinucleotide_run"] >= T["dinucleotide_repeat"]
            or r["entropy2"] < T["entropy_repeat_below"] or l >= T["lmf_partial"]
            or m < T["mappability_partial_below"]):
        return "repeat_adjacent"
    if (l == 0 and m == 1.0 and r["max_homopolymer"] <= T["homopolymer_clean_max"]
            and r["max_dinucleotide_run"] <= T["dinucleotide_clean_max"] and r["entropy2"] >= T["entropy_clean_min"]):
        return "clean_unique"
    return "unclassified"


def main(out):
    scan = json.load(open(os.path.join(SIM, "scan", "candidates.json")))
    rows = scan["candidates"]
    for r in rows:
        r["class"] = classify(r)
    idx = {(r["chrom"], r["pos"]): r for r in rows}
    pools = {(c, k): sorted(r["pos"] for r in rows if r["chrom"] == c and r["class"] == k)
             for c in ("chr20", "chr21") for k in set(CLASS_OF.values())}
    checks, chosen, implants = [], {"chr20": [], "chr21": []}, []
    for iid in sorted(CLASS_OF):
        want = CLASS_OF[iid]
        imp = {"id": iid, "class": want, "source": {}, "stats": {}}
        for chrom in ("chr20", "chr21"):
            rec = RECORDED.get(iid, {}).get(chrom)
            if rec is not None:
                r = idx.get((chrom, rec))
                got = r["class"] if r else "not a candidate (masked or N in flank)"
                ok = got == want
                checks.append({"implant": iid, "breakend": f"{chrom}:{rec}", "recorded_class": want,
                               "recomputed_class": got, "confirmed": ok})
                if ok:
                    imp[chrom], imp["source"][chrom] = rec, "recorded (Phase 6), class confirmed"
                    chosen[chrom].append(rec)
                    imp["stats"][chrom] = r
                    continue
                imp["source"][chrom] = f"recorded {rec} no longer {want} ({got}); replaced by a draw"
            seed = seed_for(f"selection:{iid}:{chrom}")
            pool = [p for p in pools[(chrom, want)]
                    if all(abs(p - q) >= MIN_SPACING for q in chosen[chrom])
                    and all(abs(p - v) >= MIN_SPACING for i2, e in RECORDED.items() if i2 > iid
                            for c2, v in e.items() if c2 == chrom)]
            if not pool:
                raise SystemExit(f"{iid} {chrom}: no {want} candidate left at {MIN_SPACING:,} bp spacing")
            p = random.Random(seed).choice(pool)
            imp[chrom] = p
            imp["source"].setdefault(chrom, "")
            imp["source"][chrom] = (imp["source"][chrom] + "; " if imp["source"][chrom] else "") + \
                f"drawn (seed {seed}, pool {len(pool)})"
            imp["stats"][chrom] = idx[(chrom, p)]
            chosen[chrom].append(p)
        implants.append(imp)
    rec = {"master_seed": MASTER_SEED, "min_spacing_bp": MIN_SPACING, "thresholds": T,
           "threshold_rule": __doc__.split("CLASS THRESHOLDS")[1].split("SELECTION")[0].strip(),
           "pool_sizes": {f"{c}:{k}": len(v) for (c, k), v in sorted(pools.items())},
           "class_counts": {c: {k: sum(1 for r in rows if r["chrom"] == c and r["class"] == k)
                                for k in ("clean_unique", "repeat_adjacent", "low_mappability", "unclassified")}
                            for c in ("chr20", "chr21")},
           "recorded_coordinate_checks": checks, "implants": implants}
    write_json(out, rec)
    for c in checks:
        print(f"recorded {c['implant']} {c['breakend']:>18s}: recorded {c['recorded_class']:16s} "
              f"recomputed {c['recomputed_class']:16s} {'CONFIRMED' if c['confirmed'] else 'DOES NOT QUALIFY'}")
    for imp in implants:
        print(f"{imp['id']} {imp['class']:16s} chr20:{imp['chr20']:>11,}  chr21:{imp['chr21']:>11,}   "
              f"[{imp['source']['chr20']} | {imp['source']['chr21']}]")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(SIM, "selection.json"))
    main(ap.parse_args().out)
