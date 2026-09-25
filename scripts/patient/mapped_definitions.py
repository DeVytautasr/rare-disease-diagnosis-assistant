#!/usr/bin/env python3
"""Which definition of "mapped reads" reproduces the Phase 0 figures?

flagstat's "primary mapped" did not equal the Phase 0 figure for either sample
(2026-09-25) while byte size, @PG count and contig md5 all matched. Before
concluding that the files differ, this counts candidate definitions for both
samples and reports every one that reproduces BOTH Phase 0 figures exactly.

Per contig class (chr1-22, chrX, chrY, chrM, chrEBV, *_alt, HLA-*, *_decoy,
*_random, chrUn_*):
  index_mapped    idxstats "mapped" from the .bai (all records: primary,
                  secondary and supplementary)
  primary_mapped  records that are mapped, not secondary, not supplementary --
                  counted by reading the records, for every class except the
                  main chromosomes, whose figure is the flagstat total minus
                  the rest (read from verify_<label>_flagstat.json)

Candidates are sums over class subsets that contain every autosome, i.e. "the
chromosomes plus some extra classes"; a hit is reported with the subset named.
Counts only; no names, no coordinates.
"""
import itertools
import json
import os
import re
import sys

import pysam

from common import PHASE0, RUN_DIR, locate


def cls(name):
    if re.fullmatch(r"chr([1-9]|1\d|2[0-2])", name):
        return "chr1-22"
    for pat, c in ((r"chrX", "chrX"), (r"chrY", "chrY"), (r"chrM", "chrM"), (r"chrEBV", "chrEBV"),
                   (r".*_alt", "alt"), (r"HLA-.*", "HLA"), (r".*_decoy", "decoy"),
                   (r".*_random", "random"), (r"chrUn_.*", "chrUn")):
        if re.fullmatch(pat, name):
            return c
    return "other"


def main():
    samples = locate()
    out = {}
    for label, s in samples.items():
        fs = json.load(open(os.path.join(RUN_DIR, f"verify_{label}_flagstat.json")))
        total_primary_mapped = fs["flagstat"]["QC-passed reads"]["primary mapped"] + \
            fs["flagstat"]["QC-failed reads"]["primary mapped"]
        idx, prim = {}, {}
        for line in pysam.idxstats(s["bam"]).splitlines():
            name, length, m, u = line.split("\t")
            if name == "*":
                continue
            c = cls(name)
            idx[c] = idx.get(c, 0) + int(m)
        with pysam.AlignmentFile(s["bam"], "rb") as f:
            for name in f.references:
                c = cls(name)
                if c in ("chr1-22",):
                    continue
                n = 0
                for r in f.fetch(name):
                    if not (r.is_unmapped or r.is_secondary or r.is_supplementary):
                        n += 1
                prim[c] = prim.get(c, 0) + n
        prim["chr1-22"] = total_primary_mapped - sum(prim.values())
        out[label] = {"index_mapped": idx, "primary_mapped": prim}

    classes = sorted(set(out["SAMPLE_A"]["index_mapped"]) | set(out["SAMPLE_B"]["index_mapped"]))
    extras = [c for c in classes if c != "chr1-22"]
    print("per class           " + "".join(f"{l + ' ' + k:>30}" for l in PHASE0
                                              for k in ("index_mapped", "primary_mapped")))
    for c in classes:
        print(f"  {c:16s}" + "".join(f"{out[l][k].get(c, 0):>30,}" for l in PHASE0
                                     for k in ("index_mapped", "primary_mapped")))
    hits = []
    for kind in ("index_mapped", "primary_mapped"):
        for r in range(len(extras) + 1):
            for sub in itertools.combinations(extras, r):
                keep = ("chr1-22",) + sub
                vals = {l: sum(out[l][kind].get(c, 0) for c in keep) for l in PHASE0}
                if all(vals[l] == PHASE0[l]["primary_mapped"] for l in PHASE0):
                    hits.append((kind, keep, vals))
    print()
    if not hits:
        print("NO candidate definition reproduces both Phase 0 figures")
    for kind, keep, vals in hits:
        print(f"REPRODUCES BOTH: {kind} summed over {' + '.join(keep)} -> "
              + ", ".join(f"{l} {v:,}" for l, v in vals.items()))
    with open(os.path.join(RUN_DIR, "mapped_definitions.json"), "w") as f:
        json.dump({"per_sample": out, "hits": [[k, list(c), v] for k, c, v in hits]}, f, indent=1)
    sys.exit(0 if hits else 1)


if __name__ == "__main__":
    main()
