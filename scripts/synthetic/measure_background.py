#!/usr/bin/env python3
"""Measure the NA12878 background that the implants are simulated to match.

    measure_background.py background [--out DIR]
        base-quality bins, fragment length, sequencing error rate
    measure_background.py profiles [--out DIR] [--threads N]
        simulate every 150 bp-capable ART profile over the same region, align it
        exactly as the implants will be aligned, and measure it the same way

DEFINITIONS (fixed before measuring)

Quality bins -- every base quality of the first 200,000 primary reads in each of
five fixed 1 Mb windows (chr20 at 5, 20 and 45 Mb; chr21 at 20 and 35 Mb).

Fragment length -- pairs flagged paired + proper, primary, not duplicate, not
QC-fail, MAPQ >= 20, mate on the same contig, counted once through the read
with TLEN > 0, TLEN <= 1500; chr20:31,000,000-41,000,000 (the q arm past the
centromere). Mean and population SD go to ART; median, MAD and the share of
proper pairs above 1500 are reported with them.

Error rate -- "isolated" mismatches per aligned base. Reads: primary, mapped,
not duplicate, not QC-fail, MAPQ >= 20. Bases: aligned read/reference pairs
(M, =, X; soft clips and indels excluded) inside chr20:35,000,000-36,000,000,
neither base N. A mismatch is keyed by (reference position, read base): seen in
exactly one read it is ISOLATED (sequencing error); in two or more, RECURRENT
(a variant, not counted as error). Rate = isolated / aligned bases.

Profiles -- HS25, HSXn, HSXt, MSv1, MSv3 (every built-in ART profile that can
produce 150 bp reads), paired, 150 bp, 30x over the region plus 1 kb each side,
fragment mean/SD from the background measurement, one fixed seed per profile.
All five are aligned in ONE bwa 0.7.15 mem -Y -K 100000000 call against the
full hs38DH index with its .alt, then measured by the SAME function as the
background. The profile closest to the background's rate is reported.
"""
import argparse
import concurrent.futures as cf
import json
import os
import statistics
import subprocess
import sys
import time

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from synth_common import (ART, BG, BWA, BWA_ARGS, REF, SAMTOOLS, SIM, tilde,  # noqa: E402
                          tool_versions, write_json)

QUAL_WINDOWS = [("chr20", 5_000_000), ("chr20", 20_000_000), ("chr20", 45_000_000),
                ("chr21", 20_000_000), ("chr21", 35_000_000)]
FRAG_REGION = ("chr20", 31_000_000, 41_000_000)
FRAG_CAP = 1500
ERR_REGION = ("chr20", 35_000_000, 36_000_000)
PROFILES = ["HS25", "HSXn", "HSXt", "MSv1", "MSv3"]
PROFILE_SEED_BASE = 20260925
PROFILE_FOLD = 30


def keep_read(r, min_mapq=20):
    return not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate
                or r.is_qcfail) and r.mapping_quality >= min_mapq


def quality_bins():
    counts = {}
    with pysam.AlignmentFile(BG) as f:
        for chrom, start in QUAL_WINDOWS:
            n = 0
            for r in f.fetch(chrom, start, start + 1_000_000):
                if r.is_secondary or r.is_supplementary or r.query_qualities is None:
                    continue
                for q in r.query_qualities:
                    counts[q] = counts.get(q, 0) + 1
                n += 1
                if n >= 200_000:
                    break
    return {"distinct": sorted(counts), "counts": {str(k): counts[k] for k in sorted(counts)},
            "bases": sum(counts.values())}


def fragment_length():
    chrom, a, b = FRAG_REGION
    tl, over = [], 0
    with pysam.AlignmentFile(BG) as f:
        for r in f.fetch(chrom, a, b):
            if not (r.is_paired and r.is_proper_pair) or not keep_read(r):
                continue
            if r.reference_id != r.next_reference_id or r.template_length <= 0:
                continue
            if r.template_length > FRAG_CAP:
                over += 1
                continue
            tl.append(r.template_length)
    med = statistics.median(tl)
    mad = statistics.median(abs(x - med) for x in tl)
    return {"region": f"{chrom}:{a:,}-{b:,}", "pairs": len(tl), "mean": round(statistics.fmean(tl), 2),
            "sd": round(statistics.pstdev(tl), 2), "median": med, "mad": mad,
            "proper_pairs_above_cap": over, "cap": FRAG_CAP}


def error_rate(bam_path, region=ERR_REGION):
    chrom, a, b = region
    with pysam.FastaFile(REF) as fa:
        ref = fa.fetch(chrom, a, b).upper()
    bases = reads = 0
    mism = {}
    with pysam.AlignmentFile(bam_path) as f:
        for r in f.fetch(chrom, a, b):
            if not keep_read(r):
                continue
            reads += 1
            seq = r.query_sequence
            for q, p in r.get_aligned_pairs(matches_only=True):
                if p < a or p >= b:
                    continue
                rb, qb = ref[p - a], seq[q]
                if rb == "N" or qb == "N":
                    continue
                bases += 1
                if qb != rb:
                    k = (p, qb)
                    mism[k] = mism.get(k, 0) + 1
    isolated = sum(1 for v in mism.values() if v == 1)
    recurrent_keys = sum(1 for v in mism.values() if v >= 2)
    recurrent_bases = sum(v for v in mism.values() if v >= 2)
    return {"bam": tilde(bam_path), "region": f"{chrom}:{a:,}-{b:,}", "reads": reads,
            "aligned_bases": bases, "isolated_mismatches": isolated,
            "recurrent_sites": recurrent_keys, "recurrent_mismatch_bases": recurrent_bases,
            "error_rate_isolated": round(isolated / bases, 6) if bases else None,
            "mismatch_rate_all": round(sum(mism.values()) / bases, 6) if bases else None}


def step_background(out):
    t0 = time.monotonic()
    with cf.ProcessPoolExecutor(max_workers=3) as ex:
        fq, ff, fe = ex.submit(quality_bins), ex.submit(fragment_length), ex.submit(error_rate, BG)
        rec = {"background": tilde(BG), "background_bytes": os.path.getsize(BG),
               "quality_bins": fq.result(), "fragment_length": ff.result(), "error": fe.result()}
    rec["wall_seconds"] = round(time.monotonic() - t0, 1)
    write_json(os.path.join(out, "background.json"), rec)
    print(json.dumps({k: v for k, v in rec.items() if k != "quality_bins"}, indent=1))
    print("quality values seen:", rec["quality_bins"]["distinct"])


def step_profiles(out, threads):
    bg = json.load(open(os.path.join(out, "background.json")))
    mean, sd = bg["fragment_length"]["mean"], bg["fragment_length"]["sd"]
    work = os.path.join(out, "profiles")
    os.makedirs(work, exist_ok=True)
    chrom, a, b = ERR_REGION
    fa_path = os.path.join(work, "region.fa")
    with pysam.FastaFile(REF) as fa, open(fa_path, "w") as o:
        o.write(f">{chrom}_{a - 1000}\n{fa.fetch(chrom, a - 1000, b + 1000)}\n")
    seeds = {}
    r1_all, r2_all = os.path.join(work, "pool_R1.fq"), os.path.join(work, "pool_R2.fq")
    with open(r1_all, "w") as o1, open(r2_all, "w") as o2:
        for i, prof in enumerate(PROFILES):
            seed = PROFILE_SEED_BASE + i
            seeds[prof] = seed
            pre = os.path.join(work, prof)
            subprocess.run([ART, "-ss", prof, "-i", fa_path, "-p", "-l", "150", "-f", str(PROFILE_FOLD),
                            "-m", str(round(mean)), "-s", str(round(sd)), "-rs", str(seed), "-na",
                            "-q", "-o", pre], check=True, stdout=subprocess.DEVNULL)
            for src, dst in ((pre + "1.fq", o1), (pre + "2.fq", o2)):
                with open(src) as f:
                    for n, line in enumerate(f):
                        dst.write(f"@{prof}_{line[1:]}" if n % 4 == 0 else line)
    sam = os.path.join(work, "pool.sam")
    t0 = time.monotonic()
    with open(sam, "w") as o, open(os.path.join(work, "bwa.stderr"), "w") as e:
        subprocess.run([BWA, *BWA_ARGS, "-t", str(threads), REF, r1_all, r2_all],
                       stdout=o, stderr=e, check=True)
    bwa_s = round(time.monotonic() - t0, 1)
    bams = {}
    with pysam.AlignmentFile(sam) as f:
        outs = {p: pysam.AlignmentFile(os.path.join(work, f"{p}.unsorted.bam"), "wb", template=f)
                for p in PROFILES}
        for r in f:
            outs[r.query_name.split("_", 1)[0]].write(r)
        for o in outs.values():
            o.close()
    for p in PROFILES:
        bams[p] = os.path.join(work, f"{p}.bam")
        subprocess.run([SAMTOOLS, "sort", "-@", "2", "-o", bams[p],
                        os.path.join(work, f"{p}.unsorted.bam")], check=True)
        subprocess.run([SAMTOOLS, "index", bams[p]], check=True)
    with cf.ProcessPoolExecutor(max_workers=len(PROFILES)) as ex:
        res = dict(zip(PROFILES, ex.map(error_rate, [bams[p] for p in PROFILES])))
    target = bg["error"]["error_rate_isolated"]
    ranked = sorted(PROFILES, key=lambda p: abs(res[p]["error_rate_isolated"] - target))
    rec = {"background_error_rate_isolated": target, "fragment_mean": round(mean),
           "fragment_sd": round(sd), "fold": PROFILE_FOLD, "seeds": seeds, "bwa_seconds": bwa_s,
           "bwa_command": " ".join([tilde(BWA), *BWA_ARGS, "-t", str(threads), tilde(REF),
                                    "pool_R1.fq", "pool_R2.fq"]),
           "profiles": res, "closest_first": ranked,
           "abs_difference": {p: round(abs(res[p]["error_rate_isolated"] - target), 6) for p in ranked},
           "versions": tool_versions()}
    write_json(os.path.join(out, "profiles.json"), rec)
    print(f"background isolated-mismatch rate {target}")
    for p in ranked:
        print(f"  {p}: {res[p]['error_rate_isolated']}  |diff| {rec['abs_difference'][p]}  "
              f"({res[p]['aligned_bases']:,} bases, {res[p]['reads']:,} reads)")
    print(f"closest profile: {ranked[0]}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["background", "profiles"])
    ap.add_argument("--out", default=os.path.join(SIM, "measure"))
    ap.add_argument("--threads", type=int, default=10)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.step == "background":
        step_background(a.out)
    else:
        step_profiles(a.out, a.threads)
