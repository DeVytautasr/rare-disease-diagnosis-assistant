#!/usr/bin/env python3
"""Scan candidate breakends on chr20 and chr21 and measure their context.

    scan_breakpoints.py stats        [--out DIR]   statistics that need no aligner
    scan_breakpoints.py mappability  [--out DIR] [--threads N]   empirical mappability
    scan_breakpoints.py summary      [--out DIR]   distributions, to choose class thresholds

CANDIDATES -- every multiple of 100 kb on chr20 and chr21 (the grid all recorded
Phase 6 coordinates lie on), kept only if the +/-5 kb flank overlaps no interval
of the unmodified delly v2.6.0 exclude template (telomeres, centromeres, chr21
heterochromatin) and contains no N. A position p means the junction falls after
base p (1-based): der(20) = chr20[1..A] + chr21[B+1..].

PER CANDIDATE
  low_mapq_fraction_200  the evidence tools' OWN function,
                         bam_tools.get_bam_stats_at_locus(bg, chrom, p-200, p+200):
                         MAPQ < 20 among primary mapped reads, duplicates counted,
                         unmapped/secondary/supplementary skipped. Class anchor.
  low_mapq_fraction_500  the same over +/-500, the window the tools' quality gate uses
  entropy2               Shannon entropy (bits, max 4) of the 2-mers in ref[p-200, p+200)
  max_homopolymer        longest single-base run in that window (bp)
  max_dinucleotide_run   longest run of one repeated 2-mer (bp), homopolymers excluded
  mappability            share of error-free 150 bp reads tiled over the window
                         (centres p-200..p+200 every 10 bp, both strands: 82 reads)
                         whose primary alignment returns to its origin (same contig,
                         start within 5 bp) with MAPQ >= 20. All candidates' tiles go
                         through ONE bwa 0.7.15 mem -Y -K 100000000 call against the
                         full hs38DH index with its .alt, as the implants will.
"""
import argparse
import collections
import json
import math
import os
import statistics
import subprocess
import sys
import time

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from synth_common import BG, BWA, BWA_ARGS, CONTIGS, REF, SIM, masks, tilde, write_json  # noqa: E402
from stage1_igv_assistant.tools import bam_tools  # noqa: E402

GRID = 100_000
FLANK = 5_000
WIN = 200
TILE_STEP = 10
TILE_LEN = 150
MAPQ_OK = 20
START_TOL = 5


def candidates():
    m = masks()
    out = []
    with pysam.FastaFile(REF) as fa:
        for chrom in CONTIGS:
            n = fa.get_reference_length(chrom)
            for p in range(GRID, n - FLANK, GRID):
                a, b = p - FLANK, p + FLANK
                if any(s < b and a < e for s, e, _ in m[chrom]):
                    continue
                if "N" in fa.fetch(chrom, a, b).upper():
                    continue
                out.append((chrom, p))
    return out


def entropy2(seq):
    c = collections.Counter(seq[i:i + 2] for i in range(len(seq) - 1))
    n = sum(c.values())
    return -sum(v / n * math.log2(v / n) for v in c.values())


def max_homopolymer(seq):
    best = cur = 0
    prev = None
    for ch in seq:
        cur = cur + 1 if ch == prev else 1
        prev = ch
        best = max(best, cur)
    return best


def max_dinuc(seq):
    best = 0
    # runs of a repeated 2-mer whose two bases differ (homopolymers are counted separately)
    for off in (0, 1):
        i = off
        while i + 2 <= len(seq):
            u = seq[i:i + 2]
            if u[0] == u[1]:
                i += 2
                continue
            j = i
            while j + 4 <= len(seq) and seq[j + 2:j + 4] == u:
                j += 2
            best = max(best, j + 2 - i)
            i = j + 2
    return best


def step_stats(out):
    t0 = time.monotonic()
    rows = []
    with pysam.FastaFile(REF) as fa:
        for chrom, p in candidates():
            s200 = bam_tools.get_bam_stats_at_locus(BG, chrom, p - WIN, p + WIN)
            s500 = bam_tools.get_bam_stats_at_locus(BG, chrom, p - 500, p + 500)
            seq = fa.fetch(chrom, p - WIN, p + WIN).upper()
            rows.append({
                "chrom": chrom, "pos": p,
                "low_mapq_fraction_200": s200["low_mapq_fraction"],
                "reads_200": s200["total_reads"], "mean_depth_200": s200["mean_depth"],
                "low_mapq_fraction_500": s500["low_mapq_fraction"],
                "entropy2": round(entropy2(seq), 4),
                "max_homopolymer": max_homopolymer(seq),
                "max_dinucleotide_run": max_dinuc(seq),
            })
    write_json(os.path.join(out, "candidates.json"),
               {"definition": __doc__, "background": tilde(BG), "grid_bp": GRID, "flank_bp": FLANK,
                "window_bp": WIN, "wall_seconds": round(time.monotonic() - t0, 1), "candidates": rows})
    print(f"{len(rows)} candidates ({sum(r['chrom'] == 'chr20' for r in rows)} chr20, "
          f"{sum(r['chrom'] == 'chr21' for r in rows)} chr21) in {time.monotonic() - t0:.0f}s")


def step_mappability(out, threads):
    path = os.path.join(out, "candidates.json")
    rec = json.load(open(path))
    rows = rec["candidates"]
    work = os.path.join(out, "mappability")
    os.makedirs(work, exist_ok=True)
    fq = os.path.join(work, "tiles.fq")
    comp = str.maketrans("ACGTN", "TGCAN")
    n_tiles = 0
    with pysam.FastaFile(REF) as fa, open(fq, "w") as o:
        for i, r in enumerate(rows):
            for centre in range(r["pos"] - WIN, r["pos"] + WIN + 1, TILE_STEP):
                s = centre - TILE_LEN // 2
                seq = fa.fetch(r["chrom"], s, s + TILE_LEN).upper()
                for strand, sq in (("f", seq), ("r", seq.translate(comp)[::-1])):
                    o.write(f"@c{i}_{s}_{strand}\n{sq}\n+\n{'I' * TILE_LEN}\n")
                    n_tiles += 1
    sam = os.path.join(work, "tiles.sam")
    t0 = time.monotonic()
    with open(sam, "w") as o, open(os.path.join(work, "bwa.stderr"), "w") as e:
        subprocess.run([BWA, *BWA_ARGS, "-t", str(threads), REF, fq], stdout=o, stderr=e, check=True)
    bwa_s = time.monotonic() - t0
    ok = collections.Counter()
    tot = collections.Counter()
    with pysam.AlignmentFile(sam) as f:
        for a in f:
            if a.is_secondary or a.is_supplementary:
                continue
            i, s, _ = a.query_name.split("_")
            i, s = int(i[1:]), int(s)
            tot[i] += 1
            if (not a.is_unmapped and a.reference_name == rows[i]["chrom"]
                    and abs(a.reference_start - s) <= START_TOL and a.mapping_quality >= MAPQ_OK):
                ok[i] += 1
    for i, r in enumerate(rows):
        r["mappability"] = round(ok[i] / tot[i], 4) if tot[i] else None
        r["tiles"] = tot[i]
    rec["mappability"] = {"tiles": n_tiles, "bwa_seconds": round(bwa_s, 1),
                          "bwa_command": " ".join([tilde(BWA), *BWA_ARGS, "-t", str(threads),
                                                   tilde(REF), "tiles.fq"])}
    write_json(path, rec)
    print(f"{n_tiles} tiles aligned in {bwa_s:.0f}s; mappability added to {len(rows)} candidates")


def step_summary(out):
    rec = json.load(open(os.path.join(out, "candidates.json")))
    rows = rec["candidates"]
    keys = ["low_mapq_fraction_200", "low_mapq_fraction_500", "mappability", "entropy2",
            "max_homopolymer", "max_dinucleotide_run", "mean_depth_200"]
    qs = [0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1]
    print(f"{len(rows)} candidates")
    print(f"{'statistic':24s}" + "".join(f"{'q' + str(q):>9s}" for q in qs))
    for k in keys:
        v = sorted(r[k] for r in rows if r.get(k) is not None)
        if not v:
            continue
        print(f"{k:24s}" + "".join(f"{v[min(len(v) - 1, int(q * (len(v) - 1)))]:>9}" for q in qs))
    for thr in (0.05, 0.1, 0.2, 0.3, 0.4, 0.5):
        print(f"  low_mapq_fraction_200 > {thr}: " + ", ".join(
            f"{c} {sum(1 for r in rows if r['chrom'] == c and r['low_mapq_fraction_200'] > thr)}"
            for c in CONTIGS))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["stats", "mappability", "summary"])
    ap.add_argument("--out", default=os.path.join(SIM, "scan"))
    ap.add_argument("--threads", type=int, default=10)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    {"stats": lambda: step_stats(a.out), "mappability": lambda: step_mappability(a.out, a.threads),
     "summary": lambda: step_summary(a.out)}[a.step]()
