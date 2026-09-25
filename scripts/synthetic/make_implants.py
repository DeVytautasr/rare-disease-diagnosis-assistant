#!/usr/bin/env python3
"""Implant heterozygous balanced reciprocal t(20;21) events into the NA12878
chr20+chr21 slice -- one implant per BAM. Rebuilt 2026-09-25 from the approved
Phase 6 design (the Phase 6 generator was never committed and is lost). This is
a NEW experiment: new random draws, an ALT-aware index, partly new breakpoints.

    make_implants.py prepare [--selection F] [--profile HS25]   per implant: removal sets,
                                                                junction reads (ART), FASTQs
    make_implants.py align   [--threads N]                      ONE bwa call for every implant
    make_implants.py build   IMPxx [IMPyy ...]                  filtered background + simulated
                                                                reads -> ~/public_data/sim/bams
    make_implants.py truth                                      ground truth for every implant

THE EVENT. For breakpoints A (chr20) and B (chr21), 1-based:
    der(20) = chr20[1..A] + chr21[B+1..]      junction J20: chr20:A | chr21:B+1
    der(21) = chr21[1..B] + chr20[A+1..]      junction J21: chr21:B | chr20:A+1
One homolog of each chromosome carries it (heterozygous); nothing is gained or
lost (balanced).

DEPTH COMPENSATION, so copy number is conserved as a balanced event requires.
At each breakpoint X, a background read pair SPANS the boundary after base X if
its fragment covers both X and X+1: for a proper pair with both primary mates
within +/-3 kb, the fragment runs from the leftmost mate start to the rightmost
mate end; otherwise any primary alignment of the pair that itself crosses the
boundary counts. Exactly floor(n/2) of the n spanning pairs are removed, drawn
with a seeded RNG from the sorted names -- the removed homolog's share. Removal
is by read name over a full-BAM stream (samtools view -N ^FILE), so mates,
secondary and supplementary records go too. Non-spanning pairs are untouched.
Simulated junction-spanning fragments are added to match the NON-duplicate
pairs removed: n_add = nondup_removed(A) + nondup_removed(B), split between the
two junctions as round(n_add/2) to J20 and the rest to J21. (Each junction
fragment covers one side of A and one side of B, so an even split keeps every
side within |removed(A) - removed(B)| / 2 of balance.) Simulated reads carry no
duplicate flag.

SIMULATION. For each junction a 4 kb template: 2 kb of reference on each side
(J20 = chr20[A-2000..A] + chr21[B+1..B+2000]; J21 = chr21[B-2000..B] +
chr20[A+1..A+2000]). ART_Illumina 2.5.8, profile HS25, paired, 150 bp,
fragment mean/SD as measured on the background, fold 60 (doubled, up to 3
times, if too few spanning fragments result), one recorded seed per template.
Fragments whose extent crosses the junction are the candidates; the required
number is drawn with a seeded RNG. ART's -sam output is kept as per-read ground
truth. Base qualities are requantized to the background's measured bins by
nearest bin, ties to the lower bin. Reads are named IMPxx_J20_nnnn /
IMPxx_J21_nnnn.

ALIGNMENT. All implants' reads in ONE bwa 0.7.15-r1140 mem -Y -K 100000000
call against the full hs38DH index with its .alt (ALT-aware, as the background
was aligned), read group @RG ID:SIM SM:NA12878 LB:NA12878 PL:illumina, then
split by the name prefix. -K fixes the batch size, so the output does not
depend on the thread count.

SEEDS. Every random draw has its own seed, derived from MASTER_SEED and the
draw's name (sha256), and every seed is recorded in the ground truth.
"""
import argparse
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
import time

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from synth_common import (ART, BAMS, BG, BG_BYTES, BWA, BWA_ARGS, READ_GROUP, REF,  # noqa: E402
                          SAMTOOLS, SIM, WORK, tilde, tool_versions, write_json)
from stage1_igv_assistant.tools import bam_tools  # noqa: E402

MASTER_SEED = 20260925
HALF = 2000             # junction template: this much reference on each side
SPAN_WINDOW = 3000      # spanning-pair search, +/- this around each breakpoint
FOLD_START = 60
MEASURE = os.path.join(SIM, "measure")
SELECTION = os.path.join(SIM, "selection.json")


def seed_for(name):
    """A 31-bit seed for one named random draw, derived from MASTER_SEED."""
    return int(hashlib.sha256(f"{MASTER_SEED}:{name}".encode()).hexdigest()[:8], 16) & 0x7FFFFFFF


def requantize(quals, bins):
    """Nearest bin; a tie goes to the lower bin."""
    out = []
    for q in quals:
        out.append(min(bins, key=lambda b: (abs(b - q), b)))
    return out


def spanning_pairs(bam, chrom, x):
    """{qname: is_duplicate} for pairs whose fragment covers 1-based x and x+1."""
    frags = {}
    for r in bam.fetch(chrom, max(0, x - SPAN_WINDOW), x + SPAN_WINDOW):
        if r.is_secondary or r.is_supplementary:
            continue
        rec = frags.setdefault(r.query_name, {"reads": [], "dup": False})
        rec["dup"] |= r.is_duplicate
        rec["reads"].append(r)
    span = {}
    for q, rec in frags.items():
        mapped = [r for r in rec["reads"] if not r.is_unmapped and r.reference_name == chrom]
        pair = [r for r in mapped if r.is_proper_pair]
        crosses = False
        if len(pair) == 2:
            s = min(r.reference_start for r in pair)
            e = max(r.reference_end for r in pair)
            crosses = s < x < e
        if not crosses:
            crosses = any(r.reference_start < x < r.reference_end for r in mapped)
        if crosses:
            span[q] = rec["dup"]
    return span


def art_truth(sam_path, junction):
    """ART's ground truth: {pair: {1: (start, end, strand), 2: ...}} in template coords."""
    pairs = {}
    with pysam.AlignmentFile(sam_path, check_sq=False) as f:
        for r in f:
            mate = 1 if r.is_read1 else 2
            pairs.setdefault(r.query_name, {})[mate] = (r.reference_start, r.reference_end,
                                                         "-" if r.is_reverse else "+")
    return pairs


def read_fastq(path):
    out = {}
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq, _, qual = f.readline().rstrip("\n"), f.readline(), f.readline().rstrip("\n")
            out[h[1:].strip().rsplit("/", 1)[0]] = (seq, qual)
    return out


def simulate_junction(imp, jname, template, n_needed, mean, sd, profile, bins, wd):
    """ART on one junction template; returns (reads for the FASTQs, truth, record)."""
    fa = os.path.join(wd, f"{jname}.fa")
    with open(fa, "w") as o:
        o.write(f">{imp}_{jname}\n{template}\n")
    seed = seed_for(f"{imp}:{jname}:art")
    fold = FOLD_START
    for attempt in range(4):
        pre = os.path.join(wd, f"{jname}.art")
        subprocess.run([ART, "-ss", profile, "-sam", "-na", "-q", "-i", fa, "-p", "-l", "150",
                        "-f", str(fold), "-m", str(round(mean)), "-s", str(round(sd)),
                        "-rs", str(seed), "-o", pre], check=True, stdout=subprocess.DEVNULL)
        truth = art_truth(pre + ".sam", HALF)
        spanning = sorted(q for q, m in truth.items() if 1 in m and 2 in m and
                          min(m[1][0], m[2][0]) < HALF < max(m[1][1], m[2][1]))
        if len(spanning) >= n_needed:
            break
        fold *= 2
    else:
        raise SystemExit(f"{imp} {jname}: only {len(spanning)} spanning fragments at fold {fold}")
    pick_seed = seed_for(f"{imp}:{jname}:pick")
    chosen = sorted(random.Random(pick_seed).sample(spanning, n_needed),
                    key=lambda q: int(q.rsplit("-", 1)[1]))
    r1, r2 = read_fastq(pre + "1.fq"), read_fastq(pre + "2.fq")
    reads, truth_rows = [], []
    for k, q in enumerate(chosen):
        name = f"{imp}_{jname}_{k:04d}"
        s1, q1 = r1[q]
        s2, q2 = r2[q]
        q1 = pysam.qualities_to_qualitystring(requantize(pysam.qualitystring_to_array(q1), bins))
        q2 = pysam.qualities_to_qualitystring(requantize(pysam.qualitystring_to_array(q2), bins))
        reads.append((name, s1, q1, s2, q2))
        m = truth[q]
        for mate in (1, 2):
            s, e, strand = m[mate]
            truth_rows.append({"read": name, "mate": mate, "art_name": q, "template_start": s,
                               "template_end": e, "strand": strand,
                               "spans_junction": s < HALF < e,
                               "bases_left": max(0, min(e, HALF) - s),
                               "bases_right": max(0, e - max(s, HALF))})
    rec = {"template_length": len(template), "junction_index": HALF, "art_seed": seed,
           "art_fold": fold, "art_attempts": attempt + 1, "spanning_fragments_available": len(spanning),
           "pick_seed": pick_seed, "fragments_added": n_needed,
           "reads_spanning_junction": sum(1 for t in truth_rows if t["spans_junction"])}
    return reads, truth_rows, rec


def step_prepare(selection, profile):
    sel = json.load(open(selection))
    bg_meas = json.load(open(os.path.join(MEASURE, "background.json")))
    mean, sd = bg_meas["fragment_length"]["mean"], bg_meas["fragment_length"]["sd"]
    bins = bg_meas["quality_bins"]["distinct"]
    if os.path.getsize(BG) != BG_BYTES:
        raise SystemExit("the background slice is not the recorded size")
    bam = pysam.AlignmentFile(BG)
    fa = pysam.FastaFile(REF)
    for imp in sel["implants"]:
        iid, A, B = imp["id"], imp["chr20"], imp["chr21"]
        wd = os.path.join(WORK, iid)
        if os.path.isdir(wd):
            shutil.rmtree(wd)
        os.makedirs(wd)
        removal = {}
        for side, chrom, x in (("A", "chr20", A), ("B", "chr21", B)):
            span = spanning_pairs(bam, chrom, x)
            names = sorted(span)
            rseed = seed_for(f"{iid}:{side}:remove")
            removed = sorted(random.Random(rseed).sample(names, len(names) // 2))
            removal[side] = {"chrom": chrom, "pos": x, "spanning_pairs": len(names),
                             "spanning_duplicates": sum(span[q] for q in names),
                             "removed": len(removed), "removed_duplicates": sum(span[q] for q in removed),
                             "removed_nondup": sum(not span[q] for q in removed), "remove_seed": rseed,
                             "_names": removed}
        n_add = removal["A"]["removed_nondup"] + removal["B"]["removed_nondup"]
        n20 = round(n_add / 2)
        n21 = n_add - n20
        t20 = fa.fetch("chr20", A - HALF, A) + fa.fetch("chr21", B, B + HALF)
        t21 = fa.fetch("chr21", B - HALF, B) + fa.fetch("chr20", A, A + HALF)
        reads20, truth20, rec20 = simulate_junction(iid, "J20", t20.upper(), n20, mean, sd, profile, bins, wd)
        reads21, truth21, rec21 = simulate_junction(iid, "J21", t21.upper(), n21, mean, sd, profile, bins, wd)
        with open(os.path.join(wd, "removed_qnames.txt"), "w") as o:
            for q in sorted(set(removal["A"]["_names"]) | set(removal["B"]["_names"])):
                o.write(q + "\n")
        with open(os.path.join(wd, "R1.fq"), "w") as o1, open(os.path.join(wd, "R2.fq"), "w") as o2:
            for name, s1, q1, s2, q2 in reads20 + reads21:
                o1.write(f"@{name}/1\n{s1}\n+\n{q1}\n")
                o2.write(f"@{name}/2\n{s2}\n+\n{q2}\n")
        for side in removal:
            removal[side].pop("_names")
        prep = {"id": iid, "class": imp["class"], "chr20": A, "chr21": B, "removal": removal,
                "fragments_added_total": n_add, "junctions": {"J20": rec20, "J21": rec21},
                "profile": profile, "fragment_mean": round(mean), "fragment_sd": round(sd),
                "quality_bins": bins, "truth_reads": truth20 + truth21}
        write_json(os.path.join(wd, "prepare.json"), prep)
        print(f"{iid} {imp['class']:16s} chr20:{A:,} chr21:{B:,} | spanning A {removal['A']['spanning_pairs']}"
              f" B {removal['B']['spanning_pairs']} | removed A {removal['A']['removed']}"
              f" ({removal['A']['removed_nondup']} non-dup) B {removal['B']['removed']}"
              f" ({removal['B']['removed_nondup']} non-dup) | added J20 {n20} J21 {n21} | "
              f"reads spanning a junction {rec20['reads_spanning_junction'] + rec21['reads_spanning_junction']}",
              flush=True)


def step_align(threads):
    ids = sorted(d for d in os.listdir(WORK) if d.startswith("IMP"))
    pool = os.path.join(WORK, "pool")
    os.makedirs(pool, exist_ok=True)
    r1, r2 = os.path.join(pool, "R1.fq"), os.path.join(pool, "R2.fq")
    with open(r1, "w") as o1, open(r2, "w") as o2:
        for i in ids:
            o1.write(open(os.path.join(WORK, i, "R1.fq")).read())
            o2.write(open(os.path.join(WORK, i, "R2.fq")).read())
    sam = os.path.join(pool, "pool.sam")
    cmd = [BWA, *BWA_ARGS, "-t", str(threads), "-R", READ_GROUP, REF, r1, r2]
    t0 = time.monotonic()
    with open(sam, "w") as o, open(os.path.join(pool, "bwa.stderr"), "w") as e:
        subprocess.run(cmd, stdout=o, stderr=e, check=True)
    secs = round(time.monotonic() - t0, 1)
    with pysam.AlignmentFile(sam) as f:
        outs = {i: pysam.AlignmentFile(os.path.join(WORK, i, "sim.unsorted.bam"), "wb", template=f)
                for i in ids}
        counts = {i: 0 for i in ids}
        for r in f:
            i = r.query_name.split("_", 1)[0]
            outs[i].write(r)
            counts[i] += 1
        for o in outs.values():
            o.close()
    for i in ids:
        d = os.path.join(WORK, i)
        subprocess.run([SAMTOOLS, "sort", "-o", os.path.join(d, "sim.bam"),
                        os.path.join(d, "sim.unsorted.bam")], check=True)
        subprocess.run([SAMTOOLS, "index", os.path.join(d, "sim.bam")], check=True)
        os.remove(os.path.join(d, "sim.unsorted.bam"))
    rec = {"command": " ".join(tilde(str(c)) for c in cmd).replace("\t", "\\t"),
           "threads": threads, "wall_seconds": secs, "records_per_implant": counts}
    write_json(os.path.join(pool, "align.json"), rec)
    print(json.dumps(rec, indent=1))


def count_records(path, *extra):
    p = subprocess.run([SAMTOOLS, "view", "-c", *extra, path], check=True, capture_output=True, text=True)
    return int(p.stdout.strip())


def step_build(ids):
    os.makedirs(BAMS, exist_ok=True)
    for iid in ids:
        d = os.path.join(WORK, iid)
        t0 = time.monotonic()
        removed = os.path.join(d, "removed_qnames.txt")
        filt = os.path.join(d, "bg_filtered.bam")
        out = os.path.join(BAMS, f"{iid}.bam")
        subprocess.run([SAMTOOLS, "view", "-@", "2", "-b", "-N", "^" + removed, "-o", filt, BG], check=True)
        subprocess.run([SAMTOOLS, "merge", "-@", "2", "-f", "-o", out, filt, os.path.join(d, "sim.bam")],
                       check=True)
        subprocess.run([SAMTOOLS, "index", out], check=True)
        n_bg = count_records(BG)
        n_removed = count_records(BG, "-N", removed)
        n_sim = count_records(os.path.join(d, "sim.bam"))
        n_out = count_records(out)
        with pysam.AlignmentFile(out) as f:
            rgs = f.header.to_dict().get("RG", [])
        rec = {"bam": tilde(out), "bytes": os.path.getsize(out), "background_records": n_bg,
               "removed_records": n_removed, "simulated_records": n_sim, "output_records": n_out,
               "record_arithmetic_holds": n_out == n_bg - n_removed + n_sim,
               "read_groups": len(rgs), "samples_in_read_groups": sorted({r.get("SM") for r in rgs}),
               "libraries_in_read_groups": sorted({r.get("LB") for r in rgs}),
               "wall_seconds": round(time.monotonic() - t0, 1)}
        os.remove(filt)
        write_json(os.path.join(d, "build.json"), rec)
        print(f"{iid}: {rec['output_records']:,} records = {n_bg:,} - {n_removed:,} + {n_sim:,} "
              f"({'holds' if rec['record_arithmetic_holds'] else 'DOES NOT HOLD'}); read groups "
              f"{rec['read_groups']}, samples {rec['samples_in_read_groups']}, "
              f"libraries {rec['libraries_in_read_groups']}; {rec['wall_seconds']}s", flush=True)
        if not rec["record_arithmetic_holds"] or rec["samples_in_read_groups"] != ["NA12878"]:
            raise SystemExit(f"{iid}: build check failed")


def locus(bam_path, chrom, x):
    s200 = bam_tools.get_bam_stats_at_locus(bam_path, chrom, x - 200, x + 200)
    s500 = bam_tools.get_bam_stats_at_locus(bam_path, chrom, x - 500, x + 500)
    return {"mean_depth_500": s500["mean_depth"], "low_mapq_fraction_500": s500["low_mapq_fraction"],
            "low_mapq_fraction_200": s200["low_mapq_fraction"], "reads_500": s500["total_reads"]}


def step_truth(out_path):
    sel = json.load(open(SELECTION))
    implants = []
    for imp in sel["implants"]:
        iid, A, B = imp["id"], imp["chr20"], imp["chr21"]
        prep = json.load(open(os.path.join(WORK, iid, "prepare.json")))
        implants.append({
            "id": iid, "class": imp["class"], "coordinate_source": imp.get("source"),
            "breakpoints": {"chr20": A, "chr21": B, "convention": "1-based; each junction lies "
                            "between the given base and the next"},
            "junctions": {
                "J20": {"on": "der(20)", "left": f"chr20:{A}", "right": f"chr21:{B + 1}",
                        "orientation": "3' end of chr20 (after base A) joined to the 5' end of "
                                       "chr21 (from base B+1)",
                        "delly_ct_chr20_first": "3to5",
                        "vcf_bnd": {f"chr20:{A}": f"N[chr21:{B + 1}[", f"chr21:{B + 1}": f"]chr20:{A}]N"}},
                "J21": {"on": "der(21)", "left": f"chr21:{B}", "right": f"chr20:{A + 1}",
                        "orientation": "3' end of chr21 (after base B) joined to the 5' end of "
                                       "chr20 (from base A+1)",
                        "delly_ct_chr20_first": "5to3",
                        "vcf_bnd": {f"chr21:{B}": f"N[chr20:{A + 1}[", f"chr20:{A + 1}": f"]chr21:{B}]N"}}},
            "zygosity": "heterozygous", "balanced": True,
            "pre_implant": {"chr20": locus(BG, "chr20", A), "chr21": locus(BG, "chr21", B)},
            "scan_statistics": imp.get("stats"),
            "compensation": {"removal": prep["removal"], "fragments_added_total": prep["fragments_added_total"],
                             "J20": prep["junctions"]["J20"], "J21": prep["junctions"]["J21"]},
            "simulation": {"profile": prep["profile"], "read_length": 150, "fragment_mean": prep["fragment_mean"],
                           "fragment_sd": prep["fragment_sd"], "quality_bins": prep["quality_bins"],
                           "quality_rule": "nearest bin, ties to the lower bin"},
            "files": {"bam": tilde(os.path.join(BAMS, f"{iid}.bam")),
                      "art_sam_ground_truth": [tilde(os.path.join(WORK, iid, f"{j}.art.sam")) for j in ("J20", "J21")]},
        })
    rec = {"experiment": "synthetic positive control, rebuilt 2026-09-25 (new draws; not a reproduction)",
           "master_seed": MASTER_SEED, "seed_derivation": "int(sha256(f'{master}:{draw name}')[:8], 16) & 0x7FFFFFFF",
           "background": {"path": tilde(BG), "bytes": BG_BYTES}, "selection": sel.get("thresholds"),
           "versions": tool_versions(), "implants": implants}
    write_json(out_path, rec)
    print(f"ground truth for {len(implants)} implants -> {tilde(out_path)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["prepare", "align", "build", "truth"])
    ap.add_argument("ids", nargs="*")
    ap.add_argument("--selection", default=SELECTION)
    ap.add_argument("--profile", default="HS25")
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--out")
    a = ap.parse_args()
    if a.step == "prepare":
        step_prepare(a.selection, a.profile)
    elif a.step == "align":
        step_align(a.threads)
    elif a.step == "build":
        step_build(a.ids)
    else:
        step_truth(a.out)
