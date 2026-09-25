#!/usr/bin/env python3
"""Attribution at a reused coordinate (Task 2f of the Phase 6 rerun): realign the
simulated reads WITHOUT the ALT index, rebuild the implant BAM, rerun delly, and
see whether ALT-aware alignment or the random draw explains a detection that
differs from Phase 6.

    realign_noalt.py align             every implant's simulated reads, pooled exactly as
                                       make_implants.py align pooled them, twice:
                                         WITH the .alt  -- a control that must reproduce the
                                                           committed pool alignment record
                                                           for record, or nothing is attributed
                                         WITHOUT it     -- the experiment
    realign_noalt.py build IMPxx ...   background minus the removed pairs + the no-ALT reads
    realign_noalt.py delly IMPxx ...   delly with the committed defaults (-q 1 -r 20)
    realign_noalt.py compare IMPxx ... detection and the reads, ALT-aware against no-ALT

WITHOUT THE .alt. bwa mem runs ALT-aware only when <prefix>.alt exists. A prefix
of symbolic links to the five index files, with no .alt beside them, is the same
index without it. Checked, not assumed: bwa marks ALT contigs with AH:* on their
@SQ lines only when it loaded an .alt, so the control's header must carry them
(3,171, as the committed pool.sam does) and the no-ALT header none.

POOLED, as the original was: bwa mem estimates the insert-size distribution per
batch, so the same reads aligned in a different batch could align differently for
a reason unrelated to the .alt. The command is the committed one (align.json)
with only the index prefix and the thread count changed; -K 100000000 makes the
output independent of the thread count, which the control verifies.

Outputs under ~/public_data/sim/noalt/ (not autodiscovered by the interface):
index/ (links), control/pool.sam, pool.sam, work/IMPxx/sim.bam, bams/IMPxx.bam,
delly/IMPxx.bcf. Records under results/.../analysis_2026-09-25/noalt/.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
from synth_common import (BAMS, BG, BWA, BWA_ARGS, DELLY_DIR, EXCL, READ_GROUP, REF, SAMTOOLS,  # noqa: E402
                          SIM, WORK)
import evidence_chain as ec  # noqa: E402
import delly_ladder as dl  # noqa: E402

NOALT = os.path.join(SIM, "noalt")
INDEX = os.path.join(NOALT, "index", os.path.basename(REF))
INDEX_EXT = (".amb", ".ann", ".bwt", ".pac", ".sa")
POOL = os.path.join(WORK, "pool")
THREADS = 8


def alignment_lines(sam):
    with open(sam) as f:
        return [l for l in f if not l.startswith("@")]


def header_counts(sam):
    sq = ah = 0
    with open(sam) as f:
        for l in f:
            if not l.startswith("@"):
                break
            if l.startswith("@SQ"):
                sq += 1
                ah += "\tAH:" in l
    return {"sq": sq, "ah": ah}


def bwa(prefix, out, log):
    cmd = [BWA, *BWA_ARGS, "-t", str(THREADS), "-R", READ_GROUP, prefix,
           os.path.join(POOL, "R1.fq"), os.path.join(POOL, "R2.fq")]
    t0 = time.monotonic()
    with open(out, "w") as o, open(log, "w") as e:
        subprocess.run(cmd, stdout=o, stderr=e, check=True)
    return {"command": " ".join(cmd).replace(os.path.expanduser("~"), "~").replace("\t", "\\t"),
            "wall_seconds": round(time.monotonic() - t0, 1)}


def step_align():
    os.makedirs(os.path.dirname(INDEX), exist_ok=True)
    for ext in INDEX_EXT:
        link = INDEX + ext
        if not os.path.islink(link):
            os.symlink(REF + ext, link)
    if os.path.exists(INDEX + ".alt"):
        raise SystemExit("the no-ALT prefix has an .alt beside it")
    os.makedirs(os.path.join(NOALT, "control"), exist_ok=True)
    ctl_sam = os.path.join(NOALT, "control", "pool.sam")
    ctl = bwa(REF, ctl_sam, os.path.join(NOALT, "control", "bwa.stderr"))
    committed = alignment_lines(os.path.join(POOL, "pool.sam"))
    control_identical = alignment_lines(ctl_sam) == committed
    sam = os.path.join(NOALT, "pool.sam")
    exp = bwa(INDEX, sam, os.path.join(NOALT, "bwa.stderr"))
    heads = {"committed": header_counts(os.path.join(POOL, "pool.sam")), "control": header_counts(ctl_sam),
             "no_alt": header_counts(sam)}
    # bwa reports "read N ALT contigs" in both modes; the number is what counts. (The
    # first version of this check counted lines naming ALT contigs, so "read 0 ALT
    # contigs" failed it -- a defect in the check, recorded in align_first_check_defective.json.)
    def alt_read(path):
        import re as _re
        m = _re.search(r"read (\d+) ALT contigs", open(path).read())
        return int(m.group(1)) if m else None
    stderr_alt = alt_read(os.path.join(NOALT, "bwa.stderr"))
    ctl_alt = alt_read(os.path.join(NOALT, "control", "bwa.stderr"))
    ok = (control_identical and heads["control"]["ah"] == heads["committed"]["ah"] > 0
          and heads["no_alt"]["ah"] == 0 and stderr_alt == 0 and ctl_alt == heads["committed"]["ah"])
    # split by implant, as make_implants.py align did
    ids = sorted(d for d in os.listdir(WORK) if d.startswith("IMP"))
    counts = {i: 0 for i in ids}
    with pysam.AlignmentFile(sam) as f:
        outs = {}
        for i in ids:
            os.makedirs(os.path.join(NOALT, "work", i), exist_ok=True)
            outs[i] = pysam.AlignmentFile(os.path.join(NOALT, "work", i, "sim.unsorted.bam"), "wb", template=f)
        for r in f:
            i = r.query_name.split("_", 1)[0]
            outs[i].write(r)
            counts[i] += 1
        for o in outs.values():
            o.close()
    for i in ids:
        d = os.path.join(NOALT, "work", i)
        subprocess.run([SAMTOOLS, "sort", "-o", os.path.join(d, "sim.bam"), os.path.join(d, "sim.unsorted.bam")],
                       check=True)
        subprocess.run([SAMTOOLS, "index", os.path.join(d, "sim.bam")], check=True)
        os.remove(os.path.join(d, "sim.unsorted.bam"))
    committed_counts = json.load(open(os.path.join(POOL, "align.json")))["records_per_implant"]
    rec = {"control_with_alt": dict(ctl, alignment_records_identical_to_committed=control_identical, alt_contigs_reported_by_bwa=ctl_alt),
           "no_alt": dict(exp, ah_lines=heads["no_alt"]["ah"], alt_contigs_reported_by_bwa=stderr_alt),
           "headers": heads, "records_per_implant_no_alt": counts,
           "records_per_implant_committed": committed_counts, "verified": ok}
    ec.write_record(os.path.join(ec.OUT, "noalt", "align.json"), rec)
    print(json.dumps(rec, indent=1))
    return 0 if ok else 1


def count_records(path, *extra):
    p = subprocess.run([SAMTOOLS, "view", "-c", *extra, path], check=True, capture_output=True, text=True)
    return int(p.stdout.strip())


def step_build(ids):
    os.makedirs(os.path.join(NOALT, "bams"), exist_ok=True)
    for iid in ids:
        t0 = time.monotonic()
        removed = os.path.join(WORK, iid, "removed_qnames.txt")
        sim = os.path.join(NOALT, "work", iid, "sim.bam")
        filt = os.path.join(NOALT, "work", iid, "bg_filtered.bam")
        out = os.path.join(NOALT, "bams", f"{iid}.bam")
        subprocess.run([SAMTOOLS, "view", "-@", "2", "-b", "-N", "^" + removed, "-o", filt, BG], check=True)
        subprocess.run([SAMTOOLS, "merge", "-@", "2", "-f", "-o", out, filt, sim], check=True)
        subprocess.run([SAMTOOLS, "index", out], check=True)
        n_bg, n_rm, n_sim, n_out = (count_records(BG), count_records(BG, "-N", removed), count_records(sim),
                                    count_records(out))
        os.remove(filt)
        rec = {"bam": out, "background_records": n_bg, "removed_records": n_rm, "simulated_records": n_sim,
               "output_records": n_out, "record_arithmetic_holds": n_out == n_bg - n_rm + n_sim,
               "wall_seconds": round(time.monotonic() - t0, 1)}
        ec.write_record(os.path.join(ec.OUT, "noalt", f"build_{iid}.json"), rec)
        print(f"{iid}: {rec}", flush=True)
        if not rec["record_arithmetic_holds"]:
            return 1
    return 0


def step_delly(ids):
    if not dl.preflight():
        return 3
    jobs = []
    for iid in ids:
        jobs.append((iid, *dl.launch(iid, os.path.join(NOALT, "bams", f"{iid}.bam"),
                                     os.path.join(NOALT, "delly"), 1, 20, os.path.join(NOALT, "logs"))))
    rc_all = 0
    for iid, p, out, log, fh in jobs:
        rc = p.wait()
        fh.close()
        s = dl.job_summary(out, log, rc)
        ec.write_record(os.path.join(ec.OUT, "noalt", f"delly_{iid}.json"), s)
        print(f"{iid}: {s}", flush=True)
        rc_all |= rc != 0
    return rc_all


def primaries(sim_bam):
    out = {}
    with pysam.AlignmentFile(sim_bam) as f:
        for r in f:
            if r.is_secondary or r.is_supplementary:
                continue
            sa = r.get_tag("SA").rstrip(";").split(";") if r.has_tag("SA") else []
            out[(r.query_name, 1 if r.is_read1 else 2)] = {
                "chrom": None if r.is_unmapped else r.reference_name,
                "pos": None if r.is_unmapped else r.reference_start + 1, "mapq": r.mapping_quality,
                "sa": [e.split(",")[0] + ":" + e.split(",")[1] + ":mq" + e.split(",")[4] for e in sa if e],
                "mate_chrom": r.next_reference_name, "pa": r.get_tag("pa") if r.has_tag("pa") else None}
    return out


def step_compare(ids):
    T = ec.truth()
    tol = ec.tolerance()
    rec = ec.recorder()
    out = {}
    for iid in ids:
        det = {}
        for name, bcf in (("alt_aware", os.path.join(DELLY_DIR, f"{iid}.bcf")),
                          ("no_alt", os.path.join(NOALT, "delly", f"{iid}.bcf"))):
            ld = ec.result(rec.call("bridge", "load_candidate_set", {"path": bcf, "label": f"{iid}:{name}"}))
            unf = ec.result(rec.call("bridge", "list_candidates", {"set_id": ld["set_id"], "limit": ec.LIMIT}))
            fun = ec.result(rec.call("bridge", "list_candidates",
                                     {"set_id": ld["set_id"], "filter_pass": True, "min_pe": 3, "min_sr": 1,
                                      "primary_only": True, "mask_path": EXCL, "limit": ec.LIMIT}))
            surv = {c["candidate_id"] for c in fun["candidates"]}
            det[name] = {}
            for jn, exp in ec.expected_junctions(T[iid]).items():
                m, _ = ec.best_match(unf["candidates"], exp, tol)
                det[name][jn] = {"called": m is not None, "survives_funnel": bool(m and m["candidate_id"] in surv),
                                 "pe": m and m["pe"], "sr": m and m["sr"], "filter": m and m["filter"]}
            det[name]["survivors"] = len(surv)
        a = primaries(os.path.join(WORK, iid, "sim.bam"))
        n = primaries(os.path.join(NOALT, "work", iid, "sim.bam"))
        keys = sorted(set(a) & set(n))
        changed = {"position_or_contig": 0, "mapq": 0, "sa": 0,
                   "primary_in_only_one_alignment": len(set(a) ^ set(n))}
        diffs = []
        for k in keys:
            x, y = a.get(k), n.get(k)
            d = {}
            if (x["chrom"], x["pos"]) != (y["chrom"], y["pos"]):
                changed["position_or_contig"] += 1
                d["alignment"] = [f"{x['chrom']}:{x['pos']}", f"{y['chrom']}:{y['pos']}"]
            if x["mapq"] != y["mapq"]:
                changed["mapq"] += 1
                d["mapq"] = [x["mapq"], y["mapq"]]
            if x["sa"] != y["sa"]:
                changed["sa"] += 1
                d["sa"] = [x["sa"], y["sa"]]
            if d:
                diffs.append({"read": f"{k[0]}/{k[1]}", **d})

        def pairs_q20(p):
            """Junction pairs (mates on the two different chromosomes) with both mates MAPQ >= 20 --
            what delly's -r 20 asks of a translocation's paired-end evidence."""
            by = {}
            for (q, m), v in p.items():
                by.setdefault(q, {})[m] = v
            inter = [v for v in by.values() if len(v) == 2 and v[1]["chrom"] and v[2]["chrom"]
                     and v[1]["chrom"] != v[2]["chrom"]]
            return {"interchromosomal_pairs": len(inter),
                    "both_mates_mapq_ge_20": sum(1 for v in inter if v[1]["mapq"] >= 20 and v[2]["mapq"] >= 20),
                    "reads_with_pa_tag": sum(1 for v in p.values() if v["pa"] is not None)}
        same = all(det["alt_aware"][j]["survives_funnel"] == det["no_alt"][j]["survives_funnel"] for j in ("J20", "J21"))
        out[iid] = {"detection": det, "detection_identical": same,
                    "reads": {"primaries": len(keys), "changed": changed, "differences": diffs,
                              "alt_aware": pairs_q20(a), "no_alt": pairs_q20(n)},
                    "reading": ("ALT-aware alignment does not explain the difference from Phase 6: without the "
                                ".alt the implant is detected exactly as with it" if same else
                                "ALT-aware alignment changes detection here")}
        print(f"{iid}: detection ALT-aware {det['alt_aware']}; no-ALT {det['no_alt']}; identical {same}; "
              f"reads changed {changed}; {out[iid]['reads']['alt_aware']} -> {out[iid]['reads']['no_alt']}", flush=True)
    ec.write_record(os.path.join(ec.OUT, "noalt", "compare.json"), {"implants": out, "calls": rec.calls})
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["align", "build", "delly", "compare"])
    ap.add_argument("ids", nargs="*")
    a = ap.parse_args()
    bad = [i for i in a.ids if i not in ec.IMPLANTS]
    if bad or (a.step != "align" and not a.ids):
        raise SystemExit(f"give implant ids from {ec.IMPLANTS}")
    sys.exit({"align": lambda: step_align(), "build": lambda: step_build(a.ids),
              "delly": lambda: step_delly(a.ids), "compare": lambda: step_compare(a.ids)}[a.step]())
