#!/usr/bin/env python3
"""Per-read adjudication of the split-read and soft-clip layers at every
implanted breakpoint, against ART's ground truth (Task 3 of the Phase 6 rerun).

    adjudicate_reads.py run IMPxx     one implant (run the twelve concurrently)
    adjudicate_reads.py merge         the twelve records + the chain records -> the figures

DEFINITIONS, fixed before running
  truly spanning read   the gate's: a simulated read whose own ART template
                        interval crosses its junction with >= 1 base on each side
  breakpoints           each implant's chr20:A and chr21:B, queried with the tools'
                        own windows (split reads and soft clips +/-200)
  the tool's count      the split_reads MCP tool at its own default, min_mapq 0,
                        through ui.ToolRecorder; again at min_mapq 20 for (c)
  per-read identity     the reads the tool counted are enumerated here with pysam
                        under the tool's own inclusion rules (primary records only;
                        MAPQ filter on the read; SA entries filtered by their own
                        mapQ; counted when one entry survives). The enumeration must
                        reproduce the tool's count and its total_reads_in_window
                        EXACTLY at every breakpoint and both thresholds, or the
                        decomposition is refused (exit 1). The same for the
                        soft-clip tool (defaults: clip >= 10 bp, MAPQ >= 20,
                        supplementary records included, as the tool does).
  a counted split read is
      genuine                simulated, truly spanning, and an SA entry that survives
                             the SA mapQ filter lies on the partner chromosome within
                             500 bp of the partner breakpoint
      spanning, SA elsewhere simulated and truly spanning, without such an entry
      simulated non-spanning simulated, not truly spanning
      background             read group not SIM (background SA reads pointing at the
                             partner breakpoint are counted separately)
  delly's SR            INFO/SR of the junction in the default-run BCF that the chain
                        matched to each implanted junction; 0 when none matched
  truth                 the number of truly spanning reads
  false-positive rate   (counted - genuine) / counted, per breakpoint and pooled by
                        class; Pearson r against the low_mapq_fraction the
                        bam_stats_at_locus tool reports over +/-500 bp there
  short overhang (b)    truly spanning reads whose PRIMARY record carries no SA tag;
                        "caught" when that record is among the reads the soft-clip
                        tool counts at the breakpoint whose window holds it. As at
                        the gate (gate_a1.py): a read whose shorter side is under
                        30 bp but whose bases could still reach -T 30 with a 19 bp
                        seed against the partner reference is flagged.
"""
import argparse
import json
import math
import os
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
from synth_common import BAMS, REF, WORK  # noqa: E402
import evidence_chain as ec  # noqa: E402
import gate_a1 as ga  # noqa: E402

OUT = os.path.join(ec.OUT, "adjudication")
WIN = 200                 # the split_reads and soft_clipped_reads defaults
STATS_WIN = 500           # the quality gate's window
SUPPORT_TOL = 500
MIN_CLIP, CLIP_MAPQ = 10, 20
COMP = str.maketrans("ACGTN", "TGCAN")


def truly_spanning(prep):
    return {(t["read"], t["mate"]): t for t in prep["truth_reads"]
            if t["spans_junction"] and t["bases_left"] >= 1 and t["bases_right"] >= 1}


def _int(s):
    try:
        return int(s)
    except (TypeError, ValueError):
        return None


def window(f, chrom, x):
    return max(0, x - WIN), min(x + WIN, f.get_reference_length(chrom))


def enumerate_split(f, chrom, x, m, pchrom, p, ts, iid):
    start, end = window(f, chrom, x)
    total = below = 0
    counted = []
    for r in f.fetch(chrom, start, end):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if r.mapping_quality < m:
            below += 1
            continue
        total += 1
        if not r.has_tag("SA"):
            continue
        kept = []
        for e in r.get_tag("SA").rstrip(";").split(";"):
            fl = e.split(",")
            if len(fl) < 2:
                continue
            if len(fl) >= 5 and _int(fl[4]) is not None and _int(fl[4]) < m:
                continue
            kept.append(fl)
        if not kept:
            continue
        sim = r.has_tag("RG") and r.get_tag("RG") == "SIM"
        key = (r.query_name, 1 if r.is_read1 else 2)
        supports = any(fl[0] == pchrom and _int(fl[1]) is not None and abs(_int(fl[1]) - p) <= SUPPORT_TOL
                       for fl in kept)
        if sim:
            if not r.query_name.startswith(iid + "_J"):
                raise SystemExit(f"a SIM read named {r.query_name} in {iid}")
            cat = ("genuine" if key in ts and supports else
                   "spanning_sa_elsewhere" if key in ts else "simulated_nonspanning")
        else:
            cat = "background_sa_to_partner" if supports else "background"
        counted.append({"read": r.query_name if sim else None, "mate": key[1], "category": cat,
                        "mapq": r.mapping_quality,
                        "junction": r.query_name.split("_")[1] if sim else None,
                        "key": list(key)})
    return {"window": [start, end], "total_reads_in_window": total, "reads_below_min_mapq": below,
            "counted": counted}


def enumerate_clips(f, chrom, x):
    start, end = window(f, chrom, x)
    total = below = 0
    counted = []
    for r in f.fetch(chrom, start, end):
        if r.is_unmapped or r.is_secondary:
            continue
        if r.mapping_quality < CLIP_MAPQ:
            below += 1
            continue
        total += 1
        c = r.cigartuples
        if not c:
            continue
        side = ("left" if c[0][0] == 4 and c[0][1] >= MIN_CLIP else
                "right" if c[-1][0] == 4 and c[-1][1] >= MIN_CLIP else None)
        if side:
            counted.append({"key": [r.query_name, 1 if r.is_read1 else 2], "supplementary": r.is_supplementary,
                            "side": side, "sim": r.has_tag("RG") and r.get_tag("RG") == "SIM"})
    return {"window": [start, end], "total_reads_in_window": total, "reads_below_min_mapq": below,
            "counted": counted}


def sim_records(iid):
    """Every simulated record, from the implant's own sim.bam (merged unchanged
    into the implant BAM): primary records by (name, mate), and every record's
    signature for the consistency check against the implant BAM."""
    prim, sig = {}, {}
    with pysam.AlignmentFile(os.path.join(WORK, iid, "sim.bam")) as f:
        for r in f:
            k = (r.query_name, 1 if r.is_read1 else 2)
            sig[(k, r.flag, r.reference_name, r.reference_start)] = (
                r.cigarstring, r.get_tag("SA") if r.has_tag("SA") else None, r.mapping_quality)
            if not (r.is_secondary or r.is_supplementary):
                prim[k] = {"unmapped": r.is_unmapped, "chrom": r.reference_name,
                           "start": r.reference_start, "end": r.reference_end, "mapq": r.mapping_quality,
                           "sa": r.get_tag("SA") if r.has_tag("SA") else None, "cigar": r.cigarstring,
                           "seq": r.query_sequence, "is_reverse": r.is_reverse}
    return prim, sig


def step_run(iid):
    imp = ec.truth()[iid]
    A, B = imp["breakpoints"]["chr20"], imp["breakpoints"]["chr21"]
    bam = os.path.join(BAMS, f"{iid}.bam")
    prep = json.load(open(os.path.join(WORK, iid, "prepare.json")))
    ts = truly_spanning(prep)
    prim, sig = sim_records(iid)
    rec = ec.recorder()
    fa = pysam.FastaFile(REF)
    art = ga.art_records(iid)
    points, problems, sim_seen, sim_mismatch = [], [], 0, 0
    with pysam.AlignmentFile(bam) as f:
        for chrom, x, pchrom, p in (("chr20", A, "chr21", B), ("chr21", B, "chr20", A)):
            c = {"split_m0": rec.call("evidence", "split_reads", {"bam_path": bam, "chromosome": chrom, "position": x}),
                 "split_m20": rec.call("evidence", "split_reads",
                                       {"bam_path": bam, "chromosome": chrom, "position": x, "min_mapq": 20}),
                 "soft_clips": rec.call("evidence", "soft_clipped_reads",
                                        {"bam_path": bam, "chromosome": chrom, "position": x}),
                 "stats_500": rec.call("evidence", "bam_stats_at_locus",
                                       {"bam_path": bam, "chromosome": chrom,
                                        "start": max(0, x - STATS_WIN), "end": x + STATS_WIN})}
            r0, r20, rc = (ec.result(c["split_m0"]), ec.result(c["split_m20"]), ec.result(c["soft_clips"]))
            e0 = enumerate_split(f, chrom, x, 0, pchrom, p, ts, iid)
            e20 = enumerate_split(f, chrom, x, 20, pchrom, p, ts, iid)
            ec_ = enumerate_clips(f, chrom, x)
            checks = {
                "split_m0_count": (len(e0["counted"]), r0["split_reads"]),
                "split_m0_total": (e0["total_reads_in_window"], r0["total_reads_in_window"]),
                "split_m0_applied": (0, r0.get("min_mapq_applied")),
                "split_m20_count": (len(e20["counted"]), r20["split_reads"]),
                "split_m20_total": (e20["total_reads_in_window"], r20["total_reads_in_window"]),
                "split_m20_applied": (20, r20.get("min_mapq_applied")),
                "clips_count": (len(ec_["counted"]), rc["soft_clipped_reads"]),
                "clips_total": (ec_["total_reads_in_window"], rc["total_reads_in_window"]),
            }
            for k, (mine, tool) in checks.items():
                if mine != tool:
                    problems.append(f"{chrom}:{x} {k}: enumerated {mine}, tool {tool}")
            # simulated records met in the implant BAM must be sim.bam's, unchanged
            for r in f.fetch(chrom, *window(f, chrom, x)):
                if r.has_tag("RG") and r.get_tag("RG") == "SIM":
                    sim_seen += 1
                    k = ((r.query_name, 1 if r.is_read1 else 2), r.flag, r.reference_name, r.reference_start)
                    want = sig.get(k)
                    got = (r.cigarstring, r.get_tag("SA") if r.has_tag("SA") else None, r.mapping_quality)
                    if want != got:
                        sim_mismatch += 1
            start, end = e0["window"]
            in_window = sorted([list(k) for k, v in prim.items()
                                if k in ts and not v["unmapped"] and v["chrom"] == chrom
                                and v["start"] < end and v["end"] > start])
            points.append({"chrom": chrom, "position": x, "partner": f"{pchrom}:{p}",
                           "calls": {k: v["id"] for k, v in c.items()},
                           "low_mapq_fraction_500": ec.result(c["stats_500"])["low_mapq_fraction"],
                           "checks": {k: {"enumerated": a, "tool": b, "equal": a == b} for k, (a, b) in checks.items()},
                           "split_m0": e0, "split_m20": e20, "soft_clips": ec_,
                           "truly_spanning_primary_in_window": in_window})
    # every truly spanning read
    reads = []
    for k, t in sorted(ts.items()):
        pr = prim.get(k)
        row = {"read": k[0], "mate": k[1], "junction": k[0].split("_")[1],
               "bases_left": t["bases_left"], "bases_right": t["bases_right"],
               "shorter_side": min(t["bases_left"], t["bases_right"])}
        if pr is None or pr["unmapped"]:
            row["primary"] = "unmapped" if pr else "not found"
            reads.append(row)
            continue
        partner = "chr21" if pr["chrom"] == "chr20" else "chr20"
        sa = [e.split(",") for e in (pr["sa"] or "").rstrip(";").split(";") if e]
        row.update({"primary": f"{pr['chrom']}:{pr['start'] + 1}", "mapq": pr["mapq"], "has_sa": bool(sa),
                    "sa_to_partner": any(e[0] == partner for e in sa)})
        where = next((pt for pt in points if pt["chrom"] == pr["chrom"]
                      and [k[0], k[1]] in pt["truly_spanning_primary_in_window"]), None)
        row["breakpoint_window"] = f"{where['chrom']}:{where['position']}" if where else None
        if where:
            row["counted_split_m0"] = any(c_["key"] == [k[0], k[1]] for c_ in where["split_m0"]["counted"])
            row["genuine_split_m0"] = any(c_["key"] == [k[0], k[1]] and c_["category"] == "genuine"
                                          for c_ in where["split_m0"]["counted"])
            row["genuine_split_m20"] = any(c_["key"] == [k[0], k[1]] and c_["category"] == "genuine"
                                           for c_ in where["split_m20"]["counted"])
            row["caught_by_soft_clip"] = any(c_["key"] == [k[0], k[1]] and not c_["supplementary"]
                                             for c_ in where["soft_clips"]["counted"])
        if not sa:
            # the gate's partner window and alignment, unchanged (gate_a1.py)
            if row["junction"] == "J20":
                w = (fa.fetch("chr21", B - ga.WIN_BEFORE, B + ga.WIN_AFTER) if partner == "chr21"
                     else fa.fetch("chr20", A - ga.WIN_AFTER, A + ga.WIN_BEFORE))
            else:
                w = (fa.fetch("chr20", A - ga.WIN_BEFORE, A + ga.WIN_AFTER) if partner == "chr20"
                     else fa.fetch("chr21", B - ga.WIN_AFTER, B + ga.WIN_BEFORE))
            w = w.upper()
            s_f, l_f = ga.sw(pr["seq"], w)
            s_r, l_r = ga.sw(pr["seq"].translate(COMP)[::-1], w)
            st = ga.side_stats(ga.art_columns(art[(t["art_name"], k[1])]))
            row.update({"sw_best_score_vs_partner": max(s_f, s_r), "longest_exact_vs_partner": max(l_f, l_r),
                        "art_shorter_side": st,
                        "alignable_to_partner_despite_side_rule":
                            row["shorter_side"] < ga.T_MIN and max(s_f, s_r) >= ga.T_MIN and max(l_f, l_r) >= ga.K_MIN})
        reads.append(row)
    record = {"implant": iid, "class": imp["class"], "breakpoints": {"chr20": A, "chr21": B},
              "truly_spanning": len(ts), "points": points, "reads": reads,
              "enumeration_problems": problems, "sim_records_checked": sim_seen,
              "sim_records_differing_from_sim_bam": sim_mismatch, "calls": rec.calls}
    ec.write_record(os.path.join(OUT, f"{iid}.json.gz"), record)
    print(f"{iid}: {len(ts)} truly spanning; enumeration problems {len(problems)}; simulated records checked "
          f"{sim_seen}, differing {sim_mismatch}", flush=True)
    for p_ in problems:
        print("  PROBLEM", p_)
    return 0 if not problems and sim_mismatch == 0 else 1


def pooled(points, key):
    counted = sum(len(p[key]["counted"]) for p in points)
    genuine = sum(sum(c["category"] == "genuine" for c in p[key]["counted"]) for p in points)
    cats = {}
    for p in points:
        for c in p[key]["counted"]:
            cats[c["category"]] = cats.get(c["category"], 0) + 1
    return {"counted": counted, "genuine": genuine, "false_positive": counted - genuine,
            "fp_rate": round((counted - genuine) / counted, 4) if counted else None, "categories": cats}


def step_merge():
    T = ec.truth()
    tol = ec.tolerance()
    recs = {i: ec.read_record(os.path.join(OUT, f"{i}.json.gz")) for i in ec.IMPLANTS}
    bad = {i: r["enumeration_problems"] for i, r in recs.items() if r["enumeration_problems"]}
    if bad or any(r["sim_records_differing_from_sim_bam"] for r in recs.values()):
        raise SystemExit(f"refusing to merge: enumeration problems {bad}")
    chains = {i: ec.read_record(os.path.join(ec.OUT, "chain", f"{i}.json.gz")) for i in ec.IMPLANTS}
    classes = {i: {"chosen": T[i]["class"], "strict": ec.implant_class(T[i], ec.T_STRICT)[0]} for i in ec.IMPLANTS}

    # (a) tool vs delly vs truth, per junction
    rows = []
    for i in ec.IMPLANTS:
        calls = {c["id"]: c for c in chains[i]["calls"]}
        unf = calls[chains[i]["unfiltered"]]["result"]["candidates"]
        final = {c["candidate_id"] for c in calls[chains[i]["funnel_cumulative"][-1]["call"]]["result"]["candidates"]}
        for jn, exp in ec.expected_junctions(T[i]).items():
            m, _ = ec.best_match(unf, exp, tol)
            truth_n = sum(1 for r in recs[i]["reads"] if r["junction"] == jn)
            genuine0 = sum(1 for r in recs[i]["reads"] if r["junction"] == jn and r.get("genuine_split_m0"))
            genuine20 = sum(1 for r in recs[i]["reads"] if r["junction"] == jn and r.get("genuine_split_m20"))
            rows.append({"implant": i, "junction": jn, "class": classes[i]["chosen"],
                         "truly_spanning": truth_n, "tool_genuine_m0": genuine0, "tool_genuine_m20": genuine20,
                         "delly_sr": (m["sr"] or 0) if m else 0, "delly_called": m is not None,
                         "delly_survives_funnel": bool(m and m["candidate_id"] in final)})
    tot_truth = sum(r["truly_spanning"] for r in rows)
    tot_g0 = sum(r["tool_genuine_m0"] for r in rows)
    tot_sr = sum(r["delly_sr"] for r in rows)
    per_imp = []
    for i in ec.IMPLANTS:
        rr = [r for r in rows if r["implant"] == i]
        raw_count = sum(len(p["split_m0"]["counted"]) for p in recs[i]["points"])
        truth_n = sum(r["truly_spanning"] for r in rr)
        g0 = sum(r["tool_genuine_m0"] for r in rr)
        sr = sum(r["delly_sr"] for r in rr)
        per_imp.append({"implant": i, "class": classes[i]["chosen"], "truly_spanning": truth_n,
                        "tool_count_m0_both_breakpoints": raw_count, "tool_genuine_m0": g0, "delly_sr": sr,
                        "closer_raw_count_vs_delly": ("tool" if abs(raw_count - truth_n) < abs(sr - truth_n) else
                                                      "delly" if abs(sr - truth_n) < abs(raw_count - truth_n) else "tie"),
                        "closer_genuine_vs_delly": ("tool" if abs(g0 - truth_n) < abs(sr - truth_n) else
                                                    "delly" if abs(sr - truth_n) < abs(g0 - truth_n) else "tie")})
    points = [dict(p, implant=i, klass=classes[i]["chosen"]) for i in ec.IMPLANTS for p in recs[i]["points"]]
    by_class = {}
    for k in ec.CLASS_ORDER:
        pts = [p for p in points if p["klass"] == k]
        by_class[k] = {"m0": pooled(pts, "split_m0"), "m20": pooled(pts, "split_m20"),
                       "breakpoints_zeroed_by_m20": sum(1 for p in pts if p["split_m0"]["counted"]
                                                        and not p["split_m20"]["counted"]),
                       "breakpoints": len(pts)}
    per_point = []
    for p in points:
        n = len(p["split_m0"]["counted"])
        g = sum(c["category"] == "genuine" for c in p["split_m0"]["counted"])
        per_point.append({"breakpoint": f"{p['implant']} {p['chrom']}:{p['position']}", "class": p["klass"],
                          "tool_count_m0": n, "genuine_m0": g, "fp_m0": n - g,
                          "fp_rate_m0": round((n - g) / n, 4) if n else None,
                          "tool_count_m20": len(p["split_m20"]["counted"]),
                          "genuine_m20": sum(c["category"] == "genuine" for c in p["split_m20"]["counted"]),
                          "low_mapq_fraction_500": p["low_mapq_fraction_500"],
                          "categories_m0": pooled([p], "split_m0")["categories"]})
    xs = [p["low_mapq_fraction_500"] for p in per_point if p["fp_rate_m0"] is not None]
    ys = [p["fp_rate_m0"] for p in per_point if p["fp_rate_m0"] is not None]
    r = ec.pearson(xs, ys)
    # The same rate with truly spanning reads whose SA names a wrong partner counted
    # as TRUE: they do cross the junction, although the tool reports their partner
    # as some other chromosome. Both readings are reported; neither is hidden.
    for p in per_point:
        n, sae = p["tool_count_m0"], p["categories_m0"].get("spanning_sa_elsewhere", 0)
        p["fp_rate_m0_spanning_counted_true"] = round((p["fp_m0"] - sae) / n, 4) if n else None
    ys2 = [p["fp_rate_m0_spanning_counted_true"] for p in per_point if p["fp_rate_m0"] is not None]
    r2 = ec.pearson(xs, ys2)
    by_class_alt = {}
    for k in ec.CLASS_ORDER:
        pts = [p for p in per_point if p["class"] == k]
        n = sum(p["tool_count_m0"] for p in pts)
        nonspan = sum(p["fp_m0"] - p["categories_m0"].get("spanning_sa_elsewhere", 0) for p in pts)
        by_class_alt[k] = {"counted": n, "not_spanning": nonspan, "rate": round(nonspan / n, 4) if n else None}
    # where the wrong-partner SA tags of truly spanning reads point, from sim.bam
    wrong = []
    for i in ec.IMPLANTS:
        sa_of = {}
        with pysam.AlignmentFile(os.path.join(WORK, i, "sim.bam")) as f:
            for x in f:
                if not (x.is_secondary or x.is_supplementary) and x.has_tag("SA"):
                    sa_of[(x.query_name, 1 if x.is_read1 else 2)] = x.get_tag("SA")
        for p in recs[i]["points"]:
            for c in p["split_m0"]["counted"]:
                if c["category"] == "spanning_sa_elsewhere":
                    ents = [e.split(",") for e in sa_of[tuple(c["key"])].rstrip(";").split(";") if e]
                    wrong.append({"read": f"{c['key'][0]}/{c['key'][1]}", "breakpoint": f"{p['chrom']}:{p['position']}",
                                  "class": classes[i]["chosen"],
                                  "sa_targets": [f"{e[0]}:{e[1]}" for e in ents],
                                  "sa_mapq": [ec_int for ec_int in (_int(e[4]) for e in ents if len(e) >= 5)]})
    wrong_mapq = {}
    for w in wrong:
        for q in w["sa_mapq"]:
            wrong_mapq[q] = wrong_mapq.get(q, 0) + 1
    # (b) short overhang
    no_sa = [dict(r_, implant=i) for i in ec.IMPLANTS for r_ in recs[i]["reads"]
             if r_.get("primary") not in ("unmapped", "not found") and not r_.get("has_sa")]
    all_ts = [r_ for i in ec.IMPLANTS for r_ in recs[i]["reads"]]
    flagged = [f"{r_['read']}/{r_['mate']} (shorter side {r_['shorter_side']} bp)" for r_ in no_sa
               if r_.get("alignable_to_partner_despite_side_rule")]
    # why an SA-less truly spanning read escapes the soft-clip layer, first reason
    # that applies, in the order the layer itself filters; "other" must stay 0
    cig = {}
    for i in ec.IMPLANTS:
        with pysam.AlignmentFile(os.path.join(WORK, i, "sim.bam")) as f:
            for x in f:
                if not (x.is_secondary or x.is_supplementary) and not x.is_unmapped:
                    c = x.cigartuples or []
                    cig[(x.query_name, 1 if x.is_read1 else 2)] = (
                        x.mapping_quality,
                        max([n for op, n in c[:1] + c[-1:] if op == 4] or [0]))
    escaped = {}
    for r_ in no_sa:
        if r_.get("caught_by_soft_clip"):
            continue
        mq, clip = cig[(r_["read"], r_["mate"])]
        why = ("primary outside both breakpoint windows" if not r_.get("breakpoint_window") else
               f"primary MAPQ < {CLIP_MAPQ}" if mq < CLIP_MAPQ else
               f"shorter side < {MIN_CLIP} bp" if r_["shorter_side"] < MIN_CLIP else
               f"longest end clip < {MIN_CLIP} bp although the shorter side is >= {MIN_CLIP} bp"
               if clip < MIN_CLIP else "other")
        escaped[why] = escaped.get(why, 0) + 1
        r_["escaped_because"] = why
    # (c) min_mapq 20
    g20 = sum(r_["tool_genuine_m20"] for r_ in rows)
    fp0 = sum(pooled([p], "split_m0")["false_positive"] for p in points)
    fp20 = sum(pooled([p], "split_m20")["false_positive"] for p in points)
    out = {
        "generated_by": "scripts/synthetic/adjudicate_reads.py merge",
        "enumeration_reproduces_every_tool_count": True,
        "a_tool_vs_delly_vs_truth": {
            "truly_spanning_reads": tot_truth,
            "tool_genuine_m0": tot_g0, "tool_genuine_m0_fraction": round(tot_g0 / tot_truth, 4),
            "delly_sr_sum": tot_sr, "delly_sr_fraction": round(tot_sr / tot_truth, 4),
            "per_junction": rows, "per_implant": per_imp,
            "closer_to_truth_raw_count": {k: sum(p["closer_raw_count_vs_delly"] == k for p in per_imp)
                                          for k in ("tool", "delly", "tie")},
            "closer_to_truth_genuine": {k: sum(p["closer_genuine_vs_delly"] == k for p in per_imp)
                                        for k in ("tool", "delly", "tie")},
            "by_class": by_class, "per_breakpoint": per_point,
            "fp_rate_pearson_r_vs_low_mapq_fraction": round(r, 3) if r is not None else None,
            "fp_rate_r_n_breakpoints": len(xs),
            "alternative_reading_wrong_partner_spanning_reads_counted_true": {
                "by_class": by_class_alt,
                "pearson_r_vs_low_mapq_fraction": round(r2, 3) if r2 is not None else None},
            "spanning_reads_with_a_wrong_partner_sa": {"reads": len(wrong), "sa_mapq_counts": wrong_mapq,
                                                       "list": wrong}},
        "b_short_overhang": {
            "truly_spanning_reads": len(all_ts),
            "primary_unmapped_or_missing": sum(1 for r_ in all_ts if r_.get("primary") in ("unmapped", "not found")),
            "without_any_sa_tag": len(no_sa),
            "without_sa_caught_by_soft_clip_layer": sum(1 for r_ in no_sa if r_.get("caught_by_soft_clip")),
            "without_sa_outside_both_windows": sum(1 for r_ in no_sa if not r_.get("breakpoint_window")),
            "without_sa_shorter_side_under_10": sum(1 for r_ in no_sa if r_["shorter_side"] < MIN_CLIP),
            "without_sa_shorter_side_27_to_29": sum(1 for r_ in no_sa if 27 <= r_["shorter_side"] <= 29),
            "flagged_alignable_despite_side_rule": flagged,
            "not_caught_first_reason": escaped,
            "with_sa_but_none_to_partner": sum(1 for r_ in all_ts if r_.get("has_sa") and not r_.get("sa_to_partner")),
            "reads": no_sa},
        "c_min_mapq_20": {
            "genuine_m0": tot_g0, "genuine_m20": g20, "genuine_lost": tot_g0 - g20,
            "false_positive_m0": fp0, "false_positive_m20": fp20, "false_positives_removed": fp0 - fp20,
            "by_class": {k: {"genuine_m0": v["m0"]["genuine"], "genuine_m20": v["m20"]["genuine"],
                             "fp_m0": v["m0"]["false_positive"], "fp_m20": v["m20"]["false_positive"],
                             "breakpoints_zeroed_by_m20": v["breakpoints_zeroed_by_m20"],
                             "breakpoints": v["breakpoints"]} for k, v in by_class.items()}},
        "strict_classes": {i: classes[i]["strict"] for i in ec.IMPLANTS},
    }
    ec.write_record(os.path.join(OUT, "summary.json"), out)
    a = out["a_tool_vs_delly_vs_truth"]
    print(f"(a) truly spanning {tot_truth}; tool genuine (min_mapq 0) {tot_g0} = {tot_g0 / tot_truth:.1%}; "
          f"delly SR {tot_sr} = {tot_sr / tot_truth:.1%}")
    print(f"    closer to truth per implant, raw count vs delly: {a['closer_to_truth_raw_count']}; "
          f"genuine vs delly: {a['closer_to_truth_genuine']}")
    for k, v in by_class.items():
        print(f"    {k}: tool count {v['m0']['counted']}, genuine {v['m0']['genuine']}, FP rate {v['m0']['fp_rate']} "
              f"{v['m0']['categories']}")
    print(f"    FP rate vs low_mapq_fraction: r = {a['fp_rate_pearson_r_vs_low_mapq_fraction']} "
          f"over {len(xs)} breakpoints with a count")
    for p in sorted(per_point, key=lambda p: -(p["fp_m0"])):
        if p["fp_m0"]:
            print(f"      {p['breakpoint']}: {p['fp_m0']}/{p['tool_count_m0']} FP, lmf {p['low_mapq_fraction_500']}")
    alt = a["alternative_reading_wrong_partner_spanning_reads_counted_true"]
    print(f"    reading 2 (spanning reads with a wrong-partner SA counted as true): {alt['by_class']}; "
          f"r = {alt['pearson_r_vs_low_mapq_fraction']}")
    w = a["spanning_reads_with_a_wrong_partner_sa"]
    print(f"    truly spanning reads whose only SA names a wrong partner: {w['reads']}; SA mapQ {w['sa_mapq_counts']}")
    b = out["b_short_overhang"]
    print(f"(b) {b['without_any_sa_tag']} of {b['truly_spanning_reads']} truly spanning reads carry no SA tag; "
          f"{b['without_sa_caught_by_soft_clip_layer']} caught by the soft-clip layer; "
          f"{b['without_sa_shorter_side_under_10']} have a shorter side under {MIN_CLIP} bp; flagged {flagged}")
    print(f"    not caught, first reason: {b['not_caught_first_reason']}")
    c = out["c_min_mapq_20"]
    print(f"(c) min_mapq 0 -> 20: genuine {c['genuine_m0']} -> {c['genuine_m20']} (lost {c['genuine_lost']}); "
          f"false positives {c['false_positive_m0']} -> {c['false_positive_m20']} (removed {c['false_positives_removed']})")
    for k, v in c["by_class"].items():
        print(f"    {k}: {v}")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["run", "merge"])
    ap.add_argument("implant", nargs="?")
    a = ap.parse_args()
    if a.step == "run":
        if a.implant not in ec.IMPLANTS:
            raise SystemExit(f"implant must be one of {ec.IMPLANTS}")
        sys.exit(step_run(a.implant))
    sys.exit(step_merge())
