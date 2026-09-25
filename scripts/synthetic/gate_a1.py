#!/usr/bin/env python3
"""REVISED criterion (a) for the implant gate, and the dose arithmetic.

    gate_a1.py [IMPxx]          (default IMP01)

Defined on 2026-09-25 AFTER the pre-registered gate failed, on the user's
decision, and recorded as such. The pre-registered result -- (a) 18/33 = 54.5%
against >= 70%, FAIL -- stands in IMP01_gate.json and is not edited. The
definition of a truly spanning read is unchanged: a read whose own ART template
interval crosses its junction with at least one base on each side.

(a1) Every truly spanning read WITHOUT an SA tag to its partner chromosome must
     be explained by bwa mem's documented limits, at bwa 0.7.15's defaults:
     minimum output score -T 30 (scores -A 1, -B 4, gap -O 6 + -E 1 per base) and
     minimum seed length -k 19. For each such read, its bases as they appear in
     the implant BAM are aligned (Smith-Waterman, affine gaps, both strands) to
     the partner chromosome's reference from 60 bp before to 250 bp after the
     junction; the best score and the longest exact match are recorded.
     Explained if the best score < 30 (nothing bwa may output) or no exact match
     reaches 19 bp (no seed). A read whose partner side could score >= 30 with a
     19 bp seed yet has no SA tag is UNEXPLAINED: it is listed and the gate FAILS.
     Reported beside it: the shorter side's length and the errors ART placed in
     it, from ART's own =/X/I/D alignment.
(a2) Consistency, not pass/fail: the observed SA fraction against
       - the geometric expectation: a 150 bp read, the junction at one of its 149
         inter-base positions with equal probability, needs a shorter side of
         >= 30 bp for an error-free supplementary to reach -T 30;
       - the same, counting ART's actual errors: over every truly spanning read of
         all twelve implants, the share whose shorter side still reaches score 30
         with a 19 bp exact run, from ART's alignment;
       - 183/304 = 60.2%, the aggregate the user quotes from the original
         experiment (no surviving record to verify it against).
(b), (c) as pre-registered, read unchanged from IMPxx's gate.json.

DOSE CHECK. From the pre-implant depth at each breakpoint (D, reads per base),
the fragment length F (mean 442.64, the background's) and the read length R:
  pairs spanning a boundary        D / (2R) * (F - 1)      (fragment starts per base
                                                            x bases each start covers)
  removed (the event's homolog)    half of that; non-duplicate from the
                                   non-duplicate depth
  added fragments                  non-dup removed at A + at B, split between J20
                                   and J21
  spanning READS per junction      added fragments x 2 E[min(R-1, F-1)] / E[F-1]
                                   under ART's fragment distribution, normal(443,
                                   105) truncated at R -- a spanning fragment is
                                   length-biased, and each of its two reads crosses
                                   the junction with probability (R-1)/(F-1)
Expected, added and observed are set side by side for every implant.
"""
import json
import math
import os
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from synth_common import BAMS, BG, REF, SAMTOOLS, SIM, WORK, write_json  # noqa: E402
from stage1_igv_assistant.tools import bam_tools  # noqa: E402

T_MIN, K_MIN = 30, 19                  # bwa 0.7.15: -T 30, -k 19
MATCH, MISMATCH, GAP_O, GAP_E = 1, 4, 6, 1
HALF = 2000                            # make_implants.py: junction at template index 2000
R = 150
WIN_BEFORE, WIN_AFTER = 60, 250
COMP = str.maketrans("ACGTN", "TGCAN")


def sw(read, ref):
    """Best local score under bwa's scoring (Gotoh affine gaps) and the longest
    exact common substring -- the two quantities -T and -k are about."""
    n, m = len(read), len(ref)
    NEG = -10 ** 9
    H = [0] * (m + 1)
    E = [NEG] * (m + 1)
    best = 0
    for i in range(1, n + 1):
        prev_diag, F, Hrow = 0, NEG, [0] * (m + 1)
        a = read[i - 1]
        for j in range(1, m + 1):
            E[j] = max(E[j] - GAP_E, H[j] - GAP_O - GAP_E)
            F = max(F - GAP_E, Hrow[j - 1] - GAP_O - GAP_E)
            s = prev_diag + (MATCH if a == ref[j - 1] else -MISMATCH)
            h = max(0, s, E[j], F)
            prev_diag = H[j]
            Hrow[j] = h
            if h > best:
                best = h
        H = Hrow
    # longest common substring
    longest, row = 0, [0] * (m + 1)
    for i in range(1, n + 1):
        new = [0] * (m + 1)
        a = read[i - 1]
        for j in range(1, m + 1):
            if a == ref[j - 1]:
                new[j] = row[j - 1] + 1
                if new[j] > longest:
                    longest = new[j]
        row = new
    return best, longest


def art_columns(rec):
    """ART's alignment of one read to its template as columns (op, template_pos)."""
    cols, t = [], rec.reference_start
    for op, n in rec.cigartuples:
        if op in (7, 8, 0):                          # = X M
            for _ in range(n):
                cols.append(("=" if op == 7 else "X" if op == 8 else "M", t))
                t += 1
        elif op == 1:                                # insertion: read only
            for _ in range(n):
                cols.append(("I", t - 1))
        elif op == 2:                                # deletion: template only
            for _ in range(n):
                cols.append(("D", t))
                t += 1
    return cols


def side_stats(cols):
    """Split columns at the junction; for the shorter side (by read bases) give its
    length, ART errors, best local score under bwa's scoring and longest '=' run."""
    left = [c for c in cols if c[1] < HALF]
    right = [c for c in cols if c[1] >= HALF]
    rb = lambda cs: sum(1 for o, _ in cs if o != "D")
    side = left if rb(left) <= rb(right) else right
    best = cur = 0
    run = longest = 0
    k = 0
    while k < len(side):
        o = side[k][0]
        if o in ("I", "D"):
            g = 1
            while k + g < len(side) and side[k + g][0] == o:
                g += 1
            cur = max(0, cur - (GAP_O + GAP_E * g))
            run = 0
            k += g
            continue
        if o == "=":
            cur += MATCH
            run += 1
        else:
            cur = max(0, cur - MISMATCH)
            run = 0
        best = max(best, cur)
        longest = max(longest, run)
        k += 1
    return {"shorter_side_read_bases": rb(side),
            "art_mismatches": sum(1 for o, _ in side if o == "X"),
            "art_inserted": sum(1 for o, _ in side if o == "I"),
            "art_deleted": sum(1 for o, _ in side if o == "D"),
            "best_score_from_art": best, "longest_exact_from_art": longest,
            "reaches_T_and_k_from_art": best >= T_MIN and longest >= K_MIN}


def art_records(iid):
    recs = {}
    for j in ("J20", "J21"):
        with pysam.AlignmentFile(os.path.join(WORK, iid, f"{j}.art.sam"), check_sq=False) as f:
            for r in f:
                recs[(r.query_name, 1 if r.is_read1 else 2)] = r
    return recs


def truly_spanning(prep):
    return [t for t in prep["truth_reads"]
            if t["spans_junction"] and t["bases_left"] >= 1 and t["bases_right"] >= 1]


def nondup_depth(bam_path, chrom, x, w=500):
    depth = 0
    with pysam.AlignmentFile(bam_path) as f:
        for r in f.fetch(chrom, x - w, x + w):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                continue
            depth += sum(1 for p in r.get_reference_positions() if x - w <= p < x + w)
    return depth / (2 * w)


def reads_per_spanning_fragment(mean, sd):
    """2 E[min(R-1, F-1)] / E[F-1] over normal(mean, sd) truncated at F >= R."""
    num = den = 0.0
    for F in range(R, int(mean + 8 * sd)):
        p = math.exp(-0.5 * ((F - mean) / sd) ** 2)
        num += p * min(R - 1, F - 1)
        den += p * (F - 1)
    return 2 * num / den


def main(iid):
    prep = json.load(open(os.path.join(WORK, iid, "prepare.json")))
    gate = json.load(open(os.path.join(WORK, iid, "gate.json")))
    A, B = prep["chr20"], prep["chr21"]
    bam = os.path.join(BAMS, f"{iid}.bam")
    sim_only = os.path.join(WORK, iid, "gate_sim_records.bam")
    if not os.path.exists(sim_only):
        subprocess.run([SAMTOOLS, "view", "-b", "-r", "SIM", "-o", sim_only, bam], check=True)
    primary = {}
    with pysam.AlignmentFile(sim_only) as f:
        for r in f:
            if not (r.is_secondary or r.is_supplementary):
                primary[(r.query_name, 1 if r.is_read1 else 2)] = r
    art = art_records(iid)
    fa = pysam.FastaFile(REF)
    partner_sa = {(p["read"], p["mate"]): p["partner_sa"] for p in gate["per_read"]}
    rows, unexplained = [], []
    for t in truly_spanning(prep):
        key = (t["read"], t["mate"])
        r = primary[key]
        st = side_stats(art_columns(art[(t["art_name"], t["mate"])]))
        junction = t["read"].split("_")[1]
        if r.is_unmapped or r.reference_name not in ("chr20", "chr21"):
            raise SystemExit(f"{t['read']}/{t['mate']}: primary not on chr20/chr21 -- (a1) cannot place it")
        pchrom = "chr21" if r.reference_name == "chr20" else "chr20"
        # J20 = chr20[..A] + chr21[B+1..]; J21 = chr21[..B] + chr20[A+1..]; the partner
        # window runs from 60 bp on the far side of the junction to 250 bp into the segment
        if junction == "J20":
            window = (fa.fetch("chr21", B - WIN_BEFORE, B + WIN_AFTER) if pchrom == "chr21"
                      else fa.fetch("chr20", A - WIN_AFTER, A + WIN_BEFORE))
        else:
            window = (fa.fetch("chr20", A - WIN_BEFORE, A + WIN_AFTER) if pchrom == "chr20"
                      else fa.fetch("chr21", B - WIN_AFTER, B + WIN_BEFORE))
        window = window.upper()
        seq = r.query_sequence
        s_f, l_f = sw(seq, window)
        s_r, l_r = sw(seq.translate(COMP)[::-1], window)
        score, seed = max(s_f, s_r), max(l_f, l_r)
        has = partner_sa.get(key, False)
        row = {"read": t["read"], "mate": t["mate"], "primary": r.reference_name, "partner": pchrom,
               "has_partner_sa": has, "bases_left": t["bases_left"], "bases_right": t["bases_right"],
               **st, "sw_best_score_vs_partner": score, "longest_exact_vs_partner": seed}
        if not has:
            reasons = []
            if st["shorter_side_read_bases"] < T_MIN:
                reasons.append(f"shorter side {st['shorter_side_read_bases']} bp < {T_MIN}: score cannot reach -T {T_MIN}")
            if score < T_MIN:
                reasons.append(f"best alignment score to the partner {score} < -T {T_MIN}")
            if seed < K_MIN:
                reasons.append(f"longest exact match to the partner {seed} bp < -k {K_MIN}")
            row["explanation"] = reasons
            # The side rule ("under 30 bp explains itself") is the user's definition and
            # decides the verdict; but chance matches across the junction can lift a
            # 27-29 bp side to 30. Such reads are flagged, not failed.
            row["alignable_to_partner_despite_side_rule"] = (
                st["shorter_side_read_bases"] < T_MIN and score >= T_MIN and seed >= K_MIN)
            if not reasons:
                unexplained.append(row)
        rows.append(row)
    n = len(rows)
    with_sa = sum(r["has_partner_sa"] for r in rows)
    # (a2) geometric expectation
    geo_ok = sum(1 for i in range(1, R) if min(i, R - i) >= T_MIN)
    # (a2) with ART's actual errors, over every implant prepared
    all_ok = all_n = 0
    for d in sorted(os.listdir(WORK)):
        p = os.path.join(WORK, d, "prepare.json")
        if not d.startswith("IMP") or not os.path.exists(p):
            continue
        pr = json.load(open(p))
        ar = art_records(d)
        for t in truly_spanning(pr):
            all_n += 1
            all_ok += side_stats(art_columns(ar[(t["art_name"], t["mate"])]))["reaches_T_and_k_from_art"]
    def binom_cdf(k, n, p):
        return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k + 1))
    a2 = {"observed": f"{with_sa}/{n}", "observed_fraction": round(with_sa / n, 4),
          "geometric": {"positions_with_shorter_side_ge_30": geo_ok, "of": R - 1,
                        "fraction": round(geo_ok / (R - 1), 4),
                        "p_observed_or_fewer": round(binom_cdf(with_sa, n, geo_ok / (R - 1)), 4)},
          "with_art_errors_all_implants": {"reads": all_n, "reaching_T30_with_k19_seed": all_ok,
                                           "fraction": round(all_ok / all_n, 4) if all_n else None},
          "original_aggregate_as_quoted": {"with_sa": 183, "of": 304, "fraction": round(183 / 304, 4),
                                           "verifiable": False,
                                           "p_observed_or_fewer": round(binom_cdf(with_sa, n, 183 / 304), 4)},
          "note": "consistency check, not pass/fail"}
    # dose check, every prepared implant
    meas = json.load(open(os.path.join(SIM, "measure", "background.json")))["fragment_length"]
    fmean = meas["mean"]
    k_reads = reads_per_spanning_fragment(round(fmean), round(meas["sd"]))
    dose = []
    for d in sorted(os.listdir(WORK)):
        p = os.path.join(WORK, d, "prepare.json")
        if not d.startswith("IMP") or not os.path.exists(p):
            continue
        pr = json.load(open(p))
        e = {"implant": d}
        exp_nd = {}
        for side, chrom, x in (("A", "chr20", pr["chr20"]), ("B", "chr21", pr["chr21"])):
            D = bam_tools.get_bam_stats_at_locus(BG, chrom, x - 500, x + 500)["mean_depth"]
            Dnd = nondup_depth(BG, chrom, x)
            span = D / (2 * R) * (fmean - 1)
            span_nd = Dnd / (2 * R) * (fmean - 1)
            exp_nd[side] = span_nd / 2
            rm = pr["removal"][side]
            e[side] = {"depth": D, "depth_nondup": round(Dnd, 2),
                       "expected_spanning_pairs": round(span, 1), "observed_spanning_pairs": rm["spanning_pairs"],
                       "expected_removed_nondup": round(span_nd / 2, 1), "removed_nondup": rm["removed_nondup"]}
        add_exp = exp_nd["A"] + exp_nd["B"]
        e["added_fragments"] = {"expected_total": round(add_exp, 1), "added_total": pr["fragments_added_total"],
                                "J20": pr["junctions"]["J20"]["fragments_added"],
                                "J21": pr["junctions"]["J21"]["fragments_added"]}
        obs = {j: sum(1 for t in truly_spanning(pr) if t["read"].split("_")[1] == j) for j in ("J20", "J21")}
        e["spanning_reads"] = {
            "expected_per_junction_from_depth": round(add_exp / 2 * k_reads, 1),
            "expected_J20_from_added": round(pr["junctions"]["J20"]["fragments_added"] * k_reads, 1),
            "expected_J21_from_added": round(pr["junctions"]["J21"]["fragments_added"] * k_reads, 1),
            "observed_J20": obs["J20"], "observed_J21": obs["J21"],
            "observed_total": obs["J20"] + obs["J21"]}
        dose.append(e)
    b_ok, c_ok = gate["b"]["holds"], gate["c"]["holds"]
    verdict = "PASS" if not unexplained and b_ok and c_ok else "FAIL"
    rec = {"implant": iid, "defined": "2026-09-25, after the pre-registered gate failed (user decision); "
                                      "the pre-registered result stands unedited in gate.json",
           "bwa_limits": {"T": T_MIN, "k": K_MIN, "A": MATCH, "B": MISMATCH, "O": GAP_O, "E": GAP_E},
           "a1": {"truly_spanning": n, "without_partner_sa": n - with_sa, "unexplained": len(unexplained),
                  "holds": not unexplained, "unexplained_reads": unexplained,
                  "flagged_alignable_despite_side_rule": [f"{r_['read']}/{r_['mate']}" for r_ in rows
                                                          if r_.get("alignable_to_partner_despite_side_rule")]},
           "a2": a2, "b": gate["b"], "c": gate["c"],
           "reads_per_spanning_fragment": round(k_reads, 4), "fragment_mean": fmean,
           "dose": dose, "verdict": verdict, "per_read": rows}
    write_json(os.path.join(WORK, iid, "gate_revised.json"), rec)
    print(f"REVISED GATE on {iid} (criterion (a) replaced after the pre-registered gate failed)")
    print(f"  (a1) {n - with_sa} of {n} truly spanning reads lack a partner SA tag; unexplained by -T {T_MIN}/-k {K_MIN}: "
          f"{len(unexplained)}  -> {'HOLDS' if not unexplained else 'FAILS'}")
    for r_ in sorted((r for r in rows if not r["has_partner_sa"]), key=lambda r: r["shorter_side_read_bases"]):
        print(f"       {r_['read']}/{r_['mate']}: shorter side {r_['shorter_side_read_bases']} bp, ART errors "
              f"{r_['art_mismatches']}X/{r_['art_inserted']}I/{r_['art_deleted']}D; best score to partner "
              f"{r_['sw_best_score_vs_partner']}, longest exact {r_['longest_exact_vs_partner']} bp -> "
              + ("; ".join(r_["explanation"]) if r_["explanation"] else "UNEXPLAINED"))
    for u in unexplained:
        print(f"  UNEXPLAINED: {u['read']}/{u['mate']}")
    flagged = [r_ for r_ in rows if r_.get("alignable_to_partner_despite_side_rule")]
    print(f"  (a1) note: SA-less reads under 30 bp whose partner side could still reach -T 30 with a 19 bp "
          f"seed through chance matches across the junction: {len(flagged)}"
          + (" -- " + ", ".join(f"{r_['read']}/{r_['mate']}" for r_ in flagged) if flagged else ""))
    g = a2["geometric"]
    print(f"  (a2) observed {with_sa}/{n} = {with_sa / n:.1%}; geometric {g['positions_with_shorter_side_ge_30']}/"
          f"{g['of']} = {g['fraction']:.1%} (P(<= observed) {g['p_observed_or_fewer']}); with ART errors, all "
          f"implants {all_ok}/{all_n} = {all_ok / all_n:.1%}; quoted aggregate 183/304 = 60.2% "
          f"(P(<= observed) {a2['original_aggregate_as_quoted']['p_observed_or_fewer']}) -- consistency only")
    print(f"  (b) {gate['b']['wrong_partner_sa_entries']} SA entries to a wrong partner -> {'HOLDS' if b_ok else 'FAILS'}")
    for chrom, dd in gate["c"]["depth"].items():
        print(f"  (c) {chrom}:{dd['position']:,} depth {dd['pre_mean_depth']} -> {dd['post_mean_depth']} "
              f"({dd['change']:+.1%}); non-duplicate {dd['pre_nondup']} -> {dd['post_nondup']} ({dd['change_nondup']:+.1%})")
    print(f"  (c) -> {'HOLDS' if c_ok else 'FAILS'}")
    print(f"DOSE (reads per spanning fragment {k_reads:.4f} at fragment mean {fmean})")
    for e in dose:
        s = e["spanning_reads"]
        print(f"  {e['implant']}: spanning pairs A exp {e['A']['expected_spanning_pairs']} obs {e['A']['observed_spanning_pairs']}, "
              f"B exp {e['B']['expected_spanning_pairs']} obs {e['B']['observed_spanning_pairs']} | added exp "
              f"{e['added_fragments']['expected_total']} got {e['added_fragments']['added_total']} | spanning reads "
              f"exp {s['expected_per_junction_from_depth']}/junction (from added: J20 {s['expected_J20_from_added']}, "
              f"J21 {s['expected_J21_from_added']}), observed J20 {s['observed_J20']} J21 {s['observed_J21']}")
    print(f"REVISED GATE {verdict}")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "IMP01")
