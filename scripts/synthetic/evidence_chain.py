#!/usr/bin/env python3
"""Phase 6 analysis, rerun on the rebuilt implants (2026-09-25): the evidence
chain, the headline figures, the threshold grid, the scoring ceiling and the
rescue.

    evidence_chain.py chain LABEL        LABEL = background | IMP01 .. IMP12 (one process each)
    evidence_chain.py rescue IMPxx       the evidence tools at an implant's TRUE breakpoints
    evidence_chain.py analyse            the figures, from the records the two above wrote

A NEW experiment, not a reproduction (results/synthetic_control_2026-09/README.md).
Figures quoted from the lost Phase 5/6 records are carried below only to be set
beside the new ones; every difference is reported, none is a target.

Every tool call goes through the MCP servers' own dispatch -- ui.ToolRecorder,
the path the interface uses -- so what is recorded is what the tools return,
including the server-layer fields (min_mapq_applied, position_provenance, the
attainable-ceiling echo). Every call and its return is written to the record,
with the home directory written as "~".

CHAIN, one process per BCF (the thirteen run concurrently)
  1 load_candidate_set
  2 list_candidates, unfiltered: every junction, for discovery
  3 the standard funnel, cumulatively, one call per step: PASS; PE >= 3;
    SR >= 1; primary contigs; outside the delly v2.6.0 exclude template
  4 the threshold grid PE >= 1,2,3,5 x SR >= 0,1,2, each with PASS, primary
    and the exclude template
  5 applicable_layers(BAM, sample_reads=20000), once, as ui.assess does
  6 every survivor of the standard funnel: get_candidate, then at BOTH
    breakends the calls ui.assess makes -- bam_stats_at_locus (+/-500),
    discordant_pairs (window 500), soft_clipped_reads, split_reads (the MCP
    tool's own default, min_mapq 0), read_depth_profile (+/-2000, focus on the
    breakend), breakpoint_evidence_summary (the applicable layers, window
    500) -- and reciprocal_breakpoint once per junction, breakend_1 as primary.
    gene_at_locus (network) and the two IGV tools are not run: no figure here
    depends on them.

RESCUE, a fresh process per implant, so get_candidate has handed out no
coordinate: load the implant's BCF and list the standard funnel, then the same
calls at the true breakpoints chr20:A and chr21:B, and reciprocal_breakpoint
between them. Every return must carry position_provenance "caller_supplied".
Then two controls: get_candidate on the first survivor, and discordant_pairs at
its breakend_1 must read "candidate_set" (the check can pass); the same call one
base further on must read "caller_supplied" (a near miss is not laundered).

MATCHING A CALL TO AN IMPLANTED JUNCTION. The bridge puts the higher contig
first, so both junctions of a t(20;21) are chr21-first:
    J20  chr20:A | chr21:B+1   ->  (chr21, B+1, chr20, A,   "5to3")
    J21  chr21:B | chr20:A+1   ->  (chr21, B,   chr20, A+1, "3to5")
A junction matches when both chromosomes and the orientation agree and both
breakends lie within vcf_tools.RECURRENCE_TOLERANCE_BP (500 bp); of several,
the closest. J20 and J21 are one base apart at both ends, so orientation is what
tells them apart. The raw BCF is read too (pysam), so a junction lost inside
the bridge is told apart from one delly never reported.

FIGURES (analyse)
  sensitivity   a junction is detected when a standard-funnel survivor matches
                it; a translocation when at least one of its junctions is (both
                is reported too); by class under the chosen class rule and under
                the stricter upper-decile rule (dinucleotide run >= 8 counts as
                repeat; an implant takes the least favourable class of its two
                breakends)
  loss          no raw BCF record -> discovery; a record but no bridge junction
                -> the bridge; otherwise the first funnel step whose survivors
                no longer include it
  precision     survivors per implant BAM that match no implanted junction, set
                against the background's survivors by candidate_id (the id is a
                hash of the normalised coordinates, so an identical call has an
                identical id in both sets) and by compare_candidate_sets at 500 bp
  localisation  each detected junction's offset from the truth at both ends
  grid          implanted junctions surviving / mean non-implanted survivors per
                implant BAM, in every cell
  ceiling       at every breakend of every detected junction: the four component
                scores, the observed discordant fraction, the composite and the
                attainable ceiling the tool returns
  rescue        the evidence at the true coordinates of every missed implant
  caller vs tool  at every background survivor, delly's SR against the
                split_reads tool's count (min_mapq 0) at both breakends, beside
                the low_mapq_fraction there
"""
import argparse
import gzip
import json
import math
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
from synth_common import BAMS, BG, DELLY_DIR, EXCL, HOME  # noqa: E402

RESULTS = os.path.join(REPO, "stage1_igv_assistant", "results", "synthetic_control_2026-09")
OUT = os.path.join(RESULTS, "analysis_2026-09-25")
TRUTH = os.path.join(RESULTS, "implants_ground_truth.json")
IMPLANTS = [f"IMP{i:02d}" for i in range(1, 13)]
LABELS = ["background"] + IMPLANTS
LIMIT = 100000
FUNNEL = [("filter_pass", True), ("min_pe", 3), ("min_sr", 1), ("primary_only", True),
          ("mask_path", EXCL)]
STEP_NAMES = ["PASS", "PE>=3", "SR>=1", "primary", "unmasked"]
GRID_PE = (1, 2, 3, 5)
GRID_SR = (0, 1, 2)
CLASS_ORDER = ["clean_unique", "repeat_adjacent", "low_mappability"]

# Quoted by the user from the lost Phase 5/6 records (thesis and abstract text).
# Carried only to be set beside the new figures.
QUOTED = {
    "sensitivity": "7/12 translocations, 14/24 junctions; clean 8/8, repeat-adjacent 6/8, low-mappability 0/8 [Phase 6]",
    "loss": "all at discovery [Phase 6]",
    "background_funnel": "894 -> 513 -> 155 -> 27 [Phase 6]",
    "added_per_detected_implant": "exactly 2 [Phase 6]",
    "localisation": "12/14 at 0 bp, 2/14 at 2 bp [Phase 6]",
    "reused_coordinates": {"IMP01": "detected", "IMP06": "missed", "IMP10": "missed"},
    "grid": "14/24 in every cell; SR>=1 cut non-implanted survivors 156 -> 28; SR>=2 identical to SR>=1 [Phase 6]",
    "ceiling": "all breakends 40.0-55.0 moderate; discordant fraction max 0.175; ceiling 57.5 at IMP01 chr20; "
               "IMP02 and IMP04 breakends reached 55 through a spurious depth contribution [Phase 6]",
    "rescue": "IMP06 chr20:3,900,000: 44 discordant pairs, 20 soft clips, 14 split reads [Phase 6]",
    "caller_vs_tool": "agreement within ~1.5 reads where low_mapq_fraction ~0.004; divergence up to "
                      "100-fold where ~0.457 [Phase 5]",
}


# ── plumbing ────────────────────────────────────────────────────────────────

def tilde_all(obj):
    """Every string with the home directory written as "~"."""
    if isinstance(obj, dict):
        return {k: tilde_all(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [tilde_all(v) for v in obj]
    if isinstance(obj, tuple):
        return [tilde_all(v) for v in obj]
    if isinstance(obj, str):
        return obj.replace(HOME, "~")
    return obj


def tracked(path):
    return subprocess.run(["git", "-C", REPO, "ls-files", "--error-unmatch", path],
                          capture_output=True).returncode == 0


def write_record(path, obj):
    """A committed record is never overwritten."""
    if os.path.exists(path) and tracked(path):
        raise SystemExit(f"refusing to overwrite a committed record: {os.path.relpath(path, REPO)}")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    data = json.dumps(tilde_all(obj), indent=1).encode()
    tmp = path + ".tmp"
    if path.endswith(".gz"):
        with gzip.GzipFile(tmp, "wb", mtime=0) as f:
            f.write(data)
    else:
        with open(tmp, "wb") as f:
            f.write(data)
    os.replace(tmp, path)


def read_record(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        return json.load(f)


def recorder():
    from stage1_igv_assistant import ui
    return ui.ToolRecorder()


def result(rec):
    r = rec.get("result")
    if isinstance(r, dict) and "error" in r:
        raise RuntimeError(f"{rec['tool']} returned an error: {r['error']}")
    return r


def bam_for(label):
    return BG if label == "background" else os.path.join(BAMS, f"{label}.bam")


def truth():
    return {i["id"]: i for i in json.load(open(TRUTH))["implants"]}


def expected_junctions(imp):
    A, B = imp["breakpoints"]["chr20"], imp["breakpoints"]["chr21"]
    return {"J20": ("chr21", B + 1, "chr20", A, "5to3"),
            "J21": ("chr21", B, "chr20", A + 1, "3to5")}


def tolerance():
    from stage1_igv_assistant.tools.vcf_tools import RECURRENCE_TOLERANCE_BP
    return RECURRENCE_TOLERANCE_BP


def matches(j, exp, tol):
    c1, p1, c2, p2, ct = exp
    return (j["chrom1"] == c1 and j["chrom2"] == c2 and j["orientation"] == ct
            and abs(j["pos1"] - p1) <= tol and abs(j["pos2"] - p2) <= tol)


def best_match(cands, exp, tol):
    hits = [j for j in cands if matches(j, exp, tol)]
    if not hits:
        return None, 0
    best = min(hits, key=lambda j: (abs(j["pos1"] - exp[1]) + abs(j["pos2"] - exp[3]), j["candidate_id"]))
    return best, len(hits)


def raw_bcf_matches(bcf, exp, tol):
    """Records in the BCF itself, read with pysam, that match an expected junction
    after the bridge's own canonicalisation. Independent of the bridge's parse,
    dedup and filters."""
    import pysam
    from stage1_igv_assistant.tools.vcf_tools import _canonicalise
    out = []
    with pysam.VariantFile(bcf) as vf:
        for r in vf:                          # every record: which contig delly puts first is not assumed
            if r.info.get("SVTYPE") != "BND":
                continue
            c1, p1, c2, p2, ct, _, _ = _canonicalise(r.chrom, r.pos, r.info.get("CHR2"),
                                                       r.info.get("POS2"), r.info.get("CT"), None, None)
            j = {"chrom1": c1, "pos1": p1, "chrom2": c2, "pos2": p2, "orientation": ct}
            if matches(j, exp, tol):
                out.append({"id": r.id, "chrom": r.chrom, "pos": r.pos, "chr2": r.info.get("CHR2"),
                            "pos2": r.info.get("POS2"), "ct": r.info.get("CT"),
                            "filter": ";".join(r.filter.keys()) or "PASS",
                            "pe": r.info.get("PE"), "sr": r.info.get("SR") if "SR" in r.info else None,
                            "precise": "PRECISE" in r.info})
    return out


# ── chain ───────────────────────────────────────────────────────────────────

def assess_calls(rec, bam, chrom, pos, applicable):
    """The calls ui.assess makes at one coordinate; returns their call ids."""
    ids = {}
    ids["stats"] = rec.call("evidence", "bam_stats_at_locus",
                            {"bam_path": bam, "chromosome": chrom,
                             "start": max(0, pos - 500), "end": pos + 500})["id"]
    ids["discordant_pairs"] = rec.call("evidence", "discordant_pairs",
                                       {"bam_path": bam, "chromosome": chrom, "position": pos,
                                        "window_bp": 500})["id"]
    ids["soft_clipped_reads"] = rec.call("evidence", "soft_clipped_reads",
                                         {"bam_path": bam, "chromosome": chrom, "position": pos})["id"]
    ids["split_reads"] = rec.call("evidence", "split_reads",
                                  {"bam_path": bam, "chromosome": chrom, "position": pos})["id"]
    ids["read_depth_profile"] = rec.call("evidence", "read_depth_profile",
                                         {"bam_path": bam, "chromosome": chrom, "start": pos - 2000,
                                          "end": pos + 2000, "focus_position": pos})["id"]
    ids["breakpoint_evidence_summary"] = rec.call(
        "evidence", "breakpoint_evidence_summary",
        {"bam_path": bam, "chromosome": chrom, "position": pos, "label": f"{chrom}:{pos}",
         "applicable_layers": applicable, "window_bp": 500})["id"]
    return ids


def step_chain(label):
    bcf = os.path.join(DELLY_DIR, f"{label}.bcf")
    bam = bam_for(label)
    rec = recorder()
    t0 = time.time()
    load = rec.call("bridge", "load_candidate_set", {"path": bcf, "label": label})
    sid = result(load)["set_id"]
    unfiltered = rec.call("bridge", "list_candidates", {"set_id": sid, "limit": LIMIT})
    params, cumulative = {"set_id": sid, "limit": LIMIT}, []
    for (k, v), name in zip(FUNNEL, STEP_NAMES):
        params = dict(params, **{k: v})
        cumulative.append({"step": name, "call": rec.call("bridge", "list_candidates", dict(params))["id"]})
    grid = {}
    for pe in GRID_PE:
        for sr in GRID_SR:
            grid[f"PE>={pe},SR>={sr}"] = rec.call(
                "bridge", "list_candidates",
                {"set_id": sid, "filter_pass": True, "min_pe": pe, "min_sr": sr, "primary_only": True,
                 "mask_path": EXCL, "limit": LIMIT})["id"]
    calls = {c["id"]: c for c in rec.calls}
    funnel = result(calls[cumulative[-1]["call"]])
    if funnel["truncated"]:
        raise SystemExit("the standard funnel's list was truncated -- raise LIMIT")
    layers = rec.call("evidence", "applicable_layers", {"bam_path": bam, "sample_reads": 20000})
    applicable = result(layers).get("applicable_layers")
    evidence = []
    for c in funnel["candidates"]:
        g = rec.call("bridge", "get_candidate", {"set_id": sid, "candidate_id": c["candidate_id"]})
        gr = result(g)
        ends = {}
        for which in ("breakend_1", "breakend_2"):
            b = gr[which]
            ends[which] = {"chromosome": b["chromosome"], "position": b["position"],
                           "calls": assess_calls(rec, bam, b["chromosome"], b["position"], applicable)}
        recip = rec.call("evidence", "reciprocal_breakpoint",
                         {"bam_path": bam,
                          "primary_chromosome": gr["breakend_1"]["chromosome"],
                          "primary_position": gr["breakend_1"]["position"],
                          "partner_chromosome": gr["breakend_2"]["chromosome"],
                          "partner_position": gr["breakend_2"]["position"]})
        evidence.append({"candidate_id": c["candidate_id"], "get_candidate": g["id"],
                         "breakends": ends, "reciprocal_breakpoint": recip["id"]})
    errors = [c["id"] for c in rec.calls if c["is_error"] or c.get("traceback")]
    tool_errors = [c["id"] for c in rec.calls if isinstance(c.get("result"), dict) and "error" in c["result"]]
    record = {"label": label, "bcf": bcf, "bam": bam, "set_id": sid,
              "load": load["id"], "unfiltered": unfiltered["id"], "funnel_cumulative": cumulative,
              "grid": grid, "applicable_layers": layers["id"], "applicable": applicable,
              "survivors": len(funnel["candidates"]), "evidence": evidence,
              "dispatch_errors": errors, "tool_errors": tool_errors,
              "n_calls": len(rec.calls), "wall_seconds": round(time.time() - t0, 1),
              "calls": rec.calls}
    write_record(os.path.join(OUT, "chain", f"{label}.json.gz"), record)
    print(f"{label}: {len(funnel['candidates'])} survivors, {len(rec.calls)} calls, "
          f"{len(errors)} dispatch errors, {len(tool_errors)} tool error returns, "
          f"{record['wall_seconds']}s", flush=True)
    return 0 if not errors else 1


# ── rescue ──────────────────────────────────────────────────────────────────

def step_rescue(iid):
    imp = truth()[iid]
    A, B = imp["breakpoints"]["chr20"], imp["breakpoints"]["chr21"]
    bam, bcf = bam_for(iid), os.path.join(DELLY_DIR, f"{iid}.bcf")
    rec = recorder()
    load = rec.call("bridge", "load_candidate_set", {"path": bcf, "label": iid})
    sid = result(load)["set_id"]
    funnel = rec.call("bridge", "list_candidates",
                      {"set_id": sid, "filter_pass": True, "min_pe": 3, "min_sr": 1, "primary_only": True,
                       "mask_path": EXCL, "limit": LIMIT})
    layers = rec.call("evidence", "applicable_layers", {"bam_path": bam, "sample_reads": 20000})
    applicable = result(layers).get("applicable_layers")
    ends = {"chr20": {"position": A, "calls": assess_calls(rec, bam, "chr20", A, applicable)},
            "chr21": {"position": B, "calls": assess_calls(rec, bam, "chr21", B, applicable)}}
    recip = rec.call("evidence", "reciprocal_breakpoint",
                     {"bam_path": bam, "primary_chromosome": "chr20", "primary_position": A,
                      "partner_chromosome": "chr21", "partner_position": B})
    prov = [{"call": c["id"], "tool": c["tool"],
             "source": ((c["result"] or {}).get("position_provenance") or {}).get("source")}
            for c in rec.calls if c["server"] == "evidence" and c["tool"] != "applicable_layers"]
    controls = {"positive": None, "negative": None}
    surv = result(funnel)["candidates"]
    if surv:
        g = result(rec.call("bridge", "get_candidate", {"set_id": sid, "candidate_id": surv[0]["candidate_id"]}))
        b = g["breakend_1"]
        pos = rec.call("evidence", "discordant_pairs",
                       {"bam_path": bam, "chromosome": b["chromosome"], "position": b["position"]})
        neg = rec.call("evidence", "discordant_pairs",
                       {"bam_path": bam, "chromosome": b["chromosome"], "position": b["position"] + 1})
        controls = {"positive": {"call": pos["id"], "at": f"{b['chromosome']}:{b['position']}",
                                 "source": result(pos)["position_provenance"]["source"],
                                 "holds": result(pos)["position_provenance"]["source"] == "candidate_set"},
                    "negative": {"call": neg["id"], "at": f"{b['chromosome']}:{b['position'] + 1}",
                                 "source": result(neg)["position_provenance"]["source"],
                                 "holds": result(neg)["position_provenance"]["source"] == "caller_supplied"}}
    ok = all(p["source"] == "caller_supplied" for p in prov)
    record = {"implant": iid, "bam": bam, "bcf": bcf, "true_breakpoints": {"chr20": A, "chr21": B},
              "ends": ends, "reciprocal_breakpoint": recip["id"],
              "provenance": prov, "all_caller_supplied": ok, "controls": controls,
              "dispatch_errors": [c["id"] for c in rec.calls if c["is_error"] or c.get("traceback")],
              "calls": rec.calls}
    write_record(os.path.join(OUT, "rescue", f"{iid}.json.gz"), record)
    print(f"{iid} rescue at chr20:{A:,} / chr21:{B:,}: provenance caller_supplied on "
          f"{sum(p['source'] == 'caller_supplied' for p in prov)}/{len(prov)} calls; controls "
          f"positive {controls['positive'] and controls['positive']['holds']}, "
          f"negative {controls['negative'] and controls['negative']['holds']}", flush=True)
    good = ok and (not surv or (controls["positive"]["holds"] and controls["negative"]["holds"]))
    return 0 if good and not record["dispatch_errors"] else 1


# ── analyse ─────────────────────────────────────────────────────────────────

T_CHOSEN = {"low_mapq_gate": 0.4, "homopolymer_repeat": 10, "dinucleotide_repeat": 10,
            "entropy_repeat_below": 3.6148, "lmf_partial": 0.05, "mappability_partial_below": 0.9,
            "homopolymer_clean_max": 8, "dinucleotide_clean_max": 8, "entropy_clean_min": 3.6897}
# select_breakpoints.py's own note: under a strict upper-decile rule a dinucleotide
# run >= 8 bp (the 90th percentile) counts as repeat; clean then allows at most 7.
T_STRICT = dict(T_CHOSEN, dinucleotide_repeat=8, dinucleotide_clean_max=7)


def classify(r, T):
    """select_breakpoints.classify with the thresholds passed in."""
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


def implant_class(imp, T):
    ends = {c: classify(imp["scan_statistics"][c], T) for c in ("chr20", "chr21")}
    rank = {k: i for i, k in enumerate(CLASS_ORDER + ["unclassified"])}
    return max(ends.values(), key=lambda k: rank[k]), ends


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(sxx * syy)


def frac(a, b):
    return f"{a}/{b}"


def step_analyse():
    tol = tolerance()
    T = truth()
    chains = {l: read_record(os.path.join(OUT, "chain", f"{l}.json.gz")) for l in LABELS}
    calls = {l: {c["id"]: c for c in chains[l]["calls"]} for l in LABELS}

    def survivors_of(label, call_id):
        return calls[label][call_id]["result"]["candidates"]

    # every record's own health
    health = {l: {"dispatch_errors": len(chains[l]["dispatch_errors"]),
                  "tool_error_returns": len(chains[l]["tool_errors"]), "calls": chains[l]["n_calls"]}
              for l in LABELS}

    # classes
    classes = {}
    for iid in IMPLANTS:
        chosen, ends_c = implant_class(T[iid], T_CHOSEN)
        strict, ends_s = implant_class(T[iid], T_STRICT)
        classes[iid] = {"recorded": T[iid]["class"], "chosen_rule": chosen, "chosen_ends": ends_c,
                        "strict_rule": strict, "strict_ends": ends_s}
    reproduces = all(classes[i]["chosen_rule"] == classes[i]["recorded"] for i in IMPLANTS)

    # per implant, per junction
    per = {}
    for iid in IMPLANTS:
        ch = chains[iid]
        unfilt = survivors_of(iid, ch["unfiltered"])
        steps = [(s["step"], {c["candidate_id"] for c in survivors_of(iid, s["call"])})
                 for s in ch["funnel_cumulative"]]
        final = steps[-1][1]
        bcf = os.path.join(DELLY_DIR, f"{iid}.bcf")
        js = {}
        for jn, exp in expected_junctions(T[iid]).items():
            raw = raw_bcf_matches(bcf, exp, tol)
            m, n_hits = best_match(unfilt, exp, tol)
            if not raw and m is None:
                lost = "discovery (no delly record)"
            elif m is None:
                lost = "bridge (a delly record, but no bridge junction)"
            else:
                lost = next((name for name, ids in steps if m["candidate_id"] not in ids), None)
            ent = {"expected": {"chrom1": exp[0], "pos1": exp[1], "chrom2": exp[2], "pos2": exp[3],
                                "orientation": exp[4]},
                   "raw_bcf_records": raw, "bridge_junction": m, "bridge_matches": n_hits,
                   "lost_at": lost, "detected": m is not None and lost is None}
            if m is not None:
                ent["offsets_bp"] = {"pos1": m["pos1"] - exp[1], "pos2": m["pos2"] - exp[3]}
            js[jn] = ent
        per[iid] = {"class": T[iid]["class"], "breakpoints": T[iid]["breakpoints"], "junctions": js,
                    "detected_junctions": sum(v["detected"] for v in js.values()),
                    "survivors": len(final)}

    def sens(ids):
        jd = sum(per[i]["detected_junctions"] for i in ids)
        any_ = sum(per[i]["detected_junctions"] >= 1 for i in ids)
        both = sum(per[i]["detected_junctions"] == 2 for i in ids)
        return {"junctions": frac(jd, 2 * len(ids)), "translocations_any_junction": frac(any_, len(ids)),
                "translocations_both_junctions": frac(both, len(ids))}

    by_class = {}
    for rule in ("chosen_rule", "strict_rule"):
        by_class[rule] = {}
        for k in CLASS_ORDER + ["unclassified"]:
            ids = [i for i in IMPLANTS if classes[i][rule] == k]
            if ids:
                by_class[rule][k] = dict(sens(ids), implants=ids)
    changed = {k: {"chosen": by_class["chosen_rule"].get(k), "strict": by_class["strict_rule"].get(k)}
               for k in CLASS_ORDER
               if (by_class["chosen_rule"].get(k) or {}).get("junctions")
               != (by_class["strict_rule"].get(k) or {}).get("junctions")}
    sensitivity = {"overall": sens(IMPLANTS), "by_class": by_class,
                   "per_class_figures_that_change_under_the_strict_rule": changed,
                   "chosen_rule_reproduces_the_recorded_classes": reproduces}

    losses = {f"{i} {jn}": per[i]["junctions"][jn]["lost_at"] for i in IMPLANTS
              for jn in ("J20", "J21") if not per[i]["junctions"][jn]["detected"]}
    loss_summary = {}
    for v in losses.values():
        loss_summary[v] = loss_summary.get(v, 0) + 1

    # background funnel
    bgc = chains["background"]
    bg_unf = survivors_of("background", bgc["unfiltered"])
    bg_steps = [(s["step"], calls["background"][s["call"]]["result"]["total_matching"])
                for s in bgc["funnel_cumulative"]]
    bg_final_call = calls["background"][bgc["funnel_cumulative"][-1]["call"]]["result"]
    background_funnel = {"junctions_after_dedup": len(bg_unf),
                         "load": calls["background"][bgc["load"]]["result"],
                         "cumulative": [{"after": n, "surviving": v} for n, v in bg_steps],
                         "filters_applied_standalone": bg_final_call["filters_applied"]}
    bg_surv = {c["candidate_id"]: c for c in bg_final_call["candidates"]}

    # precision, and compare_candidate_sets through the bridge
    rec = recorder()
    bg_load = result(rec.call("bridge", "load_candidate_set",
                              {"path": os.path.join(DELLY_DIR, "background.bcf"), "label": "background"}))
    precision = {}
    for iid in IMPLANTS:
        ch = chains[iid]
        surv = survivors_of(iid, ch["funnel_cumulative"][-1]["call"])
        imp_ids = {per[iid]["junctions"][jn]["bridge_junction"]["candidate_id"]
                   for jn in ("J20", "J21") if per[iid]["junctions"][jn]["detected"]}
        non = [c for c in surv if c["candidate_id"] not in imp_ids]
        non_ids = {c["candidate_id"] for c in non}
        ld = result(rec.call("bridge", "load_candidate_set",
                             {"path": os.path.join(DELLY_DIR, f"{iid}.bcf"), "label": iid}))
        cmp_ = result(rec.call("bridge", "compare_candidate_sets",
                               {"set_a": ld["set_id"], "set_b": bg_load["set_id"]}))
        matched_a = {p["candidate_id_a"] for p in cmp_["matched_pairs"]}
        precision[iid] = {
            "survivors": len(surv), "implanted_survivors": len(imp_ids), "non_implanted_survivors": len(non),
            "background_survivors": len(bg_surv),
            "survivors_minus_background_survivors": len(surv) - len(bg_surv),
            "non_implanted_not_among_background_survivors_by_id":
                sorted(non_ids - set(bg_surv)),
            "background_survivors_not_among_non_implanted_by_id":
                sorted(set(bg_surv) - non_ids),
            "non_implanted_with_no_background_junction_within_tolerance":
                sorted(c["candidate_id"] for c in non if c["candidate_id"] not in matched_a),
            "compare_candidate_sets": {k: cmp_[k] for k in ("total_in_a", "total_in_b", "matched_in_a",
                                                            "matched_in_b", "unmatched_in_a", "unmatched_in_b")},
        }
    non_counts = [precision[i]["non_implanted_survivors"] for i in IMPLANTS]

    # localisation
    loc = []
    for iid in IMPLANTS:
        for jn in ("J20", "J21"):
            j = per[iid]["junctions"][jn]
            if j["detected"]:
                o = j["offsets_bp"]
                loc.append({"junction": f"{iid} {jn}", "offset_pos1": o["pos1"], "offset_pos2": o["pos2"],
                            "max_abs_bp": max(abs(o["pos1"]), abs(o["pos2"])),
                            "ci_pos1": j["bridge_junction"]["ci_pos1"], "ci_pos2": j["bridge_junction"]["ci_pos2"],
                            "precise": j["bridge_junction"]["precise"]})
    loc_hist = {}
    for x in loc:
        loc_hist[x["max_abs_bp"]] = loc_hist.get(x["max_abs_bp"], 0) + 1

    # the grid
    grid = {}
    for cell in chains["background"]["grid"]:
        implanted, nons = 0, []
        for iid in IMPLANTS:
            s = survivors_of(iid, chains[iid]["grid"][cell])
            hit = set()
            for jn, exp in expected_junctions(T[iid]).items():
                m, _ = best_match(s, exp, tol)
                if m is not None:
                    hit.add(m["candidate_id"])
            implanted += len(hit)
            nons.append(len([c for c in s if c["candidate_id"] not in hit]))
        grid[cell] = {"implanted_junctions_surviving": frac(implanted, 24),
                      "mean_non_implanted_survivors_per_implant_bam": round(sum(nons) / len(nons), 2),
                      "non_implanted_per_bam": nons,
                      "background_survivors": len(survivors_of("background", chains["background"]["grid"][cell]))}

    # the scoring ceiling at every breakend of every detected junction
    ceiling = []
    for iid in IMPLANTS:
        ev = {e["candidate_id"]: e for e in chains[iid]["evidence"]}
        for jn in ("J20", "J21"):
            j = per[iid]["junctions"][jn]
            if not j["detected"]:
                continue
            e = ev[j["bridge_junction"]["candidate_id"]]
            for which in ("breakend_1", "breakend_2"):
                b = e["breakends"][which]
                s = calls[iid][b["calls"]["breakpoint_evidence_summary"]]["result"]
                ceiling.append({
                    "junction": f"{iid} {jn}", "breakend": f"{b['chromosome']}:{b['position']}",
                    "discordant_pair_score": s["discordant_pair_score"], "soft_clip_score": s["soft_clip_score"],
                    "split_read_score": s["split_read_score"], "depth_score": s["depth_score"],
                    "observed_discordant_fraction": s["discordant_pairs"]["discordant_fraction"],
                    "discordant_pairs": s["discordant_pairs"]["discordant_pairs"],
                    "max_clips_at_position": s["soft_clips"]["max_clips_at_position"],
                    "split_read_fraction": s["split_reads"]["split_read_fraction"],
                    "depth_ratio_min_to_mean": s["depth_profile"]["summary"]["depth_ratio_min_to_mean"],
                    "dip_is_at_focus": s["depth_profile"]["summary"]["dip_is_at_focus"],
                    "evidence_score": s["evidence_score"], "evidence_score_raw": s["evidence_score_raw"],
                    "evidence_strength": s["evidence_strength"], "signal_layers": s["signal_layers"],
                    "attainable_here": s.get("attainable_here"), "strong_band": s.get("strong_band"),
                    "strong_band_reachable_here": s.get("strong_band_reachable_here"),
                    "max_with_flat_depth": s.get("max_with_flat_depth"),
                    "attainable_ceiling_derivable": s.get("attainable_ceiling_derivable"),
                    "depth_contribution_at_a_balanced_event": (s["depth_score"] or 0) > 0,
                    "low_mapq_fraction": s["locus_stats"]["low_mapq_fraction"]})
    scored = [c["evidence_score"] for c in ceiling if c["evidence_score"] is not None]
    ceiling_summary = {
        "breakends": len(ceiling),
        "evidence_score_range": [min(scored), max(scored)] if scored else None,
        "strengths": {k: sum(c["evidence_strength"] == k for c in ceiling)
                      for k in sorted({c["evidence_strength"] for c in ceiling})},
        "max_observed_discordant_fraction": max((c["observed_discordant_fraction"] or 0) for c in ceiling) if ceiling else None,
        "attainable_here_range": ([min(c["attainable_here"] for c in ceiling if c["attainable_here"] is not None),
                                   max(c["attainable_here"] for c in ceiling if c["attainable_here"] is not None)]
                                  if any(c["attainable_here"] is not None for c in ceiling) else None),
        "any_strong_band_reachable": any(c["strong_band_reachable_here"] for c in ceiling),
        "breakends_with_a_depth_contribution": [f"{c['junction']} {c['breakend']}" for c in ceiling
                                                if c["depth_contribution_at_a_balanced_event"]],
    }

    # the rescue, for every missed implant
    rescue = {}
    for iid in IMPLANTS:
        path = os.path.join(OUT, "rescue", f"{iid}.json.gz")
        if not os.path.exists(path):
            continue
        rr = read_record(path)
        rc = {c["id"]: c for c in rr["calls"]}
        ends = {}
        for chrom, e in rr["ends"].items():
            cid = e["calls"]
            s = rc[cid["breakpoint_evidence_summary"]]["result"]
            ends[chrom] = {
                "position": e["position"],
                "discordant_pairs": rc[cid["discordant_pairs"]]["result"]["discordant_pairs"],
                "soft_clipped_reads": rc[cid["soft_clipped_reads"]]["result"]["soft_clipped_reads"],
                "max_clips_at_position": rc[cid["soft_clipped_reads"]]["result"]["max_clips_at_position"],
                "split_reads_min_mapq_0": rc[cid["split_reads"]]["result"]["split_reads"],
                "low_mapq_fraction": rc[cid["stats"]]["result"]["low_mapq_fraction"],
                "evidence_score": s["evidence_score"], "evidence_strength": s["evidence_strength"],
                "components": [s["discordant_pair_score"], s["soft_clip_score"], s["split_read_score"],
                               s["depth_score"]],
                "provenance": rc[cid["discordant_pairs"]]["result"]["position_provenance"]["source"]}
        rescue[iid] = {"missed": per[iid]["detected_junctions"] == 0,
                       "junctions_detected": per[iid]["detected_junctions"],
                       "all_caller_supplied": rr["all_caller_supplied"], "controls": rr["controls"],
                       "ends": ends}

    # caller versus tool at the background's survivors
    cvt = []
    for e in bgc["evidence"]:
        c = bg_surv[e["candidate_id"]]
        row = {"candidate_id": c["candidate_id"], "junction": f"{c['chrom1']}:{c['pos1']}-{c['chrom2']}:{c['pos2']}",
               "svtype": c["svtype"], "delly_sr": c["sr"], "delly_pe": c["pe"]}
        for which in ("breakend_1", "breakend_2"):
            b = e["breakends"][which]
            row[which] = {"split_reads_tool": calls["background"][b["calls"]["split_reads"]]["result"]["split_reads"],
                          "min_mapq_applied": calls["background"][b["calls"]["split_reads"]]["result"].get("min_mapq_applied"),
                          "low_mapq_fraction": calls["background"][b["calls"]["stats"]]["result"]["low_mapq_fraction"]}
        tool = max(row["breakend_1"]["split_reads_tool"], row["breakend_2"]["split_reads_tool"])
        row["tool_max_of_two_breakends"] = tool
        row["difference_tool_minus_delly"] = tool - (c["sr"] or 0)
        row["ratio_tool_over_delly"] = round(tool / c["sr"], 2) if c["sr"] else None
        row["low_mapq_fraction_max"] = max(row["breakend_1"]["low_mapq_fraction"], row["breakend_2"]["low_mapq_fraction"])
        cvt.append(row)
    cvt.sort(key=lambda r: r["low_mapq_fraction_max"])

    # the reused coordinates
    reused = {i: {"now": "detected" if per[i]["detected_junctions"] else "missed",
                  "junctions_detected": per[i]["detected_junctions"],
                  "phase6": QUOTED["reused_coordinates"][i],
                  "differs": (per[i]["detected_junctions"] > 0) != (QUOTED["reused_coordinates"][i] == "detected")}
              for i in ("IMP01", "IMP06", "IMP10")}

    out = {
        "generated_by": "scripts/synthetic/evidence_chain.py analyse",
        "tolerance_bp": tol,
        "funnel": "PASS, PE>=3, SR>=1, primary contigs, outside the delly v2.6.0 exclude template (no svtype filter)",
        "record_health": health,
        "classes": classes,
        "sensitivity": sensitivity,
        "losses": losses, "loss_summary": loss_summary,
        "background_funnel": background_funnel,
        "precision": precision,
        "precision_summary": {"background_survivors": len(bg_surv),
                              "non_implanted_survivors_per_implant_bam": non_counts,
                              "mean": round(sum(non_counts) / len(non_counts), 2)},
        "localisation": loc, "localisation_histogram_max_abs_bp": loc_hist,
        "grid": grid,
        "ceiling": ceiling, "ceiling_summary": ceiling_summary,
        "rescue": rescue,
        "caller_vs_tool_background": cvt,
        "caller_vs_tool_pearson_r_absdiff_vs_low_mapq": (
            round(pearson([r["low_mapq_fraction_max"] for r in cvt],
                          [abs(r["difference_tool_minus_delly"]) for r in cvt]), 3)
            if pearson([r["low_mapq_fraction_max"] for r in cvt],
                       [abs(r["difference_tool_minus_delly"]) for r in cvt]) is not None else None),
        "reused_coordinates": reused,
        "per_implant": per,
        "quoted_from_lost_records": QUOTED,
        "analyse_calls": rec.calls,
    }
    write_record(os.path.join(OUT, "figures.json"), out)
    # a readable digest
    print(f"record health: {sum(h['dispatch_errors'] for h in health.values())} dispatch errors, "
          f"{sum(h['tool_error_returns'] for h in health.values())} tool error returns over "
          f"{sum(h['calls'] for h in health.values())} calls")
    print(f"chosen class rule reproduces the recorded classes: {reproduces}")
    s = sensitivity["overall"]
    print(f"SENSITIVITY junctions {s['junctions']}; translocations (>=1 junction) {s['translocations_any_junction']}, "
          f"(both) {s['translocations_both_junctions']}   [{QUOTED['sensitivity']}]")
    for rule in ("chosen_rule", "strict_rule"):
        print(f"  {rule}: " + "; ".join(f"{k} {v['junctions']} junctions ({','.join(v['implants'])})"
                                         for k, v in by_class[rule].items()))
    print(f"  per-class figures that change under the strict rule: {sorted(changed) or 'none'}")
    print(f"LOSSES {loss_summary}   [{QUOTED['loss']}]")
    for k, v in losses.items():
        print(f"  {k}: {v}")
    print(f"BACKGROUND FUNNEL {len(bg_unf)} junctions -> " + " -> ".join(f"{n} {v}" for n, v in bg_steps)
          + f"   [{QUOTED['background_funnel']}]")
    print(f"PRECISION non-implanted survivors per implant BAM {non_counts} (background {len(bg_surv)})")
    for i in IMPLANTS:
        p = precision[i]
        print(f"  {i}: {p['survivors']} survivors = {p['implanted_survivors']} implanted + "
              f"{p['non_implanted_survivors']} other; vs background by id: +{len(p['non_implanted_not_among_background_survivors_by_id'])}"
              f" / -{len(p['background_survivors_not_among_non_implanted_by_id'])}")
    print(f"LOCALISATION max |offset| histogram {dict(sorted(loc_hist.items()))} over {len(loc)} detected junctions"
          f"   [{QUOTED['localisation']}]")
    print("GRID (implanted surviving / mean non-implanted per BAM / background)")
    for cell, g in grid.items():
        print(f"  {cell:14s} {g['implanted_junctions_surviving']:6s} {g['mean_non_implanted_survivors_per_implant_bam']:7} "
              f"{g['background_survivors']}")
    cs = ceiling_summary
    print(f"CEILING {cs['breakends']} breakends; score range {cs['evidence_score_range']}; strengths {cs['strengths']}; "
          f"max discordant fraction {cs['max_observed_discordant_fraction']}; attainable {cs['attainable_here_range']}; "
          f"strong reachable anywhere: {cs['any_strong_band_reachable']}; depth contributions "
          f"{cs['breakends_with_a_depth_contribution']}")
    for i, r in rescue.items():
        if r["missed"]:
            print(f"RESCUE {i}: " + "; ".join(
                f"{c}:{e['position']:,} disc {e['discordant_pairs']} clips {e['soft_clipped_reads']} "
                f"(max {e['max_clips_at_position']}) split {e['split_reads_min_mapq_0']} score {e['evidence_score']} "
                f"{e['evidence_strength']} lmf {e['low_mapq_fraction']} provenance {e['provenance']}"
                for c, e in r["ends"].items()) + f"; all caller_supplied {r['all_caller_supplied']}")
    print(f"REUSED COORDINATES {reused}")
    print(f"CALLER vs TOOL at {len(cvt)} background survivors: r(|tool-delly|, low_mapq) = "
          f"{out['caller_vs_tool_pearson_r_absdiff_vs_low_mapq']}")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["chain", "rescue", "analyse"])
    ap.add_argument("label", nargs="?")
    a = ap.parse_args()
    if a.step == "chain":
        if a.label not in LABELS:
            raise SystemExit(f"label must be one of {LABELS}")
        sys.exit(step_chain(a.label))
    if a.step == "rescue":
        if a.label not in IMPLANTS:
            raise SystemExit(f"label must be one of {IMPLANTS}")
        sys.exit(step_rescue(a.label))
    sys.exit(step_analyse())
