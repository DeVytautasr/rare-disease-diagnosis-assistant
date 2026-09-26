#!/usr/bin/env python3
"""Phase 12 Task 6(c): the candidate-set bridge as it was before its fix, reconstructed.

    prefix_bridge.py public     background and implant sets (public data): validation, merge groups
    prefix_bridge.py patient    both patient rerun BCFs: the pre-fix funnel -- counts only
    prefix_bridge.py record     both results -> the committed record

WHY A RECONSTRUCTION. The pre-fix code was never committed: cdd3823, the commit that
added the bridge, already holds the fixed code. Its commit message and the comment in
compare_candidate_sets describe the two defects exactly, and the two functions below
restore them and nothing else. They are swapped into vcf_tools inside this process
only, and restored afterwards; stage1_igv_assistant/tools/vcf_tools.py is not touched.
  prefix_dedup     vcf_tools._dedup with orientation removed from the merge key, which
                   the commit message gives as "(chrom1, chrom2, pos1 bucket) alone";
                   the tolerance test and the representative rule are unchanged
  prefix_compare   vcf_tools.compare_candidate_sets with the scan for each set_a junction
                   stopped at its first match ("broke out of the scan on first match"),
                   so a set_b junction counts as matched only when it is some set_a
                   junction's first partner; every other line unchanged

VALIDATION, fixed before running, against the documented pre-fix figures:
  background (public)   pre-fix dedup gives 840 junctions (fixed: 894), and 24
                        survivors (fixed: 27) through the Phase 6 background funnel
                        PASS, PE>=3, SR>=1, primary contigs, outside the exclude template
  IMP01 vs IMP02        two sets built on the identical background, 896 junctions each:
                        the fixed compare gives 894 matched / 2 unmatched in both
                        directions; the pre-fix compare is documented as 894/2 on the
                        scanning side and 840/56 on the other (which set was set_a is
                        not recorded, so 840 is looked for on the target side of
                        either direction; both target-side counts are kept)

BLINDING. Patient results are counts only: no candidate ID, coordinate, gene or
interpretation is printed or written. Candidate lists are counted and discarded.
Display the patient step through scripts/redact.sh anyway.

FIGURES AND THEIR DEFINITIONS (written into the record beside each figure)
  merge_group         a set of two or more junctions that the pre-fix dedup merges
                      into one (a connected component of its union-find), before the
                      representative is chosen
  del_dup_groups      pre-fix merge groups with at least one DEL and at least one DUP
                      junction among their members
  inv_both_groups     pre-fix merge groups with INV junctions of both orientations
                      (3to3 and 5to5) among their members
  bnd_after_dedup     list_candidates(svtype=BND) total_matching, over the set loaded
                      with the pre-fix dedup
  pass_bnd            the same plus filter_pass
  before_recurrence   the same plus min_pe 3, min_sr 1, primary_only and the exclude
                      template (the patient funnel of write_patient_record.py)
  final_survivors     before_recurrence survivors not matched (both breakends within
                      500 bp) by the pre-fix compare, from the one call
                      compare(SAMPLE_A, SAMPLE_B) the funnel makes: SAMPLE_A's matched
                      set is the A side of the matched pairs, SAMPLE_B's the B side
  final_survivors_two_calls  the same survivors when each sample's matched set is
                      taken from the call in which it is set_a (the scanning side):
                      compare(SAMPLE_A, SAMPLE_B) for SAMPLE_A, compare(SAMPLE_B,
                      SAMPLE_A) for SAMPLE_B. Added 2026-09-26 after the first patient
                      run, to test whether "the final candidate counts were unchanged"
                      depends on how the funnel called the comparison; the figures of
                      the first run are unchanged by it
  two_direction_diff  for each set X, X's matched count with X as set_a (the scanning
                      side) minus X's matched count with X as set_b, over the whole
                      sets. Primary definition: the pre-fix compare on the sets loaded
                      with the pre-fix dedup (the state before the fix). Also given:
                      the same on the sets loaded with the fixed dedup (the compare
                      defect alone), and each sum over the two sets
  fixed_control       every figure above recomputed with the fixed functions in the
                      same process; the patient ones must equal the committed record
                      patient_rerun_2026-09.json (9,172/9,655; 896/923; 57/63; 17/19)
"""
import contextlib
import json
import os
import sys
from collections import Counter
from dataclasses import asdict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
from common import RUN_DIR, die  # noqa: E402
from stage1_igv_assistant.tools import vcf_tools as v  # noqa: E402

EXCL = os.path.expanduser("~/reference/human.hg38.excl.tsv")
PUBLIC = os.path.expanduser("~/public_data/sim/delly")
PATIENT_BCF = os.path.join(RUN_DIR, "delly", "{label}.bcf")
LABELS = ("SAMPLE_A", "SAMPLE_B")
LIMIT = 100000
BG_FUNNEL = {"filter_pass": True, "min_pe": 3, "min_sr": 1, "primary_only": True, "mask_path": EXCL}
PT_FUNNEL = {"svtype": "BND", "filter_pass": True, "min_pe": 3, "min_sr": 1, "primary_only": True,
             "mask_path": EXCL}
RECORD = os.path.join(REPO, "stage1_igv_assistant", "results", "prefix_bridge_2026-09-26.json")
GROUPS = []  # the member lists of the last pre-fix dedup, for the merge-group counts


def definitions():
    doc = __doc__.split("FIGURES AND THEIR DEFINITIONS (written into the record beside each figure)")[1]
    out, key = {}, None
    for line in doc.splitlines():
        if line.startswith("  ") and not line.startswith("    ") and line.strip():
            key, _, rest = line.strip().partition(" ")
            out[key] = rest.strip()
        elif key and line.strip():
            out[key] += " " + line.strip()
    return {k: " ".join(x.split()) for k, x in out.items()}


def prefix_dedup(raw, tol):
    """vcf_tools._dedup with (chrom1, chrom2, pos1 bucket) as the merge key: no orientation."""
    GROUPS.clear()
    if tol <= 0 or not raw:
        return list(raw), 0
    parent = list(range(len(raw)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    buckets = {}
    for i, j in enumerate(raw):
        buckets.setdefault((j.chrom1, j.chrom2, j.pos1 // 100000), []).append(i)
    for i, j in enumerate(raw):
        for b in (j.pos1 // 100000 - 1, j.pos1 // 100000, j.pos1 // 100000 + 1):
            for k in buckets.get((j.chrom1, j.chrom2, b), ()):
                if k != i and abs(raw[k].pos1 - j.pos1) <= tol and abs(raw[k].pos2 - j.pos2) <= tol:
                    union(i, k)
    groups = {}
    for i in range(len(raw)):
        groups.setdefault(find(i), []).append(i)
    out = []
    for members in groups.values():
        if len(members) > 1:
            GROUPS.append([(raw[i].svtype, raw[i].orientation) for i in members])
        best = max(members, key=lambda i: (
            raw[i].filter == "PASS",
            (raw[i].pe or 0) + (raw[i].sr or 0),
            raw[i].sr or 0,
        ))
        j = raw[best]
        out.append(v.Junction(**{**asdict(j), "n_merged": len(members)}))
    return out, len(raw) - len(out)


def prefix_compare(set_a, set_b, tolerance_bp=v.RECURRENCE_TOLERANCE_BP):
    """vcf_tools.compare_candidate_sets with the per-junction scan stopped at the first match."""
    a, b = v._SETS.get(set_a), v._SETS.get(set_b)
    if a is None or b is None:
        return {"error": "unknown set_id(s)", "error_type": "unknown_set_id"}
    idx = {}
    for j in b["junctions"]:
        idx.setdefault((j.chrom1, j.chrom2, j.pos1 // 100000), []).append(j)
    matches, matched_a, matched_b = [], set(), set()
    for j in a["junctions"]:
        found = None
        for bkt in (j.pos1 // 100000 - 1, j.pos1 // 100000, j.pos1 // 100000 + 1):
            for k in idx.get((j.chrom1, j.chrom2, bkt), ()):
                if abs(k.pos1 - j.pos1) <= tolerance_bp and abs(k.pos2 - j.pos2) <= tolerance_bp:
                    found = k
                    break
            if found is not None:
                break
        if found is not None:
            matches.append({"candidate_id_a": j.candidate_id, "candidate_id_b": found.candidate_id})
            matched_a.add(j.candidate_id)
            matched_b.add(found.candidate_id)
    return {"set_a": set_a, "label_a": a["label"], "set_b": set_b, "label_b": b["label"],
            "total_in_a": len(a["junctions"]), "total_in_b": len(b["junctions"]),
            "matched_in_a": len(matched_a), "matched_in_b": len(matched_b),
            "unmatched_in_a": len(a["junctions"]) - len(matched_a),
            "unmatched_in_b": len(b["junctions"]) - len(matched_b),
            "matched_pairs": matches}


@contextlib.contextmanager
def state(dedup, compare):
    """dedup/compare: 'prefix' or 'fixed'. Swaps the functions in this process only."""
    saved = v._dedup, v.compare_candidate_sets
    if dedup == "prefix":
        v._dedup = prefix_dedup
    if compare == "prefix":
        v.compare_candidate_sets = prefix_compare
    v.reset_registry()
    try:
        yield
    finally:
        v._dedup, v.compare_candidate_sets = saved
        v.reset_registry()


def load(path, label):
    r = v.load_candidate_set(path, label)
    if "error" in r:
        die(f"{label}: the set did not load ({r.get('error_type')})")
    return r


def total(set_id, **kw):
    return v.list_candidates(set_id, limit=LIMIT, **kw)["total_matching"]


def merge_groups():
    comp = Counter()
    for g in GROUPS:
        comp[" + ".join(f"{s}:{o}" for s, o in sorted(Counter(g).keys()))] += 1
    svt = [{s for s, _ in g} for g in GROUPS]
    inv = [{o for s, o in g if s == "INV"} for g in GROUPS]
    return {"merge_groups": len(GROUPS),
            "del_dup_groups": sum({"DEL", "DUP"} <= s for s in svt),
            "inv_both_groups": sum({"3to3", "5to5"} <= o for o in inv),
            "group_compositions": dict(comp.most_common())}


def two_directions(set_x, set_y):
    xy, yx = v.compare_candidate_sets(set_x, set_y, 500), v.compare_candidate_sets(set_y, set_x, 500)
    d = {"x_scanning_side": xy["matched_in_a"], "x_target_side": yx["matched_in_b"],
         "y_scanning_side": yx["matched_in_a"], "y_target_side": xy["matched_in_b"],
         "total_x": xy["total_in_a"], "total_y": xy["total_in_b"]}
    d["diff_x"] = d["x_scanning_side"] - d["x_target_side"]
    d["diff_y"] = d["y_scanning_side"] - d["y_target_side"]
    d["diff_sum"] = d["diff_x"] + d["diff_y"]
    return d


def step_public():
    rec = {"data": "public: ~/public_data/sim/delly/background.bcf, IMP01.bcf, IMP02.bcf"}
    bg = os.path.join(PUBLIC, "background.bcf")
    for name, dedup in (("prefix", "prefix"), ("fixed", "fixed")):
        with state(dedup, "fixed"):
            r = load(bg, "background")
            groups = merge_groups() if dedup == "prefix" else None
            rec[f"background_{name}_dedup"] = {
                "junctions_after_dedup": r["junctions_after_dedup"],
                "records_merged_by_dedup": r["records_merged_by_dedup"],
                "survivors_phase6_funnel": total(r["set_id"], **BG_FUNNEL)}
            if groups:
                rec[f"background_{name}_dedup"].update(groups)
    pair = {}
    for name, compare in (("prefix", "prefix"), ("fixed", "fixed")):
        with state("fixed", compare):
            a = load(os.path.join(PUBLIC, "IMP01.bcf"), "IMP01")["set_id"]
            b = load(os.path.join(PUBLIC, "IMP02.bcf"), "IMP02")["set_id"]
            pair[f"{name}_compare"] = two_directions(a, b)
    rec["IMP01_vs_IMP02_fixed_dedup"] = pair
    p, f, c = rec["background_prefix_dedup"], rec["background_fixed_dedup"], pair["prefix_compare"]
    rec["validation"] = {
        "background_prefix_junctions_840": p["junctions_after_dedup"] == 840,
        "background_prefix_survivors_24": p["survivors_phase6_funnel"] == 24,
        "background_fixed_894_and_27": (f["junctions_after_dedup"], f["survivors_phase6_funnel"]) == (894, 27),
        "pair_fixed_894_both_directions": (pair["fixed_compare"]["x_scanning_side"],
                                           pair["fixed_compare"]["x_target_side"],
                                           pair["fixed_compare"]["y_scanning_side"],
                                           pair["fixed_compare"]["y_target_side"]) == (894, 894, 894, 894),
        "pair_prefix_894_scanning_and_840_target_in_a_direction":
            (c["x_scanning_side"], c["y_scanning_side"]) == (894, 894) and 840 in (c["x_target_side"],
                                                                                   c["y_target_side"])}
    os.makedirs(RUN_DIR, mode=0o700, exist_ok=True)
    with open(os.path.join(RUN_DIR, "prefix_bridge_public.json"), "w") as fh:
        json.dump(rec, fh, indent=1)
    print(json.dumps(rec, indent=1))
    return 0


def patient_funnel(dedup, compare):
    out, sets = {}, {}
    with state(dedup, compare):
        for label in LABELS:
            r = load(PATIENT_BCF.format(label=label), label)
            sets[label] = r["set_id"]
            full = v.list_candidates(r["set_id"], limit=LIMIT, **PT_FUNNEL)
            out[label] = {"junctions_after_dedup": r["junctions_after_dedup"],
                          "bnd_after_dedup": total(r["set_id"], svtype="BND"),
                          "pass_bnd": total(r["set_id"], svtype="BND", filter_pass=True),
                          "before_recurrence": full["total_matching"],
                          "_surv": {c["candidate_id"] for c in full["candidates"]}}
            if dedup == "prefix":
                out[label]["merge_groups"] = {k: x for k, x in merge_groups().items() if k != "group_compositions"}
        cmp_ = v.compare_candidate_sets(sets["SAMPLE_A"], sets["SAMPLE_B"], 500)
        matched = {"SAMPLE_A": {p["candidate_id_a"] for p in cmp_["matched_pairs"]},
                   "SAMPLE_B": {p["candidate_id_b"] for p in cmp_["matched_pairs"]}}
        ba = v.compare_candidate_sets(sets["SAMPLE_B"], sets["SAMPLE_A"], 500)
        scanning = {"SAMPLE_A": matched["SAMPLE_A"], "SAMPLE_B": {p["candidate_id_a"] for p in ba["matched_pairs"]}}
        for label in LABELS:
            surv = out[label].pop("_surv")
            out[label]["final_survivors"] = len(surv - matched[label])
            out[label]["final_survivors_two_calls"] = len(surv - scanning[label])
        out["whole_set_comparison_A_to_B"] = {k: cmp_[k] for k in (
            "total_in_a", "total_in_b", "matched_in_a", "matched_in_b", "unmatched_in_a", "unmatched_in_b")}
        out["two_direction_diff"] = two_directions(sets["SAMPLE_A"], sets["SAMPLE_B"])
    return out


def step_patient():
    for label in LABELS:
        if not os.path.exists(PATIENT_BCF.format(label=label)):
            die(f"{label}: the rerun BCF is missing")
    rec = {"prefix_dedup_prefix_compare (primary: the state before the fix)": patient_funnel("prefix", "prefix"),
           "fixed_dedup_prefix_compare (the compare defect alone)": patient_funnel("fixed", "prefix"),
           "fixed_dedup_fixed_compare (control)": patient_funnel("fixed", "fixed")}
    ctl = rec["fixed_dedup_fixed_compare (control)"]
    want = {"SAMPLE_A": (9172, 896, 57, 17), "SAMPLE_B": (9655, 923, 63, 19)}
    rec["control_equals_committed_record"] = all(
        (ctl[L]["bnd_after_dedup"], ctl[L]["pass_bnd"], ctl[L]["before_recurrence"], ctl[L]["final_survivors"])
        == want[L] for L in LABELS)
    with open(os.path.join(RUN_DIR, "prefix_bridge_patient.json"), "w") as fh:
        json.dump(rec, fh, indent=1)
    print(json.dumps(rec, indent=1))
    return 0


def step_record():
    out = {"what": "The candidate-set bridge before its fix, reconstructed; Phase 12 Task 6(c), 2026-09-26",
           "reconstruction": "the pre-fix functions are not in history (cdd3823 added the fixed code); "
                             "scripts/patient/prefix_bridge.py restores the two documented defects in a "
                             "scratch module and leaves vcf_tools.py untouched",
           "blinding": "patient figures are counts only",
           "definitions": definitions(),
           "thesis_values_compared": {"background merge groups DEL+DUP": 33,
                                      "background merge groups with both INV orientations": 13,
                                      "pre-fix patient BND after dedup": "8,716 / 9,144",
                                      "pre-fix patient PASS": "838 / 862",
                                      "pre-fix patient final survivors": "15 / 16",
                                      "pre-fix two-direction difference on patient data": 1738},
           "code": "scripts/patient/prefix_bridge.py"}
    for name in ("public", "patient"):
        p = os.path.join(RUN_DIR, f"prefix_bridge_{name}.json")
        if not os.path.exists(p):
            die(f"the {name} step has not been run yet")
        out[name] = json.load(open(p))
    if os.path.exists(RECORD):
        die("the record exists; a committed record is never overwritten")
    with open(RECORD, "w") as fh:
        json.dump(out, fh, indent=1)
    print("written: stage1_igv_assistant/results/prefix_bridge_2026-09-26.json")
    return 0


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else ""
    if step not in ("public", "patient", "record"):
        die("usage: prefix_bridge.py public | patient | record")
    sys.exit({"public": step_public, "patient": step_patient, "record": step_record}[step]())
