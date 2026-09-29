#!/usr/bin/env python3
"""
Regression test: the reported ceiling is computed over the layers the score counts.

From results/ceiling_layers_fix_2026-09-30/ (registration.json, before.json) and
FIGURE_MAP J.5. breakpoint_evidence_summary normalises evidence_score over the
layers it counts -- applicable_layers minus unassessable_layers -- but the ceiling
fields server.py's _with_ceiling adds were derived as if all four layers were
always counted. With only discordant pairs and soft clips counted, IMP01
chr20:200000 scored 65.0 at the default window while the ceiling said 57.5.

POSITIVE CONTROL: with a restricted layer list the ceiling equals the normalised
maximum over those layers (the registered values at IMP01: 65.0, 100, 43.3), and
whenever the depth component is 0 the score does not exceed attainable_here. The
pre-fix code fails these.
NEGATIVE CONTROL: with all four layers counted, in any order, and with no list,
every ceiling field is identical to the pre-fix output (score_tiers.py and
server.py at PRE_FIX, loaded from git).

Self-contained: tiers are derived from bam_tools' source as the tool derives them,
observed values are synthetic, no BAM is needed -- except section 4, one
integration check through the MCP dispatch at IMP01 chr20:200000 (public data),
which reports NOT RUN, exit 2, if the BAM is absent.

To run against another copy of the code (e.g. the pre-fix one):
    CEILING_IMPL_COMMIT=<commit> python3 stage1_igv_assistant/tests/test_ceiling_layers.py
Exit codes: 0 passed, 1 a check failed, 2 incomplete.
Run: python3 stage1_igv_assistant/tests/test_ceiling_layers.py
"""
import asyncio
import contextlib
import importlib.util
import itertools
import json
import os
import subprocess
import sys
import tempfile

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, REPO)
import stage1_igv_assistant as pkg  # noqa: E402
from stage1_igv_assistant import score_tiers as live_st  # noqa: E402

PRE_FIX = "e219d4a84c3b8d2a7c606e75fb9d726a78d977f9"
ALL4 = ["discordant_pairs", "soft_clipped_reads", "split_reads", "read_depth"]
CEILING = ["attainable_ceiling_derivable", "score_bands", "strong_band", "max_all_layers", "max_with_flat_depth",
           "attainable_here", "strong_band_reachable_here", "attainable_basis", "attainable_note"]
BAM = os.path.expanduser("~/public_data/sim/bams/IMP01.bam")
FAILURES, NOT_RUN = [], []


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def from_git(commit, relpath, name):
    src = subprocess.run(["git", "-C", REPO, "show", f"{commit}:{relpath}"], capture_output=True, text=True,
                         check=True).stdout
    path = os.path.join(tempfile.mkdtemp(), os.path.basename(relpath))
    open(path, "w").write(src)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    return spec, mod


@contextlib.contextmanager
def score_tiers_as(mod):
    """server._with_ceiling imports score_tiers at call time; point that import at `mod`."""
    old_mod, old_attr = sys.modules.get("stage1_igv_assistant.score_tiers"), getattr(pkg, "score_tiers", None)
    sys.modules["stage1_igv_assistant.score_tiers"] = mod
    pkg.score_tiers = mod
    try:
        yield
    finally:
        sys.modules["stage1_igv_assistant.score_tiers"] = old_mod
        pkg.score_tiers = old_attr


def load_pair(commit, tag):
    s_spec, st = from_git(commit, "stage1_igv_assistant/score_tiers.py", f"score_tiers_{tag}")
    s_spec.loader.exec_module(st)
    v_spec, sv = from_git(commit, "stage1_igv_assistant/server.py", f"server_{tag}")
    with score_tiers_as(st):
        v_spec.loader.exec_module(sv)
    return st, sv


# the code under test
if os.environ.get("CEILING_IMPL_COMMIT"):
    IMPL_ST, IMPL_SV = load_pair(os.environ["CEILING_IMPL_COMMIT"], "impl")
    print(f"code under test: score_tiers.py and server.py at {os.environ['CEILING_IMPL_COMMIT']}")
else:
    from stage1_igv_assistant import server as IMPL_SV  # noqa: E402
    IMPL_ST = live_st
    print("code under test: the working tree")
PRE_ST, PRE_SV = load_pair(PRE_FIX, "prefix")

TIERS = live_st.derive_tiers()
MAX = {k: TIERS[k]["max_score"] for k in ALL4}


def result(observed, applicable=None, unassessable=None):
    """A summary return with the fields _with_ceiling reads; the score replayed as the tool computes it."""
    applicable = list(applicable) if applicable is not None else list(ALL4)
    unassessable = dict(unassessable or {})
    comp = {k: live_st.score_for(TIERS[k], observed[k]) for k in ALL4}
    counted = [l for l in applicable if l not in unassessable]
    score = round(sum(comp[l] for l in counted) * (100.0 / (len(counted) * 25.0)), 1) if counted else None
    return {"applicable_layers": applicable, "unassessable_layers": unassessable, "evidence_score": score,
            "discordant_pair_score": comp["discordant_pairs"], "soft_clip_score": comp["soft_clipped_reads"],
            "split_read_score": comp["split_reads"], "depth_score": comp["read_depth"],
            "discordant_pairs": {"discordant_fraction": observed["discordant_pairs"]},
            "soft_clips": {"max_clips_at_position": observed["soft_clipped_reads"]},
            "split_reads": {"split_read_fraction": observed["split_reads"]},
            "depth_profile": {"summary": {"depth_ratio_min_to_mean": observed["read_depth"]}}}


def impl(r):
    with score_tiers_as(IMPL_ST):
        return IMPL_SV._with_ceiling(r)


def prefix(r):
    with score_tiers_as(PRE_ST):
        return PRE_SV._with_ceiling(r)


def fields(r, extra=True):
    keys = CEILING + (["ceiling_counted_layers"] if extra else [])
    return json.dumps({k: r[k] for k in keys if k in r}, sort_keys=True)


def expected(observed, counted):
    s = sum(0.0 if l == "read_depth" else MAX[l] if l != "discordant_pairs"
            else live_st.score_for(TIERS[l], observed[l]) for l in counted)
    m = sum(MAX[l] for l in counted)
    return round(100.0 * s / m, 1), round(100.0 * sum(MAX[l] for l in counted if l != "read_depth") / m, 1)


IMP01 = {"discordant_pairs": 0.157, "soft_clipped_reads": 19, "split_reads": 0.058, "read_depth": 0.779}
REGISTERED = {("discordant_pairs", "soft_clipped_reads"): (65.0, 100.0, False),
              ("soft_clipped_reads", "split_reads"): (100.0, 100.0, True),
              ("discordant_pairs", "soft_clipped_reads", "read_depth"): (43.3, 66.7, False),
              tuple(ALL4): (57.5, 75.0, False)}
GRID = [dict(zip(ALL4, v)) for v in itertools.product((0.0, 0.05, 0.157, 0.25, 0.6), (0, 2, 5, 19),
                                                      (0.0, 0.05, 0.2, 0.4), (1.0, 0.779))]
SUBSETS = [list(c) for n in (1, 2, 3) for c in itertools.combinations(ALL4, n)]


def run_tests():
    print("=" * 68)
    print("CEILING OVER THE COUNTED LAYERS")
    print("=" * 68)

    print("\n1. positive control: the registered values at IMP01's observed values")
    for layers, (att, flat, reach) in REGISTERED.items():
        r = impl(result(IMP01, applicable=list(layers)))
        check(f"{'+'.join(layers)}: attainable_here {att}", r.get("attainable_here") == att, f"got {r.get('attainable_here')}")
        check(f"{'+'.join(layers)}: max_with_flat_depth {flat}", r.get("max_with_flat_depth") == flat,
              f"got {r.get('max_with_flat_depth')}")
        check(f"{'+'.join(layers)}: strong_band_reachable_here {reach}", r.get("strong_band_reachable_here") is reach,
              f"got {r.get('strong_band_reachable_here')}")

    print("\n2. positive control: every restricted list over a grid of observed values")
    wrong, above, named, bad_basis = [], [], [], []
    for obs in GRID:
        for layers in SUBSETS:
            r = impl(result(obs, applicable=layers))
            att, flat = expected(obs, layers)
            if (r.get("attainable_here"), r.get("max_with_flat_depth")) != (att, flat):
                wrong.append((obs, layers, r.get("attainable_here"), att))
            if r["depth_score"] == 0 and r["evidence_score"] is not None and r["evidence_score"] > r["attainable_here"]:
                above.append((obs, layers, r["evidence_score"], r["attainable_here"]))
            if r.get("ceiling_counted_layers") != layers:
                named.append((layers, r.get("ceiling_counted_layers")))
            if not all(l in (r.get("attainable_basis") or "") for l in layers):
                bad_basis.append(layers)
    n = len(GRID) * len(SUBSETS)
    check(f"the ceiling is the normalised maximum over the counted layers ({n} cases)", not wrong, f"{len(wrong)} wrong, e.g. {wrong[:2]}")
    check("with the depth component 0, the score never exceeds attainable_here", not above, f"{len(above)} above, e.g. {above[:2]}")
    check("ceiling_counted_layers names the counted layers", not named, f"{len(named)} wrong, e.g. {named[:2]}")
    check("attainable_basis names the counted layers", not bad_basis, f"{len(bad_basis)} wrong, e.g. {bad_basis[:2]}")

    print("\n3. an unassessable layer leaves the denominator, and the ceiling follows")
    r = impl(result(IMP01, unassessable={"split_reads": "no reads in the window"}))
    att, flat = expected(IMP01, ["discordant_pairs", "soft_clipped_reads", "read_depth"])
    check(f"all four applicable, split reads unassessable: attainable_here {att}", r.get("attainable_here") == att,
          f"got {r.get('attainable_here')}")

    print("\n4. negative control: all four layers, and no list, give the pre-fix fields")
    diff_all, diff_order, extra_key = [], [], []
    for obs in GRID:
        base = result(obs)
        want = fields(prefix(base), extra=False)
        got = impl(base)
        if fields(got) != want:
            diff_all.append(obs)
        if "ceiling_counted_layers" in got:
            extra_key.append(obs)
        if fields(impl(result(obs, applicable=list(reversed(ALL4))))) != fields(prefix(result(obs, applicable=list(reversed(ALL4)))), extra=False):
            diff_order.append(obs)
    check(f"all four counted: every ceiling field identical to the pre-fix output ({len(GRID)} cases)", not diff_all,
          f"{len(diff_all)} differ")
    check("all four in another order: identical to the pre-fix output", not diff_order, f"{len(diff_order)} differ")
    check("no ceiling_counted_layers field when all four are counted", not extra_key, f"{len(extra_key)} carry it")

    print("\n5. integration: through the MCP dispatch at IMP01 chr20:200000, default window")
    if not os.path.exists(BAM):
        print(f"  NOT RUN: {BAM} is absent")
        NOT_RUN.append("integration")
    else:
        async def call(layers):
            args = {"bam_path": BAM, "chromosome": "chr20", "position": 200000}
            if layers is not None:
                args["applicable_layers"] = layers
            with score_tiers_as(IMPL_ST):
                r = await IMPL_SV.mcp.call_tool("breakpoint_evidence_summary", args)
            return r.structured_content
        for layers, (att, _, _) in REGISTERED.items():
            r = asyncio.run(call(list(layers)))
            check(f"dispatch {'+'.join(layers)}: attainable_here {att}, score {r.get('evidence_score')} within it",
                  r.get("attainable_here") == att and r.get("evidence_score") <= r.get("attainable_here"),
                  f"got attainable {r.get('attainable_here')}, score {r.get('evidence_score')}")
        r = asyncio.run(call(None))
        stripped = {k: v for k, v in r.items() if k not in CEILING + ["ceiling_counted_layers", "attainable_ceiling_reason"]}
        check("dispatch, no list: the ceiling fields equal the pre-fix wrapper's on the same return",
              fields(r) == fields(prefix(stripped), extra=False))

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    if NOT_RUN:
        print(f"INCOMPLETE: not run: {NOT_RUN}")
        sys.exit(2)
    print("ALL CEILING-LAYER TESTS PASSED")


if __name__ == "__main__":
    run_tests()
