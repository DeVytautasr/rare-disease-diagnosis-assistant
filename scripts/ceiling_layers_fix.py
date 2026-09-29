#!/usr/bin/env python3
"""Phase 20: the ceiling-field defect -- measurements before and after the fix. Public data only.

    ceiling_layers_fix.py before   -> results/ceiling_layers_fix_2026-09-30/before.json
    ceiling_layers_fix.py after    -> results/ceiling_layers_fix_2026-09-30/after.json

No step opens ~/patient_data. Every call goes through the MCP dispatch
(server.mcp.call_tool), as a client's would.

FIGURES AND THEIR DEFINITIONS
  counted_layers     the layers in evidence_score's denominator: the return's
                     applicable_layers minus its unassessable_layers, in that order
  ceiling_fields     attainable_ceiling_derivable, score_bands, strong_band,
                     max_all_layers, max_with_flat_depth, attainable_here,
                     strong_band_reachable_here, attainable_basis, attainable_note and,
                     after the fix, ceiling_counted_layers (present only when fewer than
                     four layers are counted)
  imp01_calls        breakpoint_evidence_summary at IMP01 chr20:200000, default window,
                     with applicable_layers = all four; [discordant_pairs,
                     soft_clipped_reads]; [soft_clipped_reads, split_reads];
                     [discordant_pairs, soft_clipped_reads, read_depth]
  phase10            every breakpoint_evidence_summary call recorded in the committed
                     Phase 10 run files (recorder_calls_unablated), per model and
                     condition: calls; calls that passed a restricted layer list (a list
                     that is not all four layers); and calls whose returned
                     evidence_score exceeded the returned attainable_here. Aggregate only
  breakends          the 40 implant breakends (32 detected, 8 missed) and 54 background
                     breakends of the synthetic control (analysis_2026-09-25: figures.json
                     and chain/background.json.gz), default window, applicable_layers =
                     all four: evidence_score, attainable_here, and the sha256 of the
                     ceiling fields serialised with sorted keys
  identity           after the fix, the breakends' ceiling-field hashes set against
                     before.json, position by position
  proof              the prepared-payload proof of the ceiling experiment
                     (phase10_rerun.py prove, run separately), its locus 1 (IMP01
                     chr20:200000): the model-visible payload (chat.shrink_for_model, home
                     replaced by ~) and the ceiling fields, set against the 2026-09-26 proof
  panel              ui.assess at IMP01 chr20:200000 (the evidence panel): score and
                     attainable ceiling
"""
import asyncio
import glob
import gzip
import hashlib
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
from stage1_igv_assistant import server as ev  # noqa: E402

HOME = os.path.expanduser("~")
OUTDIR = os.path.join(REPO, "stage1_igv_assistant", "results", "ceiling_layers_fix_2026-09-30")
AN = os.path.join(REPO, "stage1_igv_assistant", "results", "synthetic_control_2026-09", "analysis_2026-09-25")
P10 = os.path.join(REPO, "stage1_igv_assistant", "benchmark", "runs", "phase10_rerun_2026-09-25")
BAMS = os.path.join(HOME, "public_data", "sim", "bams")
BG = os.path.join(HOME, "public_data", "NA12878.chr20_chr21.bam")
ALL4 = ["discordant_pairs", "soft_clipped_reads", "split_reads", "read_depth"]
LISTS = {"all four": ALL4,
         "discordant_pairs + soft_clipped_reads": ["discordant_pairs", "soft_clipped_reads"],
         "soft_clipped_reads + split_reads": ["soft_clipped_reads", "split_reads"],
         "discordant_pairs + soft_clipped_reads + read_depth": ["discordant_pairs", "soft_clipped_reads", "read_depth"]}
EXPECTED = {"all four": 57.5, "discordant_pairs + soft_clipped_reads": 65.0,
            "soft_clipped_reads + split_reads": 100.0,
            "discordant_pairs + soft_clipped_reads + read_depth": 43.3}
CEILING = ["attainable_ceiling_derivable", "score_bands", "strong_band", "max_all_layers",
           "max_with_flat_depth", "attainable_here", "strong_band_reachable_here", "attainable_basis",
           "attainable_note", "ceiling_counted_layers"]


def definitions():
    doc = __doc__.split("FIGURES AND THEIR DEFINITIONS")[1]
    out, key = {}, None
    for line in doc.splitlines():
        if line.startswith("  ") and not line.startswith("    ") and line.strip():
            key, _, rest = line.strip().partition(" ")
            out[key] = rest.strip()
        elif key and line.strip():
            out[key] += " " + line.strip()
    return {k: " ".join(v.split()) for k, v in out.items()}


def summary(bam, chrom, pos, layers):
    async def go():
        args = {"bam_path": bam, "chromosome": chrom, "position": pos}
        if layers is not None:
            args["applicable_layers"] = layers
        r = await ev.mcp.call_tool("breakpoint_evidence_summary", args)
        return r.structured_content
    return asyncio.run(go())


def counted(r):
    return [l for l in (r.get("applicable_layers") or []) if l not in (r.get("unassessable_layers") or {})]


def ceiling_fields(r):
    return {k: r[k] for k in CEILING if k in r}


def fhash(r):
    return hashlib.sha256(json.dumps(ceiling_fields(r), sort_keys=True).encode()).hexdigest()


def imp01_calls():
    out = {}
    for name, layers in LISTS.items():
        r = summary(os.path.join(BAMS, "IMP01.bam"), "chr20", 200000, layers)
        out[name] = {"applicable_layers_passed": layers, "counted_layers": counted(r),
                     "evidence_score": r.get("evidence_score"), "evidence_strength": r.get("evidence_strength"),
                     "components": {k: r.get(k) for k in ("discordant_pair_score", "soft_clip_score",
                                                          "split_read_score", "depth_score")},
                     "ceiling_fields": ceiling_fields(r),
                     "score_exceeds_attainable": (r.get("evidence_score") is not None and r.get("attainable_here")
                                                  is not None and r["evidence_score"] > r["attainable_here"])}
    return out


def phase10():
    cells = defaultdict(lambda: {"summary_calls": 0, "restricted_layer_list": 0, "score_above_attainable": 0,
                                 "restricted_and_above": 0})
    for f in sorted(glob.glob(os.path.join(P10, "*", "*.json"))):
        try:
            r = json.load(open(f))
        except ValueError:
            continue
        if not isinstance(r, dict) or not r.get("recorder_calls_unablated"):
            continue
        key = f"{r.get('model')} {r.get('condition')}"
        for c in r["recorder_calls_unablated"]:
            if c.get("tool") != "breakpoint_evidence_summary":
                continue
            p, res = c.get("params") or {}, c.get("result") or {}
            lst = p.get("applicable_layers")
            restricted = lst is not None and sorted(lst) != sorted(ALL4)
            above = (res.get("evidence_score") is not None and res.get("attainable_here") is not None
                     and res["evidence_score"] > res["attainable_here"])
            cells[key]["summary_calls"] += 1
            cells[key]["restricted_layer_list"] += restricted
            cells[key]["score_above_attainable"] += above
            cells[key]["restricted_and_above"] += restricted and above
    tot = {k: sum(v[k] for v in cells.values()) for k in ("summary_calls", "restricted_layer_list",
                                                           "score_above_attainable", "restricted_and_above")}
    return {"per_model_condition": dict(sorted(cells.items())), "total": tot}


def positions():
    f = json.load(open(os.path.join(AN, "figures.json")))
    pos = []
    for r in f["ceiling"]:
        imp = r["junction"].split()[0]
        c, p = r["breakend"].split(":")
        pos.append(("detected", imp, os.path.join(BAMS, f"{imp}.bam"), c, int(p)))
    for imp, v in f["per_implant"].items():
        if not any(j.get("detected") for j in v["junctions"].values()):
            for c in ("chr20", "chr21"):
                pos.append(("missed", imp, os.path.join(BAMS, f"{imp}.bam"), c, int(v["breakpoints"][c])))
    d = json.load(gzip.open(os.path.join(AN, "chain", "background.json.gz")))
    seen = set()
    for e in d["evidence"]:
        for be in e["breakends"].values():
            k = (be["chromosome"], be["position"])
            if k not in seen:
                seen.add(k)
                pos.append(("background", "background", BG, be["chromosome"], be["position"]))
    return pos


def breakends():
    rows = []
    for group, label, bam, c, p in positions():
        r = summary(bam, c, p, ALL4)
        rows.append({"group": group, "label": label, "position": f"{c}:{p}", "counted_layers": len(counted(r)),
                     "evidence_score": r.get("evidence_score"), "attainable_here": r.get("attainable_here"),
                     "ceiling_sha256": fhash(r)})
    n = defaultdict(int)
    for r in rows:
        n[r["group"]] += 1
    return {"counts": dict(n), "rows": rows}


def write(name, rec):
    os.makedirs(OUTDIR, exist_ok=True)
    dest = os.path.join(OUTDIR, name)
    if os.path.exists(dest):
        sys.exit(f"STOPPED: {name} exists; a committed record is never overwritten")
    s = json.dumps(rec, indent=1).replace(HOME, "~")
    open(dest, "w").write(s)
    print("written:", os.path.relpath(dest, REPO))


def step_before():
    rec = {"what": "Ceiling fields vs counted layers, BEFORE the fix (Phase 20 Step 1)", "definitions": definitions(),
           "code": "scripts/ceiling_layers_fix.py before", "imp01_calls": imp01_calls(), "phase10": phase10(),
           "breakends": breakends()}
    write("before.json", rec)
    print(json.dumps({k: {x: v[x] for x in ("counted_layers", "evidence_score", "score_exceeds_attainable")}
                      | {"attainable_here": v["ceiling_fields"].get("attainable_here")}
                      for k, v in rec["imp01_calls"].items()}, indent=1))
    print(json.dumps(rec["phase10"], indent=1))
    print("breakends:", rec["breakends"]["counts"])
    return 0


def step_after():
    before = json.load(open(os.path.join(OUTDIR, "before.json")))
    calls = imp01_calls()
    agree = {k: {"attainable_here": v["ceiling_fields"].get("attainable_here"), "expected": EXPECTED[k],
                 "matches_registration": v["ceiling_fields"].get("attainable_here") == EXPECTED[k],
                 "score_within_ceiling": not v["score_exceeds_attainable"]} for k, v in calls.items()}
    b = breakends()
    old = {(r["group"], r["position"]): r["ceiling_sha256"] for r in before["breakends"]["rows"]}
    diff = [f"{r['group']} {r['position']}" for r in b["rows"] if old.get((r["group"], r["position"])) != r["ceiling_sha256"]]
    rec = {"what": "Ceiling fields vs counted layers, AFTER the fix (Phase 20 Step 5)", "definitions": definitions(),
           "code": "scripts/ceiling_layers_fix.py after", "imp01_calls": calls, "imp01_agreement": agree,
           "breakends": b, "identity": {"positions": len(b["rows"]), "matched_before": len(old),
                                        "ceiling_fields_changed": diff, "all_identical": not diff and len(old) == len(b["rows"])},
           "proof": proof(), "panel": panel()}
    write("after.json", rec)
    print(json.dumps(agree, indent=1))
    print("identity:", rec["identity"]["all_identical"], len(diff), "changed")
    print("proof:", json.dumps(rec["proof"], indent=1))
    print("panel:", rec["panel"])
    return 0


def proof():
    from stage1_igv_assistant import chat
    old = json.load(open(os.path.join(P10, "proof_prepared_payload_2026-09-26.json")))
    news = sorted(glob.glob(os.path.join(P10, "proof_prepared_payload_2026-09-30*.json")))
    if not news:
        return {"run": False, "reason": "no 2026-09-30 proof found; run phase10_rerun.py prove first"}
    new = json.load(open(news[-1]))
    norm = lambda x: json.loads(json.dumps(x, default=str).replace(HOME, "~"))  # noqa: E731
    a, b = norm(old["loci"][0]["with_result"]), norm(new["loci"][0]["with_result"])
    va = chat.shrink_for_model("breakpoint_evidence_summary", a)
    vb = chat.shrink_for_model("breakpoint_evidence_summary", b)
    return {"run": True, "file": os.path.basename(news[-1]), "verdict": new.get("verdict"),
            "checks": f"{sum(c['holds'] for c in new['checks'])}/{len(new['checks'])}",
            "locus": new["loci"][0]["args"],
            "model_visible_payload_identical": json.dumps(va, sort_keys=True) == json.dumps(vb, sort_keys=True),
            "ceiling_fields_identical": ceiling_fields(a) == ceiling_fields(b),
            "ceiling_keys_identical": old.get("ceiling_keys") == new.get("ceiling_keys")}


def panel():
    from stage1_igv_assistant import ui
    ui.discover_public()
    E = ui.assess("IMP01", "chr20", 200000)
    return {"score": E["summary"]["evidence_score"], "band": E["summary"]["evidence_strength"],
            "attainable_here": E["ceiling"].get("attainable_here"),
            "ceiling_counted_layers": E["ceiling"].get("counted_layers")}


if __name__ == "__main__":
    a = sys.argv[1:]
    if a == ["before"]:
        sys.exit(step_before())
    if a == ["after"]:
        sys.exit(step_after())
    sys.exit("usage: ceiling_layers_fix.py before | after")
