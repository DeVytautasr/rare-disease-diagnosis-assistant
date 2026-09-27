#!/usr/bin/env python3
"""Phase 14 Task 3: re-measure before and after the decoy-partner fix.

    decoy_fix_remeasure.py public            (a) background survivors, (b) implants, (c) the Phase 10 locus
    decoy_fix_remeasure.py patient LABEL     (d) the same 100 seeded positions as Phase 13 -- totals only
    decoy_fix_remeasure.py record            -> results/decoy_fix_2026-09-27/remeasure.json

BEFORE is stage1_igv_assistant/tools/bam_tools.py at PRE_FIX (7bb4379), loaded as a
separate module from the git object; AFTER is the committed fixed file. Both run in
this process on the same inputs, as registered in
results/decoy_fix_2026-09-27/registration.json.

FIGURES AND THEIR DEFINITIONS (written into the record beside each figure)
  summary          summarize_breakpoint_evidence at its defaults (window 500, split
                   window 200, min_mapq 20) with all four layers applicable:
                   evidence_score, evidence_strength, split_read_score and the
                   supporting_observations sentence(s) mentioning split reads
  split_mapq20     get_split_reads(window 200, min_mapq 20): split_reads before and
                   after, and after the fix decoy_only_reads and decoy_partners
  split_mapq0      the same at min_mapq 0, the MCP split_reads tool's default
  changed          a position counts as changed when its score, band, split
                   component or split-read sentence differs between before and after
  phase10_locus    IMP01 chr20:200,000 through ui._chat_exec (the path the Phase 10
                   runs used) now, set against the full return recorded in the
                   Phase 10 proof (proof_prepared_payload_2026-09-26.json, loci[0]):
                   fields differing other than inside split_reads' new keys, and
                   whether chat.shrink_for_model gives the identical model-visible
                   payload and keeps every ceiling key
  patient          the 100 positions per BAM of Phase 13 (seed '20260927:<label>',
                   uniform over chr1-22): totals of the figures above over the
                   positions, and counts of positions whose score, band or split
                   component changed. No position is kept
"""
import importlib.util
import json
import os
import random
import subprocess
import sys
import tempfile
from collections import Counter

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
from common import DEID_DIR, RUN_DIR, die  # noqa: E402
import decoy_partners as dp  # noqa: E402  (Phase 13: positions and seeds)
from stage1_igv_assistant.tools import bam_tools as after  # noqa: E402

PRE_FIX = "7bb4379"
LAYERS = dp.LAYERS
OUTDIR = os.path.join(REPO, "stage1_igv_assistant", "results", "decoy_fix_2026-09-27")
LOGDIR = os.path.expanduser("~/public_data/sim/logs/phase14_2026-09-27")


def load_before():
    src = subprocess.run(["git", "-C", REPO, "show", f"{PRE_FIX}:stage1_igv_assistant/tools/bam_tools.py"],
                         capture_output=True, text=True, check=True).stdout
    path = os.path.join(tempfile.mkdtemp(), "bam_tools_prefix.py")
    open(path, "w").write(src)
    spec = importlib.util.spec_from_file_location("bam_tools_prefix", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


before = load_before()


def definitions():
    doc = __doc__.split("FIGURES AND THEIR DEFINITIONS (written into the record beside each figure)")[1]
    out, key = {}, None
    for line in doc.splitlines():
        if line.startswith("  ") and not line.startswith("    ") and line.strip():
            key, _, rest = line.strip().partition(" ")
            out[key] = rest.strip()
        elif key and line.strip():
            out[key] += " " + line.strip()
    return {k: " ".join(v.split()) for k, v in out.items()}


def summ(mod, bam, c, p):
    s = mod.summarize_breakpoint_evidence(bam, c, p, applicable_layers=LAYERS)
    return {"score": s.get("evidence_score"), "band": s.get("evidence_strength"),
            "split_component": s.get("split_read_score"),
            "split_sentence": [o for o in s.get("supporting_observations", []) if "split read" in o]}


def measure(bam, c, p):
    b, a = summ(before, bam, c, p), summ(after, bam, c, p)
    row = {"before": b, "after": a,
           "changed": {k: b[k] != a[k] for k in ("score", "band", "split_component", "split_sentence")}}
    for q in (20, 0):
        sb = before.get_split_reads(bam, c, p, window_bp=200, min_mapq=q)
        sa = after.get_split_reads(bam, c, p, window_bp=200, min_mapq=q)
        row[f"split_mapq{q}"] = {"before": sb.get("split_reads"), "after": sa.get("split_reads"),
                                 "decoy_only_reads": sa.get("decoy_only_reads"),
                                 "decoy_partners": sa.get("decoy_partners")}
    return row


def phase10_locus():
    from stage1_igv_assistant import chat, ui
    ui.discover_public()
    now = ui._chat_exec("breakpoint_evidence_summary", {"dataset": "IMP01", "chromosome": "chr20", "position": 200000})
    proof = json.load(open(os.path.join(REPO, "stage1_igv_assistant", "benchmark", "runs", "phase10_rerun_2026-09-25",
                                        "proof_prepared_payload_2026-09-26.json")))
    rec = proof["loci"][0]["with_result"]
    ceiling = proof["ceiling_keys"]
    now_n = json.loads(json.dumps(now, default=str).replace(os.path.expanduser("~"), "~"))
    diff = sorted(k for k in set(now_n) | set(rec) if k != "split_reads" and now_n.get(k) != rec.get(k))
    split_diff = sorted(k for k in set(now_n.get("split_reads", {})) | set(rec.get("split_reads", {}))
                        if now_n.get("split_reads", {}).get(k) != rec.get("split_reads", {}).get(k))
    vis_now = chat.shrink_for_model("breakpoint_evidence_summary", now_n)
    vis_rec = chat.shrink_for_model("breakpoint_evidence_summary", rec)
    return {"fields_differing_outside_split_reads": diff,
            "split_reads_keys_differing": split_diff,
            "new_split_fields": {k: now_n.get("split_reads", {}).get(k) for k in ("decoy_partners", "decoy_only_reads")},
            "model_visible_identical": vis_now == vis_rec,
            "model_visible_chars": [len(json.dumps(vis_rec)), len(json.dumps(vis_now))],
            "ceiling_keys_in_model_visible_payload": sorted(k for k in ceiling if k in vis_now),
            "ceiling_keys_expected": sorted(ceiling)}


def step_public():
    rows = []
    for group, label, bam, c, p in dp.synthetic_positions():
        row = {"group": group, "label": label, "position": f"{c}:{p}", **measure(bam, c, p)}
        rows.append(row)
        print(group, label, row["position"], row["before"]["score"], "->", row["after"]["score"],
              row["split_mapq20"], row["changed"], flush=True)
    rec = {"positions": rows, "phase10_locus": phase10_locus()}
    agg = {}
    for g in ("detected", "missed", "background"):
        rr = [r for r in rows if r["group"] == g]
        agg[g] = {"positions": len(rr),
                  "changed": {k: sum(r["changed"][k] for r in rr) for k in ("score", "band", "split_component",
                                                                         "split_sentence")},
                  **{f"split_mapq{q}_{w}": sum(r[f"split_mapq{q}"][w] or 0 for r in rr)
                     for q in (20, 0) for w in ("before", "after", "decoy_only_reads")}}
    rec["aggregate"] = agg
    os.makedirs(LOGDIR, exist_ok=True)
    json.dump(rec, open(os.path.join(LOGDIR, "remeasure_public.json"), "w"), indent=1)
    print(json.dumps(agg, indent=1))
    print(json.dumps(rec["phase10_locus"], indent=1))
    return 0


def step_patient(label):
    bam = os.path.join(DEID_DIR, f"{label}.bam")
    if not os.path.exists(bam):
        die(f"{label}: the de-identified link is missing")
    with pysam.AlignmentFile(bam) as f:
        lengths = {c: f.get_reference_length(c) for c in dp.AUTOSOMES}
    rng = random.Random(f"{dp.SEED}:{label}")
    total = sum(lengths.values())
    agg = Counter()
    for _ in range(dp.N_POS):
        x = rng.randrange(total)
        for c in dp.AUTOSOMES:
            if x < lengths[c]:
                break
            x -= lengths[c]
        r = measure(bam, c, x)
        agg["positions"] += 1
        for k, v in r["changed"].items():
            agg[f"positions_changed_{k}"] += v
        for q in (20, 0):
            for w in ("before", "after", "decoy_only_reads"):
                agg[f"split_mapq{q}_{w}"] += r[f"split_mapq{q}"][w] or 0
        agg["band_strong_before"] += r["before"]["band"] == "strong"
        agg["band_strong_after"] += r["after"]["band"] == "strong"
    rec = {"label": label, "seed": f"{dp.SEED}:{label}", **dict(agg)}
    os.makedirs(RUN_DIR, exist_ok=True)
    json.dump(rec, open(os.path.join(RUN_DIR, f"decoy_fix_{label}.json"), "w"), indent=1)
    print(json.dumps(rec, indent=1))
    return 0


def step_record():
    pub = json.load(open(os.path.join(LOGDIR, "remeasure_public.json")))
    out = {"what": "Before/after re-measurement of the decoy-partner fix (Phase 14 Task 3, 2026-09-27)",
           "before": f"bam_tools.py at {PRE_FIX}", "after": "bam_tools.py as committed with the fix",
           "definitions": definitions(), "registration": "registration.json",
           "public": pub, "code": "scripts/patient/decoy_fix_remeasure.py",
           "blinding": "patient figures are totals over the 100 seeded positions; no position is kept"}
    for label in ("SAMPLE_A", "SAMPLE_B"):
        p = os.path.join(RUN_DIR, f"decoy_fix_{label}.json")
        if not os.path.exists(p):
            die(f"{label}: not run")
        out[label] = json.load(open(p))
    dest = os.path.join(OUTDIR, "remeasure.json")
    if os.path.exists(dest):
        die("the record exists")
    json.dump(out, open(dest, "w"), indent=1)
    print("written:", os.path.relpath(dest, REPO))
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if a == ["public"]:
        sys.exit(step_public())
    if len(a) == 2 and a[0] == "patient" and a[1] in ("SAMPLE_A", "SAMPLE_B"):
        sys.exit(step_patient(a[1]))
    if a == ["record"]:
        sys.exit(step_record())
    die("usage: decoy_fix_remeasure.py public | patient SAMPLE_A|SAMPLE_B | record")
