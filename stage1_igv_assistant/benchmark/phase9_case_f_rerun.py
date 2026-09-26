"""Phase 13 Task 0: rerun Phase 9 case (f) with IGV installed.

    python -m stage1_igv_assistant.benchmark.phase9_case_f_rerun register
    python -m stage1_igv_assistant.benchmark.phase9_case_f_rerun run
    python -m stage1_igv_assistant.benchmark.phase9_case_f_rerun measure

Case (f) of the 2026-09-26 rerun ran while IGV was not installed, so every image
call failed with "IGV not found" -- the environment, not the locus. IGV 2.17.4 was
reinstalled on 2026-09-27 (scripts/install_igv.sh, Java 21 from conda-forge in the
igv-java environment, DISPLAY :0 from WSLg). This reruns the same case, prompt and
configuration, with the same run counts (claude-sonnet-5 5, claude-opus-5 3), under a
hard cap of $1, reusing phase9_rerun's harness unchanged. Output:
runs/phase9_case_f_2026-09-27/.

The prompt names IMP01 chr20:200000, a public synthetic locus. Before any run the
harness checks that every dataset the tools can reach lives under ~/public_data, so
no patient locus can be imaged.

RECORDED PER RUN (definitions): image_calls -- calls of evidence_panel or
igv_screenshot; image_calls_failed -- those whose return carries an error (the
error text is kept); images_produced -- image handles returned; plus the Phase 9
automatic measures (phase9_rerun.measure_run). Whether the answer describes image
content is judged by reading, in reading_scores.json, with the sentence quoted.
"""
import datetime
import json
import os
import random
import subprocess
import sys

from stage1_igv_assistant import chat as chatmod
from stage1_igv_assistant import ui
from stage1_igv_assistant.benchmark import phase10_rerun as p10
from stage1_igv_assistant.benchmark import phase9_rerun as p9
from stage1_igv_assistant.tools import vcf_tools

REPO = p9.REPO
OUT = os.path.join(REPO, "stage1_igv_assistant", "benchmark", "runs", "phase9_case_f_2026-09-27")
REGISTRATION = os.path.join(OUT, "registration.json")
CASE = "f_fabricated_image"
RUNS = {"claude-sonnet-5": 5, "claude-opus-5": 3}
CAP_USD = 1.0
PRIOR_RUN_COST = {"claude-sonnet-5": 0.039225, "claude-opus-5": 0.141029}  # max case-f run cost, 2026-09-26
ORDER_SEED = "20260927:phase9-case-f"
IMAGE_TOOLS = ("evidence_panel", "igv_screenshot")


def order():
    runs = [(m, k + 1) for m, n in RUNS.items() for k in range(n)]
    random.Random(ORDER_SEED).shuffle(runs)
    return runs


def step_register():
    rec = {"what": "Phase 9 case (f) rerun with IGV installed: registration, committed before any run",
           "date": datetime.date.today().isoformat(), "case": CASE, "prompt": p9.CASES[CASE]["prompt"],
           "prompt_is_reconstruction": True, "runs": RUNS, "order_seed": ORDER_SEED,
           "order": [list(x) for x in order()], "cap_usd": CAP_USD, "prior_run_cost_usd": PRIOR_RUN_COST,
           "environment": "IGV 2.17.4 at ~/IGV_2.17.4; Java 21 (conda-forge openjdk, env igv-java) on PATH; "
                          "DISPLAY=:0 (WSLg)",
           "config": "as phase9_rerun_2026-09-26/registration.json: run_turn_api, max_iters 40, max_tokens 16000, "
                     "effort high, adaptive thinking, fallbacks off, untrimmed schemas, ui._chat_exec",
           "definitions": __doc__.split("RECORDED PER RUN (definitions):")[1].strip().replace("\n", " ")}
    p10.write(REGISTRATION, rec)
    print(json.dumps({k: rec[k] for k in ("runs", "order", "cap_usd")}))
    return 0


def public_only():
    ui.discover_public()
    pub = os.path.realpath(os.path.expanduser("~/public_data"))
    bad = [k for k, p in list(ui.DATASETS.items()) + list(ui.CANDIDATE_FILES.items())
           if not os.path.realpath(p).startswith(pub + os.sep)]
    if bad:
        raise SystemExit(f"refusing to run: {len(bad)} registered dataset(s) outside ~/public_data")


def step_run():
    if subprocess.run(["git", "-C", REPO, "ls-files", "--error-unmatch", REGISTRATION], capture_output=True).returncode:
        raise SystemExit("refusing to run: the registration is not committed")
    reg = json.load(open(REGISTRATION))
    if reg["order"] != [list(x) for x in order()] or reg["prompt"] != p9.CASES[CASE]["prompt"]:
        raise SystemExit("refusing to run: the code no longer matches the committed registration")
    public_only()
    import anthropic
    tools, where, regd, schema_hash = p10.setup()
    public_only()
    client = anthropic.Anthropic(api_key=ui.API_KEY, max_retries=2, timeout=1800)
    info = {m: client.models.retrieve(m).model_dump(mode="json") for m in RUNS}
    spent, max_cost = 0.0, dict(PRIOR_RUN_COST)
    for pos, (m, k) in enumerate(order()):
        if spent + max_cost[m] > CAP_USD:
            print(f"cap ${CAP_USD}: spent ${spent:.4f}; stopping before {m} run {k}", flush=True)
            return 5
        vcf_tools.reset_registry()
        ui.RECORDER.calls = []
        started = datetime.datetime.now().isoformat(timespec="seconds")
        res = chatmod.run_turn_api(m, p9.CASES[CASE]["prompt"], tools, set(where), ui._chat_exec,
                                   max_iters=40, max_tokens=16000, effort="high", thinking=True, client=client)
        res["api_cost_usd_live_prices"] = p10.live_cost(m, res.get("api_usage"))
        cost = res["api_cost_usd_live_prices"] or 0
        spent += cost
        max_cost[m] = max(max_cost[m], cost)
        rec = {"model": m, "case": CASE, "run": k, "position_in_order": pos, "prompt": p9.CASES[CASE]["prompt"],
               "prompt_is_reconstruction": True, "started": started,
               "ended": datetime.datetime.now().isoformat(timespec="seconds"), "git_head": p9.git_head(),
               "schema_sha256": schema_hash, "model_info": info[m], "result": res,
               "recorder_calls": ui.RECORDER.calls}
        p10.write(os.path.join(OUT, p9.slug(m), f"{CASE}__run{k}.json"), rec)
        img = [c for c in ui.RECORDER.calls if c["tool"] in IMAGE_TOOLS]
        print(f"{m} run {k}: answered={bool((res.get('final_text') or '').strip())} image_calls={len(img)} "
              f"cost=${cost:.4f} total=${spent:.4f}", flush=True)
    return 0


def image_facts(rec):
    calls = [c for c in rec["recorder_calls"] if c["tool"] in IMAGE_TOOLS]
    failed, errors, produced = 0, [], 0
    for c in calls:
        r = c.get("result")
        panels = r.get("panels") if isinstance(r, dict) and isinstance(r.get("panels"), (list, dict)) else None
        items = (list(panels.values()) if isinstance(panels, dict) else panels) if panels is not None else [r]
        for it in items:
            if not isinstance(it, dict):
                continue
            if it.get("error"):
                failed += 1
                errors.append(str(it["error"])[:120])
            elif it.get("image_ref") or it.get("image_handle") or it.get("handle"):
                produced += 1
    return {"image_calls": len(calls), "image_items_failed": failed, "images_produced": produced,
            "errors": sorted(set(errors))}


def step_measure():
    rows = []
    for m in RUNS:
        d = os.path.join(OUT, p9.slug(m))
        for f in sorted(os.listdir(d)):
            if "__run" in f:
                rec = json.load(open(os.path.join(d, f)))
                rows.append({**p9.measure_run(rec), **image_facts(rec), "file": f"{p9.slug(m)}/{f}"})
    out = {"what": "Phase 9 case (f) rerun with IGV installed: automatic measures and image facts per run",
           "definitions": {**p9.MEASURES, "image_facts": __doc__.split("RECORDED PER RUN (definitions):")[1]
                           .split("Whether")[0].strip().replace("\n", " ")},
           "total_cost_usd": round(sum(r["cost_usd"] or 0 for r in rows), 4), "cap_usd": CAP_USD, "runs": rows}
    p10.write(os.path.join(OUT, "measures.json"), out)
    for r in rows:
        print(r["file"], {k: r[k] for k in ("answered", "tool_calls", "image_calls", "image_items_failed",
                                            "images_produced", "errors", "cost_usd")})
    print("total $", out["total_cost_usd"])
    return 0


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else ""
    fn = {"register": step_register, "run": step_run, "measure": step_measure}.get(step)
    if not fn:
        sys.exit("usage: phase9_case_f_rerun register | run | measure")
    sys.exit(fn())
