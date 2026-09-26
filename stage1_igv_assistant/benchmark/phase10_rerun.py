"""
phase10_rerun.py -- the Phase 10 ceiling experiment, rerun on 2026-09-25 against
the rebuilt implants. Case (e) only.

    python -m stage1_igv_assistant.benchmark.phase10_rerun prove
    python -m stage1_igv_assistant.benchmark.phase10_rerun project
    python -m stage1_igv_assistant.benchmark.phase10_rerun run MODEL [MODEL ...] [--runs 5]
    python -m stage1_igv_assistant.benchmark.phase10_rerun extract

DESIGN. One question -- whether the balanced translocation at IMP01 chr20:200,000
is strong evidence -- put to four models (qwen2.5:7b and qwen3.5:4b through
Ollama; claude-sonnet-5 and claude-opus-5 through the Anthropic API) in two
conditions, WITHOUT and WITH the attainable-ceiling fields, 5 runs each.

ONE VARIABLE. The conditions differ only in the keys that the Phase 10 commit
(912c3c9) added to breakpoint_evidence_summary's return through
server._with_ceiling. The keys are read out of that commit's own diff (never
restated here) and cross-checked against tests/test_ceiling_echo.py. The switch
lives in THIS harness: in the WITHOUT condition it strips exactly those keys from
the recorded tool return before the model sees it. The tool is not touched. The
same commit also added position_provenance and min_mapq_applied to the same
return and a min_mapq paragraph to split_reads' description; all three are kept
in BOTH conditions, so they are not a difference between them. (The original
Phase 10 comparison set Phase 9's runs, taken before that commit, against runs
after it, so it confounded all four.) `prove` diffs the two returns at the same
locus -- as the tool returns them and as the model receives them after
chat.shrink_for_model -- and they must differ in those keys and nothing else,
with a planted difference to show the diff can fire and the other tools' returns
to show the switch leaves them alone. `run` refuses to start unless `prove`
passed on the current code.

CONFIGURATION, as the Phase 8 record (results/phase8_final_record.json): context
32768, max_iters 40, untrimmed tool descriptions, think='low' where Ollama's
/api/show declares the thinking capability (omitted otherwise; chat.run_turn
asks). The API models run through chat.run_turn_api with the committed defaults
of the interface's /api/chat route -- adaptive thinking, effort 'high',
max_tokens 16000, max_iters 40 -- and server-side fallbacks off (it never sets
them). The same system prompt (chat.SYSTEM_PROMPT), the same schemas
(ui._chat_tools, generated from list_tools), the same executor
(ui._chat_exec: labels resolved to paths inside the executor, results scrubbed)
and the same datasets (ui.discover_public). The bridge registry and the tool
recorder are reset before every run, so no run inherits another's coordinate
provenance. Conditions are interleaved in a seeded order per model.

THE PROMPT is a RECONSTRUCTION. Case (e)'s exact wording was never committed
and the Phase 8-10 run records were lost in the 2026-09-23 reinstall; the text
below is rebuilt from the recorded description of the case and is the same in
both conditions.

A model id the API no longer serves is reported and NOT substituted.
"""
from __future__ import annotations

import argparse
import ast
import datetime
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import time
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from stage1_igv_assistant import chat as chatmod  # noqa: E402
from stage1_igv_assistant import ui  # noqa: E402
from stage1_igv_assistant.tools import vcf_tools  # noqa: E402

PHASE10_COMMIT = "912c3c9"
OUT = os.path.join(REPO, "stage1_igv_assistant", "benchmark", "runs", "phase10_rerun_2026-09-25")
HOME = os.path.expanduser("~")
CASE_E = ("Is the balanced translocation at chr20:200000 in dataset IMP01 strong evidence? "
          "Use the tools to find out, and explain your answer.")
CASE_E_PROVENANCE = ("RECONSTRUCTION (2026-09-25): case (e)'s committed prompt text does not exist; rebuilt from "
                     "the recorded description 'whether the balanced translocation at IMP01 chr20:200,000 is "
                     "strong evidence'")
LOCAL = ("qwen2.5:7b", "qwen3.5:4b")
API = ("claude-sonnet-5", "claude-opus-5")
CONFIG = {"num_ctx": 32768, "max_iters": 40, "trim_descriptions": False, "think": "low",
          "api": {"max_tokens": 16000, "effort": "high", "thinking": "adaptive (display summarized)",
                  "max_iters": 40, "server_side_fallbacks": False}}
ORDER_SEED = 20260925
EXPECTED_OLLAMA = "0.32.9"
# platform.claude.com/docs/en/about-claude/pricing, read 2026-09-25. Sonnet 5's
# $2/$10 is its standard price: the increase to $3/$15 announced for 2026-09-01
# did not happen. chat.API_PRICES carried $3/$15 for it until it was corrected the
# same day; a run records both chat.api_cost and the cost at these prices.
PRICES_LIVE = {"claude-sonnet-5": {"in": 2.00, "cache_write_5m": 2.50, "cache_read": 0.20, "out": 10.00},
               "claude-opus-5": {"in": 5.00, "cache_write_5m": 6.25, "cache_read": 0.50, "out": 25.00}}
PHASE10_QUOTED = {"api_spend": "$0.96 for 10 API runs",
                  "without": "1/20 assert, 0/20 complete", "with": "19/20 assert, 17/20 complete",
                  "qwen3.5:4b": "completed it 5/5", "qwen2.5:7b": "misused the ceiling in 2 of 5"}


# ── the one variable ────────────────────────────────────────────────────────

def ceiling_keys_from_commit():
    """The keys server._with_ceiling assigns, read from the Phase 10 commit's diff."""
    diff = subprocess.run(["git", "-C", REPO, "show", PHASE10_COMMIT, "--", "stage1_igv_assistant/server.py"],
                          capture_output=True, text=True, check=True).stdout
    keys, inside = [], False
    for line in diff.splitlines():
        if line.startswith("+def _with_ceiling"):
            inside = True
            continue
        if inside and (not line.startswith("+") or line[1:].startswith(("def ", "@mcp.tool"))):
            break
        if inside:
            m = re.search(r'out\["([A-Za-z_]+)"\]\s*=', line)
            if m and m.group(1) not in keys:
                keys.append(m.group(1))
    if not keys:
        raise SystemExit(f"no keys found in _with_ceiling in {PHASE10_COMMIT}'s diff")
    return keys


def new_fields_from_test():
    src = open(os.path.join(REPO, "stage1_igv_assistant", "tests", "test_ceiling_echo.py")).read()
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "NEW_FIELDS" for t in node.targets):
            return ast.literal_eval(node.value)
    raise SystemExit("NEW_FIELDS not found in test_ceiling_echo.py")


def ablate(rec, keys):
    """WITHOUT condition: the summary's return minus the ceiling keys. Every other
    tool, and every other field of the record, passes through untouched."""
    if not isinstance(rec, dict) or rec.get("tool") != "breakpoint_evidence_summary":
        return rec
    res = rec.get("result")
    if not isinstance(res, dict):
        return rec
    out = dict(rec)
    out["result"] = {k: v for k, v in res.items() if k not in keys}
    return out


SUMMARY = "breakpoint_evidence_summary"

# ── the payload the model receives ──────────────────────────────────────────
# The model loop (chat.run_turn / run_turn_api) passes every tool result through
# chat.shrink_for_model (2,600 characters) before the model sees it, and the
# ceiling keys are long enough to decide whether the summary is truncated. So the
# switch cannot simply strip keys and leave truncation to the loop:
#   * proof_strip_before_truncation_NOT_PROVEN.json -- stripping first: at IMP01 the
#     WITH payload is truncated to its keep-list while the stripped one arrives whole;
#   * qwen2.5-7b/WITHOUT__run2.json -- the second design (hand the loop a result it
#     would shrink to the target) found no such result for that call (the target was
#     2,647 characters, so the loop re-truncated it) and logged the run as not holding.
# The switch therefore prepares the exact payload -- the WITH payload as
# chat.shrink_for_model makes it, minus the ceiling keys -- and marks it; the
# wrapper below, installed in THIS process only, hands a marked payload to the
# model unchanged and sends every unmarked one (all of WITH, every other tool)
# through chat.shrink_for_model exactly as before. The tool is not touched.
_PREPARED = "_phase10_prepared_payload"
_chat_shrink = chatmod.shrink_for_model


def _shrink_passing_prepared(name, result, budget=2600):
    if isinstance(result, dict) and result.get(_PREPARED) is True:
        return {k: v for k, v in result.items() if k != _PREPARED}
    return _chat_shrink(name, result, budget)


chatmod.shrink_for_model = _shrink_passing_prepared


def shrink_is_wrapped():
    """Both loops look shrink_for_model up in chat's namespace at call time."""
    return (chatmod.run_turn.__globals__.get("shrink_for_model") is _shrink_passing_prepared
            and chatmod.run_turn_api.__globals__.get("shrink_for_model") is _shrink_passing_prepared)


def without_result(W, keys):
    """The WITHOUT payload for one summary return W: exactly what the WITH condition
    shows the model, minus the ceiling keys, marked so the loop hands it over as is."""
    target = {k: v for k, v in _chat_shrink(SUMMARY, W).items() if k not in keys}
    return dict(target, **{_PREPARED: True}), "the WITH model-visible payload minus the keys, handed over as prepared"


def visible_diff_holds(W, R, keys):
    """What the model receives in each condition, through the loop's own call."""
    d = dict_diff(chatmod.shrink_for_model(SUMMARY, W), chatmod.shrink_for_model(SUMMARY, R))
    return set(d["only_in_first"]) <= set(keys) and not d["only_in_second"] and not d["changed"], d


def make_exec(condition, keys, log):
    """WITH: the committed path untouched. WITHOUT: every summary return replaced by
    without_result(), with the model-visible difference re-checked on EVERY call
    and logged; a call where it cannot hold is logged as such, never hidden."""
    def ex(name, args):
        rec, err = ui._chat_exec(name, args)
        if (err is not None or condition == "WITH" or not isinstance(rec, dict) or rec.get("tool") != SUMMARY
                or not isinstance(rec.get("result"), dict) or not (set(rec["result"]) & set(keys))):
            return rec, err
        W = rec["result"]
        R, how = without_result(W, keys)
        if R is None:
            R = {k: v for k, v in W.items() if k not in keys}
        holds, d = visible_diff_holds(W, R, keys)
        log.append({"call": rec.get("id"), "how": how, "holds": holds, "model_visible_diff": d})
        out = dict(rec)
        out["result"] = R
        return out, None
    return ex


def dict_diff(a, b):
    ka, kb = set(a), set(b)
    changed = sorted(k for k in ka & kb
                     if json.dumps(a[k], sort_keys=True, default=str) != json.dumps(b[k], sort_keys=True, default=str))
    return {"only_in_first": sorted(ka - kb), "only_in_second": sorted(kb - ka), "changed": changed}


# ── setup and records ───────────────────────────────────────────────────────

def tilde_all(obj):
    if isinstance(obj, dict):
        return {k: tilde_all(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [tilde_all(v) for v in obj]
    if isinstance(obj, str):
        return obj.replace(HOME, "~")
    return obj


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def setup():
    ui.discover_public()
    tools, where = ui._chat_tools(trim=False)
    reg = {"datasets": sorted(ui.DATASETS), "candidates": sorted(ui.CANDIDATE_FILES),
           "dataset_paths": {k: v.replace(HOME, "~") for k, v in sorted(ui.DATASETS.items())},
           "candidate_paths": {k: v.replace(HOME, "~") for k, v in sorted(ui.CANDIDATE_FILES.items())}}
    return tools, where, reg, sha(json.dumps(tools, sort_keys=True))


def write(path, obj):
    if os.path.exists(path) and subprocess.run(["git", "-C", REPO, "ls-files", "--error-unmatch", path],
                                               capture_output=True).returncode == 0:
        raise SystemExit(f"refusing to overwrite a committed record: {os.path.relpath(path, REPO)}")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as f:
        json.dump(tilde_all(obj), f, indent=1, default=str)
    os.replace(path + ".tmp", path)


def git_head():
    return subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()


def ollama_get(path):
    with urllib.request.urlopen(chatmod.OLLAMA + path, timeout=30) as r:
        return json.loads(r.read())


def unload_others(model):
    """One local model on the GPU at a time: chat.run_turn keeps a model resident
    for 600 s, so the previous model would otherwise still hold VRAM when the next
    one loads and change its CPU/GPU split. Returns what was unloaded."""
    gone = []
    for m in ollama_get("/api/ps").get("models", []):
        if m.get("name") != model:
            chatmod._post("/api/generate", {"model": m["name"], "keep_alive": 0}, timeout=120)
            gone.append(m["name"])
    return gone


def nvidia():
    try:
        return subprocess.run(["/usr/lib/wsl/lib/nvidia-smi", "--query-gpu=name,memory.total,memory.used",
                               "--format=csv,noheader"], capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception as e:
        return f"unavailable: {e}"


# ── prove ───────────────────────────────────────────────────────────────────

PROOF_LOCI = [
    {"dataset": "IMP01", "chromosome": "chr20", "position": 200000},
    {"dataset": "IMP01", "chromosome": "chr21", "position": 14100000, "window_bp": 200,
     "applicable_layers": ["discordant_pairs", "soft_clipped_reads", "split_reads", "read_depth"]},
    {"dataset": "IMP10", "chromosome": "chr20", "position": 25800000},
    # the arguments qwen2.5:7b passed in WITHOUT run 2, where the second design did not hold
    {"dataset": "IMP01", "chromosome": "chr20", "position": 200000, "label": "chr20:200000_balanced_translocation",
     "window_bp": 1500, "min_mapq": 20,
     "applicable_layers": ["discordant_pairs", "soft_clipped_reads", "split_reads", "read_depth"]},
]
PROOF_FILE = "proof_prepared_payload.json"          # the first, 2026-09-25 (committed)
TODAY = datetime.date.today().isoformat()


def newest_proof():
    """The most recent proof: proof_prepared_payload.json, then the dated ones
    (proof_prepared_payload_<YYYY-MM-DD>.json), which sort after it."""
    import glob as _glob
    found = sorted(_glob.glob(os.path.join(OUT, "proof_prepared_payload*.json")))
    return found[-1] if found else os.path.join(OUT, PROOF_FILE)


def step_prove():
    keys = ceiling_keys_from_commit()
    tools, where, reg, schema_hash = setup()
    checks = []

    def check(name, ok, detail=None):
        checks.append({"check": name, "holds": bool(ok), "detail": detail})
        print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  -- {detail}" if detail is not None else ""))

    print(f"ceiling keys from {PHASE10_COMMIT}'s diff: {keys}")
    check("both model loops look up the harness's shrink wrapper", shrink_is_wrapped())
    nf = new_fields_from_test()
    check("the diff's keys are the test's NEW_FIELDS plus the failure-path reason",
          set(keys) == set(nf) | {"attainable_ceiling_reason"}, sorted(set(keys) ^ (set(nf) | {"attainable_ceiling_reason"})))
    loci = []
    for args in PROOF_LOCI:
        vcf_tools.reset_registry()
        with_rec, e1 = ui._chat_exec("breakpoint_evidence_summary", dict(args))
        second, e2 = ui._chat_exec("breakpoint_evidence_summary", dict(args))
        check(f"{args['dataset']} {args['chromosome']}:{args['position']}: both calls ran", e1 is None and e2 is None)
        without_rec = ablate(second, keys)
        W, N = with_rec["result"], without_rec["result"]
        d = dict_diff(W, N)
        present = sorted(set(keys) & set(W))
        check(f"  tool return: keys only WITH = the ceiling keys present", d["only_in_first"] == present, d["only_in_first"])
        check(f"  tool return: nothing only WITHOUT", not d["only_in_second"], d["only_in_second"])
        check(f"  tool return: no shared key differs", not d["changed"], d["changed"])
        check(f"  the WITH return carries the ceiling", W.get("attainable_ceiling_derivable") is True,
              {k: W.get(k) for k in ("attainable_here", "strong_band", "strong_band_reachable_here")})
        check(f"  the switch changes no other field of the record",
              {k: v for k, v in without_rec.items() if k != "result"} == {k: v for k, v in second.items() if k != "result"})
        # what the model receives: the WITHOUT result goes through the loop's own
        # shrink_for_model, exactly as run_turn / run_turn_api apply it
        log = []
        wo_rec, _ = make_exec("WITHOUT", keys, log)("breakpoint_evidence_summary", dict(args))
        R = wo_rec["result"]
        sW = chatmod.shrink_for_model(SUMMARY, W)
        sR = chatmod.shrink_for_model(SUMMARY, R)
        d2 = dict_diff(sW, sR)
        check(f"  model-visible payload: keys only WITH are ceiling keys", set(d2["only_in_first"]) <= set(keys),
              d2["only_in_first"])
        check(f"  model-visible payload: nothing only WITHOUT", not d2["only_in_second"], d2["only_in_second"])
        check(f"  model-visible payload: no shared key differs", not d2["changed"], d2["changed"])
        check(f"  the switch's own per-call check agrees", bool(log) and all(e["holds"] for e in log),
              [e["how"] for e in log])
        # positive control: the first design (strip, then let the loop truncate) is
        # caught by this same check wherever the keys are what forces truncation
        old = chatmod.shrink_for_model(SUMMARY, N)
        d_old = dict_diff(sW, old)
        old_differs = bool(d_old["only_in_second"] or d_old["changed"] or set(d_old["only_in_first"]) - set(keys))
        forced = "_truncated_fields" in sW and "_truncated_fields" not in old
        check(f"  positive control: strip-before-truncation {'is caught' if forced else 'is identical here'}",
              old_differs == forced, {"keys_force_truncation": forced, "old_design_differs": old_differs})
        # the second design: a result the loop's own shrink turns into the target
        target = {k: v for k, v in _chat_shrink(SUMMARY, W).items() if k not in keys}
        second_ok = _chat_shrink(SUMMARY, N) == target or _chat_shrink(SUMMARY, target) == target
        edge = args.get("window_bp") == 1500
        check(f"  positive control: the second design {'fails here, as it did in the run' if edge else 'holds here'}",
              second_ok != edge, {"second_design_holds": second_ok, "target_chars": len(json.dumps(target))})
        loci.append({"args": args, "tool_return_diff": d, "model_visible_diff": d2, "switch_log": log,
                     "old_design_model_visible_diff": d_old,
                     "truncated": {"with": "_truncated_fields" in sW, "without": "_truncated_fields" in sR},
                     "model_visible_chars": {"with": len(json.dumps(sW)), "without": len(json.dumps(sR)),
                                             "old_design": len(json.dumps(old))},
                     "with_result": W, "without_tool_return": N, "without_result_handed_to_loop": R})
    # the diff can fire: planted differences
    base = loci[0]["without_tool_return"]
    W0 = loci[0]["with_result"]
    p1 = dict(base, evidence_score=(base.get("evidence_score") or 0) + 1)
    check("positive control: a changed score is caught", dict_diff(W0, p1)["changed"] == ["evidence_score"])
    p2 = dict(base, planted_extra=1)
    check("positive control: an added key is caught", dict_diff(W0, p2)["only_in_second"] == ["planted_extra"])
    p3 = {k: v for k, v in base.items() if k != "supporting_observations"}
    got = dict_diff(W0, p3)["only_in_first"]
    check("positive control: a removed non-ceiling key is caught and is not a ceiling key",
          "supporting_observations" in got and not set(got) <= set(keys), got)
    # the switch leaves every other tool alone, and no other tool returns a ceiling key
    vcf_tools.reset_registry()
    others = [("discordant_pairs", {"dataset": "IMP01", "chromosome": "chr20", "position": 200000}),
              ("soft_clipped_reads", {"dataset": "IMP01", "chromosome": "chr20", "position": 200000}),
              ("split_reads", {"dataset": "IMP01", "chromosome": "chr20", "position": 200000}),
              ("bam_stats_at_locus", {"dataset": "IMP01", "chromosome": "chr20", "start": 199500, "end": 200500}),
              ("read_depth_profile", {"dataset": "IMP01", "chromosome": "chr20", "start": 198000, "end": 202000,
                                      "focus_position": 200000}),
              ("reciprocal_breakpoint", {"dataset": "IMP01", "primary_chromosome": "chr20", "primary_position": 200000,
                                         "partner_chromosome": "chr21", "partner_position": 14100000}),
              ("applicable_layers", {"dataset": "IMP01"}),
              ("load_candidate_set", {"candidates": "IMP01", "label": "IMP01"})]
    for tool, a in others:
        r, e = ui._chat_exec(tool, a)
        before = json.dumps(r, sort_keys=True, default=str)
        check(f"negative control: {tool} passes through the switch unchanged",
              e is None and json.dumps(ablate(r, keys), sort_keys=True, default=str) == before)
        check(f"  {tool} returns no ceiling key", not (set((r or {}).get("result") or {}) & set(keys)))
    er, _ = ui._chat_exec("breakpoint_evidence_summary", {"dataset": "IMP01", "chromosome": "chrNOPE", "position": 5})
    check("negative control: an error return from the summary passes through unchanged",
          json.dumps(ablate(er, keys), sort_keys=True, default=str) == json.dumps(er, sort_keys=True, default=str))
    desc_hits = [k for k in keys if k in json.dumps(tools)]
    check("no tool description names a ceiling key (descriptions are the same list in both conditions)",
          not desc_hits, desc_hits)
    ok = all(c["holds"] for c in checks)
    rec = {"phase10_commit": PHASE10_COMMIT, "ceiling_keys": keys, "test_new_fields": nf, "git_head": git_head(),
           "schema_sha256": schema_hash, "registration": reg, "checks": checks, "loci": loci,
           "kept_in_both_conditions": ["position_provenance", "min_mapq_applied (top level and the three sub-dicts)",
                                       "split_reads' min_mapq description paragraph"],
           "verdict": "PROVEN" if ok else "NOT PROVEN"}
    out = os.path.join(OUT, PROOF_FILE if not os.path.exists(os.path.join(OUT, PROOF_FILE))
                       else f"proof_prepared_payload_{TODAY}.json")
    write(out, rec)
    print(f"proof written to {os.path.relpath(out, REPO)}")
    print(f"ablation proof: {rec['verdict']} ({sum(c['holds'] for c in checks)}/{len(checks)} checks)")
    return 0 if ok else 1


def proof_is_current():
    p = newest_proof()
    if not os.path.exists(p):
        return False, f"no {os.path.basename(p)} -- run `prove` first"
    rec = json.load(open(p))
    if rec.get("verdict") != "PROVEN":
        return False, "the proof did not pass"
    if rec.get("ceiling_keys") != ceiling_keys_from_commit():
        return False, "the ceiling keys changed since the proof"
    return True, None


# ── cost ────────────────────────────────────────────────────────────────────

def live_cost(model, usage):
    p = PRICES_LIVE.get(model)
    if not p or not usage:
        return None
    return round((usage["input_tokens"] * p["in"] + usage["cache_creation_input_tokens"] * p["cache_write_5m"]
                  + usage["cache_read_input_tokens"] * p["cache_read"] + usage["output_tokens"] * p["out"]) / 1e6, 6)


def step_project():
    tools, where, reg, schema_hash = setup()
    api_tools = chatmod.to_anthropic_tools(tools)
    body_chars = len(json.dumps(api_tools)) + len(chatmod.SYSTEM_PROMPT) + len(CASE_E)
    rec = {"schema_sha256": schema_hash, "first_request_chars": body_chars, "prices": PRICES_LIVE}
    counted = {}
    try:
        import anthropic
        c = anthropic.Anthropic(api_key=ui.API_KEY, max_retries=1, timeout=60)
        for m in API:
            counted[m] = c.messages.count_tokens(model=m, system=chatmod.SYSTEM_PROMPT, tools=api_tools,
                                                 messages=[{"role": "user", "content": CASE_E}]).input_tokens
    except Exception as e:
        counted["error"] = f"{type(e).__name__}: {str(e)[:160]}"
    rec["count_tokens_first_request"] = counted
    # The projection. Per run: the first request writes the cached prefix; each later
    # iteration reads it back and writes the growth; output is thinking + text.
    # Assumed per run (stated, not measured here): 8 iterations, 900 new input tokens
    # per iteration (a shrunk tool result is at most 2600 characters), 1500 output
    # tokens per iteration at effort 'high'.
    first = counted.get("claude-opus-5") or round(body_chars / 3.2)
    iters, grow, outp = 8, 900, 1500
    proj = {}
    for m in API:
        p = PRICES_LIVE[m]
        f_tok = counted.get(m) or first
        reads = sum(f_tok + k * grow for k in range(1, iters))
        cost = (f_tok * p["cache_write_5m"] + reads * p["cache_read"] + (iters - 1) * grow * p["cache_write_5m"]
                + iters * outp * p["out"]) / 1e6
        proj[m] = {"per_run_usd": round(cost, 3), "ten_runs_usd": round(10 * cost, 2)}
    total = sum(v["ten_runs_usd"] for v in proj.values())
    rec.update({"assumptions": {"iterations": iters, "new_input_tokens_per_iteration": grow,
                                "output_tokens_per_iteration": outp},
                "projection": proj, "projected_total_usd_20_runs": round(total, 2),
                "phase10_scaled": {"quoted": PHASE10_QUOTED["api_spend"], "twenty_runs_usd": 1.92},
                "budget_cap_usd": round(max(3 * total, 5.0), 2)})
    write(os.path.join(OUT, "projection.json"), rec)
    print(json.dumps({k: rec[k] for k in ("count_tokens_first_request", "projection", "projected_total_usd_20_runs",
                                          "phase10_scaled", "budget_cap_usd")}, indent=1))
    return 0


# ── run ─────────────────────────────────────────────────────────────────────

def slug(model):
    return model.replace(":", "-").replace("/", "-")


BILLING_WORDS = ("credit balance", "billing", "insufficient", "payment", "quota")


def is_billing_error(detail):
    d = (detail or "").lower()
    return any(w in d for w in BILLING_WORDS)


def run_model(model, runs, keys, tools, where, reg, schema_hash, budget, spent, replace=None, start_run=1):
    """replace = (condition, run number, reason): one run in that condition, written
    beside the run it replaces as <COND>__run<k>_replacement.json; the replaced run
    stays on record."""
    is_api = model.startswith("claude-")
    meta = {"model": model, "backend": "anthropic" if is_api else "ollama", "git_head": git_head(),
            "schema_sha256": schema_hash, "registration": reg, "ceiling_keys": keys, "config": CONFIG,
            "prompt": CASE_E, "prompt_provenance": CASE_E_PROVENANCE,
            "system_prompt_sha256": sha(chatmod.SYSTEM_PROMPT)}
    client = None
    if is_api:
        import anthropic
        meta["anthropic_sdk"] = anthropic.__version__
        client = anthropic.Anthropic(api_key=ui.API_KEY, max_retries=2, timeout=1800)
        try:
            info = client.models.retrieve(model)
            meta["model_info"] = info.model_dump(mode="json")
        except anthropic.NotFoundError as e:
            meta["not_served"] = f"{e.status_code}: model id not served; NOT substituted"
            write(os.path.join(OUT, slug(model), "status.json"), meta)
            print(f"{model}: NOT SERVED -- not substituted")
            return 3, spent
        except anthropic.AuthenticationError as e:
            meta["blocked"] = f"API key rejected ({e.status_code}); no run made"
            write(os.path.join(OUT, slug(model), "status.json"), meta)
            print(f"{model}: API key rejected -- no run made")
            return 4, spent
    else:
        ver = ollama_get("/api/version").get("version")
        tags = {m["name"]: m for m in ollama_get("/api/tags").get("models", [])}
        if model not in tags:
            print(f"{model}: not pulled")
            return 3, spent
        meta["unloaded_before_start"] = unload_others(model)
        meta.update({"ollama_version": ver, "ollama_version_matches_phase8": ver == EXPECTED_OLLAMA,
                     "model_digest": tags[model]["digest"], "model_details": tags[model].get("details"),
                     "capabilities": sorted(chatmod.model_capabilities(model)),
                     "think_sent": "low" if chatmod.supports_thinking(model) else "omitted (no thinking capability)",
                     "gpu_before": nvidia()})
    if replace:
        cond_r, k_r, reason = replace
        order, suffix = [cond_r], "_replacement"
        counters = {"WITHOUT": 0, "WITH": 0, cond_r: k_r - 1}
        meta.update({"replaces": f"{cond_r}__run{k_r}.json", "reason": reason, "order": order})
        write(os.path.join(OUT, slug(model), f"meta_replacement_{cond_r}_run{k_r}.json"), meta)
    elif start_run > 1:
        # an extension: runs start_run .. start_run+runs-1 in each condition, in their own
        # seeded order, beside the runs already on record
        order, suffix = ["WITHOUT"] * runs + ["WITH"] * runs, ""
        seed = f"{ORDER_SEED}:{model}:start{start_run}"
        random.Random(seed).shuffle(order)
        counters = {"WITHOUT": start_run - 1, "WITH": start_run - 1}
        meta.update({"order": order, "order_seed": seed, "runs": f"{start_run}-{start_run + runs - 1}"})
        write(os.path.join(OUT, slug(model), f"meta_runs{start_run}-{start_run + runs - 1}.json"), meta)
    else:
        order, suffix = ["WITHOUT"] * runs + ["WITH"] * runs, ""
        random.Random(f"{ORDER_SEED}:{model}").shuffle(order)
        meta["order"] = order
        meta["order_seed"] = f"{ORDER_SEED}:{model}"
        write(os.path.join(OUT, slug(model), "meta.json"), meta)
        counters = {"WITHOUT": 0, "WITH": 0}
    max_run_cost = 0.0
    for pos, cond in enumerate(order):
        # A HARD cap: a run starts only if the most expensive run so far would still fit.
        if is_api and spent + max_run_cost > budget:
            print(f"{model}: budget cap ${budget}: spent ${spent:.4f}, the costliest run so far "
                  f"${max_run_cost:.4f} would not fit; stopping")
            return 5, spent
        counters[cond] += 1
        vcf_tools.reset_registry()
        ui.RECORDER.calls = []
        ablation_log = []
        started = datetime.datetime.now().isoformat(timespec="seconds")
        if is_api:
            res = chatmod.run_turn_api(model, CASE_E, tools, set(where), make_exec(cond, keys, ablation_log),
                                       max_iters=CONFIG["api"]["max_iters"], max_tokens=CONFIG["api"]["max_tokens"],
                                       effort=CONFIG["api"]["effort"], thinking=True, client=client)
            res["api_cost_usd_live_prices"] = live_cost(model, res.get("api_usage"))
            spent += res["api_cost_usd_live_prices"] or 0
            max_run_cost = max(max_run_cost, res["api_cost_usd_live_prices"] or 0)
        else:
            res = chatmod.run_turn(model, CASE_E, tools, set(where), make_exec(cond, keys, ablation_log),
                                   num_ctx=CONFIG["num_ctx"], max_iters=CONFIG["max_iters"], think=CONFIG["think"])
        summaries = [e for e in res["events"] if e.get("type") == "tool" and e.get("name") == "breakpoint_evidence_summary"
                     and isinstance(e.get("result"), dict)]
        rec = {**{k: meta[k] for k in ("model", "backend", "schema_sha256", "ceiling_keys", "prompt",
                                       "prompt_provenance", "system_prompt_sha256", "git_head")},
               "condition": cond, "run": counters[cond], "position_in_order": pos, "started": started,
               "ended": datetime.datetime.now().isoformat(timespec="seconds"),
               "summary_calls_seen_by_model": len(summaries),
               "ceiling_keys_seen_by_model": sorted({k for e in summaries for k in e["result"] if k in keys}),
               "ablation_log": ablation_log,
               "ablation_held_on_every_call": all(e["holds"] for e in ablation_log),
               # chat's own verify_numbers checks the prose against the recorded results; this one
               # checks it against what the model was actually sent, the same way in both conditions
               "verification_against_model_visible": chatmod.verify_numbers(
                   res["final_text"], [chatmod.shrink_for_model(e["name"], e["result"]) for e in res["events"]
                                       if e.get("type") == "tool" and isinstance(e.get("result"), dict)]),
               "result": res, "recorder_calls_unablated": ui.RECORDER.calls}
        if not is_api and pos == 0:
            try:
                rec["ollama_ps_after_first_run"] = subprocess.run(["ollama", "ps"], capture_output=True, text=True,
                                                                  timeout=30).stdout
            except Exception as e:
                rec["ollama_ps_after_first_run"] = f"unavailable: {e}"
        if replace:
            rec.update({"replaces": meta["replaces"], "replacement_reason": meta["reason"]})
        write(os.path.join(OUT, slug(model), f"{cond}__run{counters[cond]}{suffix}.json"), rec)
        print(f"{model} {cond} run {counters[cond]}: answered={not res['ended_without_answer']} "
              f"tool_calls={res['n_tool_calls']} summaries={len(summaries)} ceiling_keys_seen="
              f"{len(rec['ceiling_keys_seen_by_model'])} wall={res['wall_s']}s"
              + (f" cost=${res['api_cost_usd_live_prices']}" if is_api else "")
              + (f" API_ERROR={res.get('api_error')}" if is_api and res.get("api_error") else ""), flush=True)
        if is_api and is_billing_error(res.get("api_error")):
            print(f"{model}: BILLING ERROR -- stopping, not retrying")
            return 6, spent
    return 0, spent


def step_run(models, runs, budget, replace=None, start_run=1):
    ok, why = proof_is_current()
    if not ok:
        raise SystemExit(f"refusing to run: {why}")
    keys = ceiling_keys_from_commit()
    tools, where, reg, schema_hash = setup()
    proof = json.load(open(newest_proof()))
    if proof["schema_sha256"] != schema_hash:
        raise SystemExit("the tool schemas differ from the ones the proof ran against")
    rc, spent = 0, 0.0
    for m in models:
        r, spent = run_model(m, runs, keys, tools, where, reg, schema_hash, budget, spent, replace=replace,
                             start_run=start_run)
        rc = rc or r
        if r == 6:
            break
    return rc


def step_tally():
    """Counts from scores.json (the per-run judgements, made by reading), two views:
    the scheduled runs as run, and with every invalid or unanswered scheduled run
    replaced by its recorded replacement. Every quote is checked against its trace."""
    S = json.load(open(os.path.join(OUT, "scores.json")))
    runs = S["runs"]
    for r in runs:
        text = json.load(open(os.path.join(OUT, r["file"])))["result"]["final_text"]
        if not r["deciding"].startswith("(") and r["deciding"] not in text:
            raise SystemExit(f"{r['file']}: the deciding sentence is not in the trace")
    crit = ("answered", "ceiling_seen", "C1", "C1b", "C2", "C3")

    def view(replaced):
        out = {}
        for r in runs:
            if r.get("replacement"):
                continue
            use = r
            if replaced:
                rep_file = r["file"].replace(".json", "_replacement.json")
                use = next((x for x in runs if x["file"] == rep_file), r)
            key = f"{r['model']} {r['condition']}"
            cell = out.setdefault(key, {"runs": 0, **{c: 0 for c in crit}, "files": []})
            cell["runs"] += 1
            cell["files"].append(use["file"])
            for c in crit:
                cell[c] += bool(use[c])
        for cond in ("WITHOUT", "WITH"):
            tot = {"runs": 0, **{c: 0 for c in crit}}
            for k, v in out.items():
                if k.endswith(" " + cond):
                    for c in tot:
                        tot[c] += v[c]
            out[f"local models {cond}"] = tot
        return out
    seen = [r for r in runs if r["condition"] == "WITH" and r["ceiling_seen"]]
    rec = {"scheduled_runs": view(False), "with_replacements": view(True),
           "WITH_runs_that_saw_the_ceiling": {m: {"saw": sum(1 for r in seen if r["model"] == m),
                                                  "C2": sum(1 for r in seen if r["model"] == m and r["C2"]),
                                                  "C3": sum(1 for r in seen if r["model"] == m and r["C3"])}
                                              for m in LOCAL},
           "phase10_quoted": PHASE10_QUOTED, "api_models": S.get("api_models")}
    write(os.path.join(OUT, "tally.json"), rec)
    for name in ("scheduled_runs", "with_replacements"):
        print(name)
        for k, v in rec[name].items():
            print(f"  {k:28s} " + " ".join(f"{c}={v[c]}/{v['runs']}" for c in crit))
    print("WITH runs that saw the ceiling:", rec["WITH_runs_that_saw_the_ceiling"])
    return 0


def step_keycheck():
    """ONE count_tokens call on case (e)'s first request -- the system prompt, the
    tool schemas and the prompt exactly as a run sends them. Records the outcome and
    the token count, never the key (not its text, length or hash)."""
    tools, where, reg, schema_hash = setup()
    rec = {"what": "Anthropic API key check before the API arm: one count_tokens call on case (e)'s first request",
           "date": TODAY, "model": API[0], "schema_sha256": schema_hash,
           "system_prompt_sha256": sha(chatmod.SYSTEM_PROMPT), "prompt_sha256": sha(CASE_E),
           "key_file_mode": oct(os.stat(os.path.join(REPO, ".api", "claude_api_key")).st_mode & 0o777)}
    try:
        import anthropic
        c = anthropic.Anthropic(api_key=ui.API_KEY, max_retries=0, timeout=60)
        n = c.messages.count_tokens(model=API[0], system=chatmod.SYSTEM_PROMPT,
                                    tools=chatmod.to_anthropic_tools(tools),
                                    messages=[{"role": "user", "content": CASE_E}]).input_tokens
        rec.update({"status": "ok", "input_tokens": n})
    except Exception as e:
        rec.update({"status": "failed", "error_class": type(e).__name__,
                    "http_status": getattr(e, "status_code", None)})
    write(os.path.join(OUT, f"key_check_{TODAY}.json"), rec)
    print({k: rec[k] for k in rec if k in ("status", "input_tokens", "error_class", "http_status")})
    return 0 if rec["status"] == "ok" else 1


def step_spend():
    """The API arm's spend: every API run's usage and cost, summed, beside the
    projection and the cap."""
    rows = []
    for model in API:
        d = os.path.join(OUT, slug(model))
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            if "__run" not in f or not f.endswith(".json"):
                continue
            r = json.load(open(os.path.join(d, f)))["result"]
            rows.append({"model": model, "file": f"{slug(model)}/{f}", "usage": r.get("api_usage"),
                         "cost_usd": r.get("api_cost_usd_live_prices"), "chat_api_cost_usd": r.get("api_cost_usd"),
                         "api_error": r.get("api_error"), "refusals": len(r.get("refusals") or [])})
    proj = json.load(open(os.path.join(OUT, "projection.json")))
    per_model = {}
    for m in API:
        rr = [r for r in rows if r["model"] == m]
        per_model[m] = {"runs": len(rr), "cost_usd": round(sum(r["cost_usd"] or 0 for r in rr), 4),
                        "input_tokens": sum((r["usage"] or {}).get("input_tokens", 0) for r in rr),
                        "cache_write_tokens": sum((r["usage"] or {}).get("cache_creation_input_tokens", 0) for r in rr),
                        "cache_read_tokens": sum((r["usage"] or {}).get("cache_read_input_tokens", 0) for r in rr),
                        "output_tokens": sum((r["usage"] or {}).get("output_tokens", 0) for r in rr)}
    total = round(sum(v["cost_usd"] for v in per_model.values()), 4)
    rec = {"what": "Phase 10 rerun, API arm: spend from each run's usage fields at the list prices in PRICES_LIVE",
           "definition": "cost = input x in + cache_creation x cache_write_5m + cache_read x cache_read + output x out, "
                         "per million tokens, per run, summed",
           "prices": PRICES_LIVE, "per_model": per_model, "total_usd": total,
           "projected_usd": proj["projected_total_usd_20_runs"], "cap_usd": proj["budget_cap_usd"], "runs": rows}
    write(os.path.join(OUT, f"spend_{TODAY}.json"), rec)
    print(json.dumps({k: rec[k] for k in ("per_model", "total_usd", "projected_usd", "cap_usd")}, indent=1))
    return 0


def step_extract():
    for model in LOCAL + API:
        d = os.path.join(OUT, slug(model))
        if not os.path.isdir(d):
            continue
        for f in sorted(os.listdir(d)):
            if "__run" not in f or not f.endswith(".json"):
                continue
            rec = json.load(open(os.path.join(d, f)))
            r = rec["result"]
            print(f"\n===== {model} {rec['condition']} run {rec['run']} | answered={not r['ended_without_answer']} "
                  f"| tool calls {r['n_tool_calls']} | ceiling keys seen {len(rec['ceiling_keys_seen_by_model'])}")
            print(r["final_text"])
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["prove", "project", "run", "extract", "tally", "keycheck", "spend"])
    ap.add_argument("models", nargs="*")
    ap.add_argument("--runs", type=int, default=5)
    ap.add_argument("--budget-usd", type=float, default=None)
    ap.add_argument("--replace", help="COND:K -- one run replacing <COND>__run<K>.json (with --reason)")
    ap.add_argument("--reason")
    ap.add_argument("--start-run", type=int, default=1,
                    help="an extension: number the new runs from this, in their own seeded order")
    a = ap.parse_args()
    if a.step == "prove":
        sys.exit(step_prove())
    if a.step == "project":
        sys.exit(step_project())
    if a.step == "extract":
        sys.exit(step_extract())
    if a.step == "tally":
        sys.exit(step_tally())
    if a.step == "keycheck":
        sys.exit(step_keycheck())
    if a.step == "spend":
        sys.exit(step_spend())
    bad = [m for m in a.models if m not in LOCAL + API]
    if bad or not a.models:
        raise SystemExit(f"models must be from {LOCAL + API}")
    budget = a.budget_usd
    if budget is None:
        p = os.path.join(OUT, "projection.json")
        budget = json.load(open(p))["budget_cap_usd"] if os.path.exists(p) else 5.0
    replace = None
    if a.replace:
        cond_r, k_r = a.replace.split(":")
        if cond_r not in ("WITHOUT", "WITH") or not a.reason or len(a.models) != 1:
            raise SystemExit("--replace COND:K needs one model and a --reason")
        replace = (cond_r, int(k_r), a.reason)
    sys.exit(step_run(a.models, a.runs, budget, replace=replace, start_run=a.start_run))
