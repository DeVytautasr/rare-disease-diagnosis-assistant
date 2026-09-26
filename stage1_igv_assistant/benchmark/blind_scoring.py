"""Phase 12 Task 5: blind scoring of every Phase 10 ceiling run.

    python -m stage1_igv_assistant.benchmark.blind_scoring build
    python -m stage1_igv_assistant.benchmark.blind_scoring check SCORES.json
    python -m stage1_igv_assistant.benchmark.blind_scoring unblind SCORES.json

THE ITEMS. Every run file of case (e) under runs/phase10_rerun_2026-09-25/: the 23
local runs scored unblinded on 2026-09-25 (runs 1-5 of qwen2.5:7b and qwen3.5:4b and
the three replacement runs), the 40 runs of the registered extension (runs 6-15) and
the 20 API runs (claude-sonnet-5, claude-opus-5). Originals that were replaced are
scored too; a cell uses the replacement (the registered rule).

BLINDING. `build` writes the packet outside the repository: for each item a random
ID (secrets.token_hex), the prompt and the final answer, nothing else -- no model,
condition, file name, tool call or tool return -- in an order drawn by
secrets.SystemRandom. The key (ID -> run file) goes to the session scratch area with
a random salt, and the repository receives only the commitment: sha256 of the salt,
a newline and the key serialised with sorted keys, plus sha256 of the packet. The key
and salt are committed after scoring, so anyone can check the commitment.

THE SCORER is a fresh verifier subagent that reads only the packet, which carries the
four criteria verbatim from scores.json (scores.json itself also holds the unblinded
scores, so the scorer is told not to open it). For each item it records C1, C1b, C2
and C3 and quotes, verbatim from the answer, the sentence that decides the verdict
and, where C2 or C3 holds, the sentence that decides it.

CHECKS. `check` requires every quote to occur verbatim in the item's final answer in
the packet and in the run file the key names (the trace); a quote that fails is
reported, never repaired silently.

UNBLIND. `unblind` joins scores to runs and writes blind_scores.json: the per-item
scores; the cells (original runs 1-5, runs 1-15, API runs 1-5), each with runs,
answered, runs in which the ceiling fields reached the model, and the four criteria;
the agreement with the unblinded scores of 2026-09-25 on the 23 local runs, with
every disagreement quoted from both; and the registered test -- Fisher's exact test,
two-sided, WITH vs WITHOUT, per local model, on the 15-run cells, for C3 and C2. The
earlier scores are left as they are.
"""
import hashlib
import json
import math
import os
import secrets
import sys

from stage1_igv_assistant.benchmark import phase10_rerun as p10

RUNS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs", "phase10_rerun_2026-09-25")
BLIND = os.path.join(RUNS, "blind")
PACKET_DIR = os.path.expanduser("~/public_data/sim/logs/phase12_2026-09-26/blind")
KEY_PATH = os.environ.get("BLIND_KEY_PATH", "")
MODELS = {"qwen2.5-7b": "qwen2.5:7b", "qwen3.5-4b": "qwen3.5:4b", "claude-sonnet-5": "claude-sonnet-5",
          "claude-opus-5": "claude-opus-5"}
LOCAL = ("qwen2.5-7b", "qwen3.5-4b")
CRITERIA = ("C1", "C1b", "C2", "C3")


def sha(b):
    return hashlib.sha256(b).hexdigest()


def run_files():
    out = []
    for d in MODELS:
        for f in sorted(os.listdir(os.path.join(RUNS, d))):
            if "__run" in f and f.endswith(".json"):
                out.append(f"{d}/{f}")
    return out


def parse(rel):
    d, f = rel.split("/")
    cond, rest = f[:-5].split("__run")
    num, _, repl = rest.partition("_")
    return {"model_dir": d, "model": MODELS[d], "condition": cond, "run": int(num), "replacement": repl == "replacement"}


def packet_text(items, criteria):
    lines = ["# Blind scoring packet: Phase 10 case (e)", "",
             "Score every item below on the four criteria, reading the final answer only.",
             "Criteria, verbatim from the scoring record:", ""]
    lines += [f"- **{k}**: {v}" for k, v in criteria.items()]
    lines += ["", "---", ""]
    for it in items:
        lines += [f"## Item {it['id']}", "", "### Prompt", "", it["prompt"], "", "### Final answer", "",
                  it["final_answer"] if it["final_answer"].strip() else "(empty: the run ended without an answer)",
                  "", "---", ""]
    return "\n".join(lines)


def step_build():
    if not KEY_PATH:
        sys.exit("STOPPED: set BLIND_KEY_PATH to a path outside the repository for the key")
    if os.path.exists(os.path.join(BLIND, "commitment.json")):
        sys.exit("STOPPED: a commitment exists; the packet is built once")
    criteria = json.load(open(os.path.join(RUNS, "scores.json")))["criteria"]
    rng = secrets.SystemRandom()
    files = run_files()
    rng.shuffle(files)
    ids, items, key = set(), [], {}
    for rel in files:
        while True:
            i = "B" + secrets.token_hex(4)
            if i not in ids:
                ids.add(i)
                break
        run = json.load(open(os.path.join(RUNS, rel)))
        items.append({"id": i, "prompt": run["prompt"], "final_answer": run["result"].get("final_text") or ""})
        key[i] = rel
    os.makedirs(PACKET_DIR, exist_ok=True)
    packet = {"criteria": criteria, "items": items}
    pj = json.dumps(packet, indent=1, ensure_ascii=False).encode()
    with open(os.path.join(PACKET_DIR, "packet.json"), "wb") as f:
        f.write(pj)
    with open(os.path.join(PACKET_DIR, "packet.md"), "w") as f:
        f.write(packet_text(items, criteria))
    salt = secrets.token_hex(16)
    kj = json.dumps(key, sort_keys=True).encode()
    with open(KEY_PATH, "w") as f:
        json.dump({"salt": salt, "key": key}, f, indent=1, sort_keys=True)
    os.chmod(KEY_PATH, 0o600)
    os.makedirs(BLIND, exist_ok=True)
    rec = {"what": "Commitment to the blind-scoring key, written before any item was scored (Phase 12 Task 5)",
           "items": len(items), "by_model": {d: sum(1 for r in key.values() if r.startswith(d + "/")) for d in MODELS},
           "commitment_sha256": sha(salt.encode() + b"\n" + kj),
           "commitment_method": "sha256(salt + b'\\n' + json.dumps(key, sort_keys=True).encode()); key maps item ID "
                                "to run file; salt and key are committed after scoring",
           "packet_json_sha256": sha(pj),
           "packet": "prompt and final answer per item, random IDs, order from secrets.SystemRandom; committed after "
                     "scoring beside the key",
           "code": "stage1_igv_assistant/benchmark/blind_scoring.py"}
    p10.write(os.path.join(BLIND, "commitment.json"), rec)
    print(json.dumps({k: rec[k] for k in ("items", "by_model", "commitment_sha256", "packet_json_sha256")}, indent=1))
    return 0


def load_key():
    k = json.load(open(KEY_PATH))
    com = json.load(open(os.path.join(BLIND, "commitment.json")))
    got = sha(k["salt"].encode() + b"\n" + json.dumps(k["key"], sort_keys=True).encode())
    if got != com["commitment_sha256"]:
        sys.exit("STOPPED: the key does not match the commitment")
    return k


def quotes(s):
    return [(f, s.get(f)) for f in ("deciding", "c2_c3_sentence") if s.get(f)]


def step_check(scores_path):
    key = load_key()["key"]
    packet = json.load(open(os.path.join(PACKET_DIR, "packet.json")))
    answers = {it["id"]: it["final_answer"] for it in packet["items"]}
    scores = json.load(open(scores_path))["items"]
    problems = []
    if set(scores) != set(answers):
        problems.append(f"scored IDs differ from the packet: missing {sorted(set(answers) - set(scores))}, "
                        f"extra {sorted(set(scores) - set(answers))}")
    for i, s in scores.items():
        for c in CRITERIA:
            if not isinstance(s.get(c), bool):
                problems.append(f"{i}: {c} is not true/false")
        if s.get("C3") and not s.get("C2"):
            problems.append(f"{i}: C3 without C2 (C3 is defined as C2 and both reasons)")
        trace = json.load(open(os.path.join(RUNS, key[i])))["result"].get("final_text") or ""
        for field, q in quotes(s):
            if q not in answers.get(i, ""):
                problems.append(f"{i}: {field} is not verbatim in the packet answer")
            if q not in trace:
                problems.append(f"{i}: {field} is not verbatim in the trace")
        if (s.get("C2") or s.get("C3")) and not s.get("c2_c3_sentence"):
            problems.append(f"{i}: C2/C3 true with no quoted sentence")
    print(f"{len(scores)} items; {sum(1 for i in scores for _ in quotes(scores[i]))} quotes checked; "
          f"{len(problems)} problems")
    for p in problems:
        print("  ", p)
    return 1 if problems else 0


def fisher_two_sided(a, b, c, d):
    """Table [[a, b], [c, d]]: rows WITH, WITHOUT; columns criterion met, not met."""
    r1, r2, c1 = a + b, c + d, a + c
    n = r1 + r2

    def p(x):
        return math.comb(r1, x) * math.comb(r2, c1 - x) / math.comb(n, c1)
    p0 = p(a)
    return min(1.0, sum(p(x) for x in range(max(0, c1 - r2), min(r1, c1) + 1) if p(x) <= p0 * (1 + 1e-7)))


def cell_runs(items, model_dir, cond, max_run):
    rows = {}
    for it in items:
        m = it["meta"]
        if m["model_dir"] == model_dir and m["condition"] == cond and m["run"] <= max_run:
            if m["run"] not in rows or m["replacement"]:
                rows[m["run"]] = it
    return [rows[k] for k in sorted(rows)]


def tally(rows):
    return {"runs": len(rows), "answered": sum(1 for r in rows if r["answered"]),
            "ceiling_fields_reached_model": sum(1 for r in rows if r["ceiling_seen"]),
            **{c: sum(1 for r in rows if r["blind"][c]) for c in CRITERIA}}


def step_unblind(scores_path):
    if step_check(scores_path) != 0:
        sys.exit("STOPPED: the quote check did not pass; nothing is unblinded")
    k = load_key()
    key = k["key"]
    scores = json.load(open(scores_path))
    earlier = {r["file"]: r for r in json.load(open(os.path.join(RUNS, "scores.json")))["runs"]}
    items = []
    for i, rel in sorted(key.items(), key=lambda kv: kv[1]):
        run = json.load(open(os.path.join(RUNS, rel)))
        s = scores["items"][i]
        items.append({"id": i, "file": rel, "meta": parse(rel),
                      "answered": bool((run["result"].get("final_text") or "").strip()),
                      "ceiling_seen": bool(run.get("ceiling_keys_seen_by_model")),
                      "blind": {c: s[c] for c in CRITERIA},
                      "deciding": s.get("deciding"), "c2_c3_sentence": s.get("c2_c3_sentence"),
                      "scorer_note": s.get("note")})
    cells = {"original_runs_1_5": {}, "runs_1_15": {}, "api_runs_1_5": {}}
    for d in LOCAL:
        for cond in ("WITHOUT", "WITH"):
            five = cell_runs(items, d, cond, 5)
            cells["original_runs_1_5"][f"{MODELS[d]} {cond}"] = {
                "blind": tally(five),
                "unblinded_2026_09_25": {c: sum(1 for r in five if earlier[r["file"]][c]) for c in CRITERIA}}
            cells["runs_1_15"][f"{MODELS[d]} {cond}"] = {"blind": tally(cell_runs(items, d, cond, 15))}
    for d in ("claude-sonnet-5", "claude-opus-5"):
        for cond in ("WITHOUT", "WITH"):
            cells["api_runs_1_5"][f"{d} {cond}"] = {"blind": tally(cell_runs(items, d, cond, 5))}
    tests = {}
    for d in LOCAL:
        w, wo = cell_runs(items, d, "WITH", 15), cell_runs(items, d, "WITHOUT", 15)
        for c in ("C3", "C2"):
            a, b = sum(r["blind"][c] for r in w), sum(not r["blind"][c] for r in w)
            cc, dd = sum(r["blind"][c] for r in wo), sum(not r["blind"][c] for r in wo)
            tests[f"{MODELS[d]} {c}"] = {"table [[WITH met, WITH not], [WITHOUT met, WITHOUT not]]": [[a, b], [cc, dd]],
                                         "p_two_sided": round(fisher_two_sided(a, b, cc, dd), 6)}
    agreement, disagreements = {c: [0, 0] for c in CRITERIA}, []
    for it in items:
        e = earlier.get(it["file"])
        if not e:
            continue
        for c in CRITERIA:
            agreement[c][1] += 1
            if e[c] == it["blind"][c]:
                agreement[c][0] += 1
            else:
                disagreements.append({"file": it["file"], "criterion": c, "unblinded": e[c], "blind": it["blind"][c],
                                      "unblinded_deciding": e.get("deciding"), "unblinded_note": e.get("note"),
                                      "blind_deciding": it["deciding"], "blind_c2_c3_sentence": it["c2_c3_sentence"],
                                      "blind_note": it["scorer_note"]})
    out = {"what": "Blind scores of every Phase 10 case (e) run, unblinded after scoring (Phase 12 Task 5)",
           "scorer": scores.get("scorer"), "criteria": json.load(open(os.path.join(RUNS, "scores.json")))["criteria"],
           "method": __doc__.split("THE SCORER")[0].strip(),
           "items": items, "cells": cells,
           "fisher_exact_two_sided_15_run_cells": tests,
           "agreement_with_unblinded_2026_09_25": {c: f"{a} of {n}" for c, (a, n) in agreement.items()},
           "disagreements": disagreements,
           "registration": "registration_local_extension_2026-09-26.json",
           "earlier_scores": "scores.json, unchanged"}
    p10.write(os.path.join(BLIND, "blind_scores.json"), out)
    with open(os.path.join(BLIND, "key.json"), "w") as f:
        json.dump(k, f, indent=1, sort_keys=True)
    with open(os.path.join(PACKET_DIR, "packet.json"), "rb") as f:
        pj = f.read()
    if sha(pj) != json.load(open(os.path.join(BLIND, "commitment.json")))["packet_json_sha256"]:
        sys.exit("STOPPED: the packet does not match its committed hash")
    with open(os.path.join(BLIND, "packet.json"), "wb") as f:
        f.write(pj)
    print(json.dumps({"cells": cells, "tests": tests, "agreement": out["agreement_with_unblinded_2026_09_25"],
                      "disagreements": len(disagreements)}, indent=1))
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if a == ["build"]:
        sys.exit(step_build())
    if len(a) == 2 and a[0] in ("check", "unblind"):
        sys.exit({"check": step_check, "unblind": step_unblind}[a[0]](a[1]))
    sys.exit("usage: blind_scoring.py build | check SCORES.json | unblind SCORES.json")
