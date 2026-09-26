#!/usr/bin/env python3
"""Phase 12 Task 6(d): run every test suite once per network condition, and count what held.

    test_census.py run CONDITION       CONDITION = network | no_network | annotation_only
    test_census.py suite PATH OUT      (internal) one suite, instrumented, run as __main__
    test_census.py record              the three conditions' results -> the committed record

SUITES. Every stage1_igv_assistant/tests/test_*.py, once per condition, each in a fresh
interpreter (.venv/bin/python) with the repository root as working directory and the
suite's path as argv[0], as `python stage1_igv_assistant/tests/test_x.py` runs it.
Output goes to logs outside the repository (LOGS below). Logs are not committed
because they can hold local paths; the record keeps counts and outcomes.

CONDITIONS. The suites have no switch for the network: no suite reads an environment
variable that turns the annotation service or the remote BAM on or off. The three
conditions are therefore imposed from outside, and the record says how:
  network          the machine as it is
  no_network       a new network namespace with only loopback up (unshare -rn, then
                   ip link set lo up): no DNS and no route out; servers on 127.0.0.1
                   still work
  annotation_only  an APPROXIMATION: http_proxy and https_proxy (both cases) point at
                   127.0.0.1:9, where nothing listens, with no_proxy =
                   rest.ensembl.org,localhost,127.0.0.1, so a client that honours the
                   proxy variables reaches Ensembl directly and nothing else. requests
                   (the annotation client) and libcurl inside htslib (the remote BAM)
                   both honour them; a client that ignored them would not be confined.
                   The reachability probes, run inside the condition, show what was
                   actually reachable.

FIGURES AND THEIR DEFINITIONS (written into the record beside each figure)
  assertion        one evaluation of a condition that can fail the suite: a call of the
                   suite's own module-level check() (held when its condition argument
                   is truthy), or one execution of an assert statement in the suite
                   file (held when its condition is truthy). Counted by instrumenting
                   the suite before it runs: its source is parsed and compiled here,
                   each assert's condition is passed through a counter that returns it
                   unchanged, a counter call is inserted as the first statement of
                   check(), and the result runs as __main__; nothing else changes.
                   Conditions a suite evaluates any other way are not counted: its own
                   counters (test_api_leak.py, test_server.py), or anything in a child
                   process (test_checks_can_fail.py). Each suite's static count of
                   assert statements and check() definitions is recorded beside it
  suite_outcome    the suite's exit status: 0 passed, 1 failed, 2 incomplete (a check
                   did not run), anything else an error; a suite killed at the time
                   limit is recorded as a timeout
  not_run          the lines a suite prints that begin with NOT RUN or SKIPPED,
                   counted, with the check each names (home directory shortened to ~)
  reachability     probes run inside each condition before the suites: DNS for
                   rest.ensembl.org; an HTTPS GET of Ensembl's /info/ping (status, or
                   exception type); the remote GIAB BAM's header through pysam (read, or
                   exception type), from a scratch working directory; and the local
                   Ollama endpoint 127.0.0.1:11434/api/tags (status, or exception type)
"""
import ast
import builtins
import glob
import json
import os
import re
import socket
import subprocess
import sys
import time
from collections import Counter

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(REPO, ".venv", "bin", "python")
LOGS = os.path.expanduser("~/public_data/sim/logs/phase12_2026-09-26/census")
SUITES = sorted(glob.glob(os.path.join(REPO, "stage1_igv_assistant", "tests", "test_*.py")))
CONDITIONS = ("network", "no_network", "annotation_only")
TIME_LIMIT_S = 3600
GIAB_BAM = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/"
            "NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam")
DEAD_PROXY = "http://127.0.0.1:9"
NO_PROXY = "rest.ensembl.org,localhost,127.0.0.1"
RECORD = os.path.join(REPO, "stage1_igv_assistant", "results", "test_census_2026-09-26.json")


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


def conditions_text():
    doc = __doc__.split("CONDITIONS.")[1].split("FIGURES AND THEIR DEFINITIONS")[0]
    return " ".join(doc.split())


# ── one suite, instrumented ──────────────────────────────────────────────────

class Instrument(ast.NodeTransformer):
    def __init__(self):
        self.asserts = 0
        self.checks = 0

    def visit_Assert(self, node):
        self.generic_visit(node)
        self.asserts += 1
        call = ast.Call(ast.Name("__census_assert__", ast.Load()), [node.test], [])
        node.test = ast.copy_location(call, node.test)
        return node


def instrument(tree):
    ins = Instrument()
    tree = ins.visit(tree)
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "check" and len(node.args.args) >= 2:
            ins.checks += 1
            hit = ast.Expr(ast.Call(ast.Name("__census_check__", ast.Load()),
                                    [ast.Name(node.args.args[1].arg, ast.Load())], []))
            doc = 1 if (node.body and isinstance(node.body[0], ast.Expr)
                        and isinstance(node.body[0].value, ast.Constant)
                        and isinstance(node.body[0].value.value, str)) else 0
            node.body.insert(doc, ast.copy_location(hit, node.body[min(doc, len(node.body) - 1)]))
    return ast.fix_missing_locations(tree), ins


def step_suite(path, out):
    counts = Counter()

    def truth(v):
        try:
            return bool(v)
        except Exception:
            return None

    def on_assert(v):
        t = truth(v)
        counts["assert_held" if t else ("assert_failed" if t is False else "assert_unevaluable")] += 1
        return v

    def on_check(v):
        t = truth(v)
        counts["check_held" if t else ("check_failed" if t is False else "check_unevaluable")] += 1

    tree, ins = instrument(ast.parse(open(path).read(), path))
    code = compile(tree, path, "exec")
    g = {"__name__": "__main__", "__file__": path, "__builtins__": builtins, "__package__": None,
         "__spec__": None, "__census_assert__": on_assert, "__census_check__": on_check}
    sys.argv = [path]
    sys.path[0] = os.path.dirname(path)
    status, ended = 0, "normal"
    try:
        exec(code, g)
    except SystemExit as e:
        status = e.code if isinstance(e.code, int) else (0 if e.code is None else 1)
    except AssertionError:
        status, ended = 1, "AssertionError"
    except BaseException as e:  # noqa: BLE001 -- recorded, not hidden
        status, ended = 1, type(e).__name__
    finally:
        with open(out, "w") as f:
            json.dump({"exit": status, "ended": ended, "static_asserts": ins.asserts,
                       "static_check_defs": ins.checks, **counts}, f)
    sys.stdout.flush()
    sys.exit(status if isinstance(status, int) and 0 <= status < 256 else 1)  # atexit handlers still run


# ── one condition ───────────────────────────────────────────────────────────

def probe():
    r = {}
    try:
        socket.getaddrinfo("rest.ensembl.org", 443)
        r["dns_rest_ensembl_org"] = "resolved"
    except Exception as e:  # noqa: BLE001
        r["dns_rest_ensembl_org"] = type(e).__name__
    try:
        import requests
        resp = requests.get("https://rest.ensembl.org/info/ping", headers={"Content-Type": "application/json"},
                            timeout=20)
        r["ensembl_ping"] = resp.status_code
    except Exception as e:  # noqa: BLE001
        r["ensembl_ping"] = type(e).__name__
    scratch = os.path.join(LOGS, "probe_cwd")
    os.makedirs(scratch, exist_ok=True)
    child = ("import pysam,sys\ntry:\n    f=pysam.AlignmentFile(sys.argv[1]); print('read', len(f.references))\n"
             "except Exception as e:\n    print(type(e).__name__, getattr(e,'errno',None))\n")
    try:
        p = subprocess.run([PY, "-c", child, GIAB_BAM], cwd=scratch, capture_output=True, text=True, timeout=180)
        r["giab_remote_bam_header"] = (p.stdout.strip().splitlines() or [f"exit {p.returncode}"])[-1]
    except subprocess.TimeoutExpired:
        r["giab_remote_bam_header"] = "timeout after 180 s"
    try:
        import requests
        resp = requests.get("http://127.0.0.1:11434/api/tags", timeout=5)
        r["local_ollama"] = resp.status_code
    except Exception as e:  # noqa: BLE001
        r["local_ollama"] = type(e).__name__
    return r


def not_run_lines(text):
    home = os.path.expanduser("~")
    out = []
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("NOT RUN") or s.startswith("SKIPPED"):
            out.append(re.sub(r"\s+", " ", s.replace(home, "~"))[:200])
    return out


def run_condition(condition):
    d = os.path.join(LOGS, condition)
    os.makedirs(d, exist_ok=True)
    rec = {"condition": condition, "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "proxy_env": {k: os.environ.get(k) for k in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY",
                                                       "no_proxy", "NO_PROXY") if os.environ.get(k)},
           "reachability": probe(), "suites": {}}
    print(json.dumps(rec["reachability"]), flush=True)
    for path in SUITES:
        name = os.path.basename(path)
        log, out = os.path.join(d, name + ".log"), os.path.join(d, name + ".counts.json")
        if os.path.exists(out):
            os.unlink(out)
        t0 = time.monotonic()
        with open(log, "w") as fh:
            try:
                p = subprocess.run([PY, os.path.abspath(__file__), "suite", path, out], cwd=REPO, stdout=fh,
                                   stderr=subprocess.STDOUT, timeout=TIME_LIMIT_S)
                rc = p.returncode
            except subprocess.TimeoutExpired:
                rc = "timeout"
        text = open(log, errors="replace").read()
        c = json.load(open(out)) if os.path.exists(out) else {}
        held = c.get("assert_held", 0) + c.get("check_held", 0)
        failed = c.get("assert_failed", 0) + c.get("check_failed", 0)
        rec["suites"][name] = {"exit": rc, "counts": c, "assertions_held": held, "assertions_failed": failed,
                               "not_run": not_run_lines(text), "wall_s": round(time.monotonic() - t0, 1)}
        print(f"{condition} {name}: exit={rc} held={held} failed={failed} "
              f"not_run={len(rec['suites'][name]['not_run'])} {rec['suites'][name]['wall_s']}s", flush=True)
    s = rec["suites"].values()
    rec["summary"] = {"suites": len(SUITES),
                      "passed": sum(1 for x in s if x["exit"] == 0),
                      "failed": sum(1 for x in s if x["exit"] == 1),
                      "incomplete": sum(1 for x in s if x["exit"] == 2),
                      "other": sum(1 for x in s if x["exit"] not in (0, 1, 2)),
                      "assertions_held": sum(x["assertions_held"] for x in s),
                      "assertions_failed": sum(x["assertions_failed"] for x in s),
                      "not_run_lines": sum(len(x["not_run"]) for x in s)}
    rec["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    with open(os.path.join(LOGS, f"{condition}.json"), "w") as f:
        json.dump(rec, f, indent=1)
    print(json.dumps(rec["summary"]), flush=True)
    return 0


def step_run(condition):
    if condition == "no_network":
        if os.environ.get("CENSUS_IN_NETNS") != "1":
            env = {**os.environ, "CENSUS_IN_NETNS": "1"}
            return subprocess.call(["unshare", "-rn", "sh", "-c", 'ip link set lo up && exec "$@"', "sh",
                                    PY, os.path.abspath(__file__), "run", "no_network"], env=env)
    elif condition == "annotation_only":
        if os.environ.get("https_proxy") != DEAD_PROXY:
            env = {**os.environ, "http_proxy": DEAD_PROXY, "https_proxy": DEAD_PROXY, "HTTP_PROXY": DEAD_PROXY,
                   "HTTPS_PROXY": DEAD_PROXY, "no_proxy": NO_PROXY, "NO_PROXY": NO_PROXY}
            return subprocess.call([PY, os.path.abspath(__file__), "run", "annotation_only"], env=env)
    elif condition != "network":
        print(f"unknown condition {condition}; use one of {CONDITIONS}", file=sys.stderr)
        return 2
    return run_condition(condition)


def step_record():
    out = {"what": "Every test suite run once per network condition, Phase 12 Task 6(d), 2026-09-26",
           "switches": "the suites have no environment switch for the network conditions; they were imposed "
                       "from outside as described under conditions",
           "conditions": conditions_text(),
           "definitions": definitions(),
           "thesis_values_compared": {"suites": "seventeen", "assertions no network": 306,
                                      "assertions with the annotation service": 308,
                                      "assertions with the remote BAM": 309},
           "code": "scripts/test_census.py"}
    for c in CONDITIONS:
        p = os.path.join(LOGS, f"{c}.json")
        if not os.path.exists(p):
            print(f"STOPPED: condition {c} has not been run", file=sys.stderr)
            return 2
        out[c] = json.load(open(p))
    if os.path.exists(RECORD):
        print("STOPPED: the record exists; a committed record is never overwritten", file=sys.stderr)
        return 2
    with open(RECORD, "w") as f:
        json.dump(out, f, indent=1)
    print("written: stage1_igv_assistant/results/test_census_2026-09-26.json")
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[:1] == ["suite"] and len(a) == 3:
        step_suite(a[1], a[2])
    elif a[:1] == ["run"] and len(a) == 2:
        sys.exit(step_run(a[1]))
    elif a[:1] == ["record"]:
        sys.exit(step_record())
    else:
        print("usage: test_census.py run CONDITION | suite PATH OUT | record", file=sys.stderr)
        sys.exit(2)
