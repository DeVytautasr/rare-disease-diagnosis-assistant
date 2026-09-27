#!/usr/bin/env python3
"""Run every suite in stage1_igv_assistant/tests/ once, with the test census's
instrumentation (scripts/test_census.py: an assertion is one evaluation of a
suite's check() or assert statement), and write OUT (a JSON record).

    run_suites.py OUT_JSON

Recorded per suite: exit status (0 passed, 1 failed, 2 incomplete), assertions held
and failed, and the lines beginning NOT RUN or SKIPPED. The environment is the
caller's (for the IGV checks: Java on PATH, DISPLAY set).
"""
import json
import os
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_census as tc  # noqa: E402


def main(out_path):
    logdir = os.path.join(os.path.expanduser("~/public_data/sim/logs"), "suites_" + time.strftime("%Y-%m-%d_%H%M"))
    os.makedirs(logdir, exist_ok=True)
    rec = {"what": "every suite once", "definitions": " ".join(__doc__.split()), "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "java_on_path": subprocess.run(["which", "java"], capture_output=True).returncode == 0,
           "display": os.environ.get("DISPLAY"), "suites": {}}
    for path in tc.SUITES:
        name = os.path.basename(path)
        log, cnt = os.path.join(logdir, name + ".log"), os.path.join(logdir, name + ".counts.json")
        with open(log, "w") as fh:
            p = subprocess.run([tc.PY, os.path.join(REPO, "scripts", "test_census.py"), "suite", path, cnt],
                               cwd=REPO, stdout=fh, stderr=subprocess.STDOUT, timeout=tc.TIME_LIMIT_S)
        c = json.load(open(cnt)) if os.path.exists(cnt) else {}
        rec["suites"][name] = {"exit": p.returncode,
                               "assertions_held": c.get("assert_held", 0) + c.get("check_held", 0),
                               "assertions_failed": c.get("assert_failed", 0) + c.get("check_failed", 0),
                               "not_run": tc.not_run_lines(open(log, errors="replace").read())}
        print(name, rec["suites"][name], flush=True)
    s = rec["suites"].values()
    rec["summary"] = {"suites": len(rec["suites"]), "passed": sum(x["exit"] == 0 for x in s),
                      "failed": sum(x["exit"] == 1 for x in s), "incomplete": sum(x["exit"] == 2 for x in s),
                      "assertions_held": sum(x["assertions_held"] for x in s),
                      "assertions_failed": sum(x["assertions_failed"] for x in s),
                      "not_run_lines": sum(len(x["not_run"]) for x in s)}
    rec["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    if os.path.exists(out_path):
        sys.exit("the record exists")
    json.dump(rec, open(out_path, "w"), indent=1)
    print(json.dumps(rec["summary"]))
    return 0 if rec["summary"]["failed"] == 0 and rec["summary"]["incomplete"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
