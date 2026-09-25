#!/usr/bin/env python3
"""Gather the synthetic-control records into the repository.

    write_run_record.py [--dest DIR]

Copies what the run produced under ~/public_data/sim (measurements, the scan,
the selection, per-implant preparation summaries, the gate, build and delly
records) into DIR (default stage1_igv_assistant/results/synthetic_control_2026-09)
and writes run_record.json: commands, tool versions, seeds, class thresholds,
the gate result, delly exit codes and peak RSS where delly ran, and the run's
status. It computes no new figure; it reports which steps exist and which do not.
Large data (BAM, BCF, FASTQ, ART SAM) stays under ~/public_data/sim.
"""
import argparse
import json
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from synth_common import BAMS, SIM, WORK, tool_versions, write_json  # noqa: E402

REPO = os.path.dirname(os.path.dirname(HERE))
DEST = os.path.join(REPO, "stage1_igv_assistant", "results", "synthetic_control_2026-09")


def load(path):
    return json.load(open(path)) if os.path.exists(path) else None


def delly_log(label):
    log = os.path.join(SIM, "logs", f"delly_{label}.log")
    rc = os.path.join(SIM, "logs", f"delly_{label}.rc")
    if not os.path.exists(log):
        return None
    text = open(log).read()
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", text)
    samples = re.search(r"samples in BCF: (\d+) \(([^)]*)\)", text)
    records = re.search(r"records in BCF: (\d+)", text)
    return {"exit_code": int(open(rc).read()) if os.path.exists(rc) else None,
            "peak_rss_kbytes": int(rss.group(1)) if rss else None,
            "wall_clock": wall.group(1) if wall else None,
            "samples_in_bcf": int(samples.group(1)) if samples else None,
            "sample_names": samples.group(2).split() if samples else None,
            "records_in_bcf": int(records.group(1)) if records else None}


def main(dest):
    os.makedirs(dest, exist_ok=True)
    copies = {"background_measurements.json": os.path.join(SIM, "measure", "background.json"),
              "art_profiles.json": os.path.join(SIM, "measure", "profiles.json"),
              "scan_candidates.json": os.path.join(SIM, "scan", "candidates.json"),
              "selection.json": os.path.join(SIM, "selection.json"),
              "pooled_alignment.json": os.path.join(WORK, "pool", "align.json")}
    for name, src in copies.items():
        if os.path.exists(src):
            shutil.copyfile(src, os.path.join(dest, name))
    sel = load(os.path.join(SIM, "selection.json"))
    implants = []
    for imp in (sel or {}).get("implants", []):
        iid = imp["id"]
        prep = load(os.path.join(WORK, iid, "prepare.json"))
        gate = load(os.path.join(WORK, iid, "gate.json"))
        build = load(os.path.join(WORK, iid, "build.json"))
        if gate:
            shutil.copyfile(os.path.join(WORK, iid, "gate.json"), os.path.join(dest, f"{iid}_gate.json"))
        implants.append({
            "id": iid, "class": imp["class"], "chr20": imp["chr20"], "chr21": imp["chr21"],
            "source": imp["source"],
            "prepared": prep is not None,
            "compensation": ({"removal": prep["removal"], "fragments_added_total": prep["fragments_added_total"],
                              "junctions": prep["junctions"]} if prep else None),
            "bam_built": build is not None and os.path.exists(os.path.join(BAMS, f"{iid}.bam")),
            "build": build, "gate_verdict": gate["verdict"] if gate else None,
            "delly": delly_log(iid)})
    gates = {i["id"]: i["gate_verdict"] for i in implants if i["gate_verdict"]}
    built = [i["id"] for i in implants if i["bam_built"]]
    delly_ran = {i["id"]: i["delly"] for i in implants if i["delly"]}
    bg_delly = delly_log("background")
    if gates.get("IMP01") == "FAIL":
        status = ("STOPPED at the pre-registered IMP01 gate: at least one criterion failed. Nothing "
                  "downstream of the gate was run -- no further implant BAMs, no delly on the implants "
                  "or on the background, no demo bundle.")
    elif gates.get("IMP01") == "PASS" and len(delly_ran) == 12 and bg_delly:
        status = "COMPLETE through delly (12 implants and the background); the evidence chain has not been run."
    else:
        status = "IN PROGRESS"
    rec = {"experiment": "synthetic positive control, rebuilt 2026-09-25 -- a NEW experiment "
                         "(new draws, ALT-aware index, partly new breakpoints), not a reproduction of Phase 6",
           "status": status,
           "design": "heterozygous balanced reciprocal t(20;21), one implant per BAM, in the NA12878 "
                     "chr20+chr21 slice; see scripts/synthetic/make_implants.py for every rule",
           "code": {"measure": "scripts/synthetic/measure_background.py",
                    "scan": "scripts/synthetic/scan_breakpoints.py",
                    "select": "scripts/synthetic/select_breakpoints.py",
                    "generate": "scripts/synthetic/make_implants.py",
                    "gate": "scripts/synthetic/gate_imp01.py",
                    "delly": "scripts/synthetic/run_delly_synth.sh"},
           "versions": tool_versions(),
           "seeds": {"master_seed": (sel or {}).get("master_seed"),
                     "selection_draws": {i["id"]: i["source"] for i in (sel or {}).get("implants", [])},
                     "per_implant": {i["id"]: ({"remove_A": i["compensation"]["removal"]["A"]["remove_seed"],
                                                "remove_B": i["compensation"]["removal"]["B"]["remove_seed"],
                                                "art_J20": i["compensation"]["junctions"]["J20"]["art_seed"],
                                                "pick_J20": i["compensation"]["junctions"]["J20"]["pick_seed"],
                                                "art_J21": i["compensation"]["junctions"]["J21"]["art_seed"],
                                                "pick_J21": i["compensation"]["junctions"]["J21"]["pick_seed"]}
                                               if i["compensation"] else None) for i in implants},
                     "art_profile_seeds": (load(os.path.join(SIM, "measure", "profiles.json")) or {}).get("seeds")},
           "class_thresholds": (sel or {}).get("thresholds"),
           "recorded_coordinate_checks": (sel or {}).get("recorded_coordinate_checks"),
           "gate": {i: load(os.path.join(WORK, i, "gate.json")) and
                    {k: v for k, v in load(os.path.join(WORK, i, "gate.json")).items() if k != "per_read"}
                    for i in gates},
           "bams_built": built, "delly": delly_ran, "delly_background": bg_delly,
           "implants": implants}
    write_json(os.path.join(dest, "run_record.json"), rec)
    print(f"records -> {os.path.relpath(dest, REPO)}: {sorted(os.listdir(dest))}")
    print(f"status: {status}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", default=DEST)
    main(ap.parse_args().dest)
