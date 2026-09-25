#!/usr/bin/env python3
"""Gather the synthetic-control records into the repository.

    write_run_record.py [--dest DIR]

Copies what the run produced under ~/public_data/sim (measurements, the scan,
the selection, per-implant preparation summaries, the gate and any revised gate,
build and delly records) into DIR (default
stage1_igv_assistant/results/synthetic_control_2026-09) and writes the record
under --name (a committed record is never overwritten; a later stage gets a new,
dated name): commands, tool versions, seeds, class thresholds,
the gate result, delly exit codes and peak RSS where delly ran, and the run's
status. It computes no new figure; it reports which steps exist and which do not.
Large data (BAM, BCF, FASTQ, ART SAM) stays under ~/public_data/sim.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from synth_common import BAMS, SIM, WORK, tool_versions, write_json  # noqa: E402

REPO = os.path.dirname(os.path.dirname(HERE))
DEST = os.path.join(REPO, "stage1_igv_assistant", "results", "synthetic_control_2026-09")


def tracked(path):
    return subprocess.run(["git", "-C", REPO, "ls-files", "--error-unmatch", path],
                          capture_output=True).returncode == 0


def copy_unless_committed_differs(src, dst):
    """A committed record is never changed: a copy that would alter a tracked
    file is skipped and reported instead."""
    if os.path.exists(dst) and tracked(dst):
        if open(src, "rb").read() != open(dst, "rb").read():
            print(f"NOT overwritten (committed, and the source now differs): {os.path.relpath(dst, REPO)}")
        return
    shutil.copyfile(src, dst)


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


LIMITATIONS = [
    "The implants are aligned as their own NA12878 background was (bwa 0.7.15 mem -Y -K 100000000, "
    "no -M): supplementary alignments soft-clipped and flagged 0x800. The two patient BAMs were aligned "
    "with bwa mem -M and without -Y (every bwa @PG record), so their split alignments are represented "
    "differently: hard-clipped, and marked secondary by -M (although flagstat still counts 3,694,802 and "
    "3,726,497 supplementary records in them, which -M alone would not produce). Consistency with the "
    "background was chosen over matching the clinical pipeline (user decision, 2026-09-25); the soft-clip "
    "and split-read layers may therefore behave differently on patient data than on these implants.",
    "Simulated reads come from the reference sequence and carry none of NA12878's own variants on the "
    "rearranged haplotype.",
]


def main(dest, name, revised_name):
    os.makedirs(dest, exist_ok=True)
    copies = {"background_measurements.json": os.path.join(SIM, "measure", "background.json"),
              "art_profiles.json": os.path.join(SIM, "measure", "profiles.json"),
              "scan_candidates.json": os.path.join(SIM, "scan", "candidates.json"),
              "selection.json": os.path.join(SIM, "selection.json"),
              "pooled_alignment.json": os.path.join(WORK, "pool", "align.json")}
    for fname, src in copies.items():
        if os.path.exists(src):
            copy_unless_committed_differs(src, os.path.join(dest, fname))
    sel = load(os.path.join(SIM, "selection.json"))
    implants = []
    for imp in (sel or {}).get("implants", []):
        iid = imp["id"]
        prep = load(os.path.join(WORK, iid, "prepare.json"))
        gate = load(os.path.join(WORK, iid, "gate.json"))
        build = load(os.path.join(WORK, iid, "build.json"))
        revised = load(os.path.join(WORK, iid, "gate_revised.json"))
        if gate:
            copy_unless_committed_differs(os.path.join(WORK, iid, "gate.json"),
                                          os.path.join(dest, f"{iid}_gate.json"))
        if revised:
            copy_unless_committed_differs(os.path.join(WORK, iid, "gate_revised.json"),
                                          os.path.join(dest, f"{iid}_{revised_name}"))
        implants.append({
            "id": iid, "class": imp["class"], "chr20": imp["chr20"], "chr21": imp["chr21"],
            "source": imp["source"],
            "prepared": prep is not None,
            "compensation": ({"removal": prep["removal"], "fragments_added_total": prep["fragments_added_total"],
                              "junctions": prep["junctions"]} if prep else None),
            "bam_built": build is not None and os.path.exists(os.path.join(BAMS, f"{iid}.bam")),
            "build": build, "gate_verdict": gate["verdict"] if gate else None,
            "revised_gate_verdict": revised["verdict"] if revised else None,
            "delly": delly_log(iid)})
    gates = {i["id"]: i["gate_verdict"] for i in implants if i["gate_verdict"]}
    revised = {i["id"]: i["revised_gate_verdict"] for i in implants if i["revised_gate_verdict"]}
    built = [i["id"] for i in implants if i["bam_built"]]
    delly_ran = {i["id"]: i["delly"] for i in implants if i["delly"]}
    bg_delly = delly_log("background")
    if revised.get("IMP01") == "PASS" and len(delly_ran) == 12 and bg_delly:
        status = ("COMPLETE through delly (12 implants and the background) after the REVISED gate: "
                  "criterion (a) was replaced on 2026-09-25, after the pre-registered gate failed, by "
                  "(a1)/(a2) -- see the revised gate record; the evidence chain has not been run.")
    elif revised.get("IMP01") == "PASS":
        status = "IN PROGRESS after the REVISED gate passed (criterion (a) replaced after the pre-registered failure)."
    elif revised.get("IMP01") == "FAIL":
        status = "STOPPED at the REVISED gate: an SA-less truly spanning read is not explained by -T 30 / -k 19."
    elif gates.get("IMP01") == "FAIL":
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
           "revised_gate": {i: load(os.path.join(WORK, i, "gate_revised.json")) and
                            {k: v for k, v in load(os.path.join(WORK, i, "gate_revised.json")).items()
                             if k not in ("per_read",)} for i in revised},
           "limitations": LIMITATIONS,
           "bams_built": built, "delly": delly_ran, "delly_background": bg_delly,
           "implants": implants}
    write_json(os.path.join(dest, name), rec)
    print(f"records -> {os.path.relpath(dest, REPO)}: {sorted(os.listdir(dest))}")
    print(f"status: {status}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", default=DEST)
    ap.add_argument("--name", default="run_record.json",
                    help="a committed record is never overwritten: a later stage gets a new, dated name")
    ap.add_argument("--revised-name", default="gate_revised_2026-09-25.json")
    a = ap.parse_args()
    if tracked(os.path.join(a.dest, a.name)):
        raise SystemExit(f"{a.name} is committed; write a later stage under a new --name")
    main(a.dest, a.name, a.revised_name)
