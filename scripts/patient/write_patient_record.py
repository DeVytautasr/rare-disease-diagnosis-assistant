#!/usr/bin/env python3
"""Reproduce the Phase 3 patient funnel from the rerun BCFs, and write the run
record (Tasks 1 and 6 of the 2026-09-25 recovery prompt).

    write_patient_record.py funnel     both samples, every step against the Phase 3 targets
    write_patient_record.py record     the run record committed to results/

BLINDING: counts only. No coordinate, gene, candidate or interpretation of any
call is printed or written; candidate lists returned by the tools are counted
and discarded. Nothing here reads a sample column or the BCF header's sample
name. Display the output through scripts/redact.sh anyway.

THE FUNNEL is the tools' own chain, exactly as the interface's compare step
runs it (stage1_igv_assistant/ui.py, runCompare and /api/compare):
  total records        load_candidate_set: total_records
  BND records          load_candidate_set: counts_by_svtype["BND"]
  BND after dedup      list_candidates(svtype=BND): junctions after the tools'
                       500 bp orientation-aware deduplication
  PASS BND             ... + filter_pass
  before recurrence    ... + min_pe 3, min_sr 1, both breakends on primary
                       contigs, neither in the unmodified exclude template
  after recurrence     survivors of one sample NOT matched, both breakends
                       within 500 bp, by compare_candidate_sets against the
                       other sample's whole set
Any difference from a Phase 3 count is a finding, reported as such.
"""
import json
import os
import re
import sys
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, REPO)
from common import PHASE0, RUN_DIR  # noqa: E402
from stage1_igv_assistant.tools import vcf_tools  # noqa: E402

EXCL = os.path.expanduser("~/reference/human.hg38.excl.tsv")
BCF = os.path.join(RUN_DIR, "delly", "{label}.bcf")
JOBS = os.path.join(RUN_DIR, "jobs")
RECORD = os.path.join(REPO, "stage1_igv_assistant", "results", "patient_rerun_2026-09.json")

PHASE3 = {
    "SAMPLE_A": {"total_records": 30_980, "bnd_records": 9_187, "bnd_after_dedup": 9_172,
                 "pass_bnd": 896, "before_recurrence": 57, "after_recurrence_500bp": 17,
                 "peak_rss_mib_approx": 1_184, "wall_phase3": "45:33"},
    "SAMPLE_B": {"total_records": 32_451, "bnd_records": 9_676, "bnd_after_dedup": 9_655,
                 "pass_bnd": 923, "before_recurrence": 63, "after_recurrence_500bp": 19,
                 "peak_rss_mib_approx": 1_187, "wall_phase3": "2:03:43"},
}
FILTERS = {"svtype": "BND", "filter_pass": True, "min_pe": 3, "min_sr": 1, "primary_only": True}
COUNTED = ["total_records", "bnd_records", "bnd_after_dedup", "pass_bnd", "before_recurrence",
           "after_recurrence_500bp"]


def funnel():
    vcf_tools.reset_registry()
    sets, out = {}, {}
    for label in PHASE3:
        r = vcf_tools.load_candidate_set(BCF.format(label=label), label)
        if "error" in r:
            raise SystemExit(f"{label}: {r['error']}")
        sets[label] = r["set_id"]
        n = lambda **kw: vcf_tools.list_candidates(r["set_id"], limit=1, **kw)["total_matching"]
        full = vcf_tools.list_candidates(r["set_id"], mask_path=EXCL, limit=100000, **FILTERS)
        out[label] = {
            "total_records": r["total_records"],
            "bnd_records": r["counts_by_svtype"].get("BND", 0),
            "bnd_after_dedup": n(svtype="BND"),
            "pass_bnd": n(svtype="BND", filter_pass=True),
            "before_recurrence": full["total_matching"],
            "per_step": [{k: s[k] for k in ("filter", "value", "surviving_after_this_step",
                                            "removed_cumulatively_here",
                                            "would_remove_from_unfiltered_set")}
                         for s in full["filters_applied"]],
            "dedup": {k: r[k] for k in ("junctions_before_dedup", "junctions_after_dedup",
                                        "records_merged_by_dedup", "malformed_records_skipped",
                                        "unpaired_mate_records")},
            "counts_by_filter": r["counts_by_filter"], "counts_by_svtype": r["counts_by_svtype"],
            "caller_convention": r["caller_convention"],
            "_survivor_ids": {c["candidate_id"] for c in full["candidates"]},
        }
    cmp_ = vcf_tools.compare_candidate_sets(sets["SAMPLE_A"], sets["SAMPLE_B"], tolerance_bp=500)
    matched = {"SAMPLE_A": {p["candidate_id_a"] for p in cmp_["matched_pairs"]},
               "SAMPLE_B": {p["candidate_id_b"] for p in cmp_["matched_pairs"]}}
    for label in PHASE3:
        surv = out[label].pop("_survivor_ids")
        out[label]["after_recurrence_500bp"] = len(surv - matched[label])
        out[label]["survivors_recurrent_in_other_sample"] = len(surv & matched[label])
    out["comparison_whole_sets"] = {k: cmp_[k] for k in ("total_in_a", "total_in_b", "matched_in_a",
                                                         "matched_in_b", "unmatched_in_a", "unmatched_in_b")}
    out["filters"] = {**FILTERS, "exclude_template": "~/reference/human.hg38.excl.tsv (v2.6.0, unmodified)",
                      "dedup_tolerance_bp": vcf_tools.DEDUP_TOLERANCE_BP, "recurrence_tolerance_bp": 500}
    return out


def delly_run(label):
    log = open(os.path.join(JOBS, f"delly_{label[-1]}.log")).read()
    rc_path = os.path.join(JOBS, f"delly_{label[-1]}.rc")
    marks = [(datetime.strptime(m.group(1), "%Y-%b-%d %H:%M:%S"), m.group(2).strip())
             for m in re.finditer(r"^\[(\d{4}-\w{3}-\d{2} \d{2}:\d{2}:\d{2})\] (.+)$", log, re.M)]
    ended = re.search(r"^ended: (\S+) delly exit=(\d+)", log, re.M)
    phases = []
    for i, (t, name) in enumerate(marks):
        if name.startswith("delly ") or name.startswith("/"):
            continue                                   # the echoed command line, not a phase
        nxt = marks[i + 1][0] if i + 1 < len(marks) else (
            datetime.fromisoformat(ended.group(1)).replace(tzinfo=None) if ended else None)
        phases.append({"phase": name, "started": t.strftime("%H:%M:%S"),
                       "seconds": int((nxt - t).total_seconds()) if nxt else None})
    g = lambda pat: (re.search(pat, log).group(1) if re.search(pat, log) else None)
    rss = g(r"Maximum resident set size \(kbytes\): (\d+)")
    return {"exit_code": int(open(rc_path).read()) if os.path.exists(rc_path) else None,
            "delly_exit_reported": int(ended.group(2)) if ended else None,
            "wall_clock": g(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)"),
            "peak_rss_kbytes": int(rss) if rss else None,
            "peak_rss_mib": round(int(rss) / 1024) if rss else None,
            "user_seconds": g(r"User time \(seconds\): (\S+)"),
            "cpu_percent": g(r"Percent of CPU this job got: (\S+)"),
            "phases": phases,
            "command": "OMP_NUM_THREADS=1 /usr/bin/time -v ~/tools/delly/delly sr -h 1 -g "
                       "~/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa -x "
                       f"~/reference/human.hg38.excl.tsv -o ~/patient_data/rerun_2026-09/delly/{label}.bcf "
                       f"~/patient_data/deid/{label}.bam"}


def table(res):
    rows, same = [], True
    for label in PHASE3:
        for k in COUNTED:
            got, want = res[label][k], PHASE3[label][k]
            same &= got == want
            rows.append((label, k, got, want, "MATCH" if got == want else f"DIFFERS ({got - want:+,})"))
    return rows, same


def main(step):
    res = funnel()
    rows, same = table(res)
    for label, k, got, want, v in rows:
        print(f"{label}  {k:24s} {got:>8,}   Phase 3 {want:>8,}   {v}")
    for label in PHASE3:
        print(f"{label}  per step: " + " -> ".join(f"{s['filter']} {s['surviving_after_this_step']}"
                                                  for s in res[label]["per_step"]))
    print("ALL COUNTS MATCH PHASE 3" if same else "AT LEAST ONE COUNT DIFFERS FROM PHASE 3 -- a finding")
    with open(os.path.join(RUN_DIR, "funnel.json"), "w") as f:
        json.dump(res, f, indent=1)
    if step != "record":
        return
    ver = {}
    for label in PHASE0:
        v = json.load(open(os.path.join(RUN_DIR, f"verify_{label}_header.json")))
        fs = json.load(open(os.path.join(RUN_DIR, f"verify_{label}_flagstat.json")))["flagstat"]["QC-passed reads"]
        pa = json.load(open(os.path.join(RUN_DIR, f"verify_{label}_pa.json")))
        ver[label] = {
            "bytes": v["bytes"], "phase0_bytes": PHASE0[label]["bytes"],
            "index_mapped_chr1_22_X_Y": v["idxstats_mapped_chr1_22_X_Y"],
            "phase0_primary_mapped_figure": PHASE0[label]["primary_mapped"],
            "flagstat_primary_mapped": fs["primary mapped"], "flagstat_total": fs["total"],
            "flagstat_mapped": fs["mapped"], "index_mapped_sum": v["idxstats_mapped_sum"],
            "pg_records_samtools_head": v["pg_records_samtools_head"],
            "pg_records_parsed": v["pg_records_parsed"], "phase0_pg_records": PHASE0[label]["pg_records"],
            "contig_md5": v["contig_md5"], "ah_contigs": v["ah_contigs"],
            "ah_set_equals_hs38DH_alt": v.get("ah_set_equals_hs38DH_alt_file"),
            "pa_tagged_chr20_chr21": pa["pa_tagged"], "primary_reads_chr20_chr21": pa["primary_reads_scanned"],
            "bwa_flags_in_all_command_lines": {"-M": all(f["-M"] for f in v["bwa_cl_flags"]),
                                              "-Y": any(f["-Y"] for f in v["bwa_cl_flags"]),
                                              "-j": any(f["-j (ignore .alt)"] for f in v["bwa_cl_flags"])}}
    mem = {k: v for k, v in (l.split(":", 1) for l in open("/proc/meminfo") if l.split(":")[0]
                             in ("MemTotal", "SwapTotal"))}
    rec = {"what": "patient delly rerun and Phase 0/Phase 3 reproduction, 2026-09-25",
           "blinding": "counts only -- no coordinates, genes, candidates or interpretation",
           "samples": "SAMPLE_A and SAMPLE_B, assigned by Phase 0 byte size; identifiers never leave ~/patient_data",
           "machine": {"MemTotal": mem["MemTotal"].strip(), "SwapTotal": mem["SwapTotal"].strip(),
                       "nproc": os.cpu_count(),
                       "note": "the two delly runs were concurrent with each other and with other work; "
                               "wall times are new measurements, not comparable with Phase 3"},
           "task1_verification": ver,
           "task1_finding": ("the Phase 0 'primary mapped reads' figure equals idxstats-mapped over chr1-22, "
                             "X and Y (every mapped record, secondary and supplementary included), the only one "
                             "of 1,024 candidate definitions reproducing both figures; flagstat's primary-mapped "
                             "count exceeds it by 5,163,213 and 6,329,908"),
           "delly": {label: delly_run(label) for label in PHASE3},
           "funnel": {label: {k: res[label][k] for k in COUNTED + ["survivors_recurrent_in_other_sample",
                                                                   "per_step", "dedup", "counts_by_filter",
                                                                   "counts_by_svtype", "caller_convention"]}
                      for label in PHASE3},
           "funnel_filters": res["filters"], "comparison_whole_sets": res["comparison_whole_sets"],
           "phase3_targets": PHASE3,
           "reproduction": [{"sample": l, "count": k, "rerun": g, "phase3": w, "verdict": v}
                            for l, k, g, w, v in rows],
           "all_counts_match": same,
           "code": ["scripts/patient/verify_bams.py", "scripts/patient/mapped_definitions.py",
                    "scripts/patient/run_delly.sh", "scripts/patient/write_patient_record.py"]}
    with open(RECORD, "w") as f:
        json.dump(rec, f, indent=1)
    print(f"record -> {os.path.relpath(RECORD, REPO)}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "funnel")
