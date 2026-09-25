#!/usr/bin/env python3
"""The delly MAPQ ladder (Task 4b of the Phase 6 rerun): what lowering delly's
mapping-quality floors buys at the missed implants, and what it costs on the
untouched background.

    delly_ladder.py run LABEL [LABEL ...]       every setting x label, all at once
    delly_ladder.py analyse LABEL [LABEL ...]   junctions called and surviving; background cost

SETTINGS (delly v2.6.0 sr: -q min. paired-end mapping quality, -r min. PE quality
for translocations)
    q1_r20   -q 1 -r 20   the defaults, rerun -- a determinism check: its records
                          must equal the committed run's BCF record for record
    q1_r5    -q 1 -r 5
    q0_r5    -q 0 -r 5
    q0_r0    -q 0 -r 0
Everything else exactly as scripts/synthetic/run_delly_synth.sh: sr, one thread
(-h 1, OMP_NUM_THREADS=1), the unmodified v2.6.0 exclude template, hs38DH, the
same pre-flight checks (delly sha256, template git blob, .fai md5). Output under
~/public_data/sim/ladder/<setting>/<LABEL>.bcf -- outside every directory the
interface autodiscovers -- with /usr/bin/time -v in the log beside it. Every job
starts at once, each after MemAvailable >= 3000 MiB. Exit 0 only if every job
exits 0.

analyse loads each BCF through the bridge (MCP dispatch, recorded), lists it
unfiltered and through the standard funnel, and reports for every implant label
which of its two junctions was called (evidence_chain's matching rule) and
whether it survived, and for every label the raw record count, the BND count
and the survivors.
"""
import argparse
import hashlib
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
from synth_common import BAMS, BCFTOOLS, BG, DELLY, DELLY_DIR, EXCL, REF, SIM  # noqa: E402
import evidence_chain as ec  # noqa: E402

LADDER = os.path.join(SIM, "ladder")
SETTINGS = [("q1_r20", 1, 20), ("q1_r5", 1, 5), ("q0_r5", 0, 5), ("q0_r0", 0, 0)]
MINFREE_MIB = 3000
PREFLIGHT = {"delly sha256": ("sha256", DELLY, "85ecf4d64e23672a51c71f6e3a2dfda997753266157ade79ef598c8f96ecb469"),
             "exclude template git blob (v2.6.0, unmodified)": ("blob", EXCL, "3125a61491cb1c62ce46b365897eb4cd4e9fdb66"),
             "hs38DH .fai md5": ("md5", REF + ".fai", "5ccc91e56dc4a05448dd5b9507ec6bc6")}


def digest(kind, path):
    if kind == "blob":
        return subprocess.run(["git", "hash-object", "--no-filters", path], capture_output=True,
                              text=True, check=True).stdout.strip()
    h = hashlib.sha256() if kind == "sha256" else hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def preflight():
    ok = True
    for name, (kind, path, want) in PREFLIGHT.items():
        got = digest(kind, path)
        print(f"pre-flight {'OK  ' if got == want else 'FAIL'} {name}" + ("" if got == want else f": got {got}"))
        ok &= got == want
    return ok


def mem_available_mib():
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable"):
            return int(line.split()[1]) // 1024
    return 0


def delly_cmd(bam, out, q, r):
    return [DELLY, "sr", "-h", "1", "-q", str(q), "-r", str(r), "-g", REF, "-x", EXCL, "-o", out, bam]


def launch(label, bam, outdir, q, r, logdir):
    """One delly job under /usr/bin/time -v; returns (Popen, out, log)."""
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(logdir, exist_ok=True)
    out = os.path.join(outdir, f"{label}.bcf")
    log = os.path.join(logdir, f"{os.path.basename(outdir)}__{label}.log")
    while mem_available_mib() < MINFREE_MIB:
        time.sleep(10)
    env = dict(os.environ, OMP_NUM_THREADS="1")
    fh = open(log, "w")
    fh.write("command: OMP_NUM_THREADS=1 /usr/bin/time -v " + " ".join(delly_cmd(bam, out, q, r))
             .replace(os.path.expanduser("~"), "~") + "\n")
    fh.flush()
    p = subprocess.Popen(["/usr/bin/time", "-v", *delly_cmd(bam, out, q, r)], stdout=fh, stderr=subprocess.STDOUT,
                         env=env)
    return p, out, log, fh


def job_summary(out, log, rc):
    text = open(log).read()
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    wall = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)", text)
    rec = {"exit_code": rc, "peak_rss_kbytes": int(rss.group(1)) if rss else None,
           "wall_clock": wall.group(1) if wall else None}
    if rc == 0:
        s = subprocess.run([BCFTOOLS, "query", "-l", out], capture_output=True, text=True).stdout.split()
        n = subprocess.run([BCFTOOLS, "view", "-H", out], capture_output=True, text=True).stdout.count("\n")
        rec.update({"samples_in_bcf": len(s), "sample_names": s, "records_in_bcf": n})
    return rec


def records(bcf):
    """The BCF's records, header excluded (the header carries the date and command)."""
    return subprocess.run([BCFTOOLS, "view", "-H", bcf], capture_output=True, text=True, check=True).stdout


def bam_for(label):
    return BG if label == "background" else os.path.join(BAMS, f"{label}.bam")


def step_run(labels):
    if not preflight():
        print("REFUSING to run delly: a pre-flight check failed")
        return 3
    jobs = []
    for name, q, r in SETTINGS:
        for label in labels:
            p, out, log, fh = launch(label, bam_for(label), os.path.join(LADDER, name), q, r,
                                     os.path.join(LADDER, "logs"))
            jobs.append((name, label, p, out, log, fh))
            print(f"launched {name} {label} (pid {p.pid})", flush=True)
    summary = {}
    for name, label, p, out, log, fh in jobs:
        rc = p.wait()
        fh.close()
        summary[f"{name}/{label}"] = job_summary(out, log, rc)
        print(f"{name} {label}: exit {rc}, {summary[f'{name}/{label}']}", flush=True)
    det = {}
    for label in labels:
        new = os.path.join(LADDER, "q1_r20", f"{label}.bcf")
        old = os.path.join(DELLY_DIR, f"{label}.bcf")
        if os.path.exists(new) and os.path.exists(old):
            det[label] = records(new) == records(old)
    rec = {"settings": {n: {"q": q, "r": r} for n, q, r in SETTINGS}, "labels": labels, "jobs": summary,
           "default_rerun_identical_to_committed_run": det}
    ec.write_record(os.path.join(ec.OUT, "ladder", "run.json"), rec)
    print(f"default rerun record-for-record identical to the committed run: {det}")
    return 0 if all(v["exit_code"] == 0 for v in summary.values()) else 1


def step_analyse(labels):
    T = ec.truth()
    tol = ec.tolerance()
    rec = ec.recorder()
    out = {}
    for name, q, r in [("committed_default", 1, 20)] + SETTINGS:
        for label in labels:
            bcf = (os.path.join(DELLY_DIR, f"{label}.bcf") if name == "committed_default"
                   else os.path.join(LADDER, name, f"{label}.bcf"))
            ld = ec.result(rec.call("bridge", "load_candidate_set", {"path": bcf, "label": f"{label}:{name}"}))
            unf = ec.result(rec.call("bridge", "list_candidates", {"set_id": ld["set_id"], "limit": ec.LIMIT}))
            fun = ec.result(rec.call("bridge", "list_candidates",
                                     {"set_id": ld["set_id"], "filter_pass": True, "min_pe": 3, "min_sr": 1,
                                      "primary_only": True, "mask_path": EXCL, "limit": ec.LIMIT}))
            surv = {c["candidate_id"] for c in fun["candidates"]}
            row = {"total_records": ld["total_records"], "bnd_records": ld["counts_by_svtype"].get("BND", 0),
                   "junctions_after_dedup": ld["junctions_after_dedup"], "survivors": len(surv),
                   "survivor_svtypes": {}}
            for c in fun["candidates"]:
                row["survivor_svtypes"][c["svtype"]] = row["survivor_svtypes"].get(c["svtype"], 0) + 1
            if label != "background":
                js = {}
                for jn, exp in ec.expected_junctions(T[label]).items():
                    m, _ = ec.best_match(unf["candidates"], exp, tol)
                    js[jn] = {"called": m is not None,
                              "survives_funnel": bool(m and m["candidate_id"] in surv),
                              "call": m and {k: m[k] for k in ("pos1", "pos2", "orientation", "filter", "pe", "sr",
                                                               "precise")}}
                    if m and not js[jn]["survives_funnel"]:
                        # every step that would remove it, by the bridge's own predicates
                        from stage1_igv_assistant.tools.vcf_tools import PRIMARY_CONTIGS, in_mask, load_mask
                        mask = load_mask(EXCL)
                        js[jn]["would_be_removed_by"] = [n for n, hit in (
                            ("PASS", m["filter"] != "PASS"), ("PE>=3", (m["pe"] or 0) < 3),
                            ("SR>=1", (m["sr"] or 0) < 1),
                            ("primary", not (m["chrom1"] in PRIMARY_CONTIGS and m["chrom2"] in PRIMARY_CONTIGS)),
                            ("unmasked", in_mask(mask, m["chrom1"], m["pos1"]) or in_mask(mask, m["chrom2"], m["pos2"])))
                            if hit]
                row["junctions"] = js
                row["implanted_called"] = sum(v["called"] for v in js.values())
                row["implanted_surviving"] = sum(v["survives_funnel"] for v in js.values())
                row["non_implanted_survivors"] = len(surv) - row["implanted_surviving"]
            out.setdefault(name, {})[label] = row
    summary = {}
    for name in out:
        imps = [l for l in labels if l != "background"]
        summary[name] = {
            "missed_junctions_called": sum(out[name][l]["implanted_called"] for l in imps),
            "missed_junctions_surviving": sum(out[name][l]["implanted_surviving"] for l in imps),
            "of": 2 * len(imps),
            "background_records": out[name].get("background", {}).get("total_records"),
            "background_bnd": out[name].get("background", {}).get("bnd_records"),
            "background_survivors": out[name].get("background", {}).get("survivors")}
    base = summary["committed_default"]
    for name in summary:
        if base["background_bnd"]:
            summary[name]["background_bnd_vs_default"] = round(summary[name]["background_bnd"] / base["background_bnd"], 2)
    rec_out = {"labels": labels, "per_setting": out, "summary": summary, "calls": rec.calls}
    ec.write_record(os.path.join(ec.OUT, "ladder", "analysis.json"), rec_out)
    for name, s in summary.items():
        print(f"{name:18s} missed junctions called {s['missed_junctions_called']}/{s['of']}, surviving "
              f"{s['missed_junctions_surviving']}/{s['of']}; background records {s['background_records']}, "
              f"BND {s['background_bnd']} (x{s.get('background_bnd_vs_default')}), survivors {s['background_survivors']}")
        for l in labels:
            if l != "background":
                j = out[name][l]["junctions"]
                print(f"    {l}: " + "; ".join(f"{jn} called={v['called']} surviving={v['survives_funnel']}"
                                             + (f" {v['call']}" if v["call"] else "") for jn, v in j.items()))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["run", "analyse"])
    ap.add_argument("labels", nargs="+")
    a = ap.parse_args()
    bad = [l for l in a.labels if l not in ec.LABELS]
    if bad:
        raise SystemExit(f"unknown labels {bad}")
    sys.exit(step_run(a.labels) if a.step == "run" else step_analyse(a.labels))
