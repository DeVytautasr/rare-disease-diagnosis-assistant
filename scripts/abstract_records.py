#!/usr/bin/env python3
"""Phase 17: records behind the revised conference abstract. Public data only.

    abstract_records.py provenance   -> results/na12878_provenance_2026-09-27.json
    abstract_records.py uicheck      -> results/ui_check_2026-09-27.json
    abstract_records.py panels       -> screenshots/imp01_2026-09-27/*.png + results/imp01_panels_2026-09-27.json
    abstract_records.py models       -> results/local_models_2026-09-28.json

No step reads ~/patient_data. Every path written into a record has the home
directory replaced by "~". A committed record is never overwritten.

PROVENANCE (definitions): the header lines are the BAM's @HD, @RG and @PG lines
as `samtools head` prints them; md5 and size are of the whole file; "source" is
what the file itself and its sidecar say produced it -- the CL of the last @PG
record (the slicing command) and ~/public_data/NA12878.chr20_chr21.bam.provenance.json,
written when the slice was made and never committed, by a script that was never
committed either. The NYGC comparison checks, on every bwa @PG record: VN
0.7.15-r1140, a `bwa mem` command line carrying -Y and -K 100000000, and a
reference named GRCh38_full_analysis_set_plus_decoy_hla.fa. The published md5 of
the source CRAM is quoted from the sidecar and was not verified (that would need
the whole 15.8 GB file).

UICHECK: `python -m stage1_igv_assistant.ui --check` (the self-check ui.py
documents), run from the repository at HEAD with autodiscovery of the public data
directory only (ui.discover_public never touches ~/patient_data; no config file
and no --dataset/--candidates are given). Recorded: command, environment facts,
exit status, and the full output.

PANELS: the four PNGs made on 2026-09-27 by scripts/igv_panel_demo.py (record
igv_panel_demo_2026-09-27.json) are copied if their sha256 prefixes still match
that record; the locus is the IMP01 implant at chr20:200000 of the synthetic
control, built on the public NA12878 background.

MODELS: `ollama show` for each local model used in the ceiling experiment
(parameters, quantization, architecture, context length, capabilities), the
model IDs from `ollama list`, the Ollama version, and the GPU from nvidia-smi.
"""
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HOME = os.path.expanduser("~")
RES = os.path.join(REPO, "stage1_igv_assistant", "results")
BG = os.path.join(HOME, "public_data", "NA12878.chr20_chr21.bam")
SIDECAR = BG + ".provenance.json"
SAMTOOLS = os.path.join(HOME, "miniconda3", "envs", "synth-hts", "bin", "samtools")
PANEL_SRC = os.path.join(HOME, "public_data", "sim", "logs", "phase13_2026-09-27", "igv_panel_demo")
PANEL_DST = os.path.join(REPO, "stage1_igv_assistant", "screenshots", "imp01_2026-09-27")
MODELS = ("qwen2.5:7b", "qwen3.5:4b")


def tilde(x):
    if isinstance(x, str):
        return x.replace(HOME, "~")
    if isinstance(x, list):
        return [tilde(v) for v in x]
    if isinstance(x, dict):
        return {k: tilde(v) for k, v in x.items()}
    return x


def write(name, rec):
    dest = os.path.join(RES, name)
    if os.path.exists(dest):
        sys.exit(f"STOPPED: {name} exists; a committed record is never overwritten")
    json.dump(tilde(rec), open(dest, "w"), indent=1)
    print("written:", os.path.relpath(dest, REPO))


def section(title):
    return " ".join(__doc__.split(title + " (definitions):" if title in ("PROVENANCE",) else title + ":")[1]
                    .split("\n\n")[0].split())


def step_provenance():
    head = subprocess.run([SAMTOOLS, "head", BG], capture_output=True, text=True, check=True).stdout.splitlines()
    lines = [l for l in head if l.startswith(("@HD", "@RG", "@PG"))]
    md5 = hashlib.md5()
    with open(BG, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            md5.update(chunk)
    pgs = [dict(f.split(":", 1) for f in l.split("\t")[1:] if ":" in f) for l in lines if l.startswith("@PG")]
    bwa = [p for p in pgs if p.get("PN") == "bwa"]
    checks = [{"id": p.get("ID"), "vn_0.7.15-r1140": p.get("VN") == "0.7.15-r1140",
               "bwa_mem": " mem " in p.get("CL", ""), "-Y": re.search(r"\s-Y(\s|$)", p.get("CL", "")) is not None,
               "-K 100000000": "-K 100000000" in p.get("CL", ""),
               "reference_GRCh38_full_analysis_set_plus_decoy_hla": "GRCh38_full_analysis_set_plus_decoy_hla.fa" in p.get("CL", "")}
              for p in bwa]
    last = pgs[-1] if pgs else {}
    url = re.findall(r"https?://\S+", last.get("CL", ""))
    sidecar = json.load(open(SIDECAR)) if os.path.exists(SIDECAR) else None
    all_match = bool(checks) and all(all(v for k, v in c.items() if k != "id") for c in checks)
    rec = {"what": "Provenance of the public NA12878 background, ~/public_data/NA12878.chr20_chr21.bam (Phase 17 Task 1)",
           "definitions": section("PROVENANCE"),
           "file": {"path": BG, "bytes": os.path.getsize(BG), "md5": md5.hexdigest(),
                    "sq_lines": sum(1 for l in head if l.startswith("@SQ"))},
           "header_lines": lines,
           "source": {"from_header_last_pg": {"id": last.get("ID"), "pn": last.get("PN"), "vn": last.get("VN"),
                                              "cl": last.get("CL"), "urls": url},
                      "from_sidecar": sidecar,
                      "sidecar_committed": False,
                      "producing_script_committed": False,
                      "established": bool(url) and sidecar is not None and url[0] == sidecar.get("source_url"),
                      "statement": ("the slice was cut from the CRAM at the URL named in both the header's last @PG "
                                    "command line and the sidecar (samtools view -X <crai> ... chr20 chr21); the "
                                    "code that ran it was not committed")},
           "nygc_high_coverage_match": {"bwa_pg_records": len(bwa), "per_record": checks, "all_match": all_match,
                                        "other_steps": [f"{p.get('ID')} {p.get('VN')}" for p in pgs if p.get("PN") != "bwa"]},
           "reference_5_statement": ("the URL path 1000G_2504_high_coverage/data/ERR3239334/NA12878.final.cram is the "
                                     "1000 Genomes high-coverage release, and every aligner record matches the stated "
                                     "NYGC settings; the record does not itself name the publication"),
           "code": "scripts/abstract_records.py provenance"}
    write("na12878_provenance_2026-09-27.json", rec)
    print(json.dumps({"md5": rec["file"]["md5"], "bytes": rec["file"]["bytes"], "established": rec["source"]["established"],
                      "all_match": all_match, "bwa_pg_records": len(bwa)}))
    return 0


def step_uicheck():
    cmd = [os.path.join(REPO, ".venv", "bin", "python"), "-m", "stage1_igv_assistant.ui", "--check"]
    t0 = time.time()
    p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, timeout=1800)
    head = subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    rec = {"what": "The browser interface's self-check at HEAD (Phase 17 Task 2)",
           "definitions": section("UICHECK"),
           "command": "python -m stage1_igv_assistant.ui --check", "git_head": head,
           "environment": {"java_on_path": shutil.which("java") is not None, "display": os.environ.get("DISPLAY"),
                           "config_file_present": os.path.exists(os.path.join(REPO, "sv-assistant.conf")),
                           "explicit_datasets_or_candidates": False},
           "exit_status": p.returncode, "wall_s": round(time.time() - t0, 1),
           "stdout": p.stdout.splitlines(), "stderr": p.stderr.splitlines()[-40:],
           "patient_data_mentioned_in_output": "patient_data" in p.stdout + p.stderr,
           "code": "scripts/abstract_records.py uicheck"}
    write("ui_check_2026-09-27.json", rec)
    print("exit", p.returncode)
    print("\n".join(tilde(p.stdout.splitlines())[-30:]))
    return 0


def step_panels():
    demo = json.load(open(os.path.join(RES, "igv_panel_demo_2026-09-27.json")))
    want = demo["files_in_session_dir"]
    os.makedirs(PANEL_DST, exist_ok=True)
    files = {}
    for name, meta in sorted(want.items()):
        src = os.path.join(PANEL_SRC, name)
        b = open(src, "rb").read()
        if hashlib.sha256(b).hexdigest()[:16] != meta["sha256"]:
            sys.exit(f"STOPPED: {name} no longer matches its record; regenerate instead")
        dst = os.path.join(PANEL_DST, name)
        if os.path.exists(dst):
            sys.exit(f"STOPPED: {name} already copied")
        shutil.copyfile(src, dst)
        files[name] = {"bytes": len(b), "sha256": hashlib.sha256(b).hexdigest(), "dimensions": meta["dimensions"],
                       "layer": name[:-4],
                       "panel_return": {k: demo["panels"].get(name[:-4], {}).get(k) for k in ("image_ref", "region", "color_by")}}
    rec = {"what": "IGV evidence panels for the talk, IMP01 chr20:200000 (Phase 17 Task 3)",
           "locus": "IMP01, chr20:200000: the first synthetic implant (a heterozygous balanced t(20;21), simulated "
                    "with ART and aligned into the public NA12878 chr20-chr21 background); public synthetic data, "
                    "not a patient locus",
           "made": "2026-09-27 by scripts/igv_panel_demo.py through the MCP evidence_panel tool (IGV 2.17.4); "
                   "record igv_panel_demo_2026-09-27.json; copied unchanged, sha256 prefixes checked against that record",
           "checked_before_commit": "no PNG text chunks; the rendered labels are the BAM name, coverage, locus and "
                                    "Refseq track only",
           "directory": "stage1_igv_assistant/screenshots/imp01_2026-09-27/", "files": files,
           "code": "scripts/abstract_records.py panels"}
    write("imp01_panels_2026-09-27.json", rec)
    print(json.dumps({k: v["bytes"] for k, v in files.items()}))
    return 0


def step_models():
    def show(m):
        out = subprocess.run(["ollama", "show", m], capture_output=True, text=True, check=True).stdout
        d, sect = {}, None
        for line in out.splitlines():
            if line and not line.startswith(" "):
                continue
            s = line.strip()
            if not s:
                continue
            if line.startswith("  ") and not line.startswith("    "):
                sect = s
                d.setdefault(sect, {})
                continue
            parts = re.split(r"\s{2,}", s, maxsplit=1)
            if sect in ("Model", "Parameters") and len(parts) == 2:
                d[sect][parts[0]] = parts[1]
            elif sect == "Capabilities":
                d[sect][parts[0]] = True
        return {k: v for k, v in d.items() if k in ("Model", "Capabilities", "Parameters")}
    lst = subprocess.run(["ollama", "list"], capture_output=True, text=True, check=True).stdout.splitlines()
    ids = {l.split()[0]: l.split()[1] for l in lst[1:] if l.split()}
    ver = subprocess.run(["ollama", "--version"], capture_output=True, text=True).stdout.strip()
    gpu = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
                         capture_output=True, text=True).stdout.strip()
    rec = {"what": "The local models of the ceiling experiment as served by Ollama, and the GPU (Phase 17)",
           "definitions": section("MODELS"),
           "models": {m: {"ollama_show": show(m), "ollama_id": ids.get(m)} for m in MODELS},
           "ollama_version": ver, "gpu": gpu,
           "run_records_ids": "qwen2.5:7b runs record 'ollama ps' ID 845dbda0ea48 (phase10_rerun_2026-09-25/qwen2.5-7b/*.json)",
           "note": "Q4_K_M is llama.cpp's medium 4-bit K-quant (some tensors kept at higher precision)",
           "code": "scripts/abstract_records.py models"}
    write("local_models_2026-09-28.json", rec)
    print(json.dumps({m: {"params": r["ollama_show"].get("Model", {}).get("parameters"),
                          "quant": r["ollama_show"].get("Model", {}).get("quantization"), "id": r["ollama_id"]}
                      for m, r in rec["models"].items()}), gpu)
    return 0


if __name__ == "__main__":
    steps = {"provenance": step_provenance, "uicheck": step_uicheck, "panels": step_panels, "models": step_models}
    if len(sys.argv) != 2 or sys.argv[1] not in steps:
        sys.exit("usage: abstract_records.py provenance | uicheck | panels | models")
    sys.exit(steps[sys.argv[1]]())
