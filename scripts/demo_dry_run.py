#!/usr/bin/env python3
"""Phase 18 Task 2: a dry run of the meeting sheet's demonstration. Public data only.

    demo_dry_run.py run      -> ~/public_data/sim/logs/phase18_2026-09-29/attempt_N.json (+ the UI's log)
    demo_dry_run.py record   -> stage1_igv_assistant/results/demo_dry_run_2026-09-29.json

The interface is started as the sheet's demo section describes -- `python -m
stage1_igv_assistant.ui`, repository root as working directory, no flags -- so it
autodiscovers ~/public_data and nothing else. Checked before the start and recorded:
no config file on the three paths config.py searches, no SV_* variable in the
environment, port 8765 free. The interface's environment adds only the IGV Java
runtime to PATH, DISPLAY=:0, and ANTHROPIC_BASE_URL=http://127.0.0.1:9 (a closed
local port), so an Anthropic API call, had one been attempted, could not have left
the machine; a positive control shows the SDK honours that variable. It is stopped
by its PID (SIGTERM, then SIGKILL if still running after 15 s), and any child
process still alive afterwards is stopped by its PID too.

The seven steps are walked through the interface's own HTTP routes with the bodies
its page sends. The page is fetched and checked to call those routes, but it is not
rendered and nothing is clicked: "shown" means the value the page renders from the
response, read from the response.
  1 load         GET /api/bootstrap; POST /api/load with the candidate set IMP01
  2 filters      POST /api/funnel with the page's defaults (FILTER PASS, PE >= 3,
                 SR >= 1, both ends on primary contigs, exclude template on)
  3 candidate    POST /api/candidate for the surviving candidate joining chr20:200,000
                 and chr21:14,100,001, then POST /api/assess at both breakends in
                 dataset IMP01, as the page's "evidence" button does
  4 tool return  GET /api/call?id= for the calls behind the discordant-pair count and
                 the combined score at chr20:200,000 (the page's call chips)
  5 hand entry   POST /api/assess in dataset IMP09 at chr20:33,700,000, a breakend of an
                 implant the caller missed (so no candidate carries it), set against
                 the committed rescue record (analysis_2026-09-25/figures.json)
  6 IGV          POST /api/igv for IMP01 chr20:200,000; every image fetched via /img/
  7 chat         GET /api/chat_models; POST /api/chat with one question to qwen2.5:7b
                 at the page's defaults (context 32768, thinking off)

FIGURES AND THEIR DEFINITIONS
  wall_s           seconds from sending a request to its response; startup_s from
                   process start to the first answered /api/bootstrap
  numbers_matched  the verification pass the page shows under a chat answer: every
                   number in the final answer checked against every tool return of that
                   turn (chat.verify_numbers: within 0.011, their roundings and
                   percentages, or literal text); `unsupported` lists the others
  no_patient_data  (a) the listed labels equal those discover_public's rules derive
                   from ~/public_data alone; (b) no label contains an entry of the
                   identifier list (a count; positive control: a label built from the
                   list's first entry is caught; no entry is printed or recorded);
                   (c) every dataset or candidate file named in the session's call log
                   is a listed label (positive control: an unlisted label is caught);
                   (d) every regular file held open by the interface or its children,
                   sampled after each step, lies under an allowed public or system
                   location (a count; positive control: a path under ~/patient_data is
                   caught -- a string test, nothing there is opened)
"""
import glob
import hashlib
import json
import os
import re
import signal
import socket
import struct
import subprocess
import sys
import time
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HOME = os.path.expanduser("~")
PUB = os.path.join(HOME, "public_data")
LOGDIR = os.path.join(PUB, "sim", "logs", "phase18_2026-09-29")
DEST = os.path.join(REPO, "stage1_igv_assistant", "results", "demo_dry_run_2026-09-29.json")
RES = os.path.join(REPO, "stage1_igv_assistant", "results")
URL = "http://127.0.0.1:8765"
GUARD = "http://127.0.0.1:9"
MODEL = "qwen2.5:7b"
QUESTION = ("In dataset IMP01, how strong is the evidence for a breakpoint at chr20:200000? "
            "Give the evidence score and its band, and say whether this position could reach "
            "the strong band.")
CANDIDATE = {("chr20", 200000), ("chr21", 14100001)}
HAND = ("IMP09", "chr20", 33700000)
ALLOWED = [REPO, PUB, os.path.join(HOME, "reference"), os.path.join(HOME, "IGV_2.17.4"),
           os.path.join(HOME, "igv"), os.path.join(HOME, "miniconda3"), os.path.join(HOME, ".cache"),
           "/usr", "/lib", "/lib64", "/etc", "/dev", "/proc", "/sys", "/tmp", "/var", "/run", "/opt"]


def tilde(x):
    if isinstance(x, str):
        return x.replace(HOME, "~")
    if isinstance(x, list):
        return [tilde(v) for v in x]
    if isinstance(x, dict):
        return {k: tilde(v) for k, v in x.items()}
    return x


def definitions():
    doc = __doc__.split("FIGURES AND THEIR DEFINITIONS")[1]
    out, key = {}, None
    for line in doc.splitlines():
        if line.startswith("  ") and not line.startswith("    ") and line.strip():
            key, _, rest = line.strip().partition(" ")
            out[key] = rest.strip()
        elif key and line.strip():
            out[key] += " " + line.strip()
    return {k: " ".join(v.split()) for k, v in out.items()}


def http(method, path, body=None, timeout=180):
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(URL + path, data=data, method=method,
                                 headers={"Content-Type": "application/json"} if data else {})
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=timeout) as r:
        raw = r.read()
        ctype = r.headers.get("Content-Type", "")
    dt = round(time.time() - t0, 2)
    return (json.loads(raw) if "json" in ctype else raw), dt


# ── no-patient-data checks ───────────────────────────────────────────────────
def derived_labels():
    """discover_public's rules with nothing registered by a config file, applied here
    to ~/public_data alone: BAMs by label only (its realpath check compares against
    config-registered files, of which there are none); candidate files also skipped
    when their realpath is already registered."""
    ds, cf = {}, {}
    for root in (PUB, os.path.join(PUB, "sim", "bams")):
        for f in sorted(os.listdir(root)) if os.path.isdir(root) else []:
            if f.endswith(".bam"):
                ds.setdefault(os.path.splitext(f)[0], os.path.join(root, f))
    for root in (os.path.join(PUB, "delly"), os.path.join(PUB, "sim", "delly")):
        for f in sorted(os.listdir(root)) if os.path.isdir(root) else []:
            full = os.path.join(root, f)
            if f.endswith((".bcf", ".vcf", ".vcf.gz")) and \
                    os.path.realpath(full) not in {os.path.realpath(v) for v in cf.values()}:
                cf.setdefault(f.split(".")[0], full)
    return sorted(ds), sorted(cf)


def identifier_hits(labels):
    path = os.path.join(HOME, "patient_data", ".identifier_list")
    ids = [l.strip() for l in open(path) if l.strip() and not l.startswith("#")]
    match = lambda lab: any(i.lower() in lab.lower() for i in ids)  # noqa: E731
    return {"identifiers_read": len(ids),
            "labels_containing_an_identifier": sum(match(l) for l in labels),
            "positive_control_caught": bool(ids) and match("x" + ids[0] + "y")}


def outside_allowed(path):
    return path.startswith("/") and not any(path == a or path.startswith(a + "/") for a in ALLOWED)


def process_tree(root):
    kids = {}
    for st in glob.glob("/proc/[0-9]*/stat"):
        try:
            txt = open(st).read()
        except OSError:
            continue
        pid, ppid = int(st.split("/")[2]), int(txt.rsplit(")", 1)[1].split()[1])
        kids.setdefault(ppid, []).append(pid)
    out, todo = [], [root]
    while todo:
        p = todo.pop()
        out.append(p)
        todo.extend(kids.get(p, []))
    return out


def fd_sample(root):
    n_files, n_out = 0, 0
    for pid in process_tree(root):
        for fd in glob.glob(f"/proc/{pid}/fd/*"):
            try:
                t = os.readlink(fd)
            except OSError:
                continue
            if t.startswith("/"):
                n_files += 1
                n_out += outside_allowed(t.replace(" (deleted)", ""))
    return {"open_files": n_files, "outside_allowed": n_out}


def guard_control():
    """The SDK must send to ANTHROPIC_BASE_URL: a dummy key, no retries, a closed port."""
    import anthropic
    os.environ["ANTHROPIC_BASE_URL"] = GUARD
    try:
        c = anthropic.Anthropic(api_key="dummy-not-a-key", max_retries=0, timeout=5)
        base = str(c.base_url)
        try:
            c.messages.create(model="claude-sonnet-5", max_tokens=1,
                              messages=[{"role": "user", "content": "x"}])
            outcome = "a response came back"
        except Exception as e:  # noqa: BLE001
            outcome = type(e).__name__
    finally:
        del os.environ["ANTHROPIC_BASE_URL"]
    return {"sdk_base_url": base, "call_outcome": outcome,
            "pass": base.rstrip("/") == GUARD and outcome == "APIConnectionError"}


# ── the run ──────────────────────────────────────────────────────────────────
def shown(E):
    S, C = E.get("summary") or {}, E.get("ceiling") or {}
    L = {x["key"]: x for x in E.get("layers", [])}
    pp = E.get("position_provenance") or {}
    return {"at": f"{E.get('chromosome')}:{E.get('position')}", "dataset": E.get("bam_label"),
            "error": E.get("error"), "score": S.get("evidence_score"), "band": S.get("evidence_strength"),
            "components": S.get("components"), "low_mapq_fraction": S.get("low_mapq_fraction"),
            "ceiling": {k: C.get(k) for k in ("derivable", "attainable_here", "attainable_basis", "strong_band",
                                              "max_with_flat_depth", "max_all_layers")},
            "page_says_cannot_reach_strong": (C.get("attainable_here") is not None and C.get("strong_band") is not None
                                              and C["attainable_here"] < C["strong_band"]),
            "layers": {k: {f: v.get(f) for f in ("count", "fraction", "consensus", "max_clips", "partners",
                                                 "mean_depth", "quality_limited", "min_mapq", "call")}
                       for k, v in L.items()},
            "summary_call": S.get("call"),
            "position_provenance": pp.get("source"),
            "page_banner": ("This position came from the loaded call set." if pp.get("source") == "candidate_set"
                            else "This position was entered by hand.")}


def run():
    os.makedirs(LOGDIR, exist_ok=True)
    rec = {"what": "Dry run of the meeting sheet's seven demonstration steps (Phase 18 Task 2, 2026-09-29)",
           "definitions": definitions(), "code": "scripts/demo_dry_run.py",
           "git_head": subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"],
                                      capture_output=True, text=True).stdout.strip(),
           "started": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    cfg = [os.environ.get("SV_CONFIG", ""), os.path.join(REPO, "sv-assistant.conf"),
           os.path.join(HOME, ".config", "sv-assistant", "config.ini")]
    s = socket.socket()
    port_free = s.connect_ex(("127.0.0.1", 8765)) != 0
    s.close()
    rec["preconditions"] = {"config_files_present": sum(bool(p) and os.path.exists(p) for p in cfg),
                            "sv_environment_variables": sum(k.startswith("SV_") for k in os.environ),
                            "port_8765_free": port_free,
                            "anthropic_guard": {"ANTHROPIC_BASE_URL": GUARD, "control": guard_control()}}
    pre = rec["preconditions"]
    if pre["config_files_present"] or pre["sv_environment_variables"] or not port_free \
            or not pre["anthropic_guard"]["control"]["pass"]:
        print(json.dumps(pre, indent=1))
        sys.exit("STOPPED: a precondition failed; the interface was not started")

    env = dict(os.environ)
    env.pop("ANTHROPIC_API_KEY", None)
    env["ANTHROPIC_BASE_URL"] = GUARD
    env["PATH"] = os.path.join(HOME, "miniconda3", "envs", "igv-java", "bin") + os.pathsep + env["PATH"]
    env["DISPLAY"] = ":0"
    log_path = os.path.join(LOGDIR, "ui_stdout.log")
    log = open(log_path, "w")
    t_start = time.time()
    proc = subprocess.Popen([sys.executable, "-m", "stage1_igv_assistant.ui"], cwd=REPO, env=env,
                            stdout=log, stderr=subprocess.STDOUT)
    rec["interface"] = {"command": "python -m stage1_igv_assistant.ui", "interpreter": tilde(sys.executable),
                        "pid": proc.pid, "env_added": ["PATH += ~/miniconda3/envs/igv-java/bin", "DISPLAY=:0",
                                                       f"ANTHROPIC_BASE_URL={GUARD}"]}
    steps, fds = {}, {}
    try:
        boot = None
        while time.time() - t_start < 240:
            if proc.poll() is not None:
                raise RuntimeError(f"the interface exited during startup (status {proc.returncode})")
            try:
                boot, _ = http("GET", "/api/bootstrap", timeout=10)
                break
            except OSError:
                time.sleep(0.5)
        if boot is None:
            raise RuntimeError("the interface did not answer within 240 s")
        rec["interface"]["startup_s"] = round(time.time() - t_start, 1)
        fds["startup"] = fd_sample(proc.pid)
        page, dt = http("GET", "/")
        page = page.decode() if isinstance(page, bytes) else str(page)
        rec["interface"]["page"] = {"bytes": len(page), "wall_s": dt,
                                    "routes_called_by_page": {r: r in page for r in (
                                        "/api/bootstrap", "/api/load", "/api/funnel", "/api/candidate",
                                        "/api/assess", "/api/call?id=", "/api/igv", "/api/chat_models",
                                        "/api/chat")}}

        # 1 load
        t0 = time.time()
        boot, dt_boot = http("GET", "/api/bootstrap")
        ld, dt_load = http("POST", "/api/load", {"candidates_label": "IMP01"})
        lr = ld.get("result") or {}
        steps["1_load"] = {"dataset_labels_listed": boot["datasets"], "candidate_set_labels_listed": boot["candidate_files"],
                           "page_default_dataset": boot["datasets"][0] if boot["datasets"] else None,
                           "candidate_set_used": "IMP01",
                           "load": {k: lr.get(k) for k in ("set_id", "label", "caller_convention", "total_records",
                                                           "junctions_before_dedup", "junctions_after_dedup",
                                                           "records_merged_by_dedup", "counts_by_svtype")},
                           "is_error": ld.get("is_error"), "call": ld.get("call"),
                           "limits_panel_headings": [i["h"] for i in (boot.get("limits") or {}).get("items", [])],
                           "hand_entry_note": boot.get("hand_entry_note"),
                           "wall_s": {"bootstrap": dt_boot, "load": dt_load}}
        fds["1_load"] = fd_sample(proc.pid)
        set_id = lr["set_id"]

        # 2 filters
        fn, dt = http("POST", "/api/funnel", {"set_id": set_id, "svtype": None, "filter_pass": True, "min_pe": "3",
                                              "min_sr": "1", "primary_only": True, "use_mask": True})
        fr = fn.get("result") or {}
        steps["2_filters"] = {"total_in_set": fr.get("total_in_set"),
                              "chain": [{k: s_.get(k) for k in ("filter", "value", "provenance", "surviving_after_this_step",
                                                                "removed_cumulatively_here",
                                                                "would_remove_from_unfiltered_set")}
                                        for s_ in fr.get("filters_applied", [])],
                              "total_matching": fr.get("total_matching"), "returned": fr.get("returned"),
                              "mask": fn.get("mask"), "call": fn.get("call"), "wall_s": dt}
        fds["2_filters"] = fd_sample(proc.pid)

        # 3 candidate
        cands = fr.get("candidates", [])
        pick = [c for c in cands if {(c["chrom1"], c["pos1"]), (c["chrom2"], c["pos2"])} == CANDIDATE]
        if not pick:
            raise RuntimeError("the IMP01 candidate joining chr20:200000 and chr21:14100001 did not survive")
        cr, dt_c = http("POST", "/api/candidate", {"set_id": set_id, "candidate_id": pick[0]["candidate_id"]})
        c = cr["result"]
        ends, walls = [], []
        for be in (c["breakend_1"], c["breakend_2"]):
            E, dt = http("POST", "/api/assess", {"bam_label": "IMP01", "chromosome": be["chromosome"],
                                                  "position": be["position"]})
            ends.append(E)
            walls.append(dt)
        steps["3_candidate"] = {"candidate_id": pick[0]["candidate_id"], "funnel_row": pick[0],
                                "dataset_assessed": "IMP01", "shown": [shown(E) for E in ends],
                                "wall_s": {"get_candidate": dt_c, "assess_each_breakend": walls}}
        fds["3_candidate"] = fd_sample(proc.pid)

        # 4 one number's tool return
        e20 = next(E for E in ends if E.get("chromosome") == "chr20")
        s20 = shown(e20)
        dcall, dt1 = http("GET", f"/api/call?id={s20['layers']['discordant_pairs']['call']}")
        scall, dt2 = http("GET", f"/api/call?id={s20['summary_call']}")
        steps["4_tool_return"] = {
            "discordant_pairs": {"call": dcall.get("id"), "tool": dcall.get("tool"), "params": dcall.get("params"),
                                 "shown": s20["layers"]["discordant_pairs"]["count"],
                                 "in_return": (dcall.get("result") or {}).get("discordant_pairs"),
                                 "return_fields": sorted((dcall.get("result") or {}).keys())},
            "combined_score": {"call": scall.get("id"), "tool": scall.get("tool"), "shown": s20["score"],
                               "in_return": (scall.get("result") or {}).get("evidence_score"),
                               "band_in_return": (scall.get("result") or {}).get("evidence_strength")},
            "wall_s": [dt1, dt2]}
        t4 = steps["4_tool_return"]
        t4["equal"] = (t4["discordant_pairs"]["shown"] == t4["discordant_pairs"]["in_return"]
                       and t4["combined_score"]["shown"] == t4["combined_score"]["in_return"])
        fds["4_tool_return"] = fd_sample(proc.pid)

        # 5 hand entry of a missed implant
        H, dt = http("POST", "/api/assess", {"bam_label": HAND[0], "chromosome": HAND[1], "position": HAND[2]})
        fig = json.load(open(os.path.join(RES, "synthetic_control_2026-09", "analysis_2026-09-25", "figures.json")))
        committed = fig["rescue"][HAND[0]]["ends"][HAND[1]]
        hs = shown(H)
        steps["5_hand_entry"] = {"missed_implant": HAND[0], "position": f"{HAND[1]}:{HAND[2]}", "shown": hs,
                                 "committed_rescue_record": {k: committed.get(k) for k in (
                                     "position", "discordant_pairs", "soft_clipped_reads", "split_reads_min_mapq_0",
                                     "evidence_score", "evidence_strength", "low_mapq_fraction")},
                                 "agrees_with_rescue_record": {
                                     "discordant_pairs": hs["layers"]["discordant_pairs"]["count"] == committed.get("discordant_pairs"),
                                     "soft_clipped_reads": hs["layers"]["soft_clipped_reads"]["count"] == committed.get("soft_clipped_reads"),
                                     "score": hs["score"] == committed.get("evidence_score"),
                                     "band": hs["band"] == committed.get("evidence_strength")},
                                 "wall_s": dt}
        fds["5_hand_entry"] = fd_sample(proc.pid)

        # 6 IGV panels
        ig, dt = http("POST", "/api/igv", {"bam_label": "IMP01", "chromosome": "chr20", "position": 200000},
                      timeout=900)
        committed_png = json.load(open(os.path.join(RES, "imp01_panels_2026-09-27.json")))["files"]
        by_hash = {v["sha256"]: k for k, v in committed_png.items()}
        imgs = []
        for ref in ig.get("image_refs", []):
            b, dti = http("GET", f"/img/{ref}")
            h = hashlib.sha256(b).hexdigest()
            dims = list(struct.unpack(">II", b[16:24])) if b[:8] == b"\x89PNG\r\n\x1a\n" else None
            imgs.append({"ref": ref, "bytes": len(b), "sha256": h, "dimensions": dims,
                         "identical_to_committed_panel": by_hash.get(h), "fetch_wall_s": dti})
        panels = (ig.get("result") or {}).get("panels") or {}
        steps["6_igv"] = {"images": imgs, "error": ig.get("error"), "panel_errors": ig.get("panel_errors"),
                          "panels": {k: {f: v.get(f) for f in ("image_ref", "region", "color_by", "error")}
                                     for k, v in panels.items() if isinstance(v, dict)},
                          "call": ig.get("call"), "wall_s": dt}
        fds["6_igv"] = fd_sample(proc.pid)

        # 7 chat, local model only
        cm, dt_m = http("GET", "/api/chat_models", timeout=60)
        if MODEL not in (cm.get("models") or []):
            raise RuntimeError(f"{MODEL} is not offered by the chat panel")
        ch, dt = http("POST", "/api/chat", {"model": MODEL, "message": QUESTION, "num_ctx": 32768, "think": False},
                      timeout=1500)
        tools = []
        for e in ch.get("events", []):
            if e.get("type") == "tool":
                r = e.get("result") if isinstance(e.get("result"), dict) else {}
                tools.append({"name": e.get("name"), "params": e.get("params"), "call_id": e.get("call_id"),
                              "rejected_or_error": bool(e.get("rejected") or e.get("is_error")),
                              **{k: r.get(k) for k in ("evidence_score", "evidence_strength", "attainable_here",
                                                       "window_bp", "error") if k in r}})
        v = ch.get("verification") or {}
        steps["7_chat"] = {"models_offered": cm.get("models"), "ollama_reachable": cm.get("reachable"),
                           "tool_schemas_offered": cm.get("n_tools"), "model": MODEL, "backend": "ollama (chat.run_turn)",
                           "question": QUESTION, "settings": {"num_ctx": 32768, "think": False, "max_iters": "page default (8)"},
                           "error": ch.get("error"), "final_text": ch.get("final_text"), "tool_calls": tools,
                           "numbers_matched": {"numbers_in_prose": v.get("numbers_in_prose"),
                                               "unsupported": v.get("unsupported"), "no_tool_calls": v.get("no_tool_calls"),
                                               "details": v.get("details")},
                           "ended_without_answer": ch.get("ended_without_answer"),
                           "malformed_tool_calls": len(ch.get("malformed_tool_calls") or []),
                           "text_tool_call_failures": len(ch.get("text_tool_call_failures") or []),
                           "context_truncated": len(ch.get("context_truncated") or []),
                           "iterations": ch.get("iterations"), "gen_tokens": ch.get("gen_tokens"),
                           "wall_s": {"chat_models": dt_m, "chat": dt, "reported_by_chat": ch.get("wall_s")}}
        fds["7_chat"] = fd_sample(proc.pid)

        calls, _ = http("GET", "/api/calls")
        rec["call_log"] = {"calls": len(calls["calls"]),
                           "by_tool": {t: sum(c_["tool"] == t for c_ in calls["calls"])
                                       for t in sorted({c_["tool"] for c_ in calls["calls"]})},
                           "errors": sum(bool(c_.get("is_error")) for c_ in calls["calls"])}
        listed = set(boot["datasets"]) | set(boot["candidate_files"])
        named = []
        for c_ in calls["calls"]:
            p = c_.get("params") or {}
            for k in ("bam_path", "path"):
                if isinstance(p.get(k), str):
                    named.append(p[k])
            for x in p.get("bam_paths") or []:
                named.append(x)
        labels_named = [re.fullmatch(r"<(.+)>", x).group(1) if re.fullmatch(r"<(.+)>", x) else x for x in named]
        ds, cf = derived_labels()
        rec["no_patient_data"] = {
            "a_listed_equal_derived_from_public_data": {"datasets": sorted(boot["datasets"]) == ds,
                                                        "candidate_files": sorted(boot["candidate_files"]) == cf},
            "b_identifier_list": identifier_hits(sorted(listed)),
            "c_call_log": {"files_named": len(labels_named),
                           "not_a_listed_label": sum(x not in listed for x in labels_named),
                           "positive_control_caught": "<UNLISTED>" not in listed},
            "d_open_files": {"samples": fds, "outside_allowed_total": sum(f["outside_allowed"] for f in fds.values()),
                             "allowed": tilde(ALLOWED),
                             "positive_control_caught": outside_allowed(os.path.join(HOME, "patient_data", "x"))},
        }
    except Exception as e:  # noqa: BLE001  -- recorded, then the interface is stopped all the same
        rec["failure"] = f"{type(e).__name__}: {e}"
    finally:
        tree = process_tree(proc.pid)
        stop = {"pid": proc.pid, "children_before_stop": len(tree) - 1}
        if proc.poll() is None:
            os.kill(proc.pid, signal.SIGTERM)
            try:
                proc.wait(timeout=15)
                stop["method"] = "SIGTERM"
            except subprocess.TimeoutExpired:
                os.kill(proc.pid, signal.SIGKILL)
                proc.wait(timeout=15)
                stop["method"] = "SIGKILL after SIGTERM went unanswered for 15 s"
        else:
            stop["method"] = "had already exited"
        stop["exit_status"] = proc.returncode
        left = [p for p in tree[1:] if os.path.exists(f"/proc/{p}")]
        for p in left:
            os.kill(p, signal.SIGTERM)
        stop["children_left_and_stopped_by_pid"] = len(left)
        time.sleep(1)
        s = socket.socket()
        stop["port_8765_free_after"] = s.connect_ex(("127.0.0.1", 8765)) != 0
        s.close()
        rec["interface"]["stopped"] = stop
        rec["total_wall_s"] = round(time.time() - t_start, 1)
        log.close()
    ui_log = open(log_path).read()
    rec["interface"]["log_mentions_patient_data"] = ui_log.count("patient_data")
    if rec["interface"]["log_mentions_patient_data"] == 0:
        rec["interface"]["banner"] = tilde(ui_log.splitlines()[:40])
    rec["steps"] = steps
    n = len(glob.glob(os.path.join(LOGDIR, "attempt_*.json"))) + 1
    out = os.path.join(LOGDIR, f"attempt_{n}.json")
    json.dump(tilde(rec), open(out, "w"), indent=1)
    print("written:", tilde(out))
    print(json.dumps({"failure": rec.get("failure"), "stopped": rec["interface"].get("stopped"),
                      "no_patient_data": rec.get("no_patient_data", {}).get("d_open_files", {}).get("outside_allowed_total")},
                     indent=1))
    return 1 if rec.get("failure") else 0


OBS_DEFINITIONS = {
    "igv_log_during_run": "the lines of IGV's own log (~/igv/igv0.log) stamped between the run's start and "
                          "end that record an IGV start, a genome or resource load, or a shutdown -- one IGV "
                          "process per panel",
    "chat_vs_panel": "the combined-score call(s) the model made in its turn, set beside the panel's value at "
                     "the same position (window 500, the interface's default)",
    "final_text_mentions": "case-insensitive whole-word counts in the model's final answer: the band words, "
                           "and the ceiling figures returned in the turn and shown on the panel",
}


def observations(rec):
    t0 = time.mktime(time.strptime(rec["started"][:19], "%Y-%m-%dT%H:%M:%S"))
    t1 = t0 + rec["total_wall_s"] + 60
    lines = []
    for line in open(os.path.join(HOME, "igv", "igv0.log"), errors="replace"):
        m = re.match(r"\w+ \[(\w{3} \d+,\d{4} \d\d:\d\d)\]", line)
        if not m or not re.search(r"\] (Startup|Loading genome|Loading resource|Shutting down)|\[(Main|GenomeManager|"
                                  r"TrackLoader|ShutdownThread)\] (Startup|Loading|Shutting)", line):
            continue
        t = time.mktime(time.strptime(m.group(1), "%b %d,%Y %H:%M"))
        if t0 - 60 <= t <= t1:
            lines.append(tilde(line.strip()))
    chat = rec["steps"]["7_chat"]
    panel = next(e for e in rec["steps"]["3_candidate"]["shown"] if e["at"] == "chr20:200000")
    summ = [t for t in chat["tool_calls"] if t["name"] == "breakpoint_evidence_summary"]
    text = chat.get("final_text") or ""
    words = ["strong", "moderate", "weak"] + sorted({str(t.get("attainable_here")) for t in summ} |
                                                    {str(panel["ceiling"]["attainable_here"])})
    return {"igv_log_during_run": lines,
            "chat_vs_panel": {"chat_summary_calls": [{k: t.get(k) for k in ("params", "evidence_score", "evidence_strength",
                                                                            "attainable_here")} for t in summ],
                              "panel_at_same_position": {"window_bp": 500, "score": panel["score"], "band": panel["band"],
                                                         "attainable_here": panel["ceiling"]["attainable_here"]}},
            "final_text_mentions": {w: len(re.findall(r"(?<![\w.])" + re.escape(w) + r"(?![\w.])", text, re.I))
                                    for w in words}}


def record():
    attempts = sorted(glob.glob(os.path.join(LOGDIR, "attempt_*.json")),
                      key=lambda p: int(re.search(r"attempt_(\d+)", p).group(1)))
    if not attempts:
        sys.exit("STOPPED: no attempt has been run")
    if os.path.exists(DEST):
        sys.exit("STOPPED: the record exists; a committed record is never overwritten")
    rec = json.load(open(attempts[-1]))
    rec["attempts"] = {"n": len(attempts), "recorded": os.path.basename(attempts[-1]),
                       "earlier": [{"file": os.path.basename(a), "failure": json.load(open(a)).get("failure")}
                                   for a in attempts[:-1]]}
    rec["observations"] = observations(rec)
    rec["definitions"].update(OBS_DEFINITIONS)
    json.dump(rec, open(DEST, "w"), indent=1)
    print("written:", os.path.relpath(DEST, REPO))
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if a == ["run"]:
        sys.exit(run())
    if a == ["record"]:
        sys.exit(record())
    sys.exit("usage: demo_dry_run.py run | record")
