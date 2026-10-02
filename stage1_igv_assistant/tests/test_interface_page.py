#!/usr/bin/env python3
"""
Regression tests for the interface page (Phase 24) and the routes it added or
changed: the page itself, the privacy guard in front of the cloud chat, the
cloud model list, and the comparison's per-candidate ids.

Nothing here reaches the network. Ollama is replaced by a fixed list, and the
cloud transport (chat.run_turn_api) by a scripted stand-in that calls the
route's OWN executor, so the guard, the label resolution and the recorder are
the real ones; only the model is fake. Every check that could pass vacuously
is first shown to fail on a control.

Run: python3 stage1_igv_assistant/tests/test_interface_page.py
"""
import json
import os
import re
import sys
import tempfile
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant import ui  # noqa: E402
from stage1_igv_assistant import chat as chatmod  # noqa: E402

FAILURES = []
PATH_RE = re.compile(r"(?<![\w.~:/@-])/[\w.+@%~-]+(?:/[\w.+@%~-]*)+")
# Every route the demo dry run of 2026-09-29 required the served page to call
# (scripts/demo_dry_run.py). The previous page (/classic) still calls them all.
DRY_RUN_ROUTES = ["/api/bootstrap", "/api/load", "/api/funnel", "/api/candidate", "/api/assess",
                  "/api/call?id=", "/api/igv", "/api/chat_models", "/api/chat"]
# The page since Phase 26 reviews every candidate at once and shows no score: it
# calls the review routes instead of /api/candidate and /api/assess.
PAGE_ROUTES = ["/api/bootstrap", "/api/load", "/api/funnel", "/api/review", "/api/junction", "/api/position",
               "/api/genes", "/api/call?id=", "/api/igv", "/api/chat_models", "/api/chat"]
EXTERNAL_RE = re.compile(r"""(?:src|href)\s*=\s*["']?(?:https?:)?//|@import|url\(\s*["']?(?:https?:)?//""", re.I)
# Upper-case names the page's script reads (indexed, a member taken, or interpolated),
# and the names it declares (const/let/var, a later name in the same declaration, a
# function). JSON and URL are the browser's.
_JS_READ = re.compile(r"\b([A-Z][A-Z0-9_]{2,})(?=\s*\[|\.\w|\})")
_JS_DECL = re.compile(r"\b(?:const|let|var)\s+([A-Z][A-Z0-9_]{2,})\s*=|,\s*([A-Z][A-Z0-9_]{2,})\s*=(?!=)"
                      r"|\bfunction\s+([A-Z][A-Z0-9_]{2,})\s*\(")
_JS_GLOBALS = {"JSON", "URL", "NaN"}


def undeclared_constants(script):
    declared = {n for t in _JS_DECL.findall(script) for n in t if n}
    return sorted(set(_JS_READ.findall(script)) - declared - _JS_GLOBALS)


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


# ── fixture ─────────────────────────────────────────────────────────────────
def write_bam(path):
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 100000}, {"SN": "chr2", "LN": 100000}]})
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for i, start in enumerate(range(8000, 32000, 10)):
            r = pysam.AlignedSegment(header)
            r.query_name, r.query_sequence = f"r{i}", "ACGT" * 25
            r.query_qualities = pysam.qualitystring_to_array("I" * 100)
            r.reference_id, r.reference_start, r.mapping_quality = 0, start, 60
            r.flag, r.cigar = 0x1 | 0x2, [(0, 100)]
            r.next_reference_id, r.next_reference_start, r.template_length = 0, start + 200, 300
            out.write(r)
    pysam.index(path)


def write_vcf(path, records):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n##contig=<ID=chr2,length=100000>\n"
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
                '##INFO=<ID=END,Number=1,Type=Integer,Description="end">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for rid, pos, end in records:
            f.write(f"chr1\t{pos}\t{rid}\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END={end}\n")


class Server:
    def __init__(self):
        self.srv = ThreadingHTTPServer(("127.0.0.1", 0), ui.Handler)
        threading.Thread(target=self.srv.serve_forever, daemon=True).start()
        self.base = f"http://127.0.0.1:{self.srv.server_address[1]}"
        self.bodies = []

    def raw(self, method, path, body=None):
        data = None if body is None else json.dumps(body).encode()
        req = urllib.request.Request(self.base + path, data=data, method=method,
                                     headers={"Content-Type": "application/json"} if data else {})
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                text = r.read().decode()
        except urllib.error.HTTPError as e:
            text = e.read().decode()
        self.bodies.append((f"{method} {path}", text))
        return text

    def json(self, method, path, body=None):
        try:
            return json.loads(self.raw(method, path, body))
        except ValueError:
            return {}

    def close(self):
        self.srv.shutdown()
        self.srv.server_close()


# ── a scripted cloud model ──────────────────────────────────────────────────
SCRIPT = []          # (tool name, model arguments) the fake model will request, in order
SEEN = {}            # what the fake model was offered and given


def _enums(tools):
    """Every label enum in the tool schemas the model was offered."""
    out = set()

    def walk(x):
        if isinstance(x, dict):
            if isinstance(x.get("enum"), list):
                out.update(v for v in x["enum"] if isinstance(v, str))
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
    walk(tools)
    return out


def fake_run_turn_api(model, message, tools, tool_names, exec_fn, **kw):
    SEEN["calls"] = SEEN.get("calls", 0) + 1
    SEEN["offered"] = _enums(tools)
    events = []
    for name, args in SCRIPT:
        rec, err = exec_fn(name, dict(args))
        events.append({"type": "tool", "name": name, "params": args,
                       "call_id": rec["id"] if rec else None,
                       "result": rec["result"] if rec else {"error": err},
                       "rejected": rec is None})
    return {"model": model, "backend": "anthropic", "events": events, "final_text": "",
            "verification": {"unsupported": [], "details": [], "numbers_in_prose": 0}}


def fake_run_turn(model, message, tools, tool_names, exec_fn, **kw):
    out = fake_run_turn_api(model, message, tools, tool_names, exec_fn, **kw)
    out["backend"] = "ollama"
    return out


def main():
    with tempfile.TemporaryDirectory() as d:
        bam_pub, bam_priv, bam_priv2 = (os.path.join(d, f) for f in ("pub.bam", "priv.bam", "privb.bam"))
        for b in (bam_pub, bam_priv, bam_priv2):
            write_bam(b)
        vcf_pub, vcf_priv, vcf_priv2 = (os.path.join(d, f) for f in ("pub.vcf", "priv.vcf", "privb.vcf"))
        # PUB and PRIV share one junction (within 500 bp) and differ in two; PRIV2 is a
        # second private sample, far from both, for the scope of a permission
        write_vcf(vcf_pub, [("a1", 10000, 10500), ("a2", 20000, 21000), ("a3", 30000, 30600)])
        write_vcf(vcf_priv, [("b1", 10050, 10520), ("b2", 25000, 26000), ("b3", 33000, 34000)])
        write_vcf(vcf_priv2, [("c1", 60000, 61000), ("c2", 70000, 71000)])

        saved = (dict(ui.DATASETS), dict(ui.CANDIDATE_FILES), {k: set(v) for k, v in ui.PUBLIC_LABELS.items()},
                 ui.API_KEY, ui.probe_ollama, chatmod.run_turn_api, chatmod.run_turn, ui.PAGE_FILE)
        ui.DATASETS.clear(); ui.CANDIDATE_FILES.clear()
        ui.DATASETS.update({"PUB": bam_pub, "PRIV": bam_priv, "PRIV2": bam_priv2})
        ui.CANDIDATE_FILES.update({"PUB": vcf_pub, "PRIV": vcf_priv, "PRIV2": vcf_priv2})
        ui.PUBLIC_LABELS["datasets"] = {"PUB"}
        ui.PUBLIC_LABELS["candidates"] = {"PUB"}
        ui._CHAT_TOOLS.clear()
        ui.probe_ollama = lambda timeout=1.5: ["qwen2.5:7b"]
        chatmod.run_turn_api = fake_run_turn_api
        chatmod.run_turn = fake_run_turn
        S = Server()
        try:
            run(S, d)
        finally:
            S.close()
            ui.DATASETS.clear(); ui.DATASETS.update(saved[0])
            ui.CANDIDATE_FILES.clear(); ui.CANDIDATE_FILES.update(saved[1])
            ui.PUBLIC_LABELS.update(saved[2])
            (ui.API_KEY, ui.probe_ollama, chatmod.run_turn_api, chatmod.run_turn, ui.PAGE_FILE) = saved[3:]
            ui._CHAT_TOOLS.clear()

    print("\n" + "=" * 60)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("ALL INTERFACE PAGE TESTS PASSED")


def run(S, d):
    # ── 1. the page ─────────────────────────────────────────────────────────
    print("the page")
    page = S.raw("GET", "/")
    classic = S.raw("GET", "/classic")
    check("/ serves the new page", 'id="view-cands"' in page and 'id="view-assistant"' in page)
    check("/classic serves the previous page", "local instrument" in classic and 'id="view-cands"' not in classic)
    missing = [r for r in PAGE_ROUTES if r not in page]
    check("the new page calls every route of the review view", not missing, str(missing))
    check("the previous page calls every route the demo dry run requires", not [r for r in DRY_RUN_ROUTES if r not in classic])
    check("the new page no longer asks for the combined score", "/api/assess" not in page)
    check("control: the route check fails on a page missing one route",
          [r for r in PAGE_ROUTES if r not in page.replace("/api/igv", "/api/xxx")] == ["/api/igv"])
    check("the new page carries no absolute path", not PATH_RE.findall(page.replace("/api/", "").replace("/img/", "")),
          str(PATH_RE.findall(page)[:3]))
    # The page links genes to OMIM entries (<a href>, opened only when clicked); those
    # anchors are navigation, not something the page loads, so they are set aside here.
    loaded = re.sub(r"<a\b[^>]*>", "", page)
    check("the new page loads nothing from the network (no external script, style or font)",
          not EXTERNAL_RE.search(loaded), str(EXTERNAL_RE.findall(loaded)[:3]))
    check("control: the network detector fires on an external script",
          bool(EXTERNAL_RE.search(page + '<script src="https://cdn.example/x.js"></script>')))
    # A name the script reads but no longer declares throws only when that line runs:
    # the page as patched in Phase 26 had dropped LAYER_TEXT, which the one-line
    # summary of an applicable_layers call still reads, so an answer whose model had
    # called that tool never appeared (the card kept its "working" spinner).
    script = page.split("<script>", 1)[-1]
    check("every constant the page's script reads is declared in it", not undeclared_constants(script),
          str(undeclared_constants(script)))
    check("control: ... the check fires when a declaration is missing",
          "READ_KIND" in undeclared_constants(script.replace("const READ_KIND", "const READ_KIND_GONE"))
          and "READ_KIND" not in undeclared_constants(script))
    ui.PAGE_FILE = os.path.join(d, "no-such-page.html")
    check("if the page file is missing, / falls back to the previous page", "local instrument" in S.raw("GET", "/"))
    ui.PAGE_FILE = os.path.join(os.path.dirname(os.path.abspath(ui.__file__)), "ui_page.html")

    # ── 2. which samples are private ────────────────────────────────────────
    print("\nbootstrap: test data and private data")
    b = S.json("GET", "/api/bootstrap")
    k = b.get("kinds") or {}
    check("an autodiscovered sample is reported as test data",
          k.get("datasets", {}).get("PUB") == "public" and k.get("candidates", {}).get("PUB") == "public", str(k))
    check("an explicitly registered sample is reported as private",
          k.get("datasets", {}).get("PRIV") == "private" and k.get("candidates", {}).get("PRIV") == "private", str(k))
    # a config file can declare some of its explicit labels test data (the demo bundle does)
    import configparser
    from stage1_igv_assistant import config as CFG
    saved_cfg = (CFG._CP, CFG.CONFIG_FILE, dict(ui.DATA_DIR_STATUS), dict(ui.DATASETS), dict(ui.CANDIDATE_FILES),
                 {kk: set(v) for kk, v in ui.PUBLIC_LABELS.items()})
    try:
        conf = os.path.join(d, "sv-assistant.conf")
        with open(conf, "w") as f:
            f.write(f"[datasets]\nDEMO = {ui.DATASETS['PUB']}\nOTHER = {ui.DATASETS['PRIV']}\n"
                    f"[candidates]\nDEMO = {ui.CANDIDATE_FILES['PUB']}\n[test_data]\nlabels = DEMO\n")
        cp = configparser.ConfigParser(); cp.optionxform = str; cp.read(conf)
        CFG._CP, CFG.CONFIG_FILE = cp, conf
        ui.DATA_DIR_STATUS["path"] = None                  # no autodiscovery: config entries only
        ui.DATASETS.clear(); ui.CANDIDATE_FILES.clear()
        ui.PUBLIC_LABELS["datasets"], ui.PUBLIC_LABELS["candidates"] = set(), set()
        ui.discover_public()
        check("a config label listed under [test_data] is test data",
              ui.label_kind("datasets", "DEMO") == "public" and ui.label_kind("candidates", "DEMO") == "public")
        check("a config label not listed there is private", ui.label_kind("datasets", "OTHER") == "private")
        cp.remove_section("test_data")
        ui.DATASETS.clear(); ui.CANDIDATE_FILES.clear()
        ui.PUBLIC_LABELS["datasets"], ui.PUBLIC_LABELS["candidates"] = set(), set()
        ui.discover_public()
        check("control: without [test_data] the same label is private", ui.label_kind("datasets", "DEMO") == "private")
    finally:
        CFG._CP, CFG.CONFIG_FILE = saved_cfg[0], saved_cfg[1]
        ui.DATA_DIR_STATUS.clear(); ui.DATA_DIR_STATUS.update(saved_cfg[2])
        ui.DATASETS.clear(); ui.DATASETS.update(saved_cfg[3])
        ui.CANDIDATE_FILES.clear(); ui.CANDIDATE_FILES.update(saved_cfg[4])
        ui.PUBLIC_LABELS.update(saved_cfg[5])

    # ── 3. cloud models are offered only with a key ─────────────────────────
    print("\ncloud model list")
    ui.API_KEY = None
    m = S.json("GET", "/api/chat_models")
    check("without a key no cloud model is offered", m.get("cloud_models") == [], str(m.get("cloud_models")))
    check("the local models are still listed", m.get("models") == ["qwen2.5:7b"], str(m.get("models")))
    ui.API_KEY = "test-key-not-real"
    m = S.json("GET", "/api/chat_models")
    cm = m.get("cloud_models") or []
    check("with a key the two evaluated Claude models are offered",
          sorted(cm) == sorted(ui.CLOUD_MODELS), str(cm))
    check("every offered cloud model is priced, so each answer can report its cost",
          all(x in chatmod.API_PRICES for x in cm), str(cm))
    check("a priced but unevaluated model is not offered",
          not any(x not in ui.CLOUD_MODELS for x in cm), str(cm))
    ui.API_KEY = None
    SEEN.clear()
    r = S.json("POST", "/api/chat", {"model": "claude-sonnet-5", "message": "In dataset PUB, chr1:10000?"})
    check("a cloud question without a key is refused before any model runs",
          "API key" in (r.get("error") or "") and not SEEN.get("calls"), str(r)[:200])
    ui.API_KEY = "test-key-not-real"

    # ── 4. the privacy guard ────────────────────────────────────────────────
    print("\nprivacy guard: tool arguments")
    pub_set = S.json("POST", "/api/load", {"candidates_label": "PUB"})["result"]["set_id"]
    priv_set = S.json("POST", "/api/load", {"candidates_label": "PRIV"})["result"]["set_id"]
    check("control: the guard names a private dataset label", ui._private_hits({"dataset": "PRIV"}) == ["PRIV"])
    check("control: ... a private label inside a list", ui._private_hits({"datasets": ["PUB", "PRIV"]}) == ["PRIV"])
    check("control: ... a private candidate file", ui._private_hits({"candidates": "PRIV"}) == ["PRIV"])
    check("control: ... the id of a set loaded from a private file", ui._private_hits({"set_id": priv_set}) == ["PRIV"])
    check("control: ... either side of a comparison", ui._private_hits({"set_a": pub_set, "set_b": priv_set}) == ["PRIV"])
    check("silent on test data", ui._private_hits({"dataset": "PUB", "set_id": pub_set, "candidates": "PUB"}) == [])
    # a private file loaded under a public-looking label is still private
    from stage1_igv_assistant.tools import vcf_tools
    disguised = vcf_tools.load_candidate_set(ui.CANDIDATE_FILES["PRIV"], "PUB")["set_id"]
    check("a private file loaded under a public label still counts as private, named by its sample",
          ui._private_hits({"set_id": disguised}) == ["PRIV"] and ui._set_is_private(vcf_tools._SETS[disguised]),
          str(ui._private_hits({"set_id": disguised})))
    check("control: the public file under its own label does not", not ui._set_is_private(vcf_tools._SETS[pub_set]))
    del vcf_tools._SETS[disguised]          # keep the rest of the test about PUB and PRIV only

    n_before = len(ui.RECORDER.calls)
    SCRIPT[:] = [("discordant_pairs", {"dataset": "PRIV", "chromosome": "chr1", "position": 10000}),
                 ("discordant_pairs", {"dataset": "PUB", "chromosome": "chr1", "position": 10000}),
                 ("list_candidates", {"set_id": priv_set})]
    SEEN.clear()
    r = S.json("POST", "/api/chat", {"model": "claude-sonnet-5", "message": "How strong is the evidence?"})
    priv_runs = [c for c in ui.RECORDER.calls[n_before:] if (c.get("params") or {}).get("bam_path") == ui.DATASETS["PRIV"]]
    set_runs = [c for c in ui.RECORDER.calls[n_before:] if (c.get("params") or {}).get("set_id") == priv_set]
    pub_runs = [c for c in ui.RECORDER.calls[n_before:] if (c.get("params") or {}).get("bam_path") == ui.DATASETS["PUB"]]
    blocked = {(x["tool"], tuple(x["labels"])) for x in (r.get("privacy") or {}).get("blocked", [])}
    check("without permission the model is not even offered the private label",
          "PUB" in SEEN.get("offered", set()) and "PRIV" not in SEEN.get("offered", set()), str(SEEN.get("offered")))
    check("without permission a call naming a private dataset is blocked and never runs",
          ("discordant_pairs", ("PRIV",)) in blocked and not priv_runs, str(blocked))
    check("without permission a call passing a private set id is blocked and never runs",
          ("list_candidates", ("PRIV",)) in blocked and not set_runs, str(blocked))
    check("the same question's test-data call runs", len(pub_runs) >= 1)
    check("the answer reports the privacy state", (r.get("privacy") or {}).get("private_allowed") is False
          and (r.get("privacy") or {}).get("cloud") is True, str(r.get("privacy")))
    ev = [e for e in r.get("events", []) if e.get("rejected")]
    check("the model is told why, in words", all("privacy setting" in json.dumps(e.get("result")) for e in ev) and len(ev) == 2)

    n_before = len(ui.RECORDER.calls)
    SEEN.clear()
    r = S.json("POST", "/api/chat", {"model": "claude-sonnet-5", "message": "How strong is the evidence?",
                                     "allow_private_cloud": "true"})
    check("permission must be the literal true: the string 'true' does not count",
          (r.get("privacy") or {}).get("private_allowed") is False and "PRIV" not in SEEN.get("offered", set()))
    # permission covers the private samples the question involves, and no others
    SEEN.clear()
    SCRIPT[:] = [("discordant_pairs", {"dataset": "PRIV", "chromosome": "chr1", "position": 10000}),
                 ("discordant_pairs", {"dataset": "PRIV2", "chromosome": "chr1", "position": 10000})]
    r = S.json("POST", "/api/chat", {"model": "claude-sonnet-5", "allow_private_cloud": True,
                                     "message": "In dataset PRIV, how strong is the evidence?"})
    runs_on = lambda lbl: [c for c in ui.RECORDER.calls[n_before:]
                           if (c.get("params") or {}).get("bam_path") == ui.DATASETS[lbl]]
    blocked = {(x["tool"], tuple(x["labels"])) for x in (r.get("privacy") or {}).get("blocked", [])}
    check("with permission the private sample the question names is offered, and only that one",
          "PRIV" in SEEN.get("offered", set()) and "PRIV2" not in SEEN.get("offered", set()), str(SEEN.get("offered")))
    check("... its call runs", len(runs_on("PRIV")) >= 1)
    check("... a call on another private sample is blocked and never runs",
          ("discordant_pairs", ("PRIV2",)) in blocked and not runs_on("PRIV2"), str(blocked))
    check("... and the answer records what the permission covered",
          (r.get("privacy") or {}).get("permitted") == ["PRIV"], str(r.get("privacy")))
    n_before = len(ui.RECORDER.calls)
    SEEN.clear()
    r = S.json("POST", "/api/chat", {"model": "claude-sonnet-5", "allow_private_cloud": True,
                                     "message": "How strong is the evidence?"})
    blocked = {(x["tool"], tuple(x["labels"])) for x in (r.get("privacy") or {}).get("blocked", [])}
    check("permission given for a question that involves no private sample opens none",
          "PRIV" not in SEEN.get("offered", set()) and ("discordant_pairs", ("PRIV",)) in blocked
          and not runs_on("PRIV"), str(blocked))
    SCRIPT[:] = [("discordant_pairs", {"dataset": "PRIV", "chromosome": "chr1", "position": 10000}),
                 ("discordant_pairs", {"dataset": "PUB", "chromosome": "chr1", "position": 10000}),
                 ("list_candidates", {"set_id": priv_set})]

    n_before = len(ui.RECORDER.calls)
    r = S.json("POST", "/api/chat", {"model": "qwen2.5:7b", "message": "How strong is the evidence?"})
    priv_runs = [c for c in ui.RECORDER.calls[n_before:] if (c.get("params") or {}).get("bam_path") == ui.DATASETS["PRIV"]]
    check("a model on this computer may read private data without the cloud guard",
          len(priv_runs) >= 1 and "privacy" not in r)

    # ── 4b. calls the schemas never offer (found when the page was applied; FIGURE_MAP N) ──
    print("\nprivacy guard: calls the schemas never offer")
    # (i) a file named by its path instead of a label. resolve_args passes such an
    # argument through untouched, so the label check alone never saw it.
    check("control: a label call on test data is no hit", ui._private_hits({"dataset": "PUB", "chromosome": "chr1"}) == [])
    check("a read file named by its path, not a label, counts as private",
          bool(ui._private_hits({"bam_path": ui.DATASETS["PRIV"]}))
          and bool(ui._private_hits({"bam_paths": [ui.DATASETS["PRIV"]]})))
    check("... and so does a candidate file named by its path", bool(ui._private_hits({"path": ui.CANDIDATE_FILES["PRIV"]})))
    n_before = len(ui.RECORDER.calls)
    SCRIPT[:] = [("discordant_pairs", {"bam_path": ui.DATASETS["PRIV"], "chromosome": "chr1", "position": 10000}),
                 ("load_candidate_set", {"path": ui.CANDIDATE_FILES["PRIV"], "label": "mine"})]
    r = S.json("POST", "/api/chat", {"model": "claude-sonnet-5", "message": "How strong is the evidence?"})
    ran = [c for c in ui.RECORDER.calls[n_before:]
           if ui.DATASETS["PRIV"] in json.dumps(c.get("params")) or ui.CANDIDATE_FILES["PRIV"] in json.dumps(c.get("params"))]
    check("without permission a call passing a private file's path is blocked and never runs",
          len((r.get("privacy") or {}).get("blocked", [])) == 2 and not ran,
          f"blocked={(r.get('privacy') or {}).get('blocked')} ran={len(ran)}")

    # (ii) an unknown label: the error must not tell the cloud model which private labels exist
    SCRIPT[:] = [("discordant_pairs", {"dataset": "NOPE", "chromosome": "chr1", "position": 10000}),
                 ("load_candidate_set", {"candidates": "NOPE", "label": "x"})]
    r = S.json("POST", "/api/chat", {"model": "claude-sonnet-5", "message": "How strong is the evidence?"})
    told = json.dumps([e.get("result") for e in r.get("events", [])])
    check("without permission an unknown-label error lists the test-data labels only",
          told.count("unknown") == 2 and "PUB" in told and "PRIV" not in told, told[:240])
    r = S.json("POST", "/api/chat", {"model": "qwen2.5:7b", "message": "How strong is the evidence?"})
    told = json.dumps([e.get("result") for e in r.get("events", [])])
    check("control: for a model on this computer the same error lists every label", "PRIV" in told, told[:240])

    # (iii) the label of a set is the model's own choice (in the recorded runs it never
    # equalled the registered label), so the file decides whether a set is private
    SCRIPT[:] = [("load_candidate_set", {"candidates": "PUB", "label": "PUB_calls"})]
    r = S.json("POST", "/api/chat", {"model": "claude-sonnet-5", "message": "List the candidates."})
    own = ((r.get("events") or [{}])[0].get("result") or {}).get("set_id")
    SCRIPT[:] = [("list_candidates", {"set_id": own})]
    r = S.json("POST", "/api/chat", {"model": "claude-sonnet-5", "message": "List the candidates."})
    listed = (r.get("events") or [{}])[0].get("result") or {}
    check("a test-data file loaded under a label of the model's own choosing stays usable without permission",
          own is not None and not (r.get("privacy") or {}).get("blocked") and listed.get("total_in_set") == 3,
          f"blocked={(r.get('privacy') or {}).get('blocked')} {str(listed)[:120]}")
    probe = S.json("POST", "/api/privacy_check", {"message": "Which gene lies at chr1:20000?"})
    check("... and its positions are not treated as private", probe == {"labels": [], "positions": []}, str(probe))
    mine = vcf_tools.load_candidate_set(ui.CANDIDATE_FILES["PRIV"], "PRIV_calls")["set_id"]
    check("control: a private file under a label of the model's own choosing is still private, named by its sample",
          ui._private_hits({"set_id": mine}) == ["PRIV"], str(ui._private_hits({"set_id": mine})))
    alias = vcf_tools.load_candidate_set(ui.CANDIDATE_FILES["PUB"], "PRIV")["set_id"]
    check("a set loaded under a registered private label is private whatever its file",
          ui._private_hits({"set_id": alias}) == ["PRIV"], str(ui._private_hits({"set_id": alias})))
    for sid in [k for k, s in vcf_tools._SETS.items() if s.get("label") in ("mine", "x", "PUB_calls", "PRIV_calls")] + [alias]:
        vcf_tools._SETS.pop(sid, None)   # keep the rest of the test about PUB and PRIV only

    # ── 5. the question text ────────────────────────────────────────────────
    print("\nprivacy guard: the question itself")
    SCRIPT[:] = []
    probe = S.json("POST", "/api/privacy_check", {"message": "In dataset PRIV, how strong is chr1:5,000?"})
    check("a question naming a private sample is flagged", probe.get("labels") == ["PRIV"], str(probe))
    probe = S.json("POST", "/api/privacy_check", {"message": "Is the evidence strong in PRIV."})
    check("... also when a full stop follows the name", probe.get("labels") == ["PRIV"], str(probe))
    probe = S.json("POST", "/api/privacy_check", {"message": "In dataset PRIVATE2, chr1:5000?"})
    check("a longer word that merely starts with the name is not flagged", probe.get("labels") == [], str(probe))
    probe = S.json("POST", "/api/privacy_check", {"message": "Which gene lies at chr1:33,000?"})
    check("a position from a private candidate list is flagged with its sample",
          probe.get("positions") == [{"written": "33,000", "labels": ["PRIV"]}], str(probe))
    for q, written in (("chr1:33 000", "33 000"), ("chr1:33.000", "33.000"), ("chr1:33001", "33001"),
                       ("chr1:32,400", "32,400"), ("0.033 Mb", "0.033 Mb"), ("33 kb", "33 kb"), ("33,4 kb", "33,4 kb")):
        probe = S.json("POST", "/api/privacy_check", {"message": f"Which gene lies at {q}?"})
        check(f"... also written as {q!r} (any digit grouping, within {ui.NEAR_BP} bp, or inside a Mb/kb rounding)",
              probe.get("positions") == [{"written": written, "labels": ["PRIV"]}], str(probe))
    for q in ("chr1:36,500", "0.04 Mb", "within 500 bp"):
        probe = S.json("POST", "/api/privacy_check", {"message": f"What is at {q}?"})
        check(f"control: {q!r} is not near any private position", probe.get("positions") == [], str(probe))
    probe = S.json("POST", "/api/privacy_check", {"message": "Which gene lies at chr1:20000?"})
    check("a position only in test data is not flagged", probe == {"labels": [], "positions": []}, str(probe))
    SEEN.clear()
    r = S.json("POST", "/api/chat", {"model": "claude-sonnet-5", "message": "Which gene lies at chr1:33000?"})
    check("without permission a question carrying a private position is not sent at all",
          (r.get("error") or "").startswith("Not sent") and not SEEN.get("calls"), str(r)[:200])
    r = S.json("POST", "/api/chat", {"model": "claude-sonnet-5", "message": "Which gene lies at chr1:33000?",
                                     "allow_private_cloud": True})
    check("with permission it is sent", SEEN.get("calls") == 1 and not r.get("error"), str(r)[:200])
    # a position read from a private READ file, not a candidate list
    S.json("POST", "/api/assess", {"bam_label": "PRIV", "chromosome": "chr1", "position": 18888})
    probe = S.json("POST", "/api/privacy_check", {"message": "what about chr1:18,888"})
    check("a position checked by hand in a private read file is flagged",
          probe.get("positions") == [{"written": "18,888", "labels": ["PRIV"]}], str(probe))

    # a file named by its path is refused for every model, not only for a cloud one
    print("\nfile paths")
    n_before = len(ui.RECORDER.calls)
    SCRIPT[:] = [("discordant_pairs", {"bam_path": ui.DATASETS["PUB"], "chromosome": "chr1", "position": 10000})]
    r = S.json("POST", "/api/chat", {"model": "qwen2.5:7b", "message": "How strong is the evidence?"})
    ev = (r.get("events") or [{}])[0]
    check("a model on this computer that names a file by its path is refused, and nothing runs",
          ev.get("rejected") and "file path is not accepted" in json.dumps(ev.get("result"))
          and len(ui.RECORDER.calls) == n_before, str(ev)[:200])
    SCRIPT[:] = [("discordant_pairs", {"dataset": "PUB", "chromosome": "chr1", "position": 10000})]
    r = S.json("POST", "/api/chat", {"model": "qwen2.5:7b", "message": "How strong is the evidence?"})
    check("control: the same call by label runs", len(ui.RECORDER.calls) > n_before)
    # mask_path is a path argument too: the schemas offer exclude_masked in its place,
    # but resolve_args passes a mask file the model names straight to the tool (found
    # when the Phase 25 patch was applied; FIGURE_MAP O)
    for model, where in (("qwen2.5:7b", "on this computer"), ("claude-sonnet-5", "in the cloud")):
        n_before = len(ui.RECORDER.calls)
        SCRIPT[:] = [("list_candidates", {"set_id": pub_set, "mask_path": ui.CANDIDATE_FILES["PRIV"]})]
        r = S.json("POST", "/api/chat", {"model": model, "message": "List the candidates."})
        ev = (r.get("events") or [{}])[0]
        check(f"a model {where} that names a mask file by its path is refused, and nothing runs",
              ev.get("rejected") and "file path is not accepted" in json.dumps(ev.get("result"))
              and len(ui.RECORDER.calls) == n_before, str(ev)[:200])
    SCRIPT[:] = [("list_candidates", {"set_id": pub_set, "exclude_masked": True})]
    n_before = len(ui.RECORDER.calls)
    r = S.json("POST", "/api/chat", {"model": "qwen2.5:7b", "message": "List the candidates."})
    check("control: asking for the exclude regions with exclude_masked runs",
          len(ui.RECORDER.calls) > n_before and not (r.get("events") or [{}])[0].get("rejected"))
    SCRIPT[:] = []

    # ── 5b. a private sample's list: counts only until the person asks ──────
    print("\ncounts only")
    full = S.json("POST", "/api/funnel", {"set_id": priv_set, "limit": 100})
    counts = S.json("POST", "/api/funnel", {"set_id": priv_set, "limit": 100, "counts_only": True})
    fr, cr = full.get("result") or {}, counts.get("result") or {}
    check("control: the full funnel lists the private sample's candidates", len(fr.get("candidates") or []) == 3)
    check("counts only returns no candidate, and every count unchanged",
          cr.get("candidates") == [] and cr.get("total_matching") == fr.get("total_matching") == 3
          and [x["surviving_after_this_step"] for x in cr.get("filters_applied", [])]
          == [x["surviving_after_this_step"] for x in fr.get("filters_applied", [])], str(cr)[:200])
    rec = next(c for c in ui.RECORDER.calls if c["id"] == counts.get("call"))
    check("... and the recorded call holds no position either",
          not (rec.get("result") or {}).get("candidates") and "pos1" not in json.dumps(rec.get("result")))
    check("counts only must be the literal true",
          len(((S.json("POST", "/api/funnel", {"set_id": priv_set, "limit": 100, "counts_only": "yes"})
                .get("result") or {}).get("candidates")) or []) == 3)

    # ── 6. comparison ids ───────────────────────────────────────────────────
    print("\ncomparison ids")
    c = S.json("POST", "/api/compare", {"label_a": "PUB", "label_b": "PRIV"})
    a = (c.get("survivor_effect") or {}).get("a") or {}
    fn = S.json("POST", "/api/funnel", {"set_id": pub_set, "limit": 100})
    survivors = {x["candidate_id"] for x in (fn.get("result") or {}).get("candidates", [])}
    rec, uni = set(a.get("recurrent_ids") or []), set(a.get("unique_ids") or [])
    check("recurrent and unique ids together are exactly the survivors", rec | uni == survivors and len(survivors) == 3,
          f"{sorted(rec)} {sorted(uni)} {sorted(survivors)}")
    check("no candidate is both", not rec & uni)
    check("the counts match the ids", a.get("recurrent") == len(rec) == 1 and a.get("unique") == len(uni) == 2, str(a)[:200])
    junction_at = {j.candidate_id: (j.pos1, j.pos2) for j in vcf_tools._SETS[pub_set]["junctions"]}
    check("the recurrent one is the junction both files share", [junction_at[i] for i in rec] == [(10500, 10000)],
          str([junction_at[i] for i in rec]))

    # ── 6a. the comparison's reply to the new page carries no position ──────
    # Found when the Phase 25 patch was applied (FIGURE_MAP O): with a comparison chosen,
    # the reply listed up to ten example junctions per sample, so the positions of a
    # private sample whose list was still hidden reached the page. The new page marks
    # survivors by id and never showed the examples; it now asks for none.
    print("\ncomparison examples")
    by_default = S.json("POST", "/api/compare", {"label_a": "PRIV", "label_b": "PUB"})
    for_page = S.json("POST", "/api/compare", {"label_a": "PRIV", "label_b": "PUB", "examples": False})
    da, pa = ((x.get("survivor_effect") or {}).get("a") or {} for x in (by_default, for_page))
    check("control: by default the reply lists example junctions with their positions (the previous page shows them)",
          len(da.get("unique_examples") or []) == 2 and '"pos1"' in json.dumps(by_default), str(da)[:200])
    check("asked for no examples, the reply carries no position of either sample",
          for_page.get("survivor_effect") and not re.search(r'"pos[12]?"|"position"', json.dumps(for_page)),
          str(re.findall(r'"pos[12]?": \d+', json.dumps(for_page))[:4]))
    same = lambda x: {k: v for k, v in x.items() if k not in ("unique_examples", "call")}
    check("... and the same counts and ids", bool(pa) and same(pa) == same(da) and pa.get("unique") == 2, str(pa)[:200])
    check("examples are left out only for the literal false",
          len((((S.json("POST", "/api/compare", {"label_a": "PRIV", "label_b": "PUB", "examples": "no"})
                 .get("survivor_effect") or {}).get("a") or {}).get("unique_examples")) or []) == 2)
    check("the new page asks the comparison for no examples",
          bool(re.search(r"post\('/api/compare',\s*\{[^}]*examples:\s*false", S.raw("GET", "/"))))

    # ── 6b. the startup self-test leaves nothing behind ─────────────────────
    # ui.main() runs verify_minimal() before serving. Its fixture candidate set stayed
    # in the registry under a label that is not test data, so its positions (5000,
    # 10000, 20000, 21000) were flagged as private data in every cloud question.
    print("\nstartup self-test")
    before = set(vcf_tools._SETS)
    ui.verify_minimal()
    left = sorted(str(s.get("label")) for k, s in vcf_tools._SETS.items() if k not in before)
    check("the self-test's candidate set is removed when the self-test ends", not left, str(left))
    probe = S.json("POST", "/api/privacy_check", {"message": "Which gene lies at chr1:5000?"})
    check("so a question containing one of its positions is not flagged as private data",
          probe == {"labels": [], "positions": []}, str(probe))

    # ── 7. nothing leaks a path ─────────────────────────────────────────────
    print("\nevery response body")
    needles = [d, os.path.expanduser("~")]
    leaked = []
    for name, text in S.bodies:
        if name in ("GET /", "GET /classic"):
            continue
        found = {m.group(0) for m in PATH_RE.finditer(text) if not m.group(0).startswith(("/api/", "/img/"))}
        found |= {n for n in needles if n and n in text}
        if found:
            leaked.append((name, sorted(found)[:2]))
    check(f"{len(S.bodies)} responses, none carrying an absolute path", not leaked, str(leaked[:3]))


if __name__ == "__main__":
    main()
