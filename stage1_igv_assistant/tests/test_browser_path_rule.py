#!/usr/bin/env python3
"""
Regression test: no absolute filesystem path reaches any response the browser gets.

Named for the condition that exposed it. /api/funnel sent the exclude
template's absolute path to the browser, in mask.path and inside mask.reason,
against the interface's own rule. A scan of every route then found more:
IGV's search locations in /api/igv and in the recorded evidence_panel result,
the recorded mask_path parameter in /api/calls and /api/call, and raw
tracebacks from do_POST's error branch. The fix moved scrubbing to the one
place every JSON body leaves (Handler._send) and made it reduce any absolute
path to its basename, so no route has to be judged case by case.

The real Handler is served on a free local port and every route is hit on its
success and error branches. /api/chat is deliberately not called: on a machine
with the API backend installed it could reach a paid API. The detector is shown
to fire on the raw, unscrubbed values before its silence is trusted.

Run: python3 stage1_igv_assistant/tests/test_browser_path_rule.py
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

FAILURES = []
PATH_RE = re.compile(r"(?<![\w.~:/@-])/[\w.+@%~-]+(?:/[\w.+@%~-]*)+")


def fixture(d):
    """Paired reads over chr1:8,000-12,000 and a two-record candidate VCF."""
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 100000}, {"SN": "chr2", "LN": 100000}]})
    bam = os.path.join(d, "fixture.bam")
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for i, start in enumerate(range(8000, 12000, 10)):
            r = pysam.AlignedSegment(header)
            r.query_name, r.query_sequence = f"r{i}", "ACGT" * 25
            r.query_qualities = pysam.qualitystring_to_array("I" * 100)
            r.reference_id, r.reference_start, r.mapping_quality = 0, start, 60
            r.flag, r.cigar = 0x1 | 0x2, [(0, 100)]
            r.next_reference_id, r.next_reference_start, r.template_length = 0, start + 200, 300
            out.write(r)
    pysam.index(bam)
    vcf = os.path.join(d, "fixture.vcf")
    with open(vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n##contig=<ID=chr2,length=100000>\n"
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
                '##INFO=<ID=END,Number=1,Type=Integer,Description="end">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t10000\td1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=10500\n"
                "chr1\t20000\td2\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=21000\n")
    return bam, vcf


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def leaks(text, needles):
    found = {m.group(0) for m in PATH_RE.finditer(text) if not m.group(0).startswith(("/api/", "/img/"))}
    return sorted(found | {n for n in needles if n and n in text})


def main():
    with tempfile.TemporaryDirectory() as d:
        bam, vcf = fixture(d)
        home = os.path.join(d, "home")                   # IGV's built-in search paths live here
        os.makedirs(home)
        excl = os.path.join(d, "ref", "human.hg38.excl.tsv")
        needles = [d, os.path.expanduser("~")]
        ui.DATASETS["FIX"] = bam
        ui.CANDIDATE_FILES["FIXC"] = vcf
        ui.CANDIDATE_FILES["FIXC2"] = vcf
        missing = {"key": "exclude_template", "path": excl, "found": False,
                   "source": f"config file {os.path.join(d, 'sv-assistant.conf')}",
                   "reason": f"no such file: {excl}"}

        print("detector controls")
        check("fires on the raw, unscrubbed mask state", leaks(json.dumps(missing), needles))
        check("fires on a traceback path", leaks('File "/srv/app/ui.py", line 3', needles))
        check("silent on routes, URLs, ratios and labels", not leaks(
            '"/api/load" "http://127.0.0.1:11434/api/tags" "15/28" "<FIX>"', needles))

        srv = ThreadingHTTPServer(("127.0.0.1", 0), ui.Handler)
        threading.Thread(target=srv.serve_forever, daemon=True).start()
        base = f"http://127.0.0.1:{srv.server_address[1]}"
        bodies = []

        def hit(method, path, body=None):
            data = None if body is None else json.dumps(body).encode()
            req = urllib.request.Request(base + path, data=data, method=method,
                                         headers={"Content-Type": "application/json"} if data else {})
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    text = r.read().decode()
            except urllib.error.HTTPError as e:
                text = e.read().decode()
            bodies.append((f"{method} {path}", text))
            try:
                return json.loads(text)
            except ValueError:
                return {}

        saved = (ui.MASK_STATUS, ui.MASK_PATH, os.environ.get("HOME"), os.environ.get("IGV_PATH"))
        try:
            os.environ["HOME"] = home
            os.environ.pop("IGV_PATH", None)
            for state in ("missing", "present"):
                if state == "missing":
                    ui.MASK_STATUS, ui.MASK_PATH = missing, excl
                else:
                    os.makedirs(os.path.dirname(excl), exist_ok=True)
                    with open(excl, "w") as f:
                        f.write("chr1\t0\t9000\tmask\n")
                    ui.MASK_STATUS = dict(missing, found=True, reason=None, source="$SV_EXCLUDE_TEMPLATE")
                    ui.MASK_PATH = excl
                hit("GET", "/api/bootstrap")
                set_id = (hit("POST", "/api/load", {"candidates_label": "FIXC"}).get("result") or {}).get("set_id")
                hit("POST", "/api/load", {"candidates_label": "NOPE"})                 # error branch
                fn = hit("POST", "/api/funnel", {"set_id": set_id, "use_mask": True, "limit": 5})
                mask = fn.get("mask") or {}
                check(f"[{state}] mask.path is the file name only", mask.get("path") == "human.hg38.excl.tsv",
                      str(mask.get("path")))
                if state == "missing":
                    check("[missing] mask.reason keeps its explanation, with the file name",
                          "exclude template not found" in (mask.get("reason") or "")
                          and "human.hg38.excl.tsv" in mask["reason"], str(mask.get("reason")))
                cid = ((fn.get("result") or {}).get("candidates") or [{}])[0].get("candidate_id")
                hit("POST", "/api/candidate", {"set_id": set_id, "candidate_id": cid})
                hit("POST", "/api/assess", {"bam_label": "FIX", "chromosome": "chr1", "position": 10000})
                hit("POST", "/api/assess", {"bam_label": "FIX", "chromosome": "chrNOPE", "position": 5})
                hit("POST", "/api/assess", {"bam_label": "NOPE", "chromosome": "chr1", "position": 5})
                hit("POST", "/api/igv", {"bam_label": "FIX", "chromosome": "chr1", "position": 10000})
                hit("POST", "/api/compare", {"label_a": "FIXC", "label_b": "FIXC2", "use_mask": True})
                hit("POST", "/api/nope", {})
                hit("GET", "/img/no-such-handle")
                for c in hit("GET", "/api/calls").get("calls", []):
                    hit("GET", f"/api/call?id={c['id']}")
        finally:
            ui.MASK_STATUS, ui.MASK_PATH = saved[0], saved[1]
            for k, v in (("HOME", saved[2]), ("IGV_PATH", saved[3])):
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
            srv.shutdown()
            srv.server_close()

        print("\nevery response body")
        found = [(name, leaks(text, needles)) for name, text in bodies]
        for name, got in found:
            if got:
                check(f"no absolute path in {name}", False, str(got[:3]))
        check(f"{len(bodies)} responses, none carrying an absolute path", not any(g for _, g in found))

    print("\n" + "=" * 60)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("ALL BROWSER PATH-RULE TESTS PASSED")


if __name__ == "__main__":
    main()
