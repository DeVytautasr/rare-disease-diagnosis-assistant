#!/usr/bin/env python3
"""Phase 13 Task 0: one IGV evidence panel on public data, through the MCP interface.

    igv_panel_demo.py        -> stage1_igv_assistant/results/igv_panel_demo_2026-09-27.json

Calls the evidence server's evidence_panel tool through FastMCP's in-memory client
(the same protocol path a model uses), at the IMP01 breakpoint chr20:200000 of the
public synthetic control (~/public_data/sim/bams/IMP01.bam). It refuses any BAM
outside ~/public_data: no patient locus is ever imaged.

Requirements, from the environment: IGV 2.17.4 (scripts/install_igv.sh), Java 21 on
PATH, a display (DISPLAY=:0 on WSLg). Images are written under
IGV_IMAGE_SESSION_DIR, outside the repository, and are not committed.

RECORDED (definitions): per layer the tool's success flag, image_ref and
image_dimensions as returned to the caller; whether the return carries any file path
(it must not); and, resolved on disk from the session directory, that a PNG file of
those dimensions exists (size and sha256 only -- the picture itself is not
interpreted).
"""
import asyncio
import glob
import hashlib
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
BAM = os.path.expanduser("~/public_data/sim/bams/IMP01.bam")
SESSION = os.path.expanduser("~/public_data/sim/logs/phase13_2026-09-27/igv_panel_demo")
RECORD = os.path.join(REPO, "stage1_igv_assistant", "results", "igv_panel_demo_2026-09-27.json")


def main():
    pub = os.path.realpath(os.path.expanduser("~/public_data"))
    if not os.path.realpath(BAM).startswith(pub + os.sep):
        sys.exit("refusing: the BAM is not public data")
    os.makedirs(SESSION, exist_ok=True)
    os.environ["IGV_IMAGE_SESSION_DIR"] = SESSION
    from fastmcp import Client
    from stage1_igv_assistant import server

    async def call():
        async with Client(server.mcp) as c:
            r = await c.call_tool("evidence_panel", {"bam_paths": [BAM], "chromosome": "chr20", "position": 200000})
            return r.structured_content if getattr(r, "structured_content", None) is not None else r.data

    res = asyncio.run(call())
    text = json.dumps(res, default=str)
    panels = res.get("panels") or {}
    pngs = {}
    for p in sorted(glob.glob(os.path.join(SESSION, "**", "*.png"), recursive=True)):
        b = open(p, "rb").read()
        pngs[os.path.basename(p)] = {"bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()[:16],
                                     "png_header": b[:8] == b"\x89PNG\r\n\x1a\n",
                                     "dimensions": f"{int.from_bytes(b[16:20], 'big')}x{int.from_bytes(b[20:24], 'big')}"}
    rec = {"what": "One IGV evidence panel through the MCP interface on public data (Phase 13 Task 0, 2026-09-27)",
           "definitions": __doc__.split("RECORDED (definitions):")[1].strip().replace("\n", " "),
           "locus": "IMP01 chr20:200000 (public synthetic control)",
           "environment": {"igv": "2.17.4 at ~/IGV_2.17.4", "java": "conda-forge openjdk 21 (env igv-java)",
                           "display": os.environ.get("DISPLAY")},
           "panels": {k: {f: v.get(f) for f in ("success", "image_ref", "image_dimensions", "error", "region",
                                                  "color_by")} for k, v in panels.items()},
           "return_contains_a_path": ("/home/" in text or ".png" in text or SESSION in text),
           "files_in_session_dir": pngs,
           "code": "scripts/igv_panel_demo.py"}
    rec["end_to_end"] = (bool(panels) and all(v.get("success") and v.get("image_ref") for v in panels.values())
                         and not rec["return_contains_a_path"] and len(pngs) == len(panels)
                         and all(x["png_header"] for x in pngs.values()))
    if os.path.exists(RECORD):
        sys.exit("the record exists; a committed record is never overwritten")
    json.dump(rec, open(RECORD, "w"), indent=1)
    print(json.dumps({k: rec[k] for k in ("panels", "return_contains_a_path", "files_in_session_dir", "end_to_end")},
                     indent=1))
    return 0 if rec["end_to_end"] else 1


if __name__ == "__main__":
    sys.exit(main())
