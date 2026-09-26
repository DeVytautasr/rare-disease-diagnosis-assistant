#!/usr/bin/env python3
"""Phase 13 Task 0 follow-up: which fields of the MCP evidence_panel return carry a
path. igv_panel_demo.py flagged its return as containing one, by a coarse test
("/home/", ".png" or the session directory anywhere in the text). This repeats the
same call (public IMP01 chr20:200000; refuses non-public BAMs) and lists every
string field whose value contains "/home/", "~/", ".png" or the image session
directory, with the kind of path -- the input BAM echoed back, or an image file --
and never the image. -> stage1_igv_assistant/results/igv_panel_paths_2026-09-27.json
"""
import asyncio
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from igv_panel_demo import BAM, SESSION  # noqa: E402

RECORD = os.path.join(REPO, "stage1_igv_assistant", "results", "igv_panel_paths_2026-09-27.json")


def walk(o, path=""):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from walk(v, f"{path}.{k}")
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from walk(v, f"{path}[{i}]")
    elif isinstance(o, str):
        yield path, o


def main():
    pub = os.path.realpath(os.path.expanduser("~/public_data"))
    if not os.path.realpath(BAM).startswith(pub + os.sep):
        sys.exit("refusing: the BAM is not public data")
    os.environ["IGV_IMAGE_SESSION_DIR"] = SESSION + "_paths"
    from fastmcp import Client
    from stage1_igv_assistant import server

    async def call():
        async with Client(server.mcp) as c:
            r = await c.call_tool("evidence_panel", {"bam_paths": [BAM], "chromosome": "chr20", "position": 200000})
            return r.structured_content if getattr(r, "structured_content", None) is not None else r.data

    res = asyncio.run(call())
    hits = []
    for key, v in walk(res):
        if any(t in v for t in ("/home/", "~/", ".png", SESSION)):
            kind = ("input BAM echoed" if os.path.basename(BAM) in v and ".png" not in v
                    else "image file" if ".png" in v or SESSION in v else "other")
            hits.append({"field": key, "kind": kind, "value": v.replace(os.path.expanduser("~"), "~")})
    rec = {"what": __doc__.splitlines()[0].strip(), "definitions": " ".join(__doc__.split()),
           "path_bearing_fields": hits,
           "image_paths_returned": sum(h["kind"] == "image file" for h in hits),
           "code": "scripts/igv_panel_paths.py"}
    if os.path.exists(RECORD):
        sys.exit("the record exists")
    json.dump(rec, open(RECORD, "w"), indent=1)
    print(json.dumps(rec, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
