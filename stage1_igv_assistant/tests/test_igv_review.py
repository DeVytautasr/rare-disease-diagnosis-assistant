#!/usr/bin/env python3
"""
Regression tests for the page's IGV image (Phase 27, tools/igv_review.py): the
supervisor asked for both ends side by side, only the abnormal reads, coloured by
what they say, normal reads hidden, ambiguous reads kept. The page draws that; these
tests hold IGV to the same reads, colours and words.

Sections 1-5 run on the synthetic reads of tests/review_fixture.py, with IGV
replaced by a stand-in that writes a PNG, so they need neither IGV nor a display.
Section 6 runs IGV itself on the public synthetic sample IMP01
(~/public_data/sim/bams/IMP01.bam) and checks the pixels of the image; without
IGV, a display or that file it prints NOT RUN and the suite exits 2 (incomplete).
Every check that could pass vacuously is first shown to fail on a control.

Run: python3 stage1_igv_assistant/tests/test_igv_review.py
"""
import asyncio
import os
import re
import struct
import sys
import tempfile
import zlib

import pysam

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from stage1_igv_assistant.tests import review_fixture as F  # noqa: E402
from stage1_igv_assistant.tools import gene_table as gt  # noqa: E402
from stage1_igv_assistant.tools import igv_review as IR  # noqa: E402
from stage1_igv_assistant.tools import junction_tools as jt  # noqa: E402
from stage1_igv_assistant import review_server as R  # noqa: E402

FAILURES, NOT_RUN = [], []
PAGE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ui_page.html")
# R1 of the fixture, as the page lists it: J1 then J2 (tests/review_fixture.py)
R1 = [{"chromosome_1": "chr19", "position_1": 48000000, "chromosome_2": "chr22", "position_2": 30000000,
       "orientation": "3to5"},
      {"chromosome_1": "chr19", "position_1": 48000201, "chromosome_2": "chr22", "position_2": 29999999,
       "orientation": "5to3"}]


def _png(w, h, rgb):
    """A valid w x h PNG of one colour, standard library only."""
    def chunk(kind, data):
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data))
    raw = b"".join(b"\x00" + bytes(rgb) * w for _ in range(h))
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw)) + chunk(b"IEND", b""))


PNG_1x1 = _png(1, 1, (255, 255, 255))      # what the stand-in for IGV writes


def check(label, condition, detail=""):
    ok = bool(condition)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + ("" if ok or not detail else f"  -- {detail}"))
    if not ok:
        FAILURES.append(label)
    return ok


def key(x):
    return (x["name"], x["start"], bool(x["supplementary"]))


def page_colours(text):
    """JOIN_COL, JOIN_UNCALLED and READ_KIND as the page's script declares them."""
    jc = re.search(r"const JOIN_COL = \[([^\]]*)\], JOIN_UNCALLED = '(#[0-9A-Fa-f]{6})'", text)
    kinds = dict((k, (c, l)) for k, c, l in
                 re.findall(r"(\w+):\s*\{col:'(#[0-9A-Fa-f]{6})', label:'([^']*)'\}", text))
    return ([c.strip().strip("'") for c in jc.group(1).split(",")] if jc else None,
            jc.group(2) if jc else None, kinds)


def png_pixels(path):
    """(width, height, Counter of RGB) for an 8-bit RGB or RGBA PNG, standard library only."""
    from collections import Counter
    d = open(path, "rb").read()
    i, idat, w = 8, b"", None
    while i < len(d):
        n, = struct.unpack(">I", d[i:i + 4])
        t, c = d[i + 4:i + 8], d[i + 8:i + 8 + n]
        i += 12 + n
        if t == b"IHDR":
            w, h, _, ct = struct.unpack(">IIBB", c[:10])
        elif t == b"IDAT":
            idat += c
    bpp = {2: 3, 6: 4}[ct]
    raw, stride = zlib.decompress(idat), w * bpp
    prev, p, cnt = bytearray(stride), 0, Counter()
    for _ in range(h):
        f, line = raw[p], bytearray(raw[p + 1:p + 1 + stride])
        p += 1 + stride
        for x in range(stride):
            a = line[x - bpp] if x >= bpp else 0
            b = prev[x]
            cc = prev[x - bpp] if x >= bpp else 0
            if f == 1:
                line[x] = (line[x] + a) & 255
            elif f == 2:
                line[x] = (line[x] + b) & 255
            elif f == 3:
                line[x] = (line[x] + (a + b) // 2) & 255
            elif f == 4:
                pa, pb, pc = abs(b - cc), abs(a - cc), abs(a + b - 2 * cc)
                line[x] = (line[x] + (a if pa <= pb and pa <= pc else b if pb <= pc else cc)) & 255
        cnt.update(tuple(line[x:x + 3]) for x in range(0, stride, bpp))
        prev = line
    return w, h, cnt


def near(cnt, rgb, tol=4):
    return sum(n for c, n in cnt.items() if all(abs(c[k] - rgb[k]) <= tol for k in range(3)))


def as_igv_draws(hex_colour, alpha=0.75, bg=250):
    """IGV paints a read at about 75% opacity over its (250, 250, 250) track background
    (measured on the first Phase 27 images: #0F766E came out as (74, 151, 145))."""
    h = hex_colour.lstrip("#")
    return tuple(round(bg - (bg - int(h[i:i + 2], 16)) * alpha) for i in (0, 2, 4))


def fake_igv(written):
    """A stand-in for run_batch: records the batch script and writes a PNG."""
    def run(igv_path, batch_path, png_path, timeout_sec=None):
        written["batch"] = open(batch_path).read()
        written["dir"] = os.path.dirname(batch_path)
        written["files"] = sorted(os.listdir(written["dir"]))
        with open(png_path, "wb") as fh:
            fh.write(PNG_1x1)
        return {"success": True, "shutdown": "clean_exit", "igv_severe_lines": []}
    return run


def run(d):
    # images go to this run's temporary directory, not the interface's session directory
    os.environ[IR._bt.IMAGE_SESSION_DIR_ENV] = os.path.join(d, "sessions")
    fx = F.build_all(d)
    bam = fx["bam"]
    gt.load(*fx["genes"])

    # ── 1. the reads are the page's ─────────────────────────────────────────
    print("the reads")
    page_view = jt.junction_support(bam, "chr19", 48000000, "chr22", 30000000, orientation="3to5",
                                    view=True, related=[("chr19", 48000201), ("chr22", 29999999)])["view"]
    page_keys = {key(x) for e in ("a", "b") for x in page_view[e]["reads"]}
    with pysam.AlignmentFile(bam) as B:
        sel = IR.select_reads(B, R1, also_at=["chr19:48000201", "chr22:29999999"])
        check("the reads chosen are exactly the reads the page draws", set(sel["reads"]) == page_keys,
              f"{len(sel['reads'])} vs {len(page_keys)}")
        check("... both ends, one window each, as the page's two boxes",
              [(w["chromosome"], w["start"], w["end"]) for w in sel["windows"]] ==
              [("chr19", page_view["a"]["window"]["start"], page_view["a"]["window"]["end"]),
               ("chr22", page_view["b"]["window"]["start"], page_view["b"]["window"]["end"])], sel["windows"])
        ev = os.path.join(d, "ev.bam")
        written, missing = IR.write_evidence_bam(B, sel, ev)
        normal_in_window = [a for a in B.fetch("chr19", 47999000, 48001000)
                            if not a.is_secondary and (a.query_name, a.reference_start + 1, bool(a.is_supplementary))
                            not in page_keys]
    with pysam.AlignmentFile(ev) as E:
        recs = list(E.fetch())
    ev_keys = {(a.query_name, a.reference_start + 1, bool(a.is_supplementary)) for a in recs}
    check("the evidence BAM holds every one of them, and nothing else", ev_keys == page_keys and missing == 0,
          f"{len(ev_keys)} written, {missing} not found again")
    check("control: the fixture has normal reads in the window, which the BAM leaves out",
          len(normal_in_window) > 0 and not ({(a.query_name, a.reference_start + 1, bool(a.is_supplementary))
                                              for a in normal_in_window} & ev_keys), len(normal_in_window))
    check("the BAM is sorted and indexed", os.path.exists(ev + ".bai") and
          [a.reference_start for a in recs if a.reference_name == "chr19"] ==
          sorted(a.reference_start for a in recs if a.reference_name == "chr19"))

    # ── 2. colours and words are the page's ────────────────────────────────
    print("colours and labels")
    page = open(PAGE).read()
    jc, unc, kinds = page_colours(page)
    check("the junction colours are the page's", jc == IR.JOIN_COL and unc == IR.JOIN_UNCALLED, f"{jc} {unc}")
    check("each kind of read has the page's colour and words",
          kinds == {k: v for k, v in IR.READ_KIND.items()}, sorted(set(kinds) ^ set(IR.READ_KIND)))
    check("control: a colour changed on the page is noticed",
          page_colours(page.replace("'#0F766E'", "'#0F766F'", 1))[0] != IR.JOIN_COL)
    by_view = {key(x): x for e in ("a", "b") for x in page_view[e]["reads"]}
    tag = {(a.query_name, a.reference_start + 1, bool(a.is_supplementary)): (a.get_tag("YC"), a.get_tag("XG"))
           for a in recs}
    j1 = [tag[k] for k, x in by_view.items() if x["kind"] in ("pair_partner", "split_partner")
          and x["orientation"] == "3to5" and not x["low_mapq"]]
    j2 = [tag[k] for k, x in by_view.items() if x["kind"] in ("pair_partner", "split_partner")
          and x["orientation"] == "5to3" and not x["low_mapq"]]
    check("reads joining as junction 1 carry its colour and the page's words",
          j1 and set(j1) == {("15,118,110", "1. joins the two ends as junction 1 does")}, set(j1))
    check("reads joining as junction 2 carry its colour",
          j2 and set(j2) == {("190,24,93", "2. joins the two ends as junction 2 does")}, set(j2))
    low = [tag[k][1] for k, x in by_view.items() if x["kind"] in ("pair_partner", "split_partner") and x["low_mapq"]]
    check("a joining read below the mapping-quality floor is labelled as not counted",
          low and all(t.endswith("(mapping quality below 20, not counted)") for t in low), low)
    named = IR._joins([dict(R1[0], name="der(19)"), R1[1]])
    check("a junction's name from the page is used in its label",
          IR._style({"kind": "pair_partner", "orientation": "3to5"}, named)[1] == "1. joins the two ends as der(19) does")
    check("a joining read no called junction explains has the page's 'uncalled' colour",
          IR._style({"kind": "pair_partner", "orientation": "5to5"}, named)[0] == IR.JOIN_UNCALLED)
    check("an ambiguous read is grey in IGV, not invisible white",
          IR._style({"kind": "ambiguous", "low_mapq": True}, named)[0] == IR.IGV_AMBIGUOUS_COL != "#FFFFFF")

    # ── 3. the IGV settings ─────────────────────────────────────────────────
    print("the IGV settings")
    lines = IR.batch_lines("/x", sel["windows"], "/x/out.png", genes_loaded=True)
    gotos = [l for l in lines if l.startswith("goto ")]
    check("both ends are opened at once, side by side",
          gotos == ["goto " + " ".join(f"{w['chromosome']}:{w['start']}-{w['end']}" for w in sel["windows"])], gotos)
    for cmd in ("colorBy YC_TAG", "group TAG XG", "preference SAM.DOWNSAMPLE_READS false",
                "preference SAM.SHOW_SOFT_CLIPPED true", "preference SAM.SHOW_CENTER_LINE true",
                "preference SAM.QUALITY_THRESHOLD 0", "preference SAM.FLAG_ZERO_QUALITY true",
                "preference SAM.CLIPPING_THRESHOLD 10", "preference SAM.SHOW_COV_TRACK false"):
        check(f"the batch script says: {cmd}", cmd in lines)
    load_at = next(i for i, l in enumerate(lines) if l.startswith("load "))
    check("every preference is set before the reads are loaded",
          all(i < load_at for i, l in enumerate(lines) if l.startswith("preference ")))
    check("the local genes replace IGV's RefSeq track", 'remove "Refseq All"' in lines
          and any(l.endswith(IR.GENES_FILE) for l in lines))
    check("control: without a gene table IGV's own track stays",
          'remove "Refseq All"' not in IR.batch_lines("/x", sel["windows"], "/x/o.png", genes_loaded=False))
    bed = os.path.join(d, "genes.bed")
    names = IR.write_genes_bed(sel["windows"], bed)
    rows = [l.split("\t") for l in open(bed).read().splitlines()]
    ga = next((r for r in rows if r[3].split("_(")[0] == "GENE_A"), None)
    check("a gene broken in an intron is written whole: every exon of its canonical transcript",
          ga is not None and int(ga[9]) == 4 and ga[11].split(",")[0] == "0"
          and int(ga[1]) + int(ga[11].split(",")[-1]) + int(ga[10].split(",")[-1]) == int(ga[2]), ga)
    check("the genes are the local table's at the two windows, OMIM-marked as on the page",
          {"GENE_A (OMIM)", "LNC_F", "GENE_B"} <= set(names), names)
    gt.load(None)
    check("control: with no gene table no gene track is written", IR.write_genes_bed(sel["windows"], bed) is None)
    gt.load(*fx["genes"])

    # ── 4. one call: files, cleanup, errors ────────────────────────────────
    print("one call")
    saved_run, saved_find, saved_mk = IR.run_batch, IR.find_igv, IR._tempfile.mkdtemp
    made = []
    try:
        written = {}
        IR.run_batch = fake_igv(written)
        IR.find_igv = lambda: ("/stand-in/igv.sh", None)

        def mk(**kw):
            p = saved_mk(**kw)
            made.append(p)
            return p
        IR._tempfile.mkdtemp = mk
        out = IR.igv_two_ends(bam, R1, os.path.join(d, "img", "r1.png"), also_at=["chr19:48000201", "chr22:29999999"])
        check("the call succeeds and reports the reads it drew", out.get("success") and
              out.get("reads_drawn") == len(page_keys) and out.get("reads_not_found_again") == 0, out.get("error"))
        check("... the reads per group, with the page's colours, joining reads first",
              out["groups"][0]["label"].startswith("1. ") and out["groups"][0]["colour"] == IR.JOIN_COL[0]
              and sum(g["reads"] for g in out["groups"]) == len(page_keys), out["groups"][:2])
        check("IGV was given the evidence BAM, its index, the genes and the script",
              written["files"] == sorted([IR.TRACK_FILE, IR.TRACK_FILE + ".bai", IR.GENES_FILE, "igv_batch.txt"]),
              written["files"])
        check("the temporary files, the sample's reads among them, are deleted afterwards",
              made and not os.path.exists(made[0]) and out.get("temporary_files_deleted") is True, made)
        kept = os.path.join(d, "kept")
        IR.igv_two_ends(bam, R1, os.path.join(d, "img", "r1b.png"), keep_dir=kept)
        check("control: kept on request, the same files are there", os.path.exists(os.path.join(kept, IR.TRACK_FILE)))
        hand = IR.igv_two_ends(bam, [{"chromosome_1": "chr19", "position_1": 48000000}], os.path.join(d, "img", "h.png"))
        check("a position typed by hand is one window, with no joining read",
              len(hand.get("windows") or []) == 1 and hand.get("reads_drawn", 0) > 0
              and not any(g["label"][0].isdigit() for g in hand["groups"]), hand.get("groups"))
        mcp_out = asyncio.run(R.mcp.call_tool("igv_junction_view", {"bam_path": bam, "junctions": R1})).structured_content
        check("through the review server the image is an opaque reference, no path",
              isinstance(mcp_out.get("image_ref"), str) and "screenshot_path" not in mcp_out
              and "/" not in str(mcp_out.get("windows")), sorted(mcp_out)[:8])
        bad = IR.igv_two_ends(bam, [{"chromosome_1": "chr19"}], os.path.join(d, "img", "x.png"))
        check("a junction without a position is refused before anything runs", bad.get("error_type") == "bad_parameter")
        bad2 = IR.igv_two_ends(bam, [{"chromosome_1": "chrNOPE", "position_1": 5, "chromosome_2": "chr22",
                                     "position_2": 30000000}], os.path.join(d, "img", "y.png"))
        check("a chromosome the file lacks is an error, not an empty image", "error" in bad2, bad2)
        IR.find_igv = lambda: (None, ["~/IGV_2.17.4/igv.sh"])
        none = IR.igv_two_ends(bam, R1, os.path.join(d, "img", "z.png"))
        check("without IGV the call says so", none.get("error_type") == "igv_missing", none)
    finally:
        IR.run_batch, IR.find_igv, IR._tempfile.mkdtemp = saved_run, saved_find, saved_mk

    # ── 5. the pixel check can fail ─────────────────────────────────────────
    print("the pixel check")
    w, h, cnt = png_pixels(os.path.join(d, "img", "r1.png"))
    check("control: the pixel check finds no junction colour in a plain white image",
          (w, h) == (1, 1) and near(cnt, as_igv_draws(IR.JOIN_COL[0])) == 0, cnt)

    # ── 6. IGV itself, on public data ───────────────────────────────────────
    print("IGV itself (IMP01, public)")
    imp01 = os.path.expanduser("~/public_data/sim/bams/IMP01.bam")
    igv, _ = IR.find_igv()
    if not (igv and os.environ.get("DISPLAY") and os.path.exists(imp01)):
        NOT_RUN.append("NOT RUN: IGV on IMP01 (needs IGV, a DISPLAY and ~/public_data/sim/bams/IMP01.bam)")
        print("  " + NOT_RUN[-1])
        return
    annot = os.path.expanduser("~/public_data/annotation/")
    if os.path.exists(annot + "genes_grch38.tsv.gz"):
        gt.load(annot + "genes_grch38.tsv.gz", annot + "mim2gene.txt", annot + "genes_to_disease.txt")
    imp = [{"chromosome_1": "chr20", "position_1": 200000, "chromosome_2": "chr21", "position_2": 14100001,
            "orientation": "3to5"},
           {"chromosome_1": "chr20", "position_1": 200001, "chromosome_2": "chr21", "position_2": 14100000,
            "orientation": "5to3"}]
    png = os.path.join(d, "img", "imp01.png")
    real = IR.igv_two_ends(imp01, imp, png)
    check("IGV draws the image", real.get("success") and os.path.getsize(png) > 10000, real.get("error"))
    if not real.get("success"):
        return
    w, h, cnt = png_pixels(png)
    j1, j2 = near(cnt, as_igv_draws(IR.JOIN_COL[0])), near(cnt, as_igv_draws(IR.JOIN_COL[1]))
    check("its reads are in the page's junction colours (IGV colours by the YC tag)",
          j1 > 2000 and j2 > 2000, f"junction 1: {j1} px, junction 2: {j2} px")
    check("both ends are in it: two windows", len(real["windows"]) == 2 and w > 900, (w, real["windows"]))
    check("no normal read was loaded", real["normal_reads_hidden"] > 0
          and real["reads_drawn"] == sum(g["reads"] for g in real["groups"]), real["normal_reads_hidden"])
    print(f"  ({real['seconds']} s, {w}x{h} px, {real['reads_drawn']} reads, "
          f"{real['normal_reads_hidden']} normal reads left out)")


def main():
    with tempfile.TemporaryDirectory() as d:
        run(d)
    print("\n" + "=" * 60)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    if NOT_RUN:
        print("INCOMPLETE: " + "; ".join(NOT_RUN))
        sys.exit(2)
    print("ALL IGV REVIEW TESTS PASSED")


if __name__ == "__main__":
    main()
