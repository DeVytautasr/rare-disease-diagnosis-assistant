"""
igv_review.py -- one IGV image of a rearrangement, as the supervisor asked (Phase 27).

The supervisor's review (2 October 2026) found the IGV images uninformative and
asked instead for two boxes: where the piece left and where it arrived, showing
only the abnormal reads, coloured by what they say, with the normal reads hidden
and the ambiguous ("white", mapping quality 0) reads kept. The page draws that
itself (junction_tools). This module makes IGV show the same thing, so a
geneticist who opens IGV sees the reads the page counted and drew:

  * the reads are exactly the page's: junction_tools.support_with(view=True), the
    call behind the page's read boxes, chooses them, and a small BAM holding only
    those reads (an "evidence BAM") is what IGV loads. Normal reads are absent
    from it, so IGV cannot show them; ambiguous reads are in it, and IGV draws
    those at mapping quality 0 hollow (SAM.FLAG_ZERO_QUALITY, its default).
  * every read carries the page's colour in its YC tag, and IGV colours by it
    (colorBy YC_TAG). IGV draws every read at about 75% opacity, so its colours
    are the page's, lightened.
  * every read carries the page's words for what it shows in an XG tag, and IGV
    groups by it (group TAG XG): each group is labelled, and the joining reads
    come first (their labels start with the junction's number).
  * both ends in one image: one `goto` with every end's window (IGV's multi-locus
    view), each window centred on its breakpoint, with IGV's centre line on.
  * soft-clipped bases shown, clipping of 10 or more bases flagged (the page's
    MIN_CLIP_BASES), no downsampling (every drawn read is in the image), no
    coverage track (it would be the coverage of the abnormal reads alone).
  * genes from the table on this computer, the same genes the page names, instead
    of IGV's RefSeq track.

The evidence BAM, the gene BED and the batch script are written to a temporary
directory that is deleted when IGV is done, so no copy of a sample's reads is left
behind. Only the PNG is kept. bam_tools.py (the protected evidence tools, and the
four-image panel the previous page shows) is not modified.

IGV still asks igv.org for the genome (reference sequence and ideogram) of the
regions shown, as every IGV image here always has.
"""

import os as _os
import shutil as _shutil
import signal as _signal
import subprocess as _subprocess
import tempfile as _tempfile
import time as _time

import pysam

from stage1_igv_assistant.tools import bam_tools as _bt
from stage1_igv_assistant.tools import gene_table as _genes
from stage1_igv_assistant.tools import junction_tools as _jt

# The page's colours (ui_page.html JOIN_COL, JOIN_UNCALLED, READ_KIND). A test
# reads them back from the page, so the two cannot drift apart.
JOIN_COL = ["#0F766E", "#BE185D", "#7C3AED", "#0369A1", "#4D7C0F"]
JOIN_UNCALLED = "#C2410C"
READ_KIND = {
    "pair_related":     ("#7FB8B4", "pair joining the other junction of this rearrangement"),
    "pair_other_chrom": ("#64748B", "mate on another chromosome"),
    "split_other":      ("#A3AEC0", "split read, other piece elsewhere"),
    "mate_unmapped":    ("#CA8A04", "mate not mapped"),
    "pair_orientation": ("#16A34A", "pair in the wrong orientation"),
    "pair_insert_size": ("#DC2626", "pair too far apart"),
    "pair_alt_contig":  ("#B8B2A9", "mate on an alternative contig"),
    "clipped":          ("#CFCAD8", "clipped, no other piece"),
    "ambiguous":        ("#FFFFFF", "placed equally well elsewhere"),
}
# The page draws an ambiguous read white with a grey outline. IGV has no outline
# for a read at mapping quality 1-19, so a white read would vanish: grey instead.
IGV_AMBIGUOUS_COL = "#8C8798"
PARTNER_KINDS = ("pair_partner", "split_partner", "pair_related")

MAX_WINDOWS = 6           # author judgement: more windows than this are too narrow to read
IGV_TIMEOUT_S = 600
GENOME = "hg38"
TRACK_FILE = "abnormal_reads_only.bam"
GENES_FILE = "genes_from_this_computer.bed"
IGV_CANDIDATES = ["~/IGV_2.17.4/igv.sh", "~/igv/igv.sh", "/opt/igv/igv.sh"]

# The IGV settings, each with why. Preferences come before `load`, so the track
# is built with them.
PREFERENCES = [
    ("SAM.SHOW_COV_TRACK", "false", "the coverage of the abnormal reads alone would mislead"),
    ("SAM.DOWNSAMPLE_READS", "false", "every read the page drew is in the image"),
    ("SAM.QUALITY_THRESHOLD", "0", "reads at mapping quality 0 are kept (drawn hollow)"),
    ("SAM.FLAG_ZERO_QUALITY", "true", "a read at mapping quality 0 is drawn hollow, as the page does"),
    ("SAM.FILTER_SUPPLEMENTARY_ALIGNMENTS", "false", "the pieces of split reads stay"),
    ("SAM.SHOW_SOFT_CLIPPED", "true", "the clipped part of a read at the breakpoint is shown"),
    ("SAM.FLAG_CLIPPING", "true", "a read cut at the breakpoint is marked"),
    ("SAM.CLIPPING_THRESHOLD", str(_jt.MIN_CLIP_BASES), "as the page: a read is cut when 10 or more bases are clipped"),
    ("SAM.SHOW_CENTER_LINE", "true", "each window is centred on its breakpoint"),
]


def _rgb(hex_colour):
    h = hex_colour.lstrip("#")
    return ",".join(str(int(h[i:i + 2], 16)) for i in (0, 2, 4))


def _locus(w):
    return f"{w['chromosome']}:{w['start']}-{w['end']}"


def _parse_also(also_at):
    out = []
    for item in also_at or []:
        c, _, p = str(item).replace(",", "").rpartition(":")
        if c and p.isdigit():
            out.append((c, int(p)))
    return out


def _check_junctions(junctions):
    if not isinstance(junctions, list) or not junctions:
        return "junctions must be a non-empty list"
    for i, J in enumerate(junctions):
        if not isinstance(J, dict):
            return f"junction {i + 1} is not an object"
        p1, p2 = J.get("position_1"), J.get("position_2")
        if not J.get("chromosome_1") or not isinstance(p1, int) or isinstance(p1, bool) or p1 < 1:
            return f"junction {i + 1}: chromosome_1 and a positive position_1 are required"
        if (J.get("chromosome_2") is None) != (p2 is None):
            return f"junction {i + 1}: give both chromosome_2 and position_2, or neither"
        if p2 is not None and (not isinstance(p2, int) or isinstance(p2, bool) or p2 < 1):
            return f"junction {i + 1}: position_2 must be a positive integer"
        if J.get("orientation") not in (None,) + _jt.ORIENTATIONS:
            return f"junction {i + 1}: orientation must be one of {list(_jt.ORIENTATIONS)}"
        if J.get("chromosome_2") is None and len(junctions) > 1:
            return "a position typed by hand (no second end) is drawn on its own"
    return None


def _joins(junctions):
    """One colour and one name per junction, as the page's joinsOf()."""
    n = len(junctions)
    out = []
    for i, J in enumerate(junctions):
        name = J.get("name") or ("the junction" if n == 1 else f"junction {i + 1}")
        out.append({"orientation": J.get("orientation"), "colour": JOIN_COL[i % len(JOIN_COL)],
                    "number": i + 1, "name": name})
    return out


def _style(read, joins):
    """(colour, group label) for one drawn read: the page's colour and words."""
    low = bool(read.get("low_mapq"))
    kind = read.get("kind")
    if kind in ("pair_partner", "split_partner"):
        j = next((j for j in joins if j["orientation"] and j["orientation"] == read.get("orientation")), None)
        if j:
            what = "the junction does" if j["name"] == "the junction" else f"{j['name']} does"
            colour, label = j["colour"], f"{j['number']}. joins the two ends as {what}"
        else:
            colour, label = JOIN_UNCALLED, "8. joins the two ends in a way no called junction explains"
        if low:
            label += " (mapping quality below 20, not counted)"
        return colour, label
    colour, label = READ_KIND.get(kind, ("#999999", kind or "other"))
    if kind == "pair_related":
        j = next((j for j in joins if j["orientation"] == read.get("orientation")), None)
        colour = j["colour"] if j else colour
    if kind == "ambiguous":
        colour = IGV_AMBIGUOUS_COL
    elif low:
        label += " (mapping quality below 20)"
    return colour, label


def select_reads(bam, junctions, also_at=None):
    """The reads the page draws for these junctions, each with its colour and label.

    Returns {"error"} or {"windows": [...], "reads": {(name, start, supplementary):
    (colour, label, window index)}, "normal_reads_hidden": n, "views": n}.
    As the page's readsSection(): when every junction's ends lie inside the first
    junction's two windows, that one view is drawn; otherwise one view per junction.
    A position typed by hand (one junction with no second end) is one window.
    """
    related = _parse_also(also_at)
    joins = _joins(junctions)
    J0 = junctions[0]
    if J0.get("chromosome_2") is None:
        E = _jt.scan_end(bam, J0["chromosome_1"], J0["position_1"], "", 0, related=related)
        if "error" in E:
            return E
        for x in E["view"]:
            x.pop("sides", None)
        views = [[(J0["chromosome_1"], {"window": E["window"], "reads": E["view"]})]]
        hidden = E["normal_reads_hidden"]
    else:
        def view_of(J):
            s = _jt.support_with(bam, J["chromosome_1"], J["position_1"], J["chromosome_2"], J["position_2"],
                                 view=True, orientation=J.get("orientation"), related=related)
            return s

        first = view_of(J0)
        if "error" in first:
            return first
        wa, wb = first["view"]["a"]["window"], first["view"]["b"]["window"]
        together = all(wa["start"] <= J["position_1"] <= wa["end"] and wb["start"] <= J["position_2"] <= wb["end"]
                       for J in junctions)
        results = [(J0, first)] + ([] if together else [(J, view_of(J)) for J in junctions[1:]])
        views, hidden = [], 0
        for J, s in results:
            if "error" in s:
                return s
            views.append([(J["chromosome_1"], s["view"]["a"]), (J["chromosome_2"], s["view"]["b"])])
            hidden += s["view"]["a"]["normal_reads_hidden"] + s["view"]["b"]["normal_reads_hidden"]
    windows, reads = [], {}
    for pair in views:
        for chrom, V in pair:
            w = {"chromosome": _jt.display_chrom(chrom), "start": V["window"]["start"], "end": V["window"]["end"]}
            if w not in windows:
                windows.append(w)
            wi = windows.index(w)
            for x in V["reads"]:
                key = (x["name"], x["start"], bool(x["supplementary"]))
                if key not in reads:
                    reads[key] = _style(x, joins) + (wi,)
    if len(windows) > MAX_WINDOWS:
        return {"error": f"{len(windows)} windows would be needed; at most {MAX_WINDOWS} fit in one image",
                "error_type": "too_many_windows"}
    return {"windows": windows, "reads": reads, "normal_reads_hidden": hidden, "views": len(views)}


def write_evidence_bam(bam, selection, path):
    """A sorted, indexed BAM holding only the selected reads, each tagged with its
    colour (YC) and label (XG). Returns the number of alignments written."""
    unsorted = path + ".unsorted.bam"
    want = selection["reads"]
    found = set()
    with pysam.AlignmentFile(unsorted, "wb", header=bam.header) as out:
        for w in selection["windows"]:
            contig = _bt._resolve_contig(bam, w["chromosome"])
            for a in bam.fetch(contig, w["start"] - 1, w["end"]):
                if a.is_unmapped or a.is_secondary or a.is_duplicate or a.is_qcfail:
                    continue
                key = (a.query_name, a.reference_start + 1, bool(a.is_supplementary))
                if key in want and key not in found:
                    colour, label, _ = want[key]
                    a.set_tag("YC", _rgb(colour), "Z")
                    a.set_tag("XG", label, "Z")
                    out.write(a)
                    found.add(key)
    pysam.sort("-o", path, unsorted)
    pysam.index(path)
    _os.remove(unsorted)
    return len(found), len(set(want) - found)


def write_genes_bed(windows, path):
    """Genes overlapping the windows, from the local table, whole (every exon of the
    canonical transcript), so a breakpoint in an intron shows as one. Returns the
    gene names written, or None when no gene table is set up."""
    t = _genes.table()
    if t is None:
        return None
    names, seen = [], set()
    with open(path, "w") as fh:
        for w in windows:
            for g in t.track(w["chromosome"], w["start"], w["end"]):
                if g["gene_id"] in seen:
                    continue
                seen.add(g["gene_id"])
                full = t.track(w["chromosome"], g["start"], g["end"])
                exons = sorted(next((x["exons"] for x in full if x["gene_id"] == g["gene_id"]), []))
                if not exons:
                    exons = [[g["start"], g["end"]]]
                s0, e0 = exons[0][0], exons[-1][1]
                label = g["name"] + (" (OMIM)" if g.get("disorders") or g.get("omim_phenotype") else "")
                fh.write("\t".join(str(v) for v in (
                    w["chromosome"], s0 - 1, e0, label.replace(" ", "_"), 0, g["strand"], s0 - 1, e0, "62,42,115",
                    len(exons), ",".join(str(b - a + 1) for a, b in exons),
                    ",".join(str(a - s0) for a, b in exons))) + "\n")
                names.append(label)
    return names


def batch_lines(directory, windows, png_path, genes_loaded):
    lines = ["new", f"genome {GENOME}"]
    lines += [f"preference {k} {v}" for k, v, _ in PREFERENCES]
    lines.append(f"load {_os.path.join(directory, TRACK_FILE)}")
    if genes_loaded:
        lines.append(f"load {_os.path.join(directory, GENES_FILE)}")
        lines.append('remove "Refseq All"')
    lines += ["maxPanelHeight 2000",
              "goto " + " ".join(_locus(w) for w in windows),
              "colorBy YC_TAG",
              "group TAG XG",
              f"expand {TRACK_FILE}",
              f"snapshot {_os.path.abspath(png_path)}",
              "exit"]
    return lines


def find_igv():
    env = _os.environ.get("IGV_PATH")
    for c in ([env] if env else []) + [_os.path.expanduser(c) for c in IGV_CANDIDATES]:
        if c and _os.path.exists(c):
            return c, None
    return None, ([env] if env else []) + IGV_CANDIDATES


def run_batch(igv_path, batch_path, png_path, timeout_sec=IGV_TIMEOUT_S):
    """Run IGV on a batch script and wait for the PNG (written twice the same size,
    1 s apart), then stop IGV's process group: IGV has been seen to hang on exit
    after a snapshot (bam_tools.run_igv_screenshot does the same)."""
    if _os.path.exists(png_path):
        _os.remove(png_path)
    log = _tempfile.NamedTemporaryFile(mode="w+", suffix=".igvlog", delete=False)
    proc = _subprocess.Popen([igv_path, "--batch", batch_path], stdout=log, stderr=_subprocess.STDOUT,
                             start_new_session=True)
    deadline, last, how = _time.time() + timeout_sec, None, None
    while True:
        if proc.poll() is not None:
            how = "clean_exit"
            break
        if _time.time() >= deadline:
            how = "timeout"
            break
        if _os.path.exists(png_path):
            size = _os.path.getsize(png_path)
            if size and size == last:
                how = "terminated_after_snapshot"
                break
            last = size or None
            _time.sleep(1.0)
            continue
        _time.sleep(0.5)
    if proc.poll() is None:
        _bt._signal_igv_process_group(proc, _signal.SIGTERM)
        try:
            proc.wait(timeout=5)
        except _subprocess.TimeoutExpired:
            _bt._signal_igv_process_group(proc, _signal.SIGKILL)
            try:
                proc.wait(timeout=5)
            except _subprocess.TimeoutExpired:
                pass
    log.seek(0)
    severe = [l.strip() for l in log.read().splitlines() if "SEVERE" in l][:3]
    log.close()
    _os.remove(log.name)
    ok = _os.path.exists(png_path) and _os.path.getsize(png_path) > 0
    return {"success": ok, "shutdown": how, "igv_severe_lines": severe}


def igv_two_ends(bam_path, junctions, png_path, also_at=None, timeout_sec=IGV_TIMEOUT_S, keep_dir=None):
    """Select the page's reads, write the evidence BAM and the gene track, run IGV.

    keep_dir (tests only) keeps the temporary files in that directory instead of
    deleting them."""
    problem = _check_junctions(junctions)
    if problem:
        return {"error": problem, "error_type": "bad_parameter"}
    igv, searched = find_igv()
    if igv is None:
        return {"error": "IGV not found", "searched": searched, "error_type": "igv_missing"}
    t0 = _time.time()
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": f"the read file could not be opened: {type(e).__name__}", "error_type": "bam_access"}
    tmp = keep_dir or _tempfile.mkdtemp(prefix="igv_two_ends_")
    _os.makedirs(tmp, exist_ok=True)
    try:
        sel = select_reads(bam, junctions, also_at)
        if "error" in sel:
            return sel
        written, missing = write_evidence_bam(bam, sel, _os.path.join(tmp, TRACK_FILE))
        genes = write_genes_bed(sel["windows"], _os.path.join(tmp, GENES_FILE))
        lines = batch_lines(tmp, sel["windows"], png_path, genes is not None)
        batch = _os.path.join(tmp, "igv_batch.txt")
        with open(batch, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        _os.makedirs(_os.path.dirname(_os.path.abspath(png_path)), exist_ok=True)
        run = run_batch(igv, batch, png_path, timeout_sec)
    finally:
        bam.close()
        if not keep_dir:
            _shutil.rmtree(tmp, ignore_errors=True)
    groups = {}
    for colour, label, _ in sel["reads"].values():
        g = groups.setdefault(label, {"label": label, "colour": colour, "reads": 0})
        g["reads"] += 1
    out = {
        "success": run["success"],
        "screenshot_path": _os.path.abspath(png_path) if run["success"] else None,
        "windows": [{**w, "locus": _locus(w)} for w in sel["windows"]],
        "reads_drawn": written,
        "reads_not_found_again": missing,
        "groups": [groups[k] for k in sorted(groups)],      # IGV's order: by label
        "normal_reads_hidden": sel["normal_reads_hidden"],
        "genes_shown": genes,
        "settings": [{"preference": k, "value": v, "why": why} for k, v, why in PREFERENCES]
                    + [{"command": "colorBy YC_TAG", "why": "each read in the page's colour"},
                       {"command": "group TAG XG", "why": "each group labelled with the page's words"},
                       {"command": "goto (every end at once)", "why": "both ends side by side"}],
        "igv_shutdown": run["shutdown"],
        "seconds": round(_time.time() - t0, 1),
        "temporary_files_deleted": not keep_dir,
        "note": ("Only the reads the page draws are in this image: normal reads are not loaded. Colours and "
                 "group labels are the page's; IGV draws reads at about 75% opacity, so its colours are "
                 "lighter. A read at mapping quality 0 is hollow. Genes are the local table's."),
    }
    if not run["success"]:
        out["error"] = "IGV produced no image" + (f" ({run['shutdown']})" if run["shutdown"] else "")
        out["igv_severe_lines"] = run["igv_severe_lines"]
    return out
