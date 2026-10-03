"""
review_server.py
FastMCP server for the geneticist's view of a candidate list (Phase 26):
four tools that report read counts, coordinates and genes, and no score.

    junction_evidence   the reads that join the two ends of one junction
    review_candidates   every candidate of a sample that passes the filters,
                        grouped into rearrangements, each measured from the
                        reads, sorted by the reads that support it
    genes_near          genes at a position or in a segment, from the local
                        gene table: nothing leaves the computer
    igv_junction_view   (Phase 27) one IGV image of both ends of a rearrangement
                        with only the reads the page draws, in its colours

A separate server, like candidate_server.py: the eleven evidence tools and the
four candidate-set tools are unchanged. The interface (ui.py) runs all three
in-process through one recorder; the assistant is offered the first three in
place of the combined score (igv_junction_view is the page's: its image is not
shown to a model, and it takes a minute).

Run from repo root: python -m stage1_igv_assistant.review_server
"""

import hashlib
import os
import time
from dataclasses import asdict
from typing import Optional

from fastmcp import FastMCP

from stage1_igv_assistant.tools import bam_tools as _bt
from stage1_igv_assistant.tools import gene_table as _genes
from stage1_igv_assistant.tools import igv_review as _igv
from stage1_igv_assistant.tools import junction_tools as _jt
from stage1_igv_assistant.tools import vcf_tools as _vcf

import pysam

mcp = FastMCP(
    name="SV Review",
    instructions="""
You help a clinical geneticist review structural-variant candidates.

RULES:
1. To say which candidates are strongest, call review_candidates ONCE. It
   measures every candidate that passes the filters, groups the junctions into
   rearrangements and returns them sorted by the reads that support them.
2. Describe evidence as read counts: read pairs that join the two ends, reads
   split across the junction, reads cut at the same base. There is no score.
3. Name genes only from genes_near or from the genes a tool returned.
4. Cautions (a centromere, reads pointing elsewhere, ambiguous mapping, a
   junction also called in another sample) belong in the answer.
5. Do not interpret the biology beyond what the tools returned.
""",
)

# Set by the interface at start-up (ui.main) and by the tests; never a tool argument,
# so a model cannot name a file through them.
_CONFIG = {"mask_path": None}
REVIEW_LABEL = "review"           # the label review_candidates registers a file under
RECURRENCE_TOLERANCE_BP = _vcf.RECURRENCE_TOLERANCE_BP
MAX_JUNCTIONS = 300
CAUTION_LOW_MAPQ = 0.2            # author judgement: share of ambiguously mapped reads worth a caution
CAUTION_ELSEWHERE_MIN = 5         # author judgement: pairs pointing to other places worth a caution
RECIPROCAL_READS_MIN = 3          # author judgement: reads for an uncalled reciprocal join worth mentioning


def configure(mask_path=None, gene_table=None, mim2gene=None, gene_disorders=None):
    _CONFIG["mask_path"] = mask_path if mask_path and os.path.isfile(mask_path) else None
    _genes.load(gene_table, mim2gene, gene_disorders)


# ── helpers ─────────────────────────────────────────────────────────────────
def _set_for(path, label):
    """A registered set for this file, re-using one already loaded from it."""
    ap = os.path.abspath(path)
    for sid, s in _vcf._SETS.items():
        if s.get("path") == ap:
            return sid, None
    r = _vcf.load_candidate_set(path, label)
    if "error" in r:
        return None, r
    return r["set_id"], None


def _genes_at(chrom, pos):
    t = _genes.table()
    if t is None:
        return None
    at = t.at(chrom, pos)
    out = {"in": at}
    if not at:
        out["nearest"] = t.nearest(chrom, pos)
    return out


def _gene_word(g):
    w = g["name"] + (f" ({g['where']})" if g.get("where") else "")
    if g.get("disorders") or g.get("omim_phenotype"):
        w += " [OMIM disorder]" if g.get("omim_disorders") or g.get("omim_phenotype") else " [disorder]"
    return w


def _support_for(J, support):
    """The share of a measurement that supports this junction's own join."""
    o = J.get("orientation")
    if o and o in support.get("by_orientation", {}):
        return support["by_orientation"][o]
    if o:
        return {"read_pairs": 0, "split_reads": 0, "fragments": 0}
    return {"read_pairs": support["read_pairs"], "split_reads": support["split_reads"],
            "fragments": support["fragments"]}


def _pair_index(set_junctions):
    """Every junction of the caller's file between two chromosomes, in genome order,
    keyed by its chromosome pair (built once per review)."""
    idx = {}
    for j in set_junctions:
        if _jt.norm_chrom(j.chrom1) == _jt.norm_chrom(j.chrom2):
            continue
        v = _jt._junction_view(asdict(j))
        idx.setdefault((_jt.norm_chrom(v["a"]["chromosome"]), _jt.norm_chrom(v["b"]["chromosome"])), []).append(v)
    return idx


def _partner_in_file(index, J, near_bp):
    """The reciprocal junction of a lone junction, if the caller's file holds one."""
    a, b = J["a"], J["b"]
    want = {"3to5": "5to3", "5to3": "3to5", "3to3": "5to5", "5to5": "3to3"}.get(J.get("orientation"))
    if not want:
        return None
    for v in index.get((_jt.norm_chrom(a["chromosome"]), _jt.norm_chrom(b["chromosome"])), ()):
        if (v["candidate_id"] != J["candidate_id"] and v["orientation"] == want
                and abs(v["a"]["position"] - a["position"]) <= near_bp
                and abs(v["b"]["position"] - b["position"]) <= near_bp):
            return v
    return None


# ── tools ───────────────────────────────────────────────────────────────────
@mcp.tool()
def junction_evidence(bam_path: str, chromosome_1: str, position_1: int,
                      chromosome_2: str, position_2: int, orientation: Optional[str] = None,
                      window_bp: int = _jt.WINDOW_BP, min_mapq: int = _jt.MIN_MAPQ,
                      list_reads: bool = True, include_view: bool = False,
                      also_at: Optional[list] = None) -> dict:
    """
    The reads that join two breakpoints: chromosome_1:position_1 and
    chromosome_2:position_2 (one junction, both ends measured together).

    Returns read_pairs (pairs with one read within window_bp of each end, both
    at mapping quality >= min_mapq), split_reads (reads with one aligned piece
    at each end), fragments (distinct reads behind either), the same split by
    the join they imply (by_orientation: 3to5, 5to3, 3to3, 5to5, written from
    end 1 to end 2), the base where clipped reads pile up at each end, read
    pairs at each end that point to other places, the share of ambiguously
    mapped reads, the cytogenetic band, genes at each end from the local table,
    and, when `orientation` is given, the chromosome the join makes (its two
    pieces, which one moved and which keeps the centromere). list_reads adds the
    supporting reads themselves (name, position, strand, CIGAR, mapping quality
    at both ends). also_at lists other breakpoints of the same rearrangement
    ("chr22:40050000"): a mate there is counted with the rearrangement, not as
    pointing elsewhere. There is no score.
    """
    for nm, v in (("position_1", position_1), ("position_2", position_2)):
        if not isinstance(v, int) or isinstance(v, bool) or v < 1:
            return {"error": f"{nm} must be a positive integer (got {v!r})", "error_type": "bad_parameter"}
    if orientation is not None and orientation not in _jt.ORIENTATIONS:
        return {"error": f"orientation must be one of {list(_jt.ORIENTATIONS)} (got {orientation!r})",
                "error_type": "bad_parameter"}
    related = []
    for item in also_at or []:
        c, _, p = str(item).replace(",", "").rpartition(":")
        if c and p.isdigit():
            related.append((c, int(p)))
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": f"the read file could not be opened: {type(e).__name__}", "error_type": "bam_access"}
    try:
        s = _jt.support_with(bam, chromosome_1, position_1, chromosome_2, position_2, window_bp=window_bp,
                             min_mapq=min_mapq, list_reads=list_reads, view=include_view,
                             orientation=orientation, related=related)
    finally:
        bam.close()
    if "error" in s:
        return s
    s["genes"] = {"end_1": _genes_at(chromosome_1, position_1), "end_2": _genes_at(chromosome_2, position_2)}
    if s["genes"]["end_1"] is None:
        s["genes"] = None
        s["genes_note"] = "no local gene table is set up; genes_near cannot name genes either"
    elif include_view:
        t = _genes.table()
        for e, key in (("a", "end_a"), ("b", "end_b")):
            w = s["view"][e]["window"]
            s["view"][e]["genes"] = t.track(s[key]["chromosome"], w["start"], w["end"])
    def as_junction(o):
        return {"orientation": o, "candidate_id": None,
                "a": {"chromosome": _jt.display_chrom(chromosome_1), "position": position_1},
                "b": {"chromosome": _jt.display_chrom(chromosome_2), "position": position_2}}
    if orientation:
        J = as_junction(orientation)
        s["join"] = _jt.derivative(J)
        s["supporting_this_join"] = _support_for(J, s)
    else:
        # no orientation given (a junction typed in by hand): the chromosome each
        # join the reads support would make
        s["joins"] = {o: _jt.derivative(as_junction(o)) for o in s["by_orientation"]}
    return s


@mcp.tool()
def review_candidates(path: str, bam_path: str, svtype: Optional[str] = "BND", filter_pass: bool = True,
                      min_pe: Optional[int] = 3, min_sr: Optional[int] = 1, primary_only: bool = True,
                      exclude_masked: bool = True, other_path: Optional[str] = None,
                      recurrence_tolerance_bp: int = RECURRENCE_TOLERANCE_BP,
                      max_junctions: int = MAX_JUNCTIONS) -> dict:
    """
    Review every candidate of one sample at once: which are strongest.

    Reads the caller's candidate file (path), keeps the junctions that pass the
    filters (the interface's defaults: translocations, caller PASS, read pairs
    >= 3, split reads >= 1, main chromosomes, outside the caller's exclude
    regions), groups junctions that belong together into rearrangements (the
    two junctions of a balanced translocation, the two ends of a moved
    segment), counts in the read file (bam_path) the read pairs and split reads
    that support each, names the genes at each breakpoint from the local table
    (with OMIM marks when set up), and returns the rearrangements sorted by the
    reads that support them, with cautions. other_path names a second sample's
    candidate file: rearrangements also called there are marked. Every filter
    is echoed in filters_applied. There is no score: the order is the count of
    supporting reads. summary_rows gives one line per rearrangement, in order.
    """
    for nm, v in (("min_pe", min_pe), ("min_sr", min_sr), ("max_junctions", max_junctions),
                  ("recurrence_tolerance_bp", recurrence_tolerance_bp)):
        if v is not None and (not isinstance(v, int) or isinstance(v, bool) or v < 0):
            return {"error": f"{nm} must be a non-negative integer (got {v!r})", "error_type": "bad_parameter"}
    sid, err = _set_for(path, REVIEW_LABEL)
    if err:
        return err
    mask = _CONFIG["mask_path"] if exclude_masked else None
    lc = _vcf.list_candidates(sid, svtype=svtype or None, filter_pass=bool(filter_pass), min_pe=min_pe,
                              min_sr=min_sr, primary_only=bool(primary_only), mask_path=mask,
                              limit=max(1, max_junctions or MAX_JUNCTIONS))
    if "error" in lc:
        return lc
    survivors = lc["candidates"]
    recurrent = set()
    other_note = None
    if other_path:
        osid, oerr = _set_for(other_path, REVIEW_LABEL + "-other")
        if oerr:
            return oerr
        cmp_ = _vcf.compare_candidate_sets(sid, osid, recurrence_tolerance_bp)
        if "error" in cmp_:
            return cmp_
        recurrent = {p["candidate_id_a"] for p in cmp_.get("matched_pairs", [])}
        other_note = (f"junctions whose both ends lie within {recurrence_tolerance_bp:,} bp of a "
                      f"junction in the comparison sample are marked recurrent")
    events = _jt.group_rearrangements(survivors)
    index = _pair_index(_vcf._SETS[sid]["junctions"])

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": f"the read file could not be opened: {type(e).__name__}", "error_type": "bam_access"}
    try:
        for ev in events:
            _measure(ev, bam, recurrent, index)
    finally:
        bam.close()
    events.sort(key=lambda e: (-e["support"]["fragments"], -sum((j.get("pe") or 0) + (j.get("sr") or 0)
                                                                  for j in e["junctions"])))
    rows = []
    for i, ev in enumerate(events, 1):
        ev["rank"] = i
        rows.append(_row(ev))
    return {
        "total_in_set": lc["total_in_set"],
        "junctions_reviewed": len(survivors),
        "junctions_passing_filters": lc["total_matching"],
        "truncated": lc["truncated"],
        "rearrangements_found": len(events),
        "filters_applied": lc["filters_applied"],
        "exclude_regions_available": _CONFIG["mask_path"] is not None,
        "recurrence": other_note,
        "gene_table": _genes.status(),
        "rearrangements": events,
        "summary_rows": rows,
        "thresholds_applied": [
            {"name": "window_bp", "value": _jt.WINDOW_BP, "provenance": _jt.PROVENANCE["window_bp"]},
            {"name": "min_mapq", "value": _jt.MIN_MAPQ, "provenance": _jt.PROVENANCE["min_mapq"]},
            {"name": "group_near_bp", "value": _jt.GROUP_NEAR_BP, "provenance": _jt.PROVENANCE["group_near_bp"]},
            {"name": "segment_max_bp", "value": _jt.SEGMENT_MAX_BP, "provenance": _jt.PROVENANCE["segment_max_bp"]},
            {"name": "caution_low_mapq", "value": CAUTION_LOW_MAPQ, "provenance": "author judgement"},
        ],
        "note": ("Sorted by the number of distinct reads (read pairs and split reads, both ends at "
                 "mapping quality >= 20) that support each rearrangement's own joins. A rearrangement "
                 "whose two junctions were both called is a candidate balanced translocation. Cautions "
                 "say why a count may mislead. Nothing here is a score or a diagnosis."
                 + (" Only the first max_junctions junctions were reviewed." if lc["truncated"] else "")),
    }


def _measure(ev, bam, recurrent, index):
    joins, cautions = [], []
    total = {"read_pairs": 0, "split_reads": 0, "fragments": 0}
    seen_orient = {}
    measured = []
    for J in ev["junctions"]:
        others = [(K[e]["chromosome"], K[e]["position"]) for K in ev["junctions"] if K is not J for e in ("a", "b")]
        s = _jt.support_with(bam, J["a"]["chromosome"], J["a"]["position"],
                             J["b"]["chromosome"], J["b"]["position"], orientation=J.get("orientation"),
                             related=others)
        if "error" in s:
            J["support"] = None
            J["error"] = s["error"]
            cautions.append(f"the reads at {J['a']['chromosome']}:{J['a']['position']:,} could not be read: {s['error']}")
            continue
        own = _support_for(J, s)
        J["support"] = own
        J["ends"] = {"a": s["end_a"], "b": s["end_b"]}
        J["recurrent"] = J["candidate_id"] in recurrent
        J["join"] = _jt.derivative(J)
        measured.append((J, s))
        o = J.get("orientation") or "?"
        prev = seen_orient.get(o)
        if prev is None or own["fragments"] > prev["fragments"]:
            seen_orient[o] = own
    for own in seen_orient.values():
        for k in total:
            total[k] += own[k]
    ev["support"] = total
    # cautions, each a plain sentence
    for bp in ev["breakpoints"]:
        if bp.get("region_note"):
            cautions.append(f"breakpoint {bp['region_note']}: repeats here often produce false junctions")
    for J, s in measured:
        for end in (s["end_a"], s["end_b"]):
            own_pairs = J["support"]["read_pairs"]
            if end["pairs_to_other_places"] >= CAUTION_ELSEWHERE_MIN and end["pairs_to_other_places"] >= own_pairs:
                top = ", ".join(f"{k} {v}" for k, v in list(end["other_places"].items())[:3])
                cautions.append(f"at {end['chromosome']}:{end['position']:,}, {end['pairs_to_other_places']} read pairs "
                                f"point to other places ({top}), as many as or more than to the partner")
            if end["low_mapq_fraction"] is not None and end["low_mapq_fraction"] >= CAUTION_LOW_MAPQ:
                cautions.append(f"at {end['chromosome']}:{end['position']:,}, {round(100 * end['low_mapq_fraction'])}% "
                                f"of reads map ambiguously")
    if measured and total["fragments"] == 0:
        cautions.append("no read in this file supports it at mapping quality 20 or more")
    if any(J.get("recurrent") for J in ev["junctions"]):
        cautions.append("also called in the comparison sample: an artefact or a common variant is likely")
    # a lone junction: is its reciprocal partner in the file, or in the reads?
    if ev.get("pattern") == "one junction" and ev["kind"] == "translocation" and measured:
        J, s = measured[0]
        p = _partner_in_file(index, J, _jt.GROUP_NEAR_BP)
        if p:
            ev["partner_in_file"] = {k: p[k] for k in ("candidate_id", "caller_id", "filter", "pe", "sr",
                                                         "orientation", "a", "b")}
            cautions.append(f"the reciprocal junction is in the caller's file ({p.get('caller_id') or p['candidate_id']}: "
                            f"{p.get('filter')}, read pairs {p.get('pe')}, split reads {p.get('sr')}) but did not pass the filters")
        want = {"3to5": "5to3", "5to3": "3to5", "3to3": "5to5", "5to5": "3to3"}.get(J.get("orientation"))
        rec = s["by_orientation"].get(want) if want else None
        if rec and rec["fragments"] >= RECIPROCAL_READS_MIN:
            ev["reciprocal_reads"] = {"orientation": want, **rec}
            cautions.append(f"{rec['fragments']} reads also support the reciprocal join, which the caller did not report")
    ev["cautions"] = cautions
    # genes
    t = _genes.table()
    if t is not None:
        ev["genes"] = []
        for bp in ev["breakpoints"]:
            ev["genes"].append({"chromosome": bp["chromosome"], "start": bp["start"], "end": bp["end"],
                                **(_genes_at(bp["chromosome"], bp["start"]) or {}),
                                **({"between": t.overlapping(bp["chromosome"], bp["start"], bp["end"])}
                                   if bp["end"] - bp["start"] > 1 else {})})
        if ev.get("segment"):
            seg = ev["segment"]
            ev["segment"]["genes"] = t.overlapping(seg["chromosome"], seg["start"], seg["end"])
        for J in ev["junctions"]:
            mv = (J.get("join") or {}).get("moved")
            if mv and mv.get("start") and mv.get("end"):
                o = t.overlapping(mv["chromosome"], mv["start"], mv["end"], limit=0)
                mv["genes"] = {k: o[k] for k in ("count", "protein_coding", "with_disorder")}


def _row(ev):
    bps = " ↔ ".join(f"{b['chromosome']}:{b['start']:,}" + (f"–{b['end']:,}" if b["end"] != b["start"] else "")
                     for b in ev["breakpoints"])
    name = ev.get("iscn") or ev["kind"]
    s = ev["support"]
    genes = []
    for g in ev.get("genes") or []:
        if g.get("in"):
            genes += [_gene_word(x) for x in g["in"][:3]]
    gtxt = ("breaks " + ", ".join(genes)) if genes else ("no gene at the breakpoints" if ev.get("genes") is not None
                                                          else "genes not looked up (no local gene table)")
    caut = ("; cautions: " + "; ".join(ev["cautions"])) if ev.get("cautions") else ""
    def n(k, one, many):
        return f"{k} {one if k == 1 else many}"
    return (f"{ev['rank']}. {name} {bps} - {ev['pattern']}; {n(s['read_pairs'], 'read pair', 'read pairs')}, "
            f"{n(s['split_reads'], 'split read', 'split reads')} ({n(s['fragments'], 'read', 'reads')}); {gtxt}{caut}")


@mcp.tool()
def genes_near(chromosome: str, position: int, end: Optional[int] = None) -> dict:
    """
    Genes at a position, or overlapping position..end, from the local gene
    table (nothing is sent over the network). For a position in no gene, the
    nearest gene on each side within 1 Mb. Each gene carries its OMIM gene
    number and linked disorders when those files are set up, and, at a
    breakpoint, where it falls in the gene's canonical transcript (exon or
    intron number). Also returns the cytogenetic band (GRCh38).
    """
    if not isinstance(position, int) or isinstance(position, bool) or position < 1:
        return {"error": f"position must be a positive integer (got {position!r})", "error_type": "bad_parameter"}
    t = _genes.table()
    band, _ = _jt.band_at(chromosome, position)
    out = {"chromosome": _jt.display_chrom(chromosome), "position": position, "band": band,
           "region_note": _jt.region_note(chromosome, position)}
    if t is None:
        out.update(_genes.status())
        out["error"] = "no local gene table is set up (see scripts/make_gene_table.py)"
        out["error_type"] = "no_gene_table"
        return out
    if end is not None and isinstance(end, int) and end > position:
        out["end"] = end
        out["genes"] = t.overlapping(chromosome, position, end)
    else:
        out["genes_at_position"] = t.at(chromosome, position)
        if not out["genes_at_position"]:
            out["nearest"] = t.nearest(chromosome, position)
    out["source"] = "local gene table" + (", OMIM marks" if t.has_disease_marks else "")
    return out


@mcp.tool()
def igv_junction_view(bam_path: str, junctions: list, also_at: Optional[list] = None) -> dict:
    """
    One IGV image of a rearrangement, as the supervisor's review asked: every end
    side by side (IGV's multi-locus view, each window centred on its breakpoint),
    only the abnormal reads (exactly the reads junction_evidence draws: normal
    reads are not loaded), each coloured and grouped as the interface's read
    boxes colour and label them, reads at mapping quality 0 hollow, soft clips
    shown, no downsampling, genes from the local table.

    junctions: the rearrangement's junctions in the interface's order, each
    {"chromosome_1", "position_1", "chromosome_2", "position_2", "orientation",
    "name"} (orientation and name optional); a position typed by hand is one
    junction with no second end. also_at: other breakpoints of the rearrangement
    ("chr22:40050000"), as for junction_evidence.

    Returns an opaque image_ref (no path), the windows, the reads drawn per group,
    the normal reads left out, the genes shown and the IGV settings used. The
    temporary BAM of the drawn reads is deleted when IGV is done.
    """
    if also_at is not None and not isinstance(also_at, list):
        return {"error": "also_at must be a list of 'chromosome:position'", "error_type": "bad_parameter"}
    session = _bt.image_session_dir()
    tag = hashlib.sha256(repr((junctions, also_at, time.time())).encode()).hexdigest()[:12]
    png = os.path.join(session, f"igv_two_ends_{tag}.png")
    result = _igv.igv_two_ends(bam_path, junctions, png, also_at=also_at)
    if "screenshot_path" in result or result.get("success") is not None:
        return _bt.to_handle_result(result, session)
    return result


if __name__ == "__main__":
    mcp.run()
