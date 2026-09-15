"""
vcf_tools.py
Caller-agnostic VCF/BCF parsing for SV candidate lists, plus the session
registry backing the candidate-set MCP tools.

Every function returns structured data only -- no genomic interpretation.
Two hard rules this module exists to enforce:

  1. No return value carries a sample name. VCF sample columns and the @RG
     SM: value they derive from are never read.
  2. No return value carries a filesystem path. A registered set is referred
     to by an opaque set_id, exactly as images are referred to by an opaque
     image_ref in bam_tools.

Candidate coordinates reach the assistant only as the return value of a
tool call, never from the prompt.
"""

import hashlib as _hashlib
import os as _os
import re as _re
from dataclasses import dataclass, asdict, field
from typing import Optional

import pysam


# ── Threshold inventory for this module ─────────────────────────────────────
#
# Same convention as bam_tools: a "threshold" is any numeric cutoff that
# changes what the assistant reports. Every one below is echoed back in the
# tool result that applied it, with a provenance label, because a candidate
# list quoted without its filters is not reproducible.
#
#   DEDUP_TOLERANCE_BP = 500     author judgement. Two junctions whose BOTH
#                                breakends fall within this distance are
#                                treated as one. Reused from the concordance
#                                tolerance, which was itself chosen from the
#                                measured CIPOS distribution (median ~600 bp,
#                                max ~1054 bp) -- so the number is informed by
#                                data but not derived from it.
#   RECURRENCE_TOLERANCE_BP = 500  author judgement, same figure, different
#                                purpose (cross-set matching).
#
# min_pe / min_sr are NOT defaults here. They are caller-supplied filters with
# no default value, so that a filtered list can never be produced without the
# caller having named the threshold explicitly.

DEDUP_TOLERANCE_BP = 500
RECURRENCE_TOLERANCE_BP = 500

PROVENANCE = {
    "filter_pass":      "tool-defined",       # DELLY/manta's own FILTER column
    "svtype":           "tool-defined",       # the caller's own SVTYPE
    "min_pe":           "author judgement",
    "min_sr":           "author judgement",
    "primary_only":     "reference-defined",  # contig classes of the reference
    "exclude_masked":   "reference-defined",  # the caller's shipped exclude template
    "dedup_tolerance_bp":      "author judgement",
    "recurrence_tolerance_bp": "author judgement",
}

PRIMARY_CONTIGS = frozenset(
    [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    + [str(i) for i in range(1, 23)] + ["X", "Y"]
)

# ALT breakend notation -> the same orientation vocabulary DELLY uses in
# INFO/CT, so orientation is comparable across callers.
#   t[p[  piece extending right of p joined after t     -> 3to5
#   t]p]  reverse-comp piece left of p joined after t   -> 3to3
#   ]p]t  piece extending left of p joined before t     -> 5to3
#   [p[t  reverse-comp piece right of p joined before t -> 5to5
_ALT_BND = _re.compile(r'^(?P<pre>[A-Za-z.]*)(?P<ob>[\[\]])(?P<chrom>[^\[\]:]+):(?P<pos>\d+)(?P<cb>[\[\]])(?P<post>[A-Za-z.]*)$')


def _scrub(message: str, *paths) -> str:
    """Remove filesystem paths from a library exception message.

    pysam embeds the full input path in its errors ("invalid file `/x/y.vcf`").
    Passing that straight through would put a real filename -- and therefore a
    possible patient identifier -- into a tool return. Every error message this
    module emits is scrubbed of the paths it was given, and of their basenames,
    before it leaves.
    """
    out = str(message)
    for pth in paths:
        if not pth:
            continue
        for variant in (str(pth), _os.path.abspath(str(pth)),
                        _os.path.basename(str(pth))):
            if variant:
                out = out.replace(variant, "<path withheld>")
    return out


def _info(rec, key, declared):
    """INFO field value, or None if the header does not declare the key.

    pysam raises ValueError("Invalid header") from rec.info.get() for a key the
    header never declared -- it does NOT return the default. Every INFO access
    in this module goes through here so a file that simply lacks a field is
    parsed rather than rejected.
    """
    if key not in declared:
        return None
    try:
        return rec.info.get(key)
    except (ValueError, KeyError):
        return None


def _flag(rec, key, declared):
    """True if a Flag INFO is present, False if declared-but-absent, None if
    the header does not declare it at all."""
    if key not in declared:
        return None
    try:
        return key in rec.info
    except (ValueError, KeyError):
        return None


def _contig_index(chrom: str) -> int:
    """Ordering key for canonical junction orientation. Unknown contigs sort
    last but deterministically, so cross-caller comparison stays stable."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    if c.isdigit():
        return int(c)
    return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(c, 1000 + sum(ord(x) for x in c))


def _flip_orientation(ct: Optional[str]) -> Optional[str]:
    """Orientation as seen from the other breakend. 3to3 and 5to5 are
    symmetric; 3to5 and 5to3 exchange."""
    if ct is None:
        return None
    return {"3to5": "5to3", "5to3": "3to5", "3to3": "3to3", "5to5": "5to5"}.get(ct, ct)


def parse_alt_breakend(alt: str):
    """(partner_chrom, partner_pos, orientation) from VCF breakend notation,
    or None if `alt` is not breakend notation (e.g. a symbolic <BND>)."""
    m = _ALT_BND.match(alt or "")
    if not m:
        return None
    ob, cb = m.group("ob"), m.group("cb")
    if ob != cb:
        return None
    pre, post = m.group("pre"), m.group("post")
    if pre and not post:
        ct = "3to5" if ob == "[" else "3to3"
    elif post and not pre:
        ct = "5to3" if ob == "]" else "5to5"
    else:
        return None
    return m.group("chrom"), int(m.group("pos")), ct


@dataclass
class Junction:
    """One SV junction, normalised across callers.

    Canonical ordering matches DELLY's: the locus with the HIGHER
    (contig_index, position) is placed first, so the same junction from two
    callers normalises to the same pair regardless of which end each caller
    chose to anchor on.
    """
    candidate_id: str          # opaque, stable within and across sessions for
                               # the same normalised junction
    caller_id: str             # the ID the caller assigned, for traceability
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    orientation: Optional[str] # "3to5" | "5to3" | "3to3" | "5to5" | None;
                               # DELLY also emits "NtoN" on INS records
    svtype: str
    filter: str                # "PASS" or the caller's own filter tag
    pe: Optional[int]          # None = the caller did not report the field,
                               # which is NOT the same as 0 support
    sr: Optional[int]          # None = field absent, not zero
    precise: Optional[bool]    # None = neither PRECISE nor IMPRECISE declared
    ci_pos1: Optional[list]    # [low, high] confidence interval, caller's own
    ci_pos2: Optional[list]
    n_merged: int = 1          # how many caller records this junction absorbed
                               # during deduplication; 1 = no merge


def _candidate_id(chrom1, pos1, chrom2, pos2, orientation) -> str:
    key = f"{chrom1}:{pos1}|{chrom2}:{pos2}|{orientation}"
    return "CAND_" + _hashlib.sha256(key.encode("utf-8")).hexdigest()[:12]


def _canonicalise(c1, p1, c2, p2, ct, ci1, ci2):
    """Put the higher (contig_index, position) locus first, flipping
    orientation and confidence intervals with it."""
    if (_contig_index(c1), p1) >= (_contig_index(c2), p2):
        return c1, p1, c2, p2, ct, ci1, ci2
    return c2, p2, c1, p1, _flip_orientation(ct), ci2, ci1


def _int_or_none(v):
    if v is None:
        return None
    if isinstance(v, (tuple, list)):
        v = v[0] if v else None
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def detect_convention(vcf_path: str, scan_limit: int = 5000) -> dict:
    """Which breakend convention this file uses, from header AND records.

    Returns {"convention": ..., "evidence": {...}} where convention is one of
      "chr2_pos2"    single record carrying INFO/CHR2 + INFO/POS2 (DELLY)
      "mateid"       paired records linked by INFO/MATEID (manta, VCF spec)
      "undetermined" no BND/translocation records to judge from (e.g. an
                     empty file, or one with only DEL/DUP/INV/INS)
    or {"error": ..., "error_type": ...} if BND records exist but use neither.
    """
    try:
        vf = pysam.VariantFile(vcf_path)
    except Exception as e:
        return {"error": _scrub(e, vcf_path), "error_type": "vcf_access"}

    hdr_info = set(vf.header.info.keys())
    header_mateid = "MATEID" in hdr_info
    header_chr2 = "CHR2" in hdr_info and "POS2" in hdr_info

    rec_mateid = rec_chr2 = rec_altbnd = bnd_records = 0
    scanned = 0
    try:
        for rec in vf.fetch() if vf.index is not None else vf:
            scanned += 1
            if scanned > scan_limit:
                break
            svtype = _info(rec, "SVTYPE", hdr_info)
            if svtype not in ("BND", "TRA"):
                continue
            bnd_records += 1
            if _info(rec, "MATEID", hdr_info):
                rec_mateid += 1
            if _info(rec, "CHR2", hdr_info) is not None and _info(rec, "POS2", hdr_info) is not None:
                rec_chr2 += 1
            alt = rec.alts[0] if rec.alts else ""
            if parse_alt_breakend(alt):
                rec_altbnd += 1
    except Exception as e:
        vf.close()
        return {"error": "failed reading records: " + _scrub(e, vcf_path),
                "error_type": "vcf_parse"}
    vf.close()

    evidence = {
        "header_declares_mateid": header_mateid,
        "header_declares_chr2_pos2": header_chr2,
        "records_scanned": min(scanned, scan_limit),
        "bnd_records_seen": bnd_records,
        "records_with_mateid": rec_mateid,
        "records_with_chr2_pos2": rec_chr2,
        "records_with_alt_breakend_notation": rec_altbnd,
    }

    # Records decide, not the header: a header may declare a field the file
    # never uses. MATEID wins when present because it implies paired records
    # that must be collapsed, which changes the junction count.
    if rec_mateid > 0:
        return {"convention": "mateid", "evidence": evidence}
    if rec_chr2 > 0:
        return {"convention": "chr2_pos2", "evidence": evidence}
    if bnd_records == 0:
        return {"convention": "undetermined", "evidence": evidence}
    return {
        "error": (
            f"{bnd_records} BND/TRA record(s) present but none carries "
            f"INFO/MATEID or INFO/CHR2+INFO/POS2. This file uses a breakend "
            f"convention this parser does not implement; refusing to guess."
        ),
        "error_type": "unknown_breakend_convention",
        "evidence": evidence,
    }


def parse_junctions(vcf_path: str, dedup_tolerance_bp: int = DEDUP_TOLERANCE_BP) -> dict:
    """Parse a VCF/BCF into normalised junctions.

    Returns {"junctions": [Junction...], "convention": str, "meta": {...}}
    or a flat error dict. Sample columns are never read.
    """
    det = detect_convention(vcf_path)
    if "error" in det:
        return det
    convention = det["convention"]

    try:
        vf = pysam.VariantFile(vcf_path)
    except Exception as e:
        return {"error": _scrub(e, vcf_path), "error_type": "vcf_access"}

    declared = set(vf.header.info.keys())
    total_records = 0
    by_svtype = {}
    by_filter = {}
    raw = []          # pre-dedup junctions
    mate_pool = {}    # id -> record fields, for the mateid convention
    malformed = 0

    try:
        for rec in vf:
            total_records += 1
            svtype = _info(rec, "SVTYPE", declared) or "UNKNOWN"
            by_svtype[svtype] = by_svtype.get(svtype, 0) + 1
            filt = ";".join(rec.filter.keys()) if len(rec.filter.keys()) else "PASS"
            by_filter[filt] = by_filter.get(filt, 0) + 1

            pe = _int_or_none(_info(rec, "PE", declared))
            sr = _int_or_none(_info(rec, "SR", declared))
            _p, _i = _flag(rec, "PRECISE", declared), _flag(rec, "IMPRECISE", declared)
            precise = True if _p else (False if _i else None)
            ci = _info(rec, "CIPOS", declared)
            ci = [int(ci[0]), int(ci[1])] if ci and len(ci) >= 2 else None
            ciend = _info(rec, "CIEND", declared)
            ciend = [int(ciend[0]), int(ciend[1])] if ciend and len(ciend) >= 2 else None
            alt = rec.alts[0] if rec.alts else ""

            if svtype in ("BND", "TRA"):
                if convention == "mateid":
                    mid = _info(rec, "MATEID", declared)
                    mid = mid[0] if isinstance(mid, (tuple, list)) and mid else mid
                    parsed = parse_alt_breakend(alt)
                    if not parsed:
                        malformed += 1
                        continue
                    pc, pp, ct = parsed
                    mate_pool[rec.id] = dict(
                        chrom=rec.chrom, pos=rec.pos, pchrom=pc, ppos=pp, ct=ct,
                        mate=mid, filt=filt, pe=pe, sr=sr, precise=precise,
                        ci=ci, ciend=ciend, cid=rec.id,
                    )
                    continue
                # chr2_pos2
                c2 = _info(rec, "CHR2", declared)
                p2 = _int_or_none(_info(rec, "POS2", declared))
                ct = _info(rec, "CT", declared)
                if c2 is None or p2 is None:
                    parsed = parse_alt_breakend(alt)
                    if not parsed:
                        malformed += 1
                        continue
                    c2, p2, ct = parsed
                c1, p1, c2, p2, ct, ci1, ci2 = _canonicalise(
                    rec.chrom, rec.pos, c2, p2, ct, ci, ciend)
                raw.append(Junction(
                    candidate_id=_candidate_id(c1, p1, c2, p2, ct),
                    caller_id=rec.id or "", chrom1=c1, pos1=p1, chrom2=c2, pos2=p2,
                    orientation=ct, svtype=svtype, filter=filt, pe=pe, sr=sr,
                    precise=precise, ci_pos1=ci1, ci_pos2=ci2))
            else:
                # Intra-chromosomal event: second breakend is END on the same contig.
                end = _int_or_none(_info(rec, "END", declared)) or rec.stop
                c1, p1, c2, p2, ct, ci1, ci2 = _canonicalise(
                    rec.chrom, rec.pos, rec.chrom, int(end), _info(rec, "CT", declared), ci, ciend)
                raw.append(Junction(
                    candidate_id=_candidate_id(c1, p1, c2, p2, ct),
                    caller_id=rec.id or "", chrom1=c1, pos1=p1, chrom2=c2, pos2=p2,
                    orientation=ct, svtype=svtype, filter=filt, pe=pe, sr=sr,
                    precise=precise, ci_pos1=ci1, ci_pos2=ci2))
    except Exception as e:
        vf.close()
        return {"error": "failed parsing records: " + _scrub(e, vcf_path),
                "error_type": "vcf_parse"}
    vf.close()

    unpaired_mates = 0
    if convention == "mateid":
        seen = set()
        for rid, r in mate_pool.items():
            if rid in seen:
                continue
            mate = mate_pool.get(r["mate"])
            if mate is None:
                # Mate absent from the file. Emit the junction anyway using the
                # ALT partner locus -- dropping it would silently lose a call.
                unpaired_mates += 1
                c1, p1, c2, p2, ct, ci1, ci2 = _canonicalise(
                    r["chrom"], r["pos"], r["pchrom"], r["ppos"], r["ct"], r["ci"], r["ciend"])
                seen.add(rid)
            else:
                seen.add(rid); seen.add(r["mate"])
                c1, p1, c2, p2, ct, ci1, ci2 = _canonicalise(
                    r["chrom"], r["pos"], mate["chrom"], mate["pos"], r["ct"], r["ci"], mate["ci"])
            pe = max([x for x in (r["pe"], (mate or {}).get("pe")) if x is not None], default=None)
            sr = max([x for x in (r["sr"], (mate or {}).get("sr")) if x is not None], default=None)
            raw.append(Junction(
                candidate_id=_candidate_id(c1, p1, c2, p2, ct),
                caller_id=r["cid"] or "", chrom1=c1, pos1=p1, chrom2=c2, pos2=p2,
                orientation=ct, svtype="BND", filter=r["filt"], pe=pe, sr=sr,
                precise=r["precise"], ci_pos1=ci1, ci_pos2=ci2))

    junctions, merged = _dedup(raw, dedup_tolerance_bp)
    return {
        "junctions": junctions,
        "convention": convention,
        "meta": {
            "total_records": total_records,
            "by_svtype": by_svtype,
            "by_filter": by_filter,
            "junctions_before_dedup": len(raw),
            "junctions_after_dedup": len(junctions),
            "records_merged_by_dedup": merged,
            "malformed_records_skipped": malformed,
            "unpaired_mate_records": unpaired_mates,
            "detection_evidence": det["evidence"],
        },
    }


def _dedup(raw, tol):
    """Merge junctions whose BOTH breakends fall within `tol` AND whose
    orientation is identical. Representative is PASS-preferred, then highest
    PE+SR -- the same rule used in the Phase 3.5 verification, restated here
    so the two are comparable.

    Orientation is part of the merge key because the two junctions of a
    balanced reciprocal translocation sit ~1bp apart on BOTH breakends and
    differ ONLY in orientation. Keying on coordinates alone collapsed them
    into one candidate and discarded the second's PE/SR evidence, which made
    the reciprocal structure undetectable and capped any per-junction
    sensitivity figure at half its true denominator.

    The key is safe across both parsing conventions: delly's INFO/CT and the
    ALT breakend notation decoded by parse_alt_breakend both yield the same
    vocabulary (3to5/5to3/3to3/5to5, plus delly's NtoN on INS), and
    _canonicalise flips it consistently when it reorders the two ends, so the
    same physical junction normalises to the same orientation either way.
    Verified by parsing the same translocation in both conventions and
    comparing the normalised output."""
    if tol <= 0 or not raw:
        return list(raw), 0
    parent = list(range(len(raw)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    buckets = {}
    for i, j in enumerate(raw):
        buckets.setdefault((j.chrom1, j.chrom2, j.orientation, j.pos1 // 100000), []).append(i)
    for i, j in enumerate(raw):
        for b in (j.pos1 // 100000 - 1, j.pos1 // 100000, j.pos1 // 100000 + 1):
            for k in buckets.get((j.chrom1, j.chrom2, j.orientation, b), ()):
                if k != i and abs(raw[k].pos1 - j.pos1) <= tol and abs(raw[k].pos2 - j.pos2) <= tol:
                    union(i, k)
    groups = {}
    for i in range(len(raw)):
        groups.setdefault(find(i), []).append(i)
    out = []
    for members in groups.values():
        best = max(members, key=lambda i: (
            raw[i].filter == "PASS",
            (raw[i].pe or 0) + (raw[i].sr or 0),
            raw[i].sr or 0,
        ))
        j = raw[best]
        out.append(Junction(**{**asdict(j), "n_merged": len(members)}))
    return out, len(raw) - len(out)


# ── mask handling ───────────────────────────────────────────────────────────

def load_mask(mask_path: str) -> dict:
    """Region masks from a caller exclude template (chrom, start, end, label).
    One-column whole-contig lines are recorded separately."""
    regions, whole = {}, set()
    try:
        with open(mask_path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) == 1 and f[0].strip():
                    whole.add(f[0])
                elif len(f) >= 3:
                    regions.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    except Exception as e:
        return {"error": _scrub(e, mask_path), "error_type": "mask_access"}
    return {"regions": regions, "whole_contigs": whole}


def in_mask(mask: dict, chrom: str, pos: int) -> bool:
    if chrom in mask.get("whole_contigs", ()):
        return True
    return any(a <= pos < b for a, b in mask.get("regions", {}).get(chrom, ()))


# ── session registry ────────────────────────────────────────────────────────
#
# set_id -> {"label", "path", "junctions", "convention", "meta", "by_id"}
# The path is held here and NEVER returned. Same contract as the image handle
# manifest in bam_tools: the model gets an opaque reference, not a filename.

_SETS = {}


def _set_id(path: str, label: str) -> str:
    return "CSET_" + _hashlib.sha256(
        (_os.path.abspath(path) + "|" + label).encode("utf-8")).hexdigest()[:8]


def reset_registry():
    """Drop all registered sets and recorded coordinates. For tests."""
    _SETS.clear()
    _COORD_PROVENANCE.clear()


def load_candidate_set(path: str, label: str, dedup_tolerance_bp: int = DEDUP_TOLERANCE_BP) -> dict:
    """Register a candidate VCF/BCF for this session.

    Returns aggregate counts only: no coordinates, no sample name, no path.
    """
    if not isinstance(label, str) or not label.strip():
        return {"error": "label must be a non-empty string", "error_type": "bad_parameter"}
    if not _os.path.exists(path):
        return {"error": "candidate file not found", "error_type": "vcf_access"}

    parsed = parse_junctions(path, dedup_tolerance_bp)
    if "error" in parsed:
        return parsed

    sid = _set_id(path, label)
    _SETS[sid] = {
        "label": label,
        "path": _os.path.abspath(path),      # withheld from every return value
        "junctions": parsed["junctions"],
        "convention": parsed["convention"],
        "meta": parsed["meta"],
        "by_id": {j.candidate_id: j for j in parsed["junctions"]},
    }
    m = parsed["meta"]
    return {
        "set_id": sid,
        "label": label,
        "caller_convention": parsed["convention"],
        "convention_evidence": m["detection_evidence"],
        "total_records": m["total_records"],
        "counts_by_svtype": m["by_svtype"],
        "counts_by_filter": m["by_filter"],
        "junctions_before_dedup": m["junctions_before_dedup"],
        "junctions_after_dedup": m["junctions_after_dedup"],
        "records_merged_by_dedup": m["records_merged_by_dedup"],
        "malformed_records_skipped": m["malformed_records_skipped"],
        "unpaired_mate_records": m["unpaired_mate_records"],
        "thresholds_applied": [
            {"name": "dedup_tolerance_bp", "value": dedup_tolerance_bp,
             "provenance": PROVENANCE["dedup_tolerance_bp"]},
        ],
        "note": (
            "Counts describe the whole registered set. No coordinates, sample "
            "name or file path are returned by this tool."
        ),
    }


def _apply_filters(junctions, svtype, filter_pass, min_pe, min_sr,
                   primary_only, mask):
    """Returns (surviving, per_step_records) where per_step_records lists each
    filter with cumulative survivors AND its standalone effect on the
    unfiltered set. Both are reported because a cumulative count can read as
    inert purely because an earlier step already removed the affected records."""
    steps = []
    cur = list(junctions)
    # Standalone counts are measured against the SAME svtype scope as the
    # query. Previously they were measured against the whole registered set,
    # so a BND query reported a mask figure that counted DEL/DUP/INV records
    # too -- two correct numbers measuring different things, in adjacent
    # columns. Scope now follows the query.
    scope = [j for j in junctions
             if svtype is None or j.svtype == svtype]
    total = len(scope)

    def standalone(pred, against=None):
        pool = scope if against is None else against
        return sum(1 for j in pool if not pred(j))

    specs = []
    if svtype is not None:
        specs.append(("svtype", svtype, PROVENANCE["svtype"], lambda j: j.svtype == svtype))
    if filter_pass:
        specs.append(("filter_pass", True, PROVENANCE["filter_pass"], lambda j: j.filter == "PASS"))
    if min_pe is not None:
        specs.append(("min_pe", min_pe, PROVENANCE["min_pe"], lambda j: (j.pe or 0) >= min_pe))
    if min_sr is not None:
        specs.append(("min_sr", min_sr, PROVENANCE["min_sr"], lambda j: (j.sr or 0) >= min_sr))
    if primary_only:
        specs.append(("primary_only", True, PROVENANCE["primary_only"],
                      lambda j: j.chrom1 in PRIMARY_CONTIGS and j.chrom2 in PRIMARY_CONTIGS))
    if mask is not None:
        specs.append(("exclude_masked", True, PROVENANCE["exclude_masked"],
                      lambda j: not in_mask(mask, j.chrom1, j.pos1)
                                and not in_mask(mask, j.chrom2, j.pos2)))

    for name, value, prov, pred in specs:
        before = len(cur)
        cur = [j for j in cur if pred(j)]
        # The svtype filter itself is scoped against the whole set, since it is
        # what defines the scope every later filter is measured against.
        alone = standalone(pred, against=junctions) if name == "svtype" else standalone(pred)
        steps.append({
            "filter": name,
            "value": value,
            "provenance": prov,
            "surviving_after_this_step": len(cur),
            "removed_cumulatively_here": before - len(cur),
            "would_remove_from_unfiltered_set": alone,
            "unfiltered_set_size": len(junctions) if name == "svtype" else total,
            "unfiltered_set_scope": "all svtypes" if name == "svtype" else (
                f"svtype={svtype}" if svtype else "all svtypes"),
        })
    return cur, steps


def list_candidates(set_id: str, svtype: Optional[str] = None,
                    filter_pass: bool = False,
                    min_pe: Optional[int] = None, min_sr: Optional[int] = None,
                    primary_only: bool = False, mask_path: Optional[str] = None,
                    limit: int = 50, offset: int = 0) -> dict:
    """Candidates matching the given filters, with every threshold echoed back.

    The response always reports total_in_set alongside total_matching, so a
    filtered list can never be presented as if it were the whole set.
    """
    s = _SETS.get(set_id)
    if s is None:
        return {"error": f"unknown set_id '{set_id}'. Call load_candidate_set first.",
                "error_type": "unknown_set_id"}
    if not isinstance(limit, int) or limit < 1:
        return {"error": f"limit must be a positive integer (got {limit!r})",
                "error_type": "bad_parameter"}
    if not isinstance(offset, int) or offset < 0:
        return {"error": f"offset must be a non-negative integer (got {offset!r})",
                "error_type": "bad_parameter"}

    mask = None
    if mask_path is not None:
        mask = load_mask(mask_path)
        if "error" in mask:
            return mask

    matching, steps = _apply_filters(
        s["junctions"], svtype, filter_pass, min_pe, min_sr, primary_only, mask)
    page = matching[offset:offset + limit]
    return {
        "set_id": set_id,
        "label": s["label"],
        "total_in_set": len(s["junctions"]),
        "total_matching": len(matching),
        "returned": len(page),
        "offset": offset,
        "limit": limit,
        "truncated": offset + len(page) < len(matching),
        "filters_applied": steps,
        "candidates": [asdict(j) for j in page],
        "note": (
            "total_matching is the count AFTER the filters listed in "
            "filters_applied; total_in_set is the whole registered set. "
            "Quoting a candidate count without its filters is not reproducible."
        ),
    }


def get_candidate(set_id: str, candidate_id: str) -> dict:
    """Full detail for one junction, shaped so the Stage 1 evidence tools can
    consume its two breakends directly."""
    s = _SETS.get(set_id)
    if s is None:
        return {"error": f"unknown set_id '{set_id}'. Call load_candidate_set first.",
                "error_type": "unknown_set_id"}
    j = s["by_id"].get(candidate_id)
    if j is None:
        return {"error": f"unknown candidate_id '{candidate_id}' in set '{set_id}'",
                "error_type": "unknown_candidate_id"}
    _record_coordinate(j.chrom1, j.pos1, j.candidate_id, set_id, s["label"])
    _record_coordinate(j.chrom2, j.pos2, j.candidate_id, set_id, s["label"])
    d = asdict(j)
    d.update({
        "set_id": set_id,
        "label": s["label"],
        "caller_convention": s["convention"],
        # Ready-to-use arguments for the Stage 1 evidence tools. Supplying them
        # explicitly means the assistant never has to retype a coordinate.
        "breakend_1": {"chromosome": j.chrom1, "position": j.pos1},
        "breakend_2": {"chromosome": j.chrom2, "position": j.pos2},
        "note": (
            "Pass breakend_1 / breakend_2 straight to the evidence tools. "
            "Coordinates originate from the candidate file, not from the model."
        ),
    })
    return d


def compare_candidate_sets(set_a: str, set_b: str,
                           tolerance_bp: int = RECURRENCE_TOLERANCE_BP) -> dict:
    """Junction-level recurrence between two registered sets.

    Serves two purposes: two-sample recurrence (an artifact indicator) and,
    once a second caller is added, two-caller concordance on one sample.
    """
    a, b = _SETS.get(set_a), _SETS.get(set_b)
    if a is None or b is None:
        missing = [x for x, y in ((set_a, a), (set_b, b)) if y is None]
        return {"error": f"unknown set_id(s): {missing}", "error_type": "unknown_set_id"}
    if not isinstance(tolerance_bp, int) or tolerance_bp < 0:
        return {"error": f"tolerance_bp must be a non-negative integer (got {tolerance_bp!r})",
                "error_type": "bad_parameter"}

    # Every matching pair is emitted, not just the first partner found for each
    # A junction. The previous version broke out of the scan on first match, so
    # matched_in_b counted only those B junctions that happened to be selected
    # as some A junction's first partner -- a B junction with a genuine partner
    # in A could be reported unmatched purely because a neighbour was reached
    # first. That made the comparison direction-dependent: asking "does A recur
    # in B" and "does B recur in A" gave different answers for the same pair of
    # sets. Measured on two sets built from an identical background, the true
    # relation was 894 matched / 2 unmatched in BOTH directions, while the tool
    # reported 894/2 one way and 840/56 the other; all 54 of the difference
    # were junctions with two or three genuine partners. Recurrence gates a
    # funnel step, so a direction-dependent answer there is not acceptable.
    idx = {}
    for j in b["junctions"]:
        idx.setdefault((j.chrom1, j.chrom2, j.pos1 // 100000), []).append(j)
    matches = []
    matched_a = set()
    matched_b = set()
    for j in a["junctions"]:
        for bkt in (j.pos1 // 100000 - 1, j.pos1 // 100000, j.pos1 // 100000 + 1):
            for k in idx.get((j.chrom1, j.chrom2, bkt), ()):
                if abs(k.pos1 - j.pos1) <= tolerance_bp and abs(k.pos2 - j.pos2) <= tolerance_bp:
                    matches.append({"candidate_id_a": j.candidate_id,
                                    "candidate_id_b": k.candidate_id})
                    matched_a.add(j.candidate_id)
                    matched_b.add(k.candidate_id)
    return {
        "set_a": set_a, "label_a": a["label"],
        "set_b": set_b, "label_b": b["label"],
        "total_in_a": len(a["junctions"]),
        "total_in_b": len(b["junctions"]),
        "matched_in_a": len(matched_a),
        "matched_in_b": len(matched_b),
        "unmatched_in_a": len(a["junctions"]) - len(matched_a),
        "unmatched_in_b": len(b["junctions"]) - len(matched_b),
        "matched_pairs": matches,
        "thresholds_applied": [
            {"name": "tolerance_bp", "value": tolerance_bp,
             "provenance": PROVENANCE["recurrence_tolerance_bp"]},
        ],
        "note": (
            "A match means both breakends agree within tolerance_bp. "
            "Recurrence across unrelated samples is an artifact indicator, "
            "not a finding; this tool reports counts and does not interpret them. "
            "Matching is exhaustive and symmetric: every matching pair is listed, "
            "so swapping set_a and set_b mirrors the counts rather than changing them."
        ),
    }

# ── coordinate provenance ───────────────────────────────────────────────────
#
# The eleven evidence tools take chromosome/position as free parameters, so a
# model can invent a coordinate and present the result as though it confirmed a
# candidate. That capability is deliberately KEPT -- a clinician typing a
# breakpoint from a karyotype is a primary use case. What must not happen is
# the two becoming indistinguishable in the output.
#
# get_candidate records every coordinate it hands out here. The evidence tools
# consult it and label their result accordingly. A coordinate is "candidate_set"
# only on an EXACT match to one a tool returned this session; anything else is
# "caller_supplied", including a near-miss. Conservative on purpose: the field
# exists to make an unverified coordinate visible, so it must not launder one.

_COORD_PROVENANCE = {}   # (chrom, pos) -> {"candidate_id", "set_id", "label"}


def _record_coordinate(chrom, pos, candidate_id, set_id, label):
    _COORD_PROVENANCE[(str(chrom), int(pos))] = {
        "candidate_id": candidate_id, "set_id": set_id, "label": label}


def position_provenance(chromosome, position=None, start=None, end=None) -> dict:
    """Where did this queried position come from?

    Returns source "candidate_set" with the candidate_id when the coordinate
    exactly matches one that get_candidate returned in this session, otherwise
    "caller_supplied" with candidate_id None. Never raises.
    """
    try:
        c = str(chromosome)
        hits = []
        if position is not None:
            h = _COORD_PROVENANCE.get((c, int(position)))
            if h:
                hits.append(h)
        if not hits and start is not None and end is not None:
            lo, hi = int(start), int(end)
            for (kc, kp), v in _COORD_PROVENANCE.items():
                if kc == c and lo <= kp <= hi:
                    hits.append(v); break
        if hits:
            h = hits[0]
            return {
                "source": "candidate_set",
                "candidate_id": h["candidate_id"],
                "set_id": h["set_id"],
                "set_label": h["label"],
                "note": ("This position was returned by get_candidate in this "
                         "session and originates from the SV caller's output."),
            }
    except Exception:
        pass
    return {
        "source": "caller_supplied",
        "candidate_id": None,
        "set_id": None,
        "set_label": None,
        "note": ("This position was supplied directly in the tool call. It did "
                 "NOT come from a registered candidate set, so nothing here "
                 "corroborates that the position is a real candidate. Do not "
                 "present this result as confirmation of a caller-derived "
                 "junction."),
    }
