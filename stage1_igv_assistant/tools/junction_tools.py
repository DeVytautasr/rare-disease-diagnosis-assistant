"""
junction_tools.py -- the evidence for one junction, read by read (Phase 26).

The supervisor's review of the interface (2 October 2026) asked for evidence a
clinical geneticist reads without a scoring system: how many read pairs join
the two ends of a junction, how many reads are split across it, the reads
themselves (only the abnormal ones, colour-coded, at both ends), where the
piece that moved starts and ends, and which genes it touches. This module
provides the read-level part and the chromosome-level description. It is new
rather than an edit of bam_tools.py, which stays as it was: its tools still
feed the combined score that the thesis analyses and the previous page shows.

Every figure here is a count of reads or a coordinate. No score is made.

Thresholds, each echoed in the result that applied it:

  WINDOW_BP = 1000        author judgement. A read counts for an end when it
                          starts within this distance of the breakpoint: a
                          fragment of up to ~1 kb can cross the junction
                          (short-read libraries are typically 300-600 bp).
  MIN_MAPQ = 20           the project's mapping-quality floor (every counting
                          layer in bam_tools uses 20). BOTH reads of a pair and
                          BOTH pieces of a split read must reach it to count.
                          Reads below it are still drawn, hollow, and counted
                          apart: in repeats the evidence can be in them.
  MIN_CLIP_BASES = 10     as bam_tools' soft-clip layer: a read is "cut" at a
                          base when at least this many bases are clipped there.
  CLIP_SEARCH_BP = 300    author judgement. The base where clipped reads pile
                          up is looked for within this distance of the caller's
                          position, so an unrelated pile-up elsewhere in the
                          window cannot be reported as this junction's.
  GROUP_NEAR_BP = 10000   author judgement. Two junctions joining the same two
                          chromosomes belong to one rearrangement when their
                          ends lie within this distance on both chromosomes
                          (the two junctions of a balanced translocation are
                          usually a few bp to a few kb apart).
  SEGMENT_MAX_BP = 10000000  author judgement. Two such junctions whose ends
                          are near on one chromosome but up to this far apart
                          on the other bound a segment of the other chromosome
                          (an insertion-like or unbalanced rearrangement).
  LIST_LIMIT = 100        reads listed per kind; counts are never capped.
  VIEW_LIMIT = 400        abnormal reads returned per end for the read view.

Orientation uses the vocabulary DELLY writes in INFO/CT and vcf_tools derives
from breakend notation, always written from the first end to the second:
"3" = the end's left part (the chromosome from its start up to the breakpoint)
takes part in the join, "5" = its right part (from the breakpoint to the end).
  3to5  first-left  + second-right
  5to3  second-left + first-right
  3to3  first-left  + second-left, the second piece reversed
  5to5  second-right reversed + first-right
For a read pair the side is read from the strand of each read (a read on the +
strand points right, so the part left of the breakpoint takes part in the join);
for a split read, from which end of each aligned piece lies at the breakpoint.
"""

import hashlib as _hashlib
import os as _os
import re as _re
from collections import Counter as _Counter

import pysam

from stage1_igv_assistant.tools import bam_tools as _bt

WINDOW_BP = 1000
MIN_MAPQ = 20
MIN_CLIP_BASES = 10
CLIP_SEARCH_BP = 300
GROUP_NEAR_BP = 10_000
SEGMENT_MAX_BP = 10_000_000
LIST_LIMIT = 100
VIEW_LIMIT = 400

PROVENANCE = {
    "window_bp": "author judgement",
    "min_mapq": "project floor (as bam_tools)",
    "min_clip_bases": "as bam_tools' soft-clip layer",
    "clip_search_bp": "author judgement",
    "group_near_bp": "author judgement",
    "segment_max_bp": "author judgement",
}

ORIENTATIONS = ("3to5", "5to3", "3to3", "5to5")
_FLIP = {"3to5": "5to3", "5to3": "3to5", "3to3": "3to3", "5to5": "5to5"}
_CIGAR_OP = _re.compile(r"(\d+)([MIDNSHP=X])")


# ── names and genome order ──────────────────────────────────────────────────
def norm_chrom(name) -> str:
    """Spelling-insensitive chromosome key: 'chr7', '7' -> '7'; 'chrMT' -> 'M'."""
    c = str(name or "")
    if c[:3].lower() == "chr":
        c = c[3:]
    if c.upper() in ("M", "MT"):
        return "M"
    if c.upper() in ("X", "Y"):
        return c.upper()
    return c


def display_chrom(name) -> str:
    """'7' and 'chr7' are both shown as chr7; other contigs keep their name."""
    n = norm_chrom(name)
    if n.isdigit() or n in ("X", "Y", "M"):
        return "chr" + n
    return str(name)


def chrom_rank(name) -> int:
    n = norm_chrom(name)
    if n.isdigit():
        return int(n)
    return {"X": 23, "Y": 24, "M": 25}.get(n, 1000 + sum(ord(x) for x in n))


def flip(orientation):
    return _FLIP.get(orientation, orientation)


def genome_order(c1, p1, c2, p2, orientation):
    """The two ends in genome order (lower chromosome, then lower position, first),
    with the orientation re-written from the new first end."""
    if (chrom_rank(c1), p1) <= (chrom_rank(c2), p2):
        return c1, p1, c2, p2, orientation
    return c2, p2, c1, p1, flip(orientation)


# ── cytogenetic bands (GRCh38) ──────────────────────────────────────────────
_BANDS_FILE = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                            "data", "grch38_cytobands.tsv")
_BANDS = None


def bands():
    """{norm_chrom: [(start0, end, band, stain), ...]} from the bundled GRCh38 table."""
    global _BANDS
    if _BANDS is None:
        out = {}
        try:
            with open(_BANDS_FILE) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    c, s, e, b, st = line.rstrip("\n").split("\t")[:5]
                    out.setdefault(norm_chrom(c), []).append((int(s), int(e), b, st))
        except OSError:
            out = {}
        _BANDS = out
    return _BANDS


def chrom_length(chrom):
    b = bands().get(norm_chrom(chrom))
    return b[-1][1] if b else None


def band_at(chrom, pos):
    """('p12', 'gneg') for a 1-based position, or (None, None)."""
    for s, e, b, st in bands().get(norm_chrom(chrom), ()):
        if s < pos <= e:
            return b, st
    return None, None


def centromere(chrom):
    """(start, end) of the acen bands, 1-based inclusive, or None."""
    ac = [(s + 1, e) for s, e, b, st in bands().get(norm_chrom(chrom), ()) if st == "acen"]
    if not ac:
        return None
    return min(a for a, _ in ac), max(b for _, b in ac)


def region_note(chrom, pos):
    """A plain warning when a breakpoint lies where short reads mislead most often:
    in a centromere (acen), in variable heterochromatin (gvar) or in an acrocentric
    stalk. None elsewhere."""
    b, st = band_at(chrom, pos)
    if st == "acen":
        return f"in the centromere ({display_chrom(chrom)}{b})"
    if st == "gvar":
        return f"in variable heterochromatin ({display_chrom(chrom)}{b})"
    if st == "stalk":
        return f"in an acrocentric stalk ({display_chrom(chrom)}{b})"
    return None


# ── reading one end ─────────────────────────────────────────────────────────
def _ref_len(cigar: str) -> int:
    return sum(int(n) for n, op in _CIGAR_OP.findall(cigar or "") if op in "MDN=X")


def _clips(cigar: str):
    ops = _CIGAR_OP.findall(cigar or "")
    left = int(ops[0][0]) if ops and ops[0][1] in "SH" else 0
    right = int(ops[-1][0]) if ops and ops[-1][1] in "SH" else 0
    return left, right


def _side_of_piece(start1, end1, breakpoint):
    """'3' when the piece's end lies at the breakpoint (the part left of it takes
    part in the join), '5' when its start does."""
    return "3" if abs(end1 - breakpoint) <= abs(start1 - breakpoint) else "5"


def _nearer_here(here, other, pos, partner_pos):
    """Within one chromosome, whether a read (or piece) at `here` lies nearer this
    end than its mate (or other piece) at `other` does, each measured against both
    ends. When the two ends lie less than two windows apart, both reads of a pair
    fall in both windows; each end keeps the one nearer to itself."""
    return abs(here - pos) - abs(here - partner_pos) < abs(other - pos) - abs(other - partner_pos)


def _parse_sa(tag):
    out = []
    for entry in str(tag or "").rstrip(";").split(";"):
        f = entry.split(",")
        if len(f) < 5:
            continue
        try:
            pos, mapq = int(f[1]), int(f[4])
        except ValueError:
            continue
        out.append({"chrom": f[0], "pos": pos, "strand": f[2], "cigar": f[3], "mapq": mapq})
    return out


def _read_record(r):
    ct = r.cigartuples or []
    left = ct[0][1] if ct and ct[0][0] in (4, 5) else 0
    right = ct[-1][1] if ct and ct[-1][0] in (4, 5) else 0
    return {"name": r.query_name, "start": r.reference_start + 1, "end": r.reference_end,
            "strand": "-" if r.is_reverse else "+", "cigar": r.cigarstring, "mapq": r.mapping_quality,
            "clip_left": left, "clip_right": right, "supplementary": bool(r.is_supplementary)}


def scan_end(bam, chrom, pos, partner_chrom, partner_pos, window_bp=WINDOW_BP,
             min_mapq=MIN_MAPQ, min_clip=MIN_CLIP_BASES, related=()):
    """Every abnormal read near one end, each classified.

    Returns {"error"} or a dict with:
      pairs   {name: record}  primary reads (MAPQ >= min) whose mate starts within
                              window_bp of the partner end
      split   {name: record}  alignments (MAPQ >= min) with a supplementary piece
                              (SA, MAPQ >= min) within window_bp of the partner end
                              (within one chromosome, both only for the read or
                              piece nearer this end: see _nearer_here)
      view    [record]        every abnormal read in the window, for drawing
      clips, elsewhere, counts
    """
    resolved = _bt._resolve_contig(bam, chrom)
    if resolved is None:
        return {"error": f"chromosome {chrom!r} is not in this file's header",
                "error_type": "invalid_region"}
    length = bam.get_reference_length(resolved)
    if not (1 <= pos <= length):
        return {"error": f"{display_chrom(chrom)}:{pos:,} lies outside the chromosome (length {length:,})",
                "error_type": "invalid_region"}
    here_n, there_n = norm_chrom(chrom), norm_chrom(partner_chrom)
    intra = here_n == there_n
    # other ends of the same rearrangement (the sibling junction of a segment): a
    # mate there belongs to the rearrangement, not to "other places"
    related = [(norm_chrom(c), p) for c, p in related]
    lo, hi = max(0, pos - window_bp), min(length, pos + window_bp)
    pairs, split = {}, {}
    view_partner, view_other = [], []
    kind_counts = _Counter()
    clip_starts, clip_ends = _Counter(), _Counter()
    elsewhere = _Counter()
    n_primary = n_low = n_normal = 0
    for r in bam.fetch(resolved, lo, hi):
        if r.is_unmapped or r.is_secondary or r.is_duplicate or r.is_qcfail:
            continue
        primary = not r.is_supplementary
        low = r.mapping_quality < min_mapq
        if primary:
            n_primary += 1
            n_low += low
        rec = _read_record(r)
        kinds = []
        partner = None
        sides = None          # (side here, side there) of the join a partner read implies
        # supplementary pieces (split reads)
        sa = _parse_sa(r.get_tag("SA")) if r.has_tag("SA") else []
        for e in sa:
            e_end = e["pos"] + _ref_len(e["cigar"]) - 1
            near = norm_chrom(e["chrom"]) == there_n and (
                abs(e["pos"] - partner_pos) <= window_bp or abs(e_end - partner_pos) <= window_bp)
            if near:
                kinds.append("split_partner")
                partner = {"chromosome": display_chrom(e["chrom"]), "start": e["pos"], "end": e_end,
                           "strand": e["strand"], "cigar": e["cigar"], "mapq": e["mapq"]}
                mine = not intra or _nearer_here((rec["start"] + rec["end"]) // 2, (e["pos"] + e_end) // 2,
                                                 pos, partner_pos)
                sides = ((_side_of_piece(rec["start"], rec["end"], pos), _side_of_piece(e["pos"], e_end, partner_pos))
                         if mine else
                         (_side_of_piece(e["pos"], e_end, pos), _side_of_piece(rec["start"], rec["end"], partner_pos)))
                if mine and not low and e["mapq"] >= min_mapq and r.query_name not in split:
                    split[r.query_name] = {"here": {**rec, "side": sides[0]}, "there": {**partner, "side": sides[1]}}
                break
        if sa and not kinds:
            kinds.append("split_other")
            e = sa[0]
            partner = {"chromosome": display_chrom(e["chrom"]), "start": e["pos"],
                       "strand": e["strand"], "mapq": e["mapq"]}
        # pairs (primary alignments only: a supplementary piece repeats its read's mate)
        if primary and r.is_paired:
            if r.mate_is_unmapped:
                kinds.append("mate_unmapped")
            elif r.next_reference_name is not None:
                mate_n = norm_chrom(r.next_reference_name)
                mpos = r.next_reference_start + 1
                mate = {"chromosome": display_chrom(r.next_reference_name), "start": mpos,
                        "strand": "-" if r.mate_is_reverse else "+"}
                if mate_n == there_n and abs(mpos - partner_pos) <= window_bp and (
                        not intra or not r.is_proper_pair):
                    kinds.append("pair_partner")
                    if kinds[0] == "split_other":
                        # its supplementary piece lies elsewhere, but its pair joins the two
                        # ends and is counted: it is drawn as a joining read
                        kinds.insert(0, kinds.pop())
                        partner = None
                    partner = partner or mate
                    side_r, side_m = ("3" if not r.is_reverse else "5"), ("3" if not r.mate_is_reverse else "5")
                    mine = not intra or _nearer_here(rec["start"], mpos, pos, partner_pos)
                    sides = sides or ((side_r, side_m) if mine else (side_m, side_r))
                    if mine and not low:
                        pairs[r.query_name] = {"here": {**rec, "side": side_r},
                                               "there": {**mate, "side": side_m,
                                                         "mapq": r.get_tag("MQ") if r.has_tag("MQ") else None}}
                elif mate_n != here_n and any(mate_n == rc and abs(mpos - rp) <= window_bp for rc, rp in related):
                    kinds.append("pair_related")
                    if kinds[0] == "split_other":
                        kinds.insert(0, kinds.pop())
                        partner = None
                    partner = partner or mate
                    sides = sides or ("3" if not r.is_reverse else "5", "3" if not r.mate_is_reverse else "5")
                elif _bt._primary_contig(r.next_reference_name) != _bt._primary_contig(resolved) \
                        and mate_n != here_n:
                    kinds.append("pair_other_chrom")
                    partner = partner or mate
                    if not low:
                        elsewhere[display_chrom(r.next_reference_name)] += 1
                elif r.next_reference_name != resolved:
                    # an alt haplotype or unlocalised scaffold of this same chromosome:
                    # not a join between chromosomes (see bam_tools._primary_contig)
                    kinds.append("pair_alt_contig")
                    partner = partner or mate
                elif not r.is_proper_pair:
                    same_strand = r.is_reverse == r.mate_is_reverse
                    if same_strand:
                        kinds.append("pair_orientation")
                    else:
                        left_read = r.reference_start <= r.next_reference_start
                        # FR is normal: the left read on +, the right read on -
                        outward = (left_read and r.is_reverse) or (not left_read and not r.is_reverse)
                        kinds.append("pair_orientation" if outward else "pair_insert_size")
                    partner = partner or mate
        # reads cut at a base: a clip on the left means the aligned part starts at
        # the breakpoint (the right part joins, side 5); on the right, it ends there (side 3)
        if not low and rec["clip_left"] >= min_clip and abs(rec["start"] - pos) <= CLIP_SEARCH_BP:
            clip_starts[rec["start"]] += 1
        if not low and rec["clip_right"] >= min_clip and abs(rec["end"] - pos) <= CLIP_SEARCH_BP:
            clip_ends[rec["end"]] += 1
        if not kinds and (rec["clip_left"] >= min_clip or rec["clip_right"] >= min_clip):
            kinds.append("clipped")
        if not kinds and low:
            # a read placed equally well elsewhere (IGV draws these white): shown,
            # because at a breakpoint inside a repeat the evidence can be in them
            kinds.append("ambiguous")
        if not kinds:
            n_normal += 1
            continue
        kind_counts[kinds[0]] += 1
        item = {**rec, "kind": kinds[0], "kinds": kinds, "low_mapq": bool(low),
                "partner": partner, "sides": sides}
        if kinds[0] in ("pair_partner", "split_partner", "pair_related"):
            if len(view_partner) < VIEW_LIMIT:
                view_partner.append(item)
        elif len(view_other) < 10 * VIEW_LIMIT:
            view_other.append(item)
    # Every read that joins the two ends is drawn (up to VIEW_LIMIT); the other
    # abnormal reads fill the rest, taken evenly across the window so a dense
    # repeat at one edge cannot crowd out the rest of it.
    room = max(0, VIEW_LIMIT - len(view_partner))
    if len(view_other) > room:
        step = len(view_other) / room if room else 0
        view_other = [view_other[int(i * step)] for i in range(room)] if room else []
    view = sorted(view_partner + view_other, key=lambda x: (x["start"], x["end"]))
    shown = len(view)
    return {"pairs": pairs, "split": split, "view": view,
            "clips": {"5": dict(clip_starts), "3": dict(clip_ends)},
            "elsewhere": dict(elsewhere.most_common()), "primary_reads": n_primary,
            "low_mapq_reads": n_low, "normal_reads_hidden": n_normal,
            "abnormal_by_kind": dict(kind_counts), "view_truncated": shown < sum(kind_counts.values()),
            "window": {"start": lo + 1, "end": hi}}


# ── one junction, both ends ─────────────────────────────────────────────────
def _brief(rec):
    keep = ("start", "end", "strand", "cigar", "mapq", "chromosome", "side")
    return {k: rec.get(k) for k in keep if rec.get(k) is not None}


def junction_support(bam_path, chrom_a, pos_a, chrom_b, pos_b, window_bp=WINDOW_BP,
                     min_mapq=MIN_MAPQ, list_reads=True, view=False, orientation=None, related=()):
    """The reads that support a junction between chrom_a:pos_a and chrom_b:pos_b.

    Counts read pairs with one read at each end and split reads with one piece
    at each end, both at MAPQ >= min_mapq on both sides, and splits them by the
    orientation they imply (written from end A to end B). A pair or read seen
    from both ends is counted once.
    """
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": f"the read file could not be opened: {type(e).__name__}", "error_type": "bam_access"}
    try:
        return support_with(bam, chrom_a, pos_a, chrom_b, pos_b, window_bp, min_mapq, list_reads, view,
                            orientation, related)
    finally:
        bam.close()


def _pile(clips, side, pos):
    """The base where most reads are cut on the given side (both sides when side is
    None), nearest the breakpoint on a tie: (position, reads)."""
    pool = _Counter()
    for sd in ((side,) if side in ("3", "5") else ("3", "5")):
        pool.update(clips.get(sd) or {})
    if not pool:
        return None, 0
    best = max(pool.items(), key=lambda kv: (kv[1], -abs(kv[0] - pos)))
    return best[0], best[1]


def support_with(bam, chrom_a, pos_a, chrom_b, pos_b, window_bp=WINDOW_BP, min_mapq=MIN_MAPQ,
                 list_reads=False, view=False, orientation=None, related=()):
    """junction_support on a read file that is already open (review of many junctions).
    With the junction's orientation, the clipped-read pile-up at each end is looked
    for on the side that joins (the other junction of a reciprocal pair, a few bp
    away, has its own pile on the other side)."""
    for name, v in (("position_a", pos_a), ("position_b", pos_b), ("window_bp", window_bp),
                    ("min_mapq", min_mapq)):
        if not isinstance(v, int) or isinstance(v, bool) or v < 0:
            return {"error": f"{name} must be a non-negative integer (got {v!r})",
                    "error_type": "bad_parameter"}
    A = scan_end(bam, chrom_a, pos_a, chrom_b, pos_b, window_bp, min_mapq, related=related)
    if "error" in A:
        return {**A, "end": "a"}
    B = scan_end(bam, chrom_b, pos_b, chrom_a, pos_a, window_bp, min_mapq, related=related)
    if "error" in B:
        return {**B, "end": "b"}

    # A pair counts when BOTH reads reached the floor: each was seen, with its
    # own mapping quality, from its own end.
    pair_names = sorted(set(A["pairs"]) & set(B["pairs"]))
    one_sided = len((set(A["pairs"]) | set(B["pairs"])) - set(pair_names))
    pairs = []
    for n in pair_names:
        a, b = A["pairs"][n]["here"], B["pairs"][n]["here"]
        pairs.append({"name": n, "orientation": f"{a['side']}to{b['side']}",
                      "a": _brief({**a, "chromosome": display_chrom(chrom_a)}),
                      "b": _brief({**b, "chromosome": display_chrom(chrom_b)})})
    split = {}
    for n, s in A["split"].items():
        split[n] = {"name": n, "orientation": f"{s['here']['side']}to{s['there']['side']}",
                    "a": _brief({**s["here"], "chromosome": display_chrom(chrom_a)}), "b": _brief(s["there"])}
    for n, s in B["split"].items():
        if n not in split:
            split[n] = {"name": n, "orientation": f"{s['there']['side']}to{s['here']['side']}",
                        "a": _brief(s["there"]), "b": _brief({**s["here"], "chromosome": display_chrom(chrom_b)})}
    split = [split[k] for k in sorted(split)]

    by_or = {}
    for o in ORIENTATIONS:
        pn = {p["name"] for p in pairs if p["orientation"] == o}
        sn = {s["name"] for s in split if s["orientation"] == o}
        if pn or sn:
            by_or[o] = {"read_pairs": len(pn), "split_reads": len(sn), "fragments": len(pn | sn)}
    fragments = len({p["name"] for p in pairs} | {s["name"] for s in split})

    sides = (orientation[0], orientation[-1]) if orientation in ORIENTATIONS else (None, None)

    def end_summary(E, chrom, pos, side):
        b, st = band_at(chrom, pos)
        cp, cn = _pile(E["clips"], side, pos)
        return {"chromosome": display_chrom(chrom), "position": pos, "band": b,
                "region_note": region_note(chrom, pos),
                "clipped_at": cp, "clipped_reads": cn,
                "pairs_to_other_places": sum(E["elsewhere"].values()),
                "other_places": dict(list(E["elsewhere"].items())[:8]),
                "low_mapq_fraction": round(E["low_mapq_reads"] / E["primary_reads"], 3) if E["primary_reads"] else None,
                "reads_in_window": E["primary_reads"], "window": E["window"]}

    out = {
        "end_a": end_summary(A, chrom_a, pos_a, sides[0]),
        "end_b": end_summary(B, chrom_b, pos_b, sides[1]),
        "read_pairs": len(pairs),
        "split_reads": len(split),
        "fragments": fragments,
        "by_orientation": by_or,
        "pairs_with_one_read_below_min_mapq": one_sided,
        "thresholds_applied": [
            {"name": "window_bp", "value": window_bp, "provenance": PROVENANCE["window_bp"]},
            {"name": "min_mapq", "value": min_mapq, "provenance": PROVENANCE["min_mapq"]},
            {"name": "min_clip_bases", "value": MIN_CLIP_BASES, "provenance": PROVENANCE["min_clip_bases"]},
        ],
        "note": ("read_pairs: pairs with one read within window_bp of each end, both at "
                 "mapping quality >= min_mapq. split_reads: reads with one aligned piece "
                 "within window_bp of each end, both pieces at >= min_mapq. fragments: "
                 "distinct reads behind either count. by_orientation splits them by the join "
                 "they imply, written from end_a to end_b (3 = the part left of the "
                 "breakpoint, 5 = the part right of it)."),
    }
    if list_reads:
        out["reads"] = {"pairs": pairs[:LIST_LIMIT], "split": split[:LIST_LIMIT],
                        "listed_limit": LIST_LIMIT,
                        "truncated": len(pairs) > LIST_LIMIT or len(split) > LIST_LIMIT}
    if view:
        # the join each partner read implies, written from end A to end B like by_orientation
        for x in A["view"]:
            sd = x.pop("sides", None)
            x["orientation"] = f"{sd[0]}to{sd[1]}" if sd else None
        for x in B["view"]:
            sd = x.pop("sides", None)
            x["orientation"] = f"{sd[1]}to{sd[0]}" if sd else None
        out["view"] = {e: {"reads": E["view"], "normal_reads_hidden": E["normal_reads_hidden"],
                           "abnormal_by_kind": E["abnormal_by_kind"], "truncated": E["view_truncated"],
                           "window": E["window"]}
                       for e, E in (("a", A), ("b", B))}
    return out


# ── junctions -> rearrangements ─────────────────────────────────────────────
def _junction_view(j):
    """A caller junction (vcf_tools.Junction as a dict) in genome order."""
    c1, p1, c2, p2, o = genome_order(j["chrom1"], j["pos1"], j["chrom2"], j["pos2"], j.get("orientation"))
    return {"candidate_id": j["candidate_id"], "caller_id": j.get("caller_id"), "svtype": j.get("svtype"),
            "filter": j.get("filter"), "pe": j.get("pe"), "sr": j.get("sr"), "precise": j.get("precise"),
            "orientation": o, "a": {"chromosome": display_chrom(c1), "position": p1},
            "b": {"chromosome": display_chrom(c2), "position": p2}}


def _complementary(o1, o2):
    return {o1, o2} in ({"3to5", "5to3"}, {"3to3", "5to5"})


def group_rearrangements(junctions, near_bp=GROUP_NEAR_BP, segment_max_bp=SEGMENT_MAX_BP):
    """Group a list of junctions into rearrangements.

    Junctions between two different chromosomes are linked when they join the
    same two chromosomes and their ends lie within near_bp on both (the two
    junctions of a balanced translocation), or within near_bp on one and up to
    segment_max_bp apart on the other (the two ends of a segment). Linked
    junctions form one rearrangement. A junction within one chromosome
    (deletion, duplication, inversion, insertion) is a rearrangement of its own.
    """
    J = [_junction_view(j) for j in junctions]
    parent = list(range(len(J)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    by_pair = {}
    for i, j in enumerate(J):
        ca, cb = norm_chrom(j["a"]["chromosome"]), norm_chrom(j["b"]["chromosome"])
        if ca != cb:
            by_pair.setdefault((ca, cb), []).append(i)
    for idx in by_pair.values():
        idx.sort(key=lambda i: J[i]["a"]["position"])
        for x in range(len(idx)):
            for y in range(x + 1, len(idx)):
                i, k = idx[x], idx[y]
                da = abs(J[i]["a"]["position"] - J[k]["a"]["position"])
                db = abs(J[i]["b"]["position"] - J[k]["b"]["position"])
                if da > segment_max_bp:
                    break
                if (da <= near_bp and db <= segment_max_bp) or (db <= near_bp and da <= segment_max_bp):
                    ri, rk = find(i), find(k)
                    if ri != rk:
                        parent[rk] = ri
    groups = {}
    for i in range(len(J)):
        groups.setdefault(find(i), []).append(J[i])
    out = []
    for members in groups.values():
        members.sort(key=lambda j: (j["a"]["position"], j["b"]["position"]))
        out.append(_describe(members, near_bp))
    return out


def _span(positions):
    return {"start": min(positions), "end": max(positions)}


def _describe(members, near_bp):
    first = members[0]
    ca, cb = first["a"]["chromosome"], first["b"]["chromosome"]
    inter = norm_chrom(ca) != norm_chrom(cb)
    pa = [m["a"]["position"] for m in members]
    pb = [m["b"]["position"] for m in members]
    rid = "REARR_" + _hashlib.sha256("|".join(sorted(m["candidate_id"] for m in members)).encode()).hexdigest()[:10]
    ev = {"rearrangement_id": rid, "junctions": members,
          "chromosomes": [ca, cb] if inter else [ca],
          "breakpoints": [], "segment": None}
    if not inter:
        svt = first.get("svtype") or "SV"
        kind = {"DEL": "deletion", "DUP": "duplication", "INV": "inversion", "INS": "insertion"}.get(svt, "rearrangement")
        lo, hi = min(pa + pb), max(pa + pb)
        ev.update(kind=kind, pattern="within one chromosome",
                  segment={"chromosome": ca, "start": lo, "end": hi, "length": hi - lo + 1},
                  breakpoints=[{"chromosome": ca, "start": lo, "end": lo}, {"chromosome": ca, "start": hi, "end": hi}])
        return _with_bands(ev)
    ev["kind"] = "translocation"
    sa, sb = _span(pa), _span(pb)
    ev["breakpoints"] = [{"chromosome": ca, **sa}, {"chromosome": cb, **sb}]
    n = len(members)
    if n == 1:
        ev["pattern"] = "one junction"
    elif n == 2 and sa["end"] - sa["start"] <= near_bp and sb["end"] - sb["start"] <= near_bp:
        ev["pattern"] = ("both junctions" if _complementary(members[0]["orientation"], members[1]["orientation"])
                         else "two junctions, same direction")
    elif n == 2:
        # near on one chromosome, apart on the other: the far one bounds a segment
        seg_chrom, seg = (cb, sb) if sb["end"] - sb["start"] > near_bp else (ca, sa)
        ev["pattern"] = "segment"
        ev["segment"] = {"chromosome": seg_chrom, "start": seg["start"], "end": seg["end"],
                         "length": seg["end"] - seg["start"] + 1}
    else:
        ev["pattern"] = f"complex ({n} junctions)"
    return _with_bands(ev)


def _with_bands(ev):
    for bp in ev["breakpoints"]:
        b, _ = band_at(bp["chromosome"], bp["start"])
        bp["band"] = b
        bp["region_note"] = region_note(bp["chromosome"], bp["start"])
    if ev["kind"] == "translocation":
        (x, y) = ev["breakpoints"]
        cx, cy = norm_chrom(x["chromosome"]), norm_chrom(y["chromosome"])
        if x.get("band") and y.get("band"):
            # ISCN lists a sex chromosome first, then autosomes in numerical order;
            # genome order already does that except X/Y, which sort last there
            first, second = (y, x) if cy in ("X", "Y") and cx not in ("X", "Y") else (x, y)
            ev["iscn"] = (f"t({norm_chrom(first['chromosome'])};{norm_chrom(second['chromosome'])})"
                          f"({first['band']};{second['band']})")
    return ev


# ── what each junction makes of the chromosomes ─────────────────────────────
def _piece(chrom, part, pos, reversed_=False):
    L = chrom_length(chrom)
    start, end = (1, pos) if part == "left" else (pos, L)
    cen = centromere(chrom)
    has_cen = None
    if cen and start is not None and end is not None:
        mid = (cen[0] + cen[1]) // 2
        has_cen = start <= mid <= end
    return {"chromosome": display_chrom(chrom), "part": part, "start": start, "end": end,
            "length": (end - start + 1) if (start is not None and end is not None) else None,
            "reversed": reversed_, "has_centromere": has_cen}


def derivative(junction):
    """The chromosome a junction makes, read from one telomere to the other.

    Pieces in order; a piece flagged reversed runs against its reference
    direction. `moved` is the piece without a centromere (the one that left its
    chromosome), `receiving` the piece that keeps one. When both or neither
    piece holds a centromere, both keys are None and `centromeres` says how many.
    """
    o = junction.get("orientation")
    a, b = junction["a"], junction["b"]
    if o == "3to5":
        pieces = [_piece(a["chromosome"], "left", a["position"]), _piece(b["chromosome"], "right", b["position"])]
    elif o == "5to3":
        pieces = [_piece(b["chromosome"], "left", b["position"]), _piece(a["chromosome"], "right", a["position"])]
    elif o == "3to3":
        pieces = [_piece(a["chromosome"], "left", a["position"]),
                  _piece(b["chromosome"], "left", b["position"], reversed_=True)]
    elif o == "5to5":
        pieces = [_piece(b["chromosome"], "right", b["position"], reversed_=True),
                  _piece(a["chromosome"], "right", a["position"])]
    else:
        return {"candidate_id": junction.get("candidate_id"), "orientation": o, "pieces": [],
                "note": "the caller reported no orientation, so the join cannot be drawn"}
    cens = [p for p in pieces if p["has_centromere"]]
    moved = receiving = None
    if len(cens) == 1:
        receiving = cens[0]
        moved = pieces[1] if pieces[0] is receiving else pieces[0]
    name = f"der({norm_chrom(receiving['chromosome'])})" if receiving else None
    return {"candidate_id": junction.get("candidate_id"), "orientation": o, "name": name,
            "pieces": pieces, "moved": moved, "receiving": receiving, "centromeres": len(cens)}
