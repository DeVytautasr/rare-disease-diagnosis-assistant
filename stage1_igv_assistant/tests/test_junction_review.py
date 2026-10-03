#!/usr/bin/env python3
"""
Regression tests for the geneticist's view (Phase 26): junction-level read
counts, the grouping of junctions into rearrangements, the chromosomes a join
makes, cytogenetic bands, the local gene table with OMIM marks, the review of
every candidate, and how the interface offers all this to the assistant.

Built on synthetic reads (tests/review_fixture.py) whose every supporting read
was placed on purpose, so each expected count is the planted count, not a
number read back from the code under test. Every check that could pass
vacuously is first shown to fail on a control.

Run: python3 stage1_igv_assistant/tests/test_junction_review.py
"""
import asyncio
import json
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant.tests import review_fixture as F  # noqa: E402
from stage1_igv_assistant.tools import junction_tools as jt  # noqa: E402
from stage1_igv_assistant.tools import gene_table as gt  # noqa: E402
from stage1_igv_assistant.tools import vcf_tools  # noqa: E402
from stage1_igv_assistant import review_server as R  # noqa: E402

FAILURES = []


def check(label, condition, detail=""):
    ok = bool(condition)
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + ("" if ok or not detail else f"  -- {detail}"))
    if not ok:
        FAILURES.append(label)
    return ok


def run(d):
    fx = F.build_all(d)
    bam = fx["bam"]

    # ── 1. one junction, both ends ──────────────────────────────────────────
    print("reads at one junction")
    s = jt.junction_support(bam, "chr19", 48000000, "chr22", 30000000, orientation="3to5")
    o = s.get("by_orientation", {})
    check("the J1 join counts exactly the planted read pairs and split reads",
          o.get("3to5") == {"read_pairs": F.J1_PAIRS, "split_reads": F.J1_SPLIT, "fragments": F.J1_PAIRS + F.J1_SPLIT},
          str(o.get("3to5")))
    check("the reciprocal join's reads, a few bp away, are counted under their own orientation",
          o.get("5to3") == {"read_pairs": F.J2_PAIRS, "split_reads": F.J2_SPLIT, "fragments": F.J2_PAIRS + F.J2_SPLIT},
          str(o.get("5to3")))
    check("a pair seen from both ends is counted once",
          s.get("read_pairs") == F.J1_PAIRS + F.J2_PAIRS, str(s.get("read_pairs")))
    check("a pair whose mate maps ambiguously is not counted, and is reported apart",
          s.get("pairs_with_one_read_below_min_mapq") == F.J1_ONE_SIDED,
          str(s.get("pairs_with_one_read_below_min_mapq")))
    strict = jt.junction_support(bam, "chr19", 48000000, "chr22", 30000000, min_mapq=61)
    check("control: above the reads' mapping quality nothing is counted",
          strict.get("read_pairs") == 0 and strict.get("split_reads") == 0, str(strict.get("read_pairs")))
    check("the clipped reads pile up at the J1 breakpoints, on the joining side",
          (s["end_a"]["clipped_at"], s["end_a"]["clipped_reads"], s["end_b"]["clipped_at"], s["end_b"]["clipped_reads"])
          == (48000000, F.J1_SPLIT, 30000000, F.J1_SPLIT),
          str((s["end_a"]["clipped_at"], s["end_a"]["clipped_reads"], s["end_b"]["clipped_at"], s["end_b"]["clipped_reads"])))
    s2 = jt.junction_support(bam, "chr19", 48000201, "chr22", 29999999, orientation="5to3")
    check("... and for J2, 201 bp away, on its own side",
          (s2["end_a"]["clipped_at"], s2["end_a"]["clipped_reads"], s2["end_b"]["clipped_at"], s2["end_b"]["clipped_reads"])
          == (48000201, F.J2_SPLIT, 29999999, F.J2_SPLIT),
          str((s2["end_a"]["clipped_at"], s2["end_b"]["clipped_at"])))
    noside = jt.junction_support(bam, "chr19", 48000201, "chr22", 29999999)
    check("control: without the orientation the larger pile (J1's) is taken",
          noside["end_a"]["clipped_at"] == 48000000, str(noside["end_a"]["clipped_at"]))
    r2 = jt.junction_support(bam, "chr4", 60000000, "chr19", 10000000)
    check("a 5to5 join (both right parts) is read from the strands and clip sides",
          (r2["by_orientation"].get("5to5") or {}).get("fragments") == F.R2_PAIRS + F.R2_SPLIT
          and (r2["by_orientation"].get("3to3") or {}).get("read_pairs") == F.R2_RECIP_PAIRS,
          json.dumps(r2["by_orientation"]))
    flipped = jt.junction_support(bam, "chr19", 10000000, "chr4", 60000000)
    check("asking with the ends swapped gives the same counts, written from the other end",
          flipped["by_orientation"].get("5to5") == r2["by_orientation"].get("5to5")
          and flipped["by_orientation"].get("3to3") == r2["by_orientation"].get("3to3"),
          json.dumps(flipped["by_orientation"]))
    lst = s.get("reads") or {}
    check("the supporting reads are listed with both ends",
          len(lst.get("pairs", [])) == F.J1_PAIRS + F.J2_PAIRS
          and all(p["a"]["chromosome"] == "chr19" and p["b"]["chromosome"] == "chr22" for p in lst["pairs"])
          and len(lst.get("split", [])) == F.J1_SPLIT + F.J2_SPLIT, str(len(lst.get("pairs", []))))
    bad = jt.junction_support(bam, "chr19", 48000000, "chrZZ", 1)
    check("an end on a chromosome the file lacks is an error, not zero reads",
          bad.get("error_type") == "invalid_region", str(bad)[:120])
    out = jt.junction_support(bam, "chr19", 99_000_000, "chr22", 30000000)
    check("a position past the chromosome's end is an error", out.get("error_type") == "invalid_region", str(out)[:120])

    # ── 2. the abnormal reads drawn at each end ─────────────────────────────
    print("\nthe read view")
    v = jt.junction_support(bam, "chr19", 48000000, "chr22", 30000000, orientation="3to5", view=True)["view"]
    kinds = {x["kind"] for e in "ab" for x in v[e]["reads"]}
    check("normal reads are hidden and counted", v["a"]["normal_reads_hidden"] > 100, str(v["a"]["normal_reads_hidden"]))
    check("every read drawn carries a kind (none is a normal read)",
          all(x.get("kind") for e in "ab" for x in v[e]["reads"]) and kinds <= set(
              ["pair_partner", "split_partner", "pair_related", "pair_other_chrom", "split_other", "mate_unmapped",
               "pair_orientation", "pair_insert_size", "pair_alt_contig", "clipped", "ambiguous"]), str(kinds))
    joining = {x["name"] for x in v["a"]["reads"] if x["kind"] in ("pair_partner", "split_partner")}
    check("every read joining the two ends at chr19 is drawn",
          {f"j1p{i}" for i in range(F.J1_PAIRS)} <= joining and {f"j1s{i}" for i in range(F.J1_SPLIT)} <= joining,
          str(sorted(joining))[:200])
    check("each joining read says which join it supports, written from end A to end B",
          {x["orientation"] for x in v["a"]["reads"] if x["name"].startswith("j1p")} == {"3to5"}
          and {x["orientation"] for x in v["b"]["reads"] if x["name"].startswith("j2p")} == {"5to3"},
          str({x["orientation"] for x in v["b"]["reads"] if x["name"].startswith("j2p")}))
    c4 = jt.junction_support(bam, "chr4", 50000000, "chr19", 26000000, view=True)
    hollow = [x for x in c4["view"]["a"]["reads"] if x["low_mapq"]]
    check("reads below the quality floor are drawn (flagged) but not counted",
          len(hollow) >= 4 and c4["read_pairs"] == 0, f"{len(hollow)} {c4['read_pairs']}")
    check("ambiguously placed normal reads are drawn too (IGV's white reads)",
          sum(1 for x in c4["view"]["b"]["reads"] if x["kind"] == "ambiguous") == F.R4_AMBIGUOUS,
          str(c4["view"]["b"]["abnormal_by_kind"]))
    check("read pairs pointing to other chromosomes are counted per chromosome",
          c4["end_a"]["pairs_to_other_places"] == F.R4_ELSEWHERE and set(c4["end_a"]["other_places"]) == {"chr1", "chr7", "chr16"},
          str(c4["end_a"]["other_places"]))
    seg = jt.junction_support(bam, "chr4", 100000000, "chr22", 40000000, view=True)
    seg_rel = jt.junction_support(bam, "chr4", 100000000, "chr22", 40000000, view=True,
                                  related=[("chr22", 40050000), ("chr4", 100000001)])
    check("mates at the other end of the same rearrangement are not 'other places'",
          seg["end_a"]["pairs_to_other_places"] == F.R3B_PAIRS and seg_rel["end_a"]["pairs_to_other_places"] == 0,
          f"{seg['end_a']['pairs_to_other_places']} {seg_rel['end_a']['pairs_to_other_places']}")
    check("... and a related position on this end's own chromosome does not swallow its normal reads",
          seg_rel["view"]["a"]["normal_reads_hidden"] == seg["view"]["a"]["normal_reads_hidden"],
          f"{seg_rel['view']['a']['normal_reads_hidden']} {seg['view']['a']['normal_reads_hidden']}")

    # ── 3. junctions into rearrangements ────────────────────────────────────
    print("\ngrouping")
    J = lambda cid, c1, p1, c2, p2, o: {"candidate_id": cid, "chrom1": c1, "pos1": p1, "chrom2": c2, "pos2": p2,
                                         "orientation": o, "svtype": "BND", "filter": "PASS", "pe": 5, "sr": 2}
    g = jt.group_rearrangements([J("a", "chr22", 30000000, "chr19", 48000000, "5to3"),
                                 J("b", "chr22", 29999999, "chr19", 48000201, "3to5")])
    check("two complementary junctions a few bp apart are one balanced translocation",
          len(g) == 1 and g[0]["pattern"] == "both junctions" and g[0]["chromosomes"] == ["chr19", "chr22"],
          json.dumps([(x["pattern"], x["chromosomes"]) for x in g]))
    check("the junctions are written in genome order, orientation flipped with them",
          sorted(j["orientation"] for j in g[0]["junctions"]) == ["3to5", "5to3"]
          and g[0]["junctions"][0]["a"]["chromosome"] == "chr19", json.dumps(g[0]["junctions"][0]))
    g = jt.group_rearrangements([J("a", "chr22", 30000000, "chr19", 48000000, "3to5"),
                                 J("b", "chr22", 30000600, "chr19", 48000700, "3to5")])
    check("two junctions in the same direction are flagged, not called reciprocal",
          len(g) == 1 and g[0]["pattern"] == "two junctions, same direction", str([x["pattern"] for x in g]))
    g = jt.group_rearrangements([J("a", "chr22", 30000000, "chr19", 48000000, "3to5"),
                                 J("b", "chr22", 30020000, "chr19", 48030000, "5to3")])
    check("control: junctions 20-30 kb apart on both chromosomes stay separate", len(g) == 2, str(len(g)))
    g = jt.group_rearrangements([J("a", "chr4", 100000000, "chr22", 40000000, "3to5"),
                                 J("b", "chr4", 100000001, "chr22", 40050000, "5to3")])
    check("near on one chromosome, 50 kb apart on the other: a segment, bounded by the far ends",
          len(g) == 1 and g[0]["pattern"] == "segment" and g[0]["segment"] == {
              "chromosome": "chr22", "start": 40000000, "end": 40050000, "length": 50001}, json.dumps(g[0].get("segment")))
    g = jt.group_rearrangements([J("a", "chr4", 100000000, "chr22", 40000000, "3to5"),
                                 J("b", "chr4", 100000001, "chr22", 40050000, "5to3"),
                                 J("c", "chr4", 100000500, "chr22", 40070000, "3to3")])
    check("three linked junctions are reported as complex", g[0]["pattern"] == "complex (3 junctions)", g[0]["pattern"])
    g = jt.group_rearrangements([{"candidate_id": "d", "chrom1": "chr19", "pos1": 50010000, "chrom2": "chr19",
                                  "pos2": 50000000, "orientation": "3to5", "svtype": "DEL"}])
    check("a deletion is its own rearrangement, its segment the deleted stretch",
          g[0]["kind"] == "deletion" and g[0]["segment"]["start"] == 50000000 and g[0]["segment"]["end"] == 50010000,
          json.dumps(g[0]["segment"]))

    # ── 4. bands, ISCN-like names, derivative chromosomes ───────────────────
    print("\nbands and derivatives")
    check("bands come from the bundled GRCh38 table", jt.band_at("chr19", 48000000) == ("q13.33", "gneg")
          and jt.chrom_length("chr19") == 58617616, str(jt.band_at("chr19", 48000000)))
    check("a centromere is recognised", jt.region_note("chr4", 50000000) == "in the centromere (chr4p11)",
          str(jt.region_note("chr4", 50000000)))
    check("control: a band outside centromeres and heterochromatin draws no note", jt.region_note("chr19", 48000000) is None)
    g = jt.group_rearrangements([J("a", "chr22", 30000000, "chr19", 48000000, "5to3")])
    check("the t(...)(...) name is built from the two bands", g[0].get("iscn") == "t(19;22)(q13.33;q12.2)", g[0].get("iscn"))
    g = jt.group_rearrangements([J("a", "chrX", 60000000, "chr5", 1000000, "3to5")])
    check("a sex chromosome is named first, as ISCN writes it", (g[0].get("iscn") or "").startswith("t(X;5)"),
          g[0].get("iscn"))
    dv = jt.derivative({"orientation": "3to5", "a": {"chromosome": "chr19", "position": 48000000},
                       "b": {"chromosome": "chr22", "position": 30000000}})
    check("3to5: chr19's centric part receives chr22's distal long arm",
          dv["name"] == "der(19)" and (dv["moved"]["chromosome"], dv["moved"]["start"], dv["moved"]["end"]) == ("chr22", 30000000, 50818468),
          json.dumps({k: dv[k] for k in ("name", "moved")}))
    dv = jt.derivative({"orientation": "5to5", "a": {"chromosome": "chr4", "position": 60000000},
                       "b": {"chromosome": "chr19", "position": 10000000}})
    check("5to5: both right parts, the second reversed; the centric one names the chromosome",
          dv["name"] == "der(19)" and dv["pieces"][0]["reversed"] and dv["moved"]["chromosome"] == "chr4",
          json.dumps({k: dv[k] for k in ("name", "pieces")})[:300])
    dv = jt.derivative({"orientation": "5to5", "a": {"chromosome": "chr4", "position": 60000000},
                       "b": {"chromosome": "chr19", "position": 40000000}})
    check("a join whose pieces both lack a centromere is reported as such (no der name)",
          dv["centromeres"] == 0 and dv["name"] is None and dv["moved"] is None, str(dv["centromeres"]))

    # ── 5. the gene table ───────────────────────────────────────────────────
    print("\ngenes")
    t = gt.GeneTable(*fx["genes"][:1], mim2gene_path=fx["genes"][1], disease_path=fx["genes"][2])
    a = t.at("chr19", 48000000)
    check("a breakpoint in GENE_A (+ strand) falls in intron 2 of 3",
          a and a[0]["name"] == "GENE_A" and a[0]["where"] == "intron 2 of 3", str(a[:1]))
    b = t.at("chr22", 30000000)
    check("on the - strand introns are numbered in the direction of transcription",
          b and b[0]["name"] == "GENE_B" and b[0]["where"] == "intron 2 of 3", str(b[:1]))
    e = t.at("chr22", 29950500)
    check("... and exons too: the lowest exon of a - strand gene is its last",
          e and e[0]["where"] == "exon 4 of 4", str(e[:1]))
    lnc = t.at("chr19", 48000100)
    check("protein-coding genes are listed before others", [x["name"] for x in lnc] == ["GENE_A", "LNC_F"],
          str([x["name"] for x in lnc]))
    check("OMIM marks: gene number and disorders", a[0].get("omim_gene") == 600001 and "OMIM:610002" in a[0].get("disorders", []),
          json.dumps(a[0]))
    check("control: a gene without a disorder carries none", "disorders" not in b[0], json.dumps(b[0]))
    n = t.nearest("chr19", 10000000)
    check("nearest gene on each side within 1 Mb", (n.get("right") or {}).get("name") == "GENE_E"
          and n["right"]["distance_bp"] == 100000 and "left" not in n, json.dumps(n))
    ov = t.overlapping("chr22", 40000000, 40050000)
    check("genes in a segment, with the count of disorder genes", ov["count"] == 1 and ov["with_disorder"] == 1,
          json.dumps(ov)[:200])
    gt.load(os.path.join(d, "no-such-table.tsv.gz"))
    check("without a table the status says so, and nothing is named", gt.status()["available"] is False and gt.table() is None)

    # ── 6. reviewing every candidate ────────────────────────────────────────
    print("\nreview_candidates")
    R.configure(mask_path=fx["mask"], gene_table=fx["genes"][0], mim2gene=fx["genes"][1], gene_disorders=fx["genes"][2])
    vcf_tools.reset_registry()
    r = R.review_candidates(fx["calls"], bam, other_path=fx["other"])
    evs = r.get("rearrangements") or []
    check("7 junctions pass the default filters and form 5 rearrangements",
          r.get("junctions_passing_filters") == 7 and r.get("rearrangements_found") == 5,
          f"{r.get('junctions_passing_filters')} {r.get('rearrangements_found')}")
    order = [(e["breakpoints"][0]["chromosome"], e["breakpoints"][0]["start"], e["support"]["fragments"]) for e in evs]
    check("sorted by supporting reads; a tie goes to the caller's own counts",
          order == [("chr19", 48000000, 31), ("chr4", 60000000, 15), ("chr4", 100000000, 12),
                    ("chr4", 50000000, 0), ("chr19", 45000000, 0)], str(order))
    check("every filter is echoed", [f["filter"] for f in r["filters_applied"]] ==
          ["svtype", "filter_pass", "min_pe", "min_sr", "primary_only", "exclude_masked"], str(r["filters_applied"])[:200])
    c = {e["breakpoints"][0]["start"]: e["cautions"] for e in evs}
    check("a lone junction whose reciprocal is in the file but filtered out says so",
          any("reciprocal junction is in the caller's file" in x for x in c[60000000]), str(c[60000000]))
    check("... and that the reads support the reciprocal join the caller did not report",
          any(f"{F.R2_RECIP_PAIRS} reads also support the reciprocal join" in x for x in c[60000000]), str(c[60000000]))
    check("a centromeric junction with reads pointing elsewhere carries both cautions",
          sum("centromere" in x for x in c[50000000]) == 2 and any("point to other places" in x for x in c[50000000]),
          str(c[50000000]))
    check("a junction also called in the other sample is marked",
          any("also called in the comparison sample" in x for x in c[45000000]), str(c[45000000]))
    check("control: the balanced translocation carries no caution", c[48000000] == [], str(c[48000000]))
    vcf_tools.reset_registry()
    r0 = R.review_candidates(fx["calls"], bam)
    c0 = {e["breakpoints"][0]["start"]: e["cautions"] for e in r0["rearrangements"]}
    check("control: without a comparison sample nothing is marked recurrent",
          not any("comparison sample" in x for x in c0[45000000]), str(c0[45000000]))
    first = evs[0]
    check("the balanced translocation's genes, with OMIM marks, are named at both breakpoints",
          [g["in"][0]["name"] for g in first["genes"]] == ["GENE_A", "GENE_B"]
          and first["genes"][0]["in"][0].get("omim_gene") == 600001, json.dumps(first["genes"])[:300])
    check("each junction knows the chromosome it makes and the piece that moved",
          {j["join"]["name"] for j in first["junctions"]} == {"der(19)", "der(22)"}
          and all(j["join"]["moved"]["genes"]["count"] >= 1 for j in first["junctions"]),
          json.dumps([j["join"]["name"] for j in first["junctions"]]))
    rows = r.get("summary_rows") or []
    check("one summary row per rearrangement, in order, with no score",
          len(rows) == 5 and rows[0].startswith("1. t(19;22)") and not any("score" in x.lower() for x in rows), rows[:1])
    def keys(x):
        if isinstance(x, dict):
            for k, v in x.items():
                yield k
                yield from keys(v)
        elif isinstance(x, list):
            for v in x:
                yield from keys(v)
    check("the result carries no score field anywhere", not [k for k in keys(r) if "score" in k.lower()],
          str([k for k in keys(r) if "score" in k.lower()][:5]))
    R.configure(mask_path=None, gene_table=fx["genes"][0])
    vcf_tools.reset_registry()
    rm = R.review_candidates(fx["calls"], bam)
    check("without a configured exclude list the step is absent and the result says so",
          rm["exclude_regions_available"] is False and "exclude_masked" not in [f["filter"] for f in rm["filters_applied"]])
    R.configure(mask_path=fx["mask"], gene_table=fx["genes"][0], mim2gene=fx["genes"][1], gene_disorders=fx["genes"][2])
    gn = R.genes_near("chr22", 40000000, end=40050000)
    check("genes_near lists the genes of a segment from the local table", gn["genes"]["count"] == 1
          and gn["genes"]["genes"][0]["name"] == "GENE_D", json.dumps(gn)[:200])
    gt.load(None)
    gx = R.genes_near("chr19", 48000000)
    check("genes_near without a table is an error that says why, with the band still given",
          gx.get("error_type") == "no_gene_table" and gx.get("band") == "q13.33", json.dumps(gx)[:200])
    R.configure(mask_path=fx["mask"], gene_table=fx["genes"][0], mim2gene=fx["genes"][1], gene_disorders=fx["genes"][2])
    je = R.junction_evidence(bam, "chr19", 48000000, "chr22", 30000000, orientation="3to5", include_view=True)
    check("junction_evidence gives the join, the genes at both ends and the gene track of each window",
          je["join"]["name"] == "der(19)" and je["genes"]["end_1"]["in"][0]["name"] == "GENE_A"
          and any(gg["name"] == "GENE_A" for gg in je["view"]["a"]["genes"]), json.dumps(je["join"])[:160])
    jn = R.junction_evidence(bam, "chr19", 48000000, "chr22", 30000000)
    check("without an orientation it describes every join the reads support",
          set(jn.get("joins", {})) == {"3to5", "5to3"}, str(list(jn.get("joins", {}))))
    check("an orientation outside the vocabulary is refused",
          R.junction_evidence(bam, "chr19", 1, "chr22", 2, orientation="sideways").get("error_type") == "bad_parameter")

    # ── 7. what the interface offers the assistant ──────────────────────────
    print("\nthe assistant's tools")
    from stage1_igv_assistant import ui
    from stage1_igv_assistant import chat as chatmod
    saved = (dict(ui.DATASETS), dict(ui.CANDIDATE_FILES), {k: set(v) for k, v in ui.PUBLIC_LABELS.items()})
    try:
        ui.DATASETS.clear(); ui.CANDIDATE_FILES.clear()
        ui.DATASETS.update({"PUB": bam, "PRIV": bam})
        ui.CANDIDATE_FILES.update({"PUB": fx["calls"], "PRIV": fx["other"]})
        ui.PUBLIC_LABELS["datasets"] = {"PUB"}
        ui.PUBLIC_LABELS["candidates"] = {"PUB"}
        ui._CHAT_TOOLS.clear()
        tools, where = ui._chat_tools()
        names = {t["function"]["name"] for t in tools}
        check("the assistant is offered the three review tools", {"junction_evidence", "review_candidates", "genes_near"} <= names)
        check("... and not the combined score", "breakpoint_evidence_summary" not in names)
        check("... nor the internet gene lookup while a local gene table is set up", "gene_at_locus" not in names)
        # Phase 27: the IGV image is the page's; a model would get an opaque reference only
        check("... nor the IGV image tool", "igv_junction_view" not in names)
        check("control: ... which the review server does provide",
              "igv_junction_view" in {t.name for t in asyncio.run(R.mcp.list_tools())})
        gt.load(None); ui._CHAT_TOOLS.clear()
        names0 = {t["function"]["name"] for t in ui._chat_tools()[0]}
        check("control: without a local table the internet lookup is offered again", "gene_at_locus" in names0)
        R.configure(mask_path=fx["mask"], gene_table=fx["genes"][0], mim2gene=fx["genes"][1], gene_disorders=fx["genes"][2])
        ui._CHAT_TOOLS.clear()
        res, err = ui._chat_exec("breakpoint_evidence_summary", {"dataset": "PUB", "chromosome": "chr19", "position": 48000000})
        check("a call to the hidden scoring tool is refused before it runs", res is None and "not offered" in (err or ""), err)
        rv = next(t for t in tools if t["function"]["name"] == "review_candidates")["function"]["parameters"]["properties"]
        check("review_candidates takes labels, never paths", {"candidates", "dataset", "other_candidates"} <= set(rv)
              and not {"path", "bam_path", "other_path"} & set(rv), str(sorted(rv)))
        check("a private comparison sample is caught by the cloud guard",
              ui._private_hits({"candidates": "PUB", "dataset": "PUB", "other_candidates": "PRIV"}) == ["PRIV"])
        check("control: test data alone is no hit",
              ui._private_hits({"candidates": "PUB", "dataset": "PUB", "other_candidates": "PUB"}) == [])
        rec, err = ui._chat_exec("review_candidates", {"candidates": "PUB", "dataset": "PUB"})
        shown = chatmod.shrink_for_model("review_candidates", rec["result"])
        check("the model is shown one line per rearrangement, not the drawing data",
              shown.get("summary_rows") and "rearrangements" not in shown and len(json.dumps(shown)) < 12000,
              str(len(json.dumps(shown))))
        rec, err = ui._chat_exec("junction_evidence", {"dataset": "PUB", "chromosome_1": "chr19", "position_1": 48000000,
                                                       "chromosome_2": "chr22", "position_2": 30000000,
                                                       "include_view": True})
        shown = chatmod.shrink_for_model("junction_evidence", rec["result"])
        check("... and the counts of a junction without its read lists",
              shown.get("read_pairs") == F.J1_PAIRS + F.J2_PAIRS and "reads" not in shown and "view" not in shown,
              str(sorted(shown)))
        check("the assistant's rules say there is no score and how to rank",
              "no evidence score" in ui.assistant_prompt() and "review_candidates ONCE" in ui.assistant_prompt())
        check("control: the shared system prompt (named by the recorded runs) is unchanged by this",
              "review_candidates" not in chatmod.SYSTEM_PROMPT)
        page = open(ui.PAGE_FILE).read()
        check("the page shows no score: none of the old gauge's words is left",
              not any(w in page for w in ("combined score", "reachable here", "evidence_score", "gauge-score")))
        check("control: the previous page (/classic) still shows the score", "combined score" in ui.PAGE.lower()
              or "evidence_score" in ui.PAGE)
    finally:
        ui.DATASETS.clear(); ui.DATASETS.update(saved[0])
        ui.CANDIDATE_FILES.clear(); ui.CANDIDATE_FILES.update(saved[1])
        ui.PUBLIC_LABELS.update(saved[2])
        ui._CHAT_TOOLS.clear()

    # ── 8. an event within one chromosome, its ends closer than the windows ──
    # Found when the Phase 26 patch was applied (FIGURE_MAP P.2, Q3): with both
    # ends on one chromosome and less than two windows apart, both reads of a pair
    # lie in both windows, and each end kept whichever read came last, so a pair
    # across an 800 bp deletion was reported as joining 5to5 with the same read at
    # both ends. Each end must keep the read (or split piece) nearer to itself.
    print("\nends closer than the two windows")
    small = os.path.join(d, "small_event.bam")
    h = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                         "SQ": [{"SN": c, "LN": n} for c, n in F.LEN.items()]})
    out = []
    # a deletion joining the left part of chr7 at 100,000,000 to its right part at
    # 100,000,800: one pair across it, and one read split at the two ends
    F._pair(h, out, "sd_pair", "chr7", 99_999_800, False, "chr7", 100_000_850, True)
    F._split_left_right(h, out, "sd_split", "chr7", 100_000_000, "chr7", 100_000_800, 70)
    # the same pair across a 10 kb deletion, whose windows do not overlap
    F._pair(h, out, "ld_pair", "chr7", 119_999_800, False, "chr7", 120_010_050, True)
    # a pair joining chr16:50,000,000 to chr22:45,500,000 (3to5) whose chr16 read also
    # has a supplementary piece on chr1, away from both ends
    out.append(F._seg(h, "sa_pair", "chr16", 49_999_800, [(0, 150)], F.P | F.R1F | F.MREV, "chr22", 45_500_100,
                      0, 60, {"MQ": 60, "SA": "chr1,5000000,+,60S90M,60,0;"}))
    out.append(F._seg(h, "sa_pair", "chr22", 45_500_100, [(0, 150)], F.P | F.R2F | F.REV, "chr16", 49_999_800,
                      0, 60, {"MQ": 60}))
    out.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(small, "wb", header=h) as fh:
        for r in out:
            fh.write(r)
    pysam.index(small)
    sd = jt.junction_support(small, "chr7", 100_000_000, "chr7", 100_000_800, view=True)
    check("a deletion whose ends lie 800 bp apart: its pair and its split read count once each, joined 3to5",
          sd.get("by_orientation") == {"3to5": {"read_pairs": 1, "split_reads": 1, "fragments": 2}},
          json.dumps(sd.get("by_orientation")))
    lp = ((sd.get("reads") or {}).get("pairs") or [{}])[0]
    check("... each end lists its own read: the + read at the first end, the - read at the second",
          ((lp.get("a") or {}).get("start"), (lp.get("a") or {}).get("strand"),
           (lp.get("b") or {}).get("start"), (lp.get("b") or {}).get("strand")) == (99_999_800, "+", 100_000_850, "-"),
          json.dumps(lp)[:200])
    drawn = {x["orientation"] for e in "ab" for x in sd["view"][e]["reads"] if x["kind"] in ("pair_partner", "split_partner")}
    check("... and every joining read drawn at either end carries that join", drawn == {"3to5"}, str(drawn))
    ld = jt.junction_support(small, "chr7", 120_000_000, "chr7", 120_010_000)
    check("a deletion whose ends lie 10 kb apart is counted as before: one pair, joined 3to5",
          ld.get("by_orientation") == {"3to5": {"read_pairs": 1, "split_reads": 0, "fragments": 1}},
          json.dumps(ld.get("by_orientation")))
    # Found with Q3 (FIGURE_MAP P.2, Q4): a read was drawn by the first thing it showed,
    # so a joining pair whose read also had a supplementary piece elsewhere was counted
    # but drawn as "split read, other piece elsewhere", among the reads that are thinned.
    sa = jt.junction_support(small, "chr16", 50_000_000, "chr22", 45_500_000, orientation="3to5", view=True)
    sa_drawn = [(x["kind"], x["orientation"]) for x in sa["view"]["a"]["reads"] if x["name"] == "sa_pair"]
    check("a joining pair whose read also has a piece elsewhere is counted, and drawn as joining the two ends",
          sa.get("read_pairs") == 1 and sa_drawn == [("pair_partner", "3to5")], f"{sa.get('read_pairs')} {sa_drawn}")


def main():
    with tempfile.TemporaryDirectory() as d:
        run(d)
    print("\n" + "=" * 60)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("ALL JUNCTION REVIEW TESTS PASSED")


if __name__ == "__main__":
    main()
