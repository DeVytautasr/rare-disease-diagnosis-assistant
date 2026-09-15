"""
test_vcf_tools.py
Covers the caller-agnostic VCF parsing layer and the four candidate-set
tools, against synthetic fixtures only -- no patient data.

Fixtures: a DELLY-style file (CHR2/POS2), a MATEID-paired file (manta/VCF
spec), a malformed file whose BND records use neither convention, a file with
a junction inside a masked region, a file with records that must deduplicate,
and an empty file.

Same style as the other suites here: a check() accumulator and sys.exit(1).
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from stage1_igv_assistant.tools import vcf_tools as V

FAILURES = []


def check(label, condition, detail=""):
    if condition:
        print(f"  PASS  {label}")
    else:
        print(f"  FAIL  {label}" + (f" -- {detail}" if detail else ""))
        FAILURES.append(label)


HDR_COMMON = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chrUn_KI270302v1,length=2274>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support">
##INFO=<ID=SR,Number=1,Type=Integer,Description="Split-read support">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="CI around POS">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="CI around END">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise">
##FILTER=<ID=LowQual,Description="low quality">
"""
HDR_DELLY = HDR_COMMON + """##INFO=<ID=CHR2,Number=1,Type=String,Description="Partner chromosome">
##INFO=<ID=POS2,Number=1,Type=Integer,Description="Partner position">
##INFO=<ID=CT,Number=1,Type=String,Description="Connection type">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
HDR_MATE = HDR_COMMON + """##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
HDR_NEITHER = HDR_COMMON + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def write(tmp, name, text):
    p = os.path.join(tmp, name)
    with open(p, "w") as fh:
        fh.write(text)
    return p


def run_tests():
    print("=" * 68)
    print("VCF PARSING LAYER AND CANDIDATE TOOLS")
    print("=" * 68)
    tmp = tempfile.mkdtemp(prefix="vcf_tools_test_")

    # ── fixtures ────────────────────────────────────────────────────────────
    # DELLY style: CHROM > CHR2 already (canonical), plus one intra-chrom DEL.
    delly = write(tmp, "delly.vcf", HDR_DELLY +
        "chr2\t3000000\tBND00000001\tN\tN[chr1:1000000[\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1000000;CT=3to5;PE=9;SR=4;PRECISE;CIPOS=-5,5;CIEND=-5,5\n"
        "chr17\t5000000\tBND00000002\tN\tN[chr1:2000000[\t.\tLowQual\t"
        "SVTYPE=BND;CHR2=chr1;POS2=2000000;CT=3to5;PE=2;SR=0;IMPRECISE;CIPOS=-300,300;CIEND=-300,300\n"
        "chr1\t7000000\tDEL00000001\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=7005000;PE=15;SR=6;PRECISE\n")

    # MATEID style: the SAME junction as delly record 1, written as two records
    # anchored the other way round.
    mate = write(tmp, "mate.vcf", HDR_MATE +
        "chr1\t1000000\tMANTA_A\tN\t]chr2:3000000]N\t.\tPASS\t"
        "SVTYPE=BND;MATEID=MANTA_B;PE=9;SR=4;PRECISE;CIPOS=-5,5\n"
        "chr2\t3000000\tMANTA_B\tN\tN[chr1:1000000[\t.\tPASS\t"
        "SVTYPE=BND;MATEID=MANTA_A;PE=9;SR=4;PRECISE;CIPOS=-5,5\n")

    # BND records carrying neither convention -> must fail loudly.
    bad = write(tmp, "bad.vcf", HDR_NEITHER +
        "chr2\t3000000\tX1\tN\t<BND>\t.\tPASS\tSVTYPE=BND;PE=5;SR=1\n")

    # Junction with one breakend inside a masked region.
    masked = write(tmp, "masked.vcf", HDR_DELLY +
        "chr2\t3000000\tM1\tN\tN[chr1:500[\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=500;CT=3to5;PE=8;SR=3;PRECISE\n"
        "chr2\t9000000\tM2\tN\tN[chr1:8000000[\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=8000000;CT=3to5;PE=8;SR=3;PRECISE\n")
    maskfile = write(tmp, "mask.tsv", "chr1\t0\t10000\ttelomere\nchrUn_KI270302v1\n")

    # Three records describing one junction within 500 bp, plus one far away.
    dup = write(tmp, "dup.vcf", HDR_DELLY +
        "chr2\t3000000\tD1\tN\tN[chr1:1000000[\t.\tLowQual\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1000000;CT=3to5;PE=3;SR=0\n"
        "chr2\t3000200\tD2\tN\tN[chr1:1000150[\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1000150;CT=3to5;PE=9;SR=5\n"
        "chr2\t3000400\tD3\tN\tN[chr1:1000300[\t.\tLowQual\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1000300;CT=3to5;PE=4;SR=1\n"
        "chr2\t8000000\tD4\tN\tN[chr1:6000000[\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=6000000;CT=3to5;PE=7;SR=2\n")

    empty = write(tmp, "empty.vcf", HDR_DELLY)

    # ── convention detection ────────────────────────────────────────────────
    print("\nconvention detection")
    d = V.detect_convention(delly)
    check("DELLY file detected as chr2_pos2", d.get("convention") == "chr2_pos2", str(d))
    m = V.detect_convention(mate)
    check("MATEID file detected as mateid", m.get("convention") == "mateid", str(m))
    b = V.detect_convention(bad)
    check("unknown convention fails loudly, not silently",
          "error" in b and b.get("error_type") == "unknown_breakend_convention", str(b))
    check("the failure message names the missing fields",
          "MATEID" in b.get("error", "") and "CHR2" in b.get("error", ""), b.get("error", ""))
    e = V.detect_convention(empty)
    check("empty file is 'undetermined', not an error",
          e.get("convention") == "undetermined", str(e))
    check("detection reports record-level evidence, not just header",
          d.get("evidence", {}).get("records_with_chr2_pos2", 0) == 2, str(d.get("evidence")))

    # ── parsing and canonical ordering ──────────────────────────────────────
    print("\nparsing and canonical ordering")
    pd = V.parse_junctions(delly)
    check("DELLY parse succeeds", "error" not in pd, str(pd)[:200])
    check("DELLY yields 3 junctions (2 BND + 1 DEL)", len(pd["junctions"]) == 3,
          str(len(pd["junctions"])))
    bnd = [j for j in pd["junctions"] if j.svtype == "BND"]
    check("canonical ordering places the higher contig first",
          all(V._contig_index(j.chrom1) >= V._contig_index(j.chrom2) for j in bnd))
    check("DELLY records are already canonical, so ordering is a no-op",
          all(j.chrom1 in ("chr2", "chr17") for j in bnd))

    pm = V.parse_junctions(mate)
    check("MATEID parse succeeds", "error" not in pm, str(pm)[:200])
    check("two MATEID records collapse to ONE junction", len(pm["junctions"]) == 1,
          f"got {len(pm['junctions'])}")

    # the cross-caller equivalence this whole layer exists for
    j_delly = [j for j in bnd if j.chrom2 == "chr1" and j.pos2 == 1000000][0]
    j_mate = pm["junctions"][0]
    check("same junction from both conventions normalises identically",
          (j_delly.chrom1, j_delly.pos1, j_delly.chrom2, j_delly.pos2)
          == (j_mate.chrom1, j_mate.pos1, j_mate.chrom2, j_mate.pos2),
          f"{j_delly.chrom1}:{j_delly.pos1}-{j_delly.chrom2}:{j_delly.pos2} vs "
          f"{j_mate.chrom1}:{j_mate.pos1}-{j_mate.chrom2}:{j_mate.pos2}")
    check("and therefore gets the SAME candidate_id across callers",
          j_delly.candidate_id == j_mate.candidate_id,
          f"{j_delly.candidate_id} vs {j_mate.candidate_id}")

    check("PE/SR absent is None, not 0",
          V.parse_junctions(write(tmp, "nope.vcf", HDR_DELLY +
              "chr2\t100000\tN1\tN\tN[chr1:5000[\t.\tPASS\tSVTYPE=BND;CHR2=chr1;POS2=5000;CT=3to5\n"
          ))["junctions"][0].pe is None)

    pe_ = V.parse_junctions(empty)
    check("empty file parses to zero junctions without error",
          "error" not in pe_ and len(pe_["junctions"]) == 0, str(pe_)[:200])

    # ── deduplication ───────────────────────────────────────────────────────
    print("\ndeduplication")
    pdup = V.parse_junctions(dup, dedup_tolerance_bp=500)
    check("three near-identical records collapse to one, far one kept",
          len(pdup["junctions"]) == 2, f"got {len(pdup['junctions'])}")
    merged = [j for j in pdup["junctions"] if j.n_merged > 1]
    check("the merged junction records how many it absorbed",
          len(merged) == 1 and merged[0].n_merged == 3,
          str([j.n_merged for j in pdup["junctions"]]))
    check("representative is the PASS record, not a LowQual one",
          merged and merged[0].filter == "PASS", str(merged[0].filter) if merged else "none")
    check("representative keeps the best support",
          merged and merged[0].pe == 9 and merged[0].sr == 5,
          f"pe={merged[0].pe} sr={merged[0].sr}" if merged else "none")
    nodup = V.parse_junctions(dup, dedup_tolerance_bp=0)
    check("dedup_tolerance_bp=0 disables merging", len(nodup["junctions"]) == 4,
          f"got {len(nodup['junctions'])}")

    # ── masking ─────────────────────────────────────────────────────────────
    print("\nmasking")
    mk = V.load_mask(maskfile)
    check("mask file loads", "error" not in mk, str(mk)[:150])
    check("region mask matches inside", V.in_mask(mk, "chr1", 500))
    check("region mask does not match outside", not V.in_mask(mk, "chr1", 50000))
    check("whole-contig mask matches", V.in_mask(mk, "chrUn_KI270302v1", 1))

    # ── the four tools ──────────────────────────────────────────────────────
    print("\nload_candidate_set")
    V.reset_registry()
    r = V.load_candidate_set(delly, "FIXTURE_DELLY")
    check("load succeeds and returns a set_id", "error" not in r and r["set_id"].startswith("CSET_"), str(r)[:200])
    check("reports the detected convention", r["caller_convention"] == "chr2_pos2")
    check("reports counts by SVTYPE", r["counts_by_svtype"] == {"BND": 2, "DEL": 1}, str(r["counts_by_svtype"]))
    check("reports counts by FILTER", r["counts_by_filter"] == {"PASS": 2, "LowQual": 1}, str(r["counts_by_filter"]))
    check("echoes the dedup threshold with provenance",
          any(t["name"] == "dedup_tolerance_bp" and t["provenance"] == "author judgement"
              for t in r["thresholds_applied"]), str(r["thresholds_applied"]))
    blob = str(r)
    check("returns NO file path", tmp not in blob and "delly.vcf" not in blob)
    check("returns no coordinate-shaped field",
          not any(k in r for k in ("chrom1", "pos1", "candidates", "junctions")))
    check("missing file is an error, not a crash",
          V.load_candidate_set(os.path.join(tmp, "definitely_absent.vcf"), "X").get("error_type") == "vcf_access")
    check("empty label rejected",
          V.load_candidate_set(delly, "").get("error_type") == "bad_parameter")
    check("malformed file surfaces the parse error through load",
          V.load_candidate_set(bad, "BAD").get("error_type") == "unknown_breakend_convention")

    sid = r["set_id"]
    check("set_id is stable for the same path+label",
          V.load_candidate_set(delly, "FIXTURE_DELLY")["set_id"] == sid)

    print("\nlist_candidates")
    lst = V.list_candidates(sid)
    check("unfiltered list returns everything", lst["total_matching"] == 3 and lst["total_in_set"] == 3)
    f1 = V.list_candidates(sid, filter_pass=True, min_pe=3, min_sr=1)
    check("filters reduce the set", f1["total_matching"] == 2, str(f1["total_matching"]))
    check("total_in_set still reports the WHOLE set alongside the filtered count",
          f1["total_in_set"] == 3, str(f1["total_in_set"]))
    check("every applied filter is echoed with a value",
          {s["filter"] for s in f1["filters_applied"]} == {"filter_pass", "min_pe", "min_sr"},
          str([s["filter"] for s in f1["filters_applied"]]))
    check("each filter carries a provenance label",
          all(s["provenance"] in ("tool-defined", "reference-defined", "author judgement")
              for s in f1["filters_applied"]))
    check("author-judgement thresholds are labelled as such",
          all(s["provenance"] == "author judgement"
              for s in f1["filters_applied"] if s["filter"] in ("min_pe", "min_sr")))
    check("each step reports its standalone effect on the unfiltered set",
          all("would_remove_from_unfiltered_set" in s for s in f1["filters_applied"]))
    check("unknown set_id is an error",
          V.list_candidates("CSET_nope").get("error_type") == "unknown_set_id")
    check("bad limit rejected", V.list_candidates(sid, limit=0).get("error_type") == "bad_parameter")
    check("bad offset rejected", V.list_candidates(sid, offset=-1).get("error_type") == "bad_parameter")
    pg = V.list_candidates(sid, limit=1, offset=0)
    check("pagination reports truncation", pg["returned"] == 1 and pg["truncated"] is True)
    pg2 = V.list_candidates(sid, limit=1, offset=2)
    check("last page is not marked truncated", pg2["truncated"] is False)
    check("list returns no file path", tmp not in str(lst))

    # masking through the tool
    V.reset_registry()
    rm = V.load_candidate_set(masked, "FIXTURE_MASKED")
    mres = V.list_candidates(rm["set_id"], mask_path=maskfile)
    check("mask filter removes the junction with a masked breakend",
          mres["total_matching"] == 1 and mres["total_in_set"] == 2,
          f"{mres['total_matching']}/{mres['total_in_set']}")
    check("mask filter is labelled reference-defined",
          [s for s in mres["filters_applied"] if s["filter"] == "exclude_masked"][0]["provenance"]
          == "reference-defined")
    check("bad mask path is an error",
          V.list_candidates(rm["set_id"], mask_path=os.path.join(tmp, "nope.tsv")).get("error_type")
          == "mask_access")

    print("\nget_candidate")
    V.reset_registry()
    r = V.load_candidate_set(delly, "FIXTURE_DELLY"); sid = r["set_id"]
    some = V.list_candidates(sid)["candidates"][0]
    g = V.get_candidate(sid, some["candidate_id"])
    check("returns the junction", "error" not in g and g["candidate_id"] == some["candidate_id"])
    check("provides breakend_1/breakend_2 shaped for the evidence tools",
          set(g["breakend_1"]) == {"chromosome", "position"} and
          set(g["breakend_2"]) == {"chromosome", "position"}, str(g.get("breakend_1")))
    check("breakends match the normalised coordinates",
          g["breakend_1"]["position"] == g["pos1"] and g["breakend_2"]["position"] == g["pos2"])
    check("unknown candidate_id is an error",
          V.get_candidate(sid, "CAND_nope").get("error_type") == "unknown_candidate_id")
    check("get_candidate returns no file path", tmp not in str(g))
    check("candidate_id is stable across calls",
          V.get_candidate(sid, some["candidate_id"])["candidate_id"] == some["candidate_id"])


    # ── error messages must not carry the input path ────────────────────────
    # A corrupt file makes pysam raise "invalid file `<full path>`". Passing
    # that message through put a real filename -- and so a possible patient
    # identifier -- into a tool return. Regression test named for that.
    print("\nerror messages carry no filesystem path")
    leaky = write(tmp, "SECRETNAME-corrupt.vcf", "this is not a vcf at all\n")
    r_bad = V.load_candidate_set(leaky, "X")
    check("corrupt file is an error, not a crash", r_bad.get("error_type") == "vcf_access", str(r_bad))
    check("the error does NOT echo the file path",
          "SECRETNAME" not in str(r_bad) and tmp not in str(r_bad), str(r_bad.get("error")))
    check("the error still says something useful",
          "withheld" in str(r_bad.get("error")) or "VCF" in str(r_bad.get("error")),
          str(r_bad.get("error")))
    r_mask = V.list_candidates(sid, mask_path=os.path.join(tmp, "SECRETMASK.tsv"))
    check("mask error does NOT echo the mask path",
          "SECRETMASK" not in str(r_mask), str(r_mask.get("error")))


    # ── coordinate provenance ───────────────────────────────────────────────
    # The evidence tools take chromosome/position as free parameters, so an
    # invented coordinate and a caller-derived one produced identical-looking
    # output. get_candidate now records what it hands out; the tools label the
    # difference. A NEAR-MISS must not be laundered into "candidate_set".
    print("\ncoordinate provenance")
    V.reset_registry()
    rp = V.load_candidate_set(delly, "PROV")
    cidp = V.list_candidates(rp["set_id"])["candidates"][0]["candidate_id"]
    jp = V.get_candidate(rp["set_id"], cidp)
    check("before get_candidate a coordinate is caller_supplied",
          V.position_provenance("chr9", position=12345)["source"] == "caller_supplied")
    pv = V.position_provenance(jp["chrom1"], position=jp["pos1"])
    check("after get_candidate the exact coordinate is candidate_set",
          pv["source"] == "candidate_set" and pv["candidate_id"] == cidp, str(pv))
    check("an off-by-one coordinate is NOT laundered",
          V.position_provenance(jp["chrom1"], position=jp["pos1"] + 1)["source"]
          == "caller_supplied")
    _cs = V.position_provenance("chr9", position=1)
    check("caller_supplied carries candidate_id None",
          _cs["candidate_id"] is None and _cs["set_id"] is None, str(_cs))
    check("caller_supplied note warns against presenting it as confirmation",
          "confirmation" in _cs["note"].lower(), _cs["note"])
    check("reset_registry clears recorded coordinates",
          (V.reset_registry() or V.position_provenance(jp["chrom1"], position=jp["pos1"])["source"])
          == "caller_supplied")

    # ── standalone counts follow the svtype scope of the query ──────────────
    # Previously measured against the whole registered set, so a BND query
    # reported a figure that counted DEL/DUP/INV records too.
    print("\nstandalone filter counts respect svtype scope")
    V.reset_registry()
    rs = V.load_candidate_set(delly, "SCOPE")
    allq = V.list_candidates(rs["set_id"], filter_pass=True)
    bndq = V.list_candidates(rs["set_id"], svtype="BND", filter_pass=True)
    a_step = [x for x in allq["filters_applied"] if x["filter"] == "filter_pass"][0]
    b_step = [x for x in bndq["filters_applied"] if x["filter"] == "filter_pass"][0]
    check("unscoped query measures against the whole set",
          a_step["unfiltered_set_size"] == 3, str(a_step))
    check("BND-scoped query measures against BND records only",
          b_step["unfiltered_set_size"] == 2, str(b_step))
    check("the two therefore give different standalone counts, correctly",
          a_step["would_remove_from_unfiltered_set"] == 1 and
          b_step["would_remove_from_unfiltered_set"] == 1,
          f"all={a_step['would_remove_from_unfiltered_set']} bnd={b_step['would_remove_from_unfiltered_set']}")
    check("the scope is stated in the response",
          b_step["unfiltered_set_scope"] == "svtype=BND", str(b_step.get("unfiltered_set_scope")))
    sv_step = [x for x in bndq["filters_applied"] if x["filter"] == "svtype"][0]
    check("the svtype filter itself is measured against the whole set",
          sv_step["unfiltered_set_size"] == 3 and sv_step["unfiltered_set_scope"] == "all svtypes",
          str(sv_step))

    print("\ncompare_candidate_sets")
    V.reset_registry()
    ra = V.load_candidate_set(delly, "SET_ONE")
    rb = V.load_candidate_set(mate, "SET_TWO")
    cmp_ = V.compare_candidate_sets(ra["set_id"], rb["set_id"], tolerance_bp=500)
    check("the shared junction matches across the two conventions",
          cmp_["matched_in_a"] == 1 and cmp_["matched_in_b"] == 1,
          f"a={cmp_['matched_in_a']} b={cmp_['matched_in_b']}")
    check("unmatched counts are reported too",
          cmp_["unmatched_in_a"] == 2 and cmp_["unmatched_in_b"] == 0,
          f"{cmp_['unmatched_in_a']}/{cmp_['unmatched_in_b']}")
    check("matched pairs carry both candidate_ids",
          cmp_["matched_pairs"] and set(cmp_["matched_pairs"][0]) == {"candidate_id_a", "candidate_id_b"})
    check("tolerance echoed with provenance",
          cmp_["thresholds_applied"][0]["name"] == "tolerance_bp" and
          cmp_["thresholds_applied"][0]["provenance"] == "author judgement")
    check("tolerance 0 still matches an exact-coordinate junction",
          V.compare_candidate_sets(ra["set_id"], rb["set_id"], 0)["matched_in_a"] == 1)
    check("unknown set_id is an error",
          V.compare_candidate_sets("CSET_x", rb["set_id"]).get("error_type") == "unknown_set_id")
    check("negative tolerance rejected",
          V.compare_candidate_sets(ra["set_id"], rb["set_id"], -1).get("error_type") == "bad_parameter")

    # ── dedup orientation key (Phase 6 Task A) ──────────────────────────────
    # The two junctions of a balanced reciprocal translocation sit ~1bp apart
    # on BOTH breakends and differ ONLY in orientation. Before the fix they
    # merged into one candidate and the second's PE/SR was discarded.
    # t(1;2) at chr1:1000000 / chr2:3000000:
    #   der(1) = chr1[..1000000] + chr2[3000001..]  -> canonical chr2:3000001 / chr1:1000000, 5to3
    #   der(2) = chr2[..3000000] + chr1[1000001..]  -> canonical chr2:3000000 / chr1:1000001, 3to5
    recip_delly = write(tmp, "recip_delly.vcf", HDR_DELLY +
        "chr1\t1000000\tR1\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr2;POS2=3000001;CT=3to5;PE=20;SR=10\n"
        "chr2\t3000000\tR2\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1000001;CT=3to5;PE=18;SR=9\n")
    rd = V.parse_junctions(recip_delly)
    check("reciprocal junctions stay distinct (chr2_pos2)",
          len(rd["junctions"]) == 2, f'got {len(rd["junctions"])}')
    check("reciprocal dedup merges nothing (chr2_pos2)",
          rd["meta"]["records_merged_by_dedup"] == 0)
    check("both reciprocal orientations preserved (chr2_pos2)",
          sorted(j.orientation for j in rd["junctions"]) == ["3to5", "5to3"],
          str(sorted(j.orientation for j in rd["junctions"])))
    check("second junction's evidence not discarded (chr2_pos2)",
          sorted((j.pe, j.sr) for j in rd["junctions"]) == [(18, 9), (20, 10)])

    # Same physical event in MATEID convention must normalise identically.
    recip_mate = write(tmp, "recip_mate.vcf", HDR_MATE +
        "chr1\t1000000\tP1a\tN\tN[chr2:3000001[\t.\tPASS\tSVTYPE=BND;MATEID=P1b;PE=20;SR=10\n"
        "chr2\t3000001\tP1b\tN\t]chr1:1000000]N\t.\tPASS\tSVTYPE=BND;MATEID=P1a;PE=20;SR=10\n"
        "chr2\t3000000\tP2a\tN\tN[chr1:1000001[\t.\tPASS\tSVTYPE=BND;MATEID=P2b;PE=18;SR=9\n"
        "chr1\t1000001\tP2b\tN\t]chr2:3000000]N\t.\tPASS\tSVTYPE=BND;MATEID=P2a;PE=18;SR=9\n")
    rm = V.parse_junctions(recip_mate)
    check("reciprocal junctions stay distinct (mateid)",
          len(rm["junctions"]) == 2, f'got {len(rm["junctions"])}')

    def norm(res):
        return sorted((j.chrom1, j.pos1, j.chrom2, j.pos2, j.orientation)
                      for j in res["junctions"])
    check("mateid and chr2_pos2 normalise the same translocation identically",
          norm(rd) == norm(rm), f"{norm(rd)} vs {norm(rm)}")

    # Mate encounter order must not change the normalised result.
    recip_mate_rev = write(tmp, "recip_mate_rev.vcf", HDR_MATE +
        "chr2\t3000001\tQ1b\tN\t]chr1:1000000]N\t.\tPASS\tSVTYPE=BND;MATEID=Q1a;PE=20;SR=10\n"
        "chr1\t1000000\tQ1a\tN\tN[chr2:3000001[\t.\tPASS\tSVTYPE=BND;MATEID=Q1b;PE=20;SR=10\n"
        "chr1\t1000001\tQ2b\tN\t]chr2:3000000]N\t.\tPASS\tSVTYPE=BND;MATEID=Q2a;PE=18;SR=9\n"
        "chr2\t3000000\tQ2a\tN\tN[chr1:1000001[\t.\tPASS\tSVTYPE=BND;MATEID=Q2b;PE=18;SR=9\n")
    check("mate encounter order does not change normalisation",
          norm(V.parse_junctions(recip_mate_rev)) == norm(rm))

    # Genuine near-duplicates with the SAME orientation must still merge.
    neardup = write(tmp, "neardup_same_ct.vcf", HDR_DELLY +
        "chr2\t3000000\tN1\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1000000;CT=3to5;PE=9;SR=4\n"
        "chr2\t3000200\tN2\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1000150;CT=3to5;PE=7;SR=2\n")
    nd = V.parse_junctions(neardup)
    check("near-duplicates with identical orientation still merge",
          len(nd["junctions"]) == 1, f'got {len(nd["junctions"])}')
    check("near-duplicate merge still records n_merged=2",
          nd["junctions"][0].n_merged == 2)

    # Orientation flip under canonical reordering: the SAME junction written
    # from either end must normalise to one orientation and still merge.
    flip = write(tmp, "flip_canonical.vcf", HDR_DELLY +
        "chr2\t3000000\tF1\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1000000;CT=3to5;PE=9;SR=4\n"
        "chr1\t1000000\tF2\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr2;POS2=3000000;CT=5to3;PE=8;SR=3\n")
    fl = V.parse_junctions(flip)
    check("junction written from either end merges after orientation flip",
          len(fl["junctions"]) == 1, f'got {len(fl["junctions"])}')
    check("flip normalises to the canonical orientation",
          fl["junctions"][0].orientation == "3to5",
          str(fl["junctions"][0].orientation))

    # Negative control: same coordinates, genuinely different orientation
    # (written from the far end WITHOUT the matching flip) must NOT merge.
    nomerge = write(tmp, "flip_negative.vcf", HDR_DELLY +
        "chr2\t3000000\tG1\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1000000;CT=3to5;PE=9;SR=4\n"
        "chr1\t1000000\tG2\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr2;POS2=3000000;CT=3to5;PE=8;SR=3\n")
    nm = V.parse_junctions(nomerge)
    check("genuinely different orientations at same coords do NOT merge",
          len(nm["junctions"]) == 2, f'got {len(nm["junctions"])}')

    # ── compare_candidate_sets symmetry (Phase 7) ───────────────────────────
    # The scan used to stop at the first partner found for each set_a junction,
    # so matched_in_b counted only B junctions that happened to be someone's
    # first match. A B junction with a genuine partner could be reported
    # unmatched because a neighbour was reached first, making "does A recur in
    # B" and "does B recur in A" disagree. Recurrence gates a funnel step, so
    # the answer must not depend on argument order.
    # Junctions are spaced 800bp apart so the 500bp dedup keeps them distinct,
    # then compared at 1000bp so each query junction has more than one partner --
    # which is the only situation the first-match break ever got wrong.
    sym_a = write(tmp, "sym_a.vcf", HDR_DELLY +
        "chr2\t3000000\tSA1\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1000000;CT=3to5;PE=9;SR=4\n"
        "chr2\t3000800\tSA2\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1000800;CT=3to5;PE=8;SR=3\n")
    sym_b = write(tmp, "sym_b.vcf", HDR_DELLY +
        "chr2\t3000400\tSB1\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1000400;CT=3to5;PE=7;SR=2\n"
        "chr2\t3001200\tSB2\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=1001200;CT=3to5;PE=6;SR=1\n"
        # no partner in sym_a at all
        "chr17\t500000\tSB3\tN\t<BND>\t.\tPASS\t"
        "SVTYPE=BND;CHR2=chr1;POS2=9000000;CT=3to5;PE=5;SR=1\n")
    V.reset_registry()
    ra_ = V.load_candidate_set(sym_a, "SYMA")
    rb_ = V.load_candidate_set(sym_b, "SYMB")
    check("symmetry fixture survives dedup (2 and 3 junctions)",
          (ra_["junctions_after_dedup"], rb_["junctions_after_dedup"]) == (2, 3),
          f'{ra_["junctions_after_dedup"]},{rb_["junctions_after_dedup"]}')
    fwd = V.compare_candidate_sets(ra_["set_id"], rb_["set_id"], 1000)
    rev = V.compare_candidate_sets(rb_["set_id"], ra_["set_id"], 1000)
    check("compare is symmetric: matched counts mirror when the order swaps",
          (fwd["matched_in_a"], fwd["matched_in_b"]) == (rev["matched_in_b"], rev["matched_in_a"]),
          f'fwd=({fwd["matched_in_a"]},{fwd["matched_in_b"]}) rev=({rev["matched_in_a"]},{rev["matched_in_b"]})')
    check("compare is symmetric: unmatched counts mirror too",
          (fwd["unmatched_in_a"], fwd["unmatched_in_b"]) == (rev["unmatched_in_b"], rev["unmatched_in_a"]))
    check("a junction with two partners is not left uncounted on the b side",
          fwd["matched_in_b"] == 2, f'matched_in_b={fwd["matched_in_b"]}')
    check("every matching pair is emitted, not just the first per query junction",
          len(fwd["matched_pairs"]) == 3, f'pairs={len(fwd["matched_pairs"])}')
    check("a junction with genuinely no partner is still reported unmatched",
          fwd["unmatched_in_b"] == 1, f'unmatched_in_b={fwd["unmatched_in_b"]}')
    check("pair ids reference real candidates in both sets",
          all(p["candidate_id_a"] and p["candidate_id_b"] for p in fwd["matched_pairs"]))

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL VCF TOOL TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
