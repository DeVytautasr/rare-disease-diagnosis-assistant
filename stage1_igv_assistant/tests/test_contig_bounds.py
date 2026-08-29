#!/usr/bin/env python3
"""
Regression tests for windows that overrun the end of a contig.

From results/REAL_PATIENT_DATA_VALIDATION.md findings 5 and 14. They share
one chokepoint (_clamp_or_error) because they are the same missing check seen
at two distances past the edge.

Finding 5 — a window overrunning the contig end was neither clamped nor
flagged, and the FULL requested width was reported. chrM:16000-17500 on a
16,569 bp contig returned the same reads as 16000-16569 but mean_depth
1,011.06 instead of 2,665.36: the denominator included 931 bases that do not
exist. Depth profiles emitted whole bins past the end (9 of 20 at chrM:16400,
20 of 20 at a past-end locus, 9,918 of 10,000 for chrM:0-2,000,000). Those
empty bins then scored a maximum depth component with a "possible deletion"
sentence attached.

Finding 14 — a coordinate entirely past the end returned a clean zero,
indistinguishable from a genuine coverage gap. pysam does not raise for an
out-of-range fetch, so nothing downstream noticed, and a typo'd coordinate
looked exactly like real data.

The quantitative assertion is the important one: mean depth over (0, L) must
equal mean depth over (0, L + slack). If it does not, the denominator is
still counting bases the reference does not have.

Run: python3 stage1_igv_assistant/tests/test_contig_bounds.py
"""
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant.tools.bam_tools import (  # noqa: E402
    get_bam_stats_at_locus,
    get_read_depth_profile,
    count_discordant_pairs,
)

FAILURES = []

SHORT_LEN = 16569           # a chrM-sized contig, where this actually bit
HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chrM", "LN": SHORT_LEN}, {"SN": "chr1", "LN": 248956422}],
})


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def _even_coverage(path):
    """Uniform coverage across the whole short contig, so any dilution shows."""
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        n = 0
        for pos in range(0, SHORT_LEN - 100, 20):
            r = pysam.AlignedSegment(HEADER)
            r.query_name = f"m{n}"
            r.query_sequence = "A" * 100
            r.flag = 0
            r.reference_id = 0
            r.reference_start = pos
            r.mapping_quality = 60
            r.cigar = [(0, 100)]
            r.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(r)
            n += 1
    pysam.index(path)


def run_tests():
    print("=" * 68)
    print("CONTIG BOUNDS REGRESSION SUITE")
    print("=" * 68)

    tmp = tempfile.mkdtemp()
    bam = os.path.join(tmp, "short.bam")
    _even_coverage(bam)

    # ── Finding 5: the denominator must not include absent bases ─────────
    print("\ndepth is not diluted by bases the contig does not have")
    exact = get_bam_stats_at_locus(bam, "chrM", 0, SHORT_LEN)
    over = get_bam_stats_at_locus(bam, "chrM", 0, SHORT_LEN + 10000)
    check("same read count either way",
          exact["total_reads"] == over["total_reads"],
          f"{exact['total_reads']} vs {over['total_reads']}")
    check("mean_depth identical whether or not the request overran",
          exact["mean_depth"] == over["mean_depth"],
          f"exact={exact['mean_depth']} over={over['mean_depth']} "
          f"(dilution ratio would be {SHORT_LEN / (SHORT_LEN + 10000):.3f})")

    print("\nthe reported span is the clamped one, not the requested one")
    check("end is clamped to the contig length",
          over["end"] == SHORT_LEN, f"got {over['end']}")
    check("clamped_to_contig is flagged",
          over["clamped_to_contig"] is True, f"got {over['clamped_to_contig']}")
    check("contig_length is reported",
          over["contig_length"] == SHORT_LEN, f"got {over['contig_length']}")
    check("an in-bounds window is not flagged as clamped",
          exact["clamped_to_contig"] is False, f"got {exact['clamped_to_contig']}")

    print("\ndepth profiles emit no bins past the contig end")
    prof = get_read_depth_profile(bam, "chrM", SHORT_LEN - 1000, SHORT_LEN + 4000,
                                  window_size=200)
    wins = prof.get("windows", [])
    past = [w for w in wins if w.get("start", 0) >= SHORT_LEN]
    check("profile returned bins", len(wins) > 0, f"got {len(wins)}")
    check("no bin starts past the contig end",
          not past, f"{len(past)} of {len(wins)} bins past {SHORT_LEN}")
    check("profile flags the clamp",
          prof["summary"]["clamped_to_contig"] is True,
          f"got {prof['summary']['clamped_to_contig']}")

    # ── Finding 14: past the end is an error, not a coverage gap ─────────
    print("\na coordinate past the contig end is out_of_range, not empty data")
    past_end = get_bam_stats_at_locus(bam, "chrM", SHORT_LEN + 10, SHORT_LEN + 500)
    check("error_type is out_of_range",
          past_end.get("error_type") == "out_of_range",
          f"got {past_end.get('error_type')!r}")
    check("the error names the contig length",
          past_end.get("contig_length") == SHORT_LEN,
          f"got {past_end.get('contig_length')}")
    check("it does not masquerade as a zero-coverage result",
          "total_reads" not in past_end, f"keys={sorted(past_end)}")

    disc_past = count_discordant_pairs(bam, "chrM", SHORT_LEN + 5000)
    check("windowed tools also report out_of_range",
          disc_past.get("error_type") == "out_of_range",
          f"got {disc_past.get('error_type')!r}")

    print("\nan unknown contig is still a different error")
    unknown = get_bam_stats_at_locus(bam, "chrNOPE", 100, 200)
    check("unknown contig stays invalid_region",
          unknown.get("error_type") == "invalid_region",
          f"got {unknown.get('error_type')!r}")

    print("\na genuine coverage gap is still reported as data, not an error")
    gap = get_bam_stats_at_locus(bam, "chr1", 1000, 2000)
    check("in-bounds empty region returns a result, not out_of_range",
          "error_type" not in gap and gap["total_reads"] == 0,
          f"got {gap.get('error_type')!r} reads={gap.get('total_reads')}")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL CONTIG BOUNDS TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
