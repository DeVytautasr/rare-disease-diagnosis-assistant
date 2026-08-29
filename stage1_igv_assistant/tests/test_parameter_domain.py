#!/usr/bin/env python3
"""
Regression tests for the three parameter-domain crashes.

From results/REAL_PATIENT_DATA_VALIDATION.md. None of them fired at a real
coordinate — every one is reachable only through a malformed argument, which
is exactly why they survived: the suites passed and real data never hit them.

  1. get_read_depth_profile(window_size=0)
       ValueError: range() arg 3 must not be zero
  2. get_read_depth_profile(window_size=-100)
       IndexError: list index out of range
       Data-dependent: over a zero-coverage region the identical call
       returned a clean result, so it only crashed where reads existed.
  3. get_bam_stats_at_locus(start="...", end="...")
       TypeError: '<' not supported between instances of 'str' and 'int'
       Raised out of _validate_range itself — the one function whose whole
       job is to turn bad input into a structured error dict.

Every other malformed input in this module returns a structured error dict.
These three escaped that contract and raised out of the tool, which for an
MCP server means a transport-level failure rather than a result the caller
can read and act on.

Run: python3 stage1_igv_assistant/tests/test_parameter_domain.py
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
    summarize_breakpoint_evidence,
)

FAILURES = []

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 1000000}],
})


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def _covered_bam(path):
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        for n, pos in enumerate(range(10000, 20000, 25)):
            r = pysam.AlignedSegment(HEADER)
            r.query_name = f"r{n}"
            r.query_sequence = "A" * 100
            r.flag = 0
            r.reference_id = 0
            r.reference_start = pos
            r.mapping_quality = 60
            r.cigar = [(0, 100)]
            r.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(r)
    pysam.index(path)


def expect_error_dict(label, fn, expected_type="invalid_parameters"):
    """A bad argument must produce a dict, never an exception."""
    try:
        result = fn()
    except Exception as exc:
        check(label, False, f"raised {type(exc).__name__}: {exc}")
        return
    check(f"{label}: returns a dict",
          isinstance(result, dict), f"got {type(result).__name__}")
    check(f"{label}: error_type is {expected_type}",
          result.get("error_type") == expected_type,
          f"got {result.get('error_type')!r} — {result.get('error', '')[:70]}")


def run_tests():
    print("=" * 68)
    print("PARAMETER DOMAIN REGRESSION SUITE")
    print("=" * 68)

    tmp = tempfile.mkdtemp()
    bam = os.path.join(tmp, "cov.bam")
    _covered_bam(bam)

    print("\ncrash 1 & 2: degenerate window_size")
    expect_error_dict("window_size=0",
                      lambda: get_read_depth_profile(bam, "chr1", 10000, 12000, window_size=0))
    expect_error_dict("window_size=-100",
                      lambda: get_read_depth_profile(bam, "chr1", 10000, 12000, window_size=-100))
    # The negative case only crashed where reads existed; assert the covered
    # region specifically, not just an empty one that never reproduced it.
    expect_error_dict("window_size=-1 over a covered region",
                      lambda: get_read_depth_profile(bam, "chr1", 10000, 20000, window_size=-1))
    expect_error_dict("window_size as a float",
                      lambda: get_read_depth_profile(bam, "chr1", 10000, 12000, window_size=1.5))

    print("\ncrash 3: non-integer coordinates")
    expect_error_dict("string start/end",
                      lambda: get_bam_stats_at_locus(bam, "chr1", "10000", "12000"))
    expect_error_dict("float start",
                      lambda: get_bam_stats_at_locus(bam, "chr1", 10000.5, 12000))
    expect_error_dict("None end",
                      lambda: get_bam_stats_at_locus(bam, "chr1", 10000, None))
    expect_error_dict("string position on a windowed tool",
                      lambda: count_discordant_pairs(bam, "chr1", "15000"))
    # Found by this suite rather than by the validation run: the windowed
    # tools compute position - window_bp before _validate_range sees either
    # value, so the subtraction raised first. summarize does its own
    # arithmetic and needed the same guard.
    expect_error_dict("string position on summarize",
                      lambda: summarize_breakpoint_evidence(bam, "chr1", "15000"))
    expect_error_dict("None window_bp on summarize",
                      lambda: summarize_breakpoint_evidence(bam, "chr1", 15000, window_bp=None))
    expect_error_dict("float window_bp",
                      lambda: count_discordant_pairs(bam, "chr1", 15000, window_bp=2.5))

    print("\nvalid parameters still work (the guards must not over-reject)")
    ok = get_read_depth_profile(bam, "chr1", 10000, 12000, window_size=200)
    check("a normal depth profile still succeeds",
          "error_type" not in ok and len(ok.get("windows", [])) > 0,
          f"got {ok.get('error_type')!r}")
    ok2 = get_bam_stats_at_locus(bam, "chr1", 10000, 12000)
    check("a normal stats call still succeeds",
          "error_type" not in ok2 and ok2["total_reads"] > 0,
          f"got {ok2.get('error_type')!r}")
    check("window_size=1 is legal, not rejected as degenerate",
          "error_type" not in get_read_depth_profile(bam, "chr1", 10000, 10500, window_size=1))

    print("\nthe pre-existing structured errors are unchanged")
    expect_error_dict("negative start",
                      lambda: get_bam_stats_at_locus(bam, "chr1", -100, -50))
    expect_error_dict("start > end",
                      lambda: get_bam_stats_at_locus(bam, "chr1", 5000, 1000))
    expect_error_dict("zero-width region",
                      lambda: get_bam_stats_at_locus(bam, "chr1", 5000, 5000))
    expect_error_dict("unknown contig",
                      lambda: get_bam_stats_at_locus(bam, "chrNOPE", 100, 200),
                      expected_type="invalid_region")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL PARAMETER DOMAIN TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
