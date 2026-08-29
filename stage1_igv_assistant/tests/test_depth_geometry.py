#!/usr/bin/env python3
"""
Regression test: the depth layer's geometry is explicit and reproducible.

From results/REAL_PATIENT_DATA_VALIDATION.md finding 10.

summarize_breakpoint_evidence hardcoded its depth geometry at +/-2000 with
200 bp bins, silently ignoring window_bp. The documented tool order has the
operator call get_read_depth_profile standalone first, so at any other
geometry the two disagreed: 4 of 10 validation positions gave opposite
likely_deletion verdicts, and dip_is_at_focus conflicted at 3 of 10, in both
directions. Nothing in either return said the standalone profile was not the
one that fed the score.

The cause is that depth_ratio_min_to_mean is strongly window-size dependent
while DEPTH_RATIO_DELETION_THRESHOLD is a fixed 0.7 calibrated at
window_size=200. Holding the span fixed, one real locus gave 0.523 at
window_size=100, 0.585 at 200, 0.830 at 500 and 0.868 at 1000 — crossing the
threshold between 200 and 500.

The geometry is now a parameter, defaulted to the old values so behaviour is
unchanged, and echoed in the result. The property that matters is
reproducibility: a caller who reads depth_window_bp/depth_window_size can
re-derive the same ratio from get_read_depth_profile.

Run: python3 stage1_igv_assistant/tests/test_depth_geometry.py
"""
import os
import re
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant.tools.bam_tools import (  # noqa: E402
    get_read_depth_profile,
    summarize_breakpoint_evidence,
    DEPTH_RATIO_DELETION_THRESHOLD,
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


def _bam_with_dip(path):
    """Coverage across chr1:8000-12000 with a gap at the focus position."""
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        n = 0
        for pos in range(8000, 12000, 25):
            if 9800 <= pos < 10200:
                continue
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
            n += 1
    pysam.index(path)


def run_tests():
    print("=" * 68)
    print("DEPTH GEOMETRY REGRESSION SUITE")
    print("=" * 68)

    tmp = tempfile.mkdtemp()
    bam = os.path.join(tmp, "dip.bam")
    _bam_with_dip(bam)

    print("\nthe defaults preserve the previous hardcoded behaviour")
    d = summarize_breakpoint_evidence(bam, "chr1", 10000)
    check("depth_window_bp defaults to 2000",
          d["depth_window_bp"] == 2000, f"got {d['depth_window_bp']}")
    check("depth_window_size defaults to 200",
          d["depth_window_size"] == 200, f"got {d['depth_window_size']}")

    print("\nthe geometry is honoured and echoed when overridden")
    o = summarize_breakpoint_evidence(bam, "chr1", 10000,
                                      depth_window_bp=1000, depth_window_size=100)
    check("override is applied", o["depth_window_bp"] == 1000
          and o["depth_window_size"] == 100,
          f"got {o['depth_window_bp']}/{o['depth_window_size']}")

    print("\nREPRODUCIBILITY: a caller can re-derive the ratio from the echoed geometry")
    for wbp, wsz in [(2000, 200), (1000, 100), (3000, 500)]:
        summ = summarize_breakpoint_evidence(bam, "chr1", 10000,
                                             depth_window_bp=wbp, depth_window_size=wsz)
        standalone = get_read_depth_profile(
            bam, "chr1", max(0, 10000 - summ["depth_window_bp"]),
            10000 + summ["depth_window_bp"],
            window_size=summ["depth_window_size"], focus_position=10000)
        inner = summ["depth_profile"]["summary"]
        check(f"geometry {wbp}/{wsz}: standalone ratio matches summarize's",
              standalone["summary"]["depth_ratio_min_to_mean"]
              == inner["depth_ratio_min_to_mean"],
              f"standalone={standalone['summary']['depth_ratio_min_to_mean']} "
              f"summarize={inner['depth_ratio_min_to_mean']}")
        check(f"geometry {wbp}/{wsz}: likely_deletion agrees",
              standalone["summary"]["likely_deletion"] == inner["likely_deletion"],
              f"standalone={standalone['summary']['likely_deletion']} "
              f"summarize={inner['likely_deletion']}")
        check(f"geometry {wbp}/{wsz}: the echoed geometry is what was used",
              summ["depth_profile"]["window_size"] == wsz,
              f"echoed={wsz} actual={summ['depth_profile']['window_size']}")

    print("\nwindow-size dependence is real (why the threshold needs its geometry)")
    ratios = {}
    for wsz in (100, 200, 500, 1000):
        prof = get_read_depth_profile(bam, "chr1", 8000, 12000,
                                      window_size=wsz, focus_position=10000)
        ratios[wsz] = prof["summary"]["depth_ratio_min_to_mean"]
    print(f"    same span, varying bin size: {ratios}")
    check("the ratio genuinely changes with bin size",
          len(set(ratios.values())) > 1, f"got {ratios}")

    print("\nthe coupling is recorded with the threshold constant")
    src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "tools", "bam_tools.py"), encoding="utf-8").read()
    # Locate the ASSIGNMENT, not the first mention: the name also appears in
    # a provenance comment far earlier in the file, and indexing on that read
    # the wrong preamble entirely.
    m = re.search(r"^DEPTH_RATIO_DELETION_THRESHOLD = 0\.7", src, re.MULTILINE)
    check("the threshold assignment was located", m is not None, "not found")
    preamble = src[max(0, m.start() - 1200):m.start()] if m else ""
    check("a window-size coupling note precedes the constant",
          "WINDOW-SIZE COUPLING" in preamble,
          "no coupling note found above the constant")
    check("the note names the calibration window size",
          "window_size=200" in preamble, "calibration bin size not stated")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL DEPTH GEOMETRY TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
