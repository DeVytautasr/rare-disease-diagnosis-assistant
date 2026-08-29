#!/usr/bin/env python3
"""
Regression test: one read must not fire the bottom scoring tier.

From results/REAL_PATIENT_DATA_VALIDATION.md finding 4, first half.

The discordant and split layers awarded 7.5/25 for ANY non-zero fraction:

    elif disc_fraction > 0:
        discordant_pair_score = 7.5

At 30x WGS a +/-500 bp window holds roughly 250 reads, and one read with an
interchromosomal mate is ordinary background. Measured over 42 arbitrary
control positions chosen only for ordinary coverage, 32 (76%) came back
"weak" rather than "none" — every discordant firing sitting at the 7.5 floor
off 1-3 reads with no dominant partner. An assistant reading
evidence_strength and skipping the prose would call three quarters of the
genome "weak evidence".

MIN_ABSOLUTE_SUPPORT = 3 is a judgement call, not a calibrated figure: the
smallest count that cannot be a single read plus one duplicate of it. It is
deliberately an ABSOLUTE count, not a fraction — a fraction cannot
distinguish 1-in-250 from 1-in-4.

The reads are still reported in supporting_observations. They just no longer
move the score.

Run: python3 stage1_igv_assistant/tests/test_minimum_support.py
"""
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant.tools.bam_tools import (  # noqa: E402
    summarize_breakpoint_evidence,
    MIN_ABSOLUTE_SUPPORT,
)

FAILURES = []

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 1000000}, {"SN": "chr8", "LN": 900000}],
})


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def _bam(path, n_disc=0, n_split=0, n_background=250):
    """Background concordant coverage plus a controlled number of signal reads."""
    recs = [(9800 + i * 2, f"n{i}", 0, 10000, None) for i in range(n_background)]
    recs += [(10000 + i, f"d{i}", 1, 500000, None) for i in range(n_disc)]
    recs += [(10000 + i, f"s{i}", 0, 10000, "chr8,500000,+,80M20S,60,0;")
             for i in range(n_split)]
    recs.sort(key=lambda x: x[0])
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        for pos, name, mate_ref, mate_pos, sa in recs:
            r = pysam.AlignedSegment(HEADER)
            r.query_name = name
            r.query_sequence = "A" * 100
            r.flag = 0x1 if (mate_ref == 1 or sa) else (0x1 | 0x2)
            r.reference_id = 0
            r.reference_start = pos
            r.mapping_quality = 60
            r.cigar = [(4, 20), (0, 80)] if sa else [(0, 100)]
            r.next_reference_id = mate_ref
            r.next_reference_start = mate_pos
            r.query_qualities = pysam.qualitystring_to_array("I" * 100)
            if sa:
                r.set_tag("SA", sa)
            bam.write(r)
    pysam.index(path)
    return path


def run_tests():
    print("=" * 68)
    print("MINIMUM SUPPORT REGRESSION SUITE")
    print("=" * 68)
    tmp = tempfile.mkdtemp()

    print(f"\nthe threshold is exposed (MIN_ABSOLUTE_SUPPORT = {MIN_ABSOLUTE_SUPPORT})")
    check("it is an integer read count, not a fraction",
          isinstance(MIN_ABSOLUTE_SUPPORT, int) and MIN_ABSOLUTE_SUPPORT == 3,
          f"got {MIN_ABSOLUTE_SUPPORT!r}")

    print("\ndiscordant: below the bar scores nothing, at the bar scores 7.5")
    for n, expect_score in [(0, 0.0), (1, 0.0), (2, 0.0), (3, 7.5), (5, 7.5)]:
        path = _bam(os.path.join(tmp, f"d{n}.bam"), n_disc=n)
        res = summarize_breakpoint_evidence(path, "chr1", 10000)
        check(f"{n} discordant read(s) -> score {expect_score}",
              res["discordant_pair_score"] == expect_score,
              f"got {res['discordant_pair_score']}")

    print("\nsplit: same rule")
    for n, expect_score in [(1, 0.0), (2, 0.0), (3, 7.5)]:
        path = _bam(os.path.join(tmp, f"s{n}.bam"), n_split=n)
        res = summarize_breakpoint_evidence(path, "chr1", 10000)
        check(f"{n} split read(s) -> score {expect_score}",
              res["split_read_score"] == expect_score,
              f"got {res['split_read_score']}")

    print("\nthe headline strength changes, which is the point")
    one = summarize_breakpoint_evidence(_bam(os.path.join(tmp, "one.bam"), n_disc=1),
                                        "chr1", 10000)
    three = summarize_breakpoint_evidence(_bam(os.path.join(tmp, "three.bam"), n_disc=3),
                                          "chr1", 10000)
    check("1 background read reads 'none', not 'weak'",
          one["evidence_strength"] == "none", f"got {one['evidence_strength']!r}")
    check("3 reads still reads 'weak'",
          three["evidence_strength"] == "weak", f"got {three['evidence_strength']!r}")

    print("\nsub-threshold reads are still REPORTED, just not scored")
    obs = " ".join(one["supporting_observations"])
    check("the single discordant read still appears in the observations",
          "discordant" in obs.lower(), f"got {obs[:200]!r}")
    check("the underlying count is still returned",
          one["discordant_pairs"]["discordant_pairs"] == 1
          if isinstance(one.get("discordant_pairs"), dict) else True,
          "count not exposed")

    print("\nthe upper tiers are untouched")
    many = summarize_breakpoint_evidence(
        _bam(os.path.join(tmp, "many.bam"), n_disc=200, n_background=200),
        "chr1", 10000)
    check("a genuinely high discordant fraction still scores above the floor",
          many["discordant_pair_score"] is not None and many["discordant_pair_score"] > 7.5,
          f"got {many['discordant_pair_score']}")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL MINIMUM SUPPORT TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
