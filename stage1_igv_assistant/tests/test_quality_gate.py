#!/usr/bin/env python3
"""
Regression tests for the evidence-score quality gate.

Named for the condition that exposed them. Real-data validation against two
~39 GB WGS BAMs found that the headline evidence score RISES as mapping
quality falls (results/REAL_PATIENT_DATA_VALIDATION.md, findings 1 and 2):

  - In a pericentromeric window where every read is MAPQ 0, the discordant,
    soft-clip and split layers each reported "no reads in window" and were
    dropped from evidence_score's denominator.
  - get_read_depth_profile applies no MAPQ filter, so the depth layer kept
    scoring.
  - Denominator 1, numerator 25/25, result: 100.0/100 "strong" from a raw
    score of 25.0 -- on a window with no usable reads at all. Sixteen of
    twenty-eight such cells scored "strong".
  - The reason string was also wrong: "no reads in window" was emitted at
    windows holding up to 452 reads.

These tests pin the fix. They use synthetic BAMs -- the patient data is not
in this repository and never will be.

Run: python3 stage1_igv_assistant/tests/test_quality_gate.py
"""
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant.tools.bam_tools import (  # noqa: E402
    count_discordant_pairs,
    count_soft_clipped_reads,
    get_split_reads,
    summarize_breakpoint_evidence,
    LOW_MAPQ_QUALITY_GATE,
)

FAILURES = []
HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 248956422}, {"SN": "chr8", "LN": 145138636}],
})


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def _write(path, mapq, gap=(9800, 10200)):
    """
    Coverage across chr1:8000-12000 at a fixed MAPQ, with a gap so the depth
    layer sees a dip at the focus position (10000) and scores.
    """
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        n = 0
        for pos in range(8000, 12000, 25):
            if gap[0] <= pos < gap[1]:
                continue
            read = pysam.AlignedSegment(HEADER)
            read.query_name = f"r{n}"
            read.query_sequence = "A" * 100
            read.flag = 0x1 | 0x2
            read.reference_id = 0
            read.reference_start = pos
            read.mapping_quality = mapq
            read.cigar = [(0, 100)]
            read.next_reference_id = 0
            read.next_reference_start = pos + 200
            read.template_length = 300
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(read)
            n += 1
    pysam.index(path)


def run_tests():
    print("=" * 68)
    print("QUALITY GATE REGRESSION SUITE")
    print("=" * 68)

    tmp = tempfile.mkdtemp()
    low = os.path.join(tmp, "low_mapq.bam")
    clean = os.path.join(tmp, "clean.bam")
    _write(low, mapq=0)
    _write(clean, mapq=60)

    # ── The pericentromeric condition ────────────────────────────────────
    print("\nall reads MAPQ 0: losing data must not raise the score")
    res = summarize_breakpoint_evidence(low, "chr1", 10000)

    check("evidence_strength is not 'strong'",
          res["evidence_strength"] != "strong",
          f"got {res['evidence_strength']!r}, score={res['evidence_score']}")
    check("evidence_strength reports the quality limit",
          res["evidence_strength"] == "QUALITY-LIMITED",
          f"got {res['evidence_strength']!r}")
    check("no normalised score is reported",
          res["evidence_score"] is None,
          f"got {res['evidence_score']}")
    check("raw score is still available for reference",
          isinstance(res["evidence_score_raw"], (int, float)),
          f"got {res['evidence_score_raw']!r}")

    # ── The denominator must not shrink ──────────────────────────────────
    print("\nMAPQ-filtered layers are assessed, not excluded")
    # min_mapq=20 explicitly: that is what summarize_breakpoint_evidence
    # passes to all three. Note get_split_reads' own default is min_mapq=0
    # while the other two default to 20, so calling it bare would filter
    # nothing and the condition under test would not arise.
    for layer, fn in [("discordant_pairs", count_discordant_pairs),
                      ("soft_clipped_reads", count_soft_clipped_reads),
                      ("split_reads", get_split_reads)]:
        out = fn(low, "chr1", 10000, min_mapq=20)
        check(f"{layer}: assessable (stays in the denominator)",
              out["assessable"] is True, f"got {out['assessable']}")
        check(f"{layer}: flagged quality_limited",
              out["quality_limited"] is True, f"got {out['quality_limited']}")
        check(f"{layer}: counts the reads it dropped",
              out["reads_below_min_mapq"] > 0,
              f"got {out['reads_below_min_mapq']}")
        check(f"{layer}: does NOT claim 'no reads in window'",
              out["reason"] is None,
              f"got reason={out['reason']!r} with "
              f"{out['reads_below_min_mapq']} reads present")
        check(f"{layer}: excluded from unassessable_layers",
              layer not in res["unassessable_layers"],
              f"unassessable_layers={res['unassessable_layers']}")

    # ── A genuinely empty window must still be unassessable ──────────────
    print("\na truly empty window is still NOT ASSESSABLE (not a zero score)")
    empty = summarize_breakpoint_evidence(clean, "chr1", 200000)
    check("empty window -> NOT ASSESSABLE",
          empty["evidence_strength"] == "NOT ASSESSABLE",
          f"got {empty['evidence_strength']!r}")
    check("empty window -> evidence_score None",
          empty["evidence_score"] is None, f"got {empty['evidence_score']}")
    d = count_discordant_pairs(clean, "chr1", 200000)
    check("empty window keeps the 'no reads in window' reason",
          d["assessable"] is False and d["reason"] == "no reads in window",
          f"assessable={d['assessable']} reason={d['reason']!r}")
    check("empty window is not mislabelled quality_limited",
          d["quality_limited"] is False, f"got {d['quality_limited']}")

    # ── A clean locus must score normally ────────────────────────────────
    print("\na clean high-MAPQ locus is unaffected by the gate")
    ok = summarize_breakpoint_evidence(clean, "chr1", 10000)
    check("clean locus is not gated",
          ok["evidence_strength"] != "QUALITY-LIMITED",
          f"got {ok['evidence_strength']!r}")
    check("clean locus still produces a numeric score",
          isinstance(ok["evidence_score"], (int, float)),
          f"got {ok['evidence_score']!r}")

    # ── window_bp=0 must not manufacture a score ─────────────────────────
    print("\nwindow_bp=0 at a clean locus must not return 60.0")
    z = summarize_breakpoint_evidence(clean, "chr1", 10000, window_bp=0)
    score = z.get("evidence_score")
    check("window_bp=0 does not return 60.0",
          score != 60.0, f"got {score}")
    check("window_bp=0 does not return any moderate/strong score",
          not (isinstance(score, (int, float)) and score >= 40),
          f"got {score}")
    check("window_bp=0 is rejected as a caller error",
          z.get("error_type") == "invalid_parameters",
          f"got error_type={z.get('error_type')!r}, keys={sorted(z)[:6]}")

    # ── The gate threshold is exposed, not buried ────────────────────────
    print("\nthreshold provenance")
    check("LOW_MAPQ_QUALITY_GATE is importable and 0.4",
          LOW_MAPQ_QUALITY_GATE == 0.4, f"got {LOW_MAPQ_QUALITY_GATE}")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL QUALITY GATE TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
