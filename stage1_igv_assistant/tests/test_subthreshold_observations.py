#!/usr/bin/env python3
"""
Regression tests for the two defects found by the first LLM session run
against real patient data on the post-fix tools
(results/LLM_SESSION_5_PATIENT_DATA_qwen.md).

DEFECT A — self-contradicting observations.
    supporting_observations ended with a blanket denial:

        "No discordant pairs, soft-clipping, split reads, or depth changes
         detected near this position."

    fired on evidence_score == 0 alone, so it printed directly beneath
    "1 discordant pair (mate on chr12)" in the same list. Two contradicting
    sentences, no way for a reader to tell which is authoritative.

    MIN_ABSOLUTE_SUPPORT amplified this — a locus with 1-2 supporting reads
    scores 0 while its counts stay non-zero — but did not cause it. On the
    42-locus control grid the contradiction fires on 7/42 and 8/42 loci in
    the pre-MIN_ABSOLUTE_SUPPORT code as well, via scattered soft clips and
    off-position depth dips, which were always counted but unscored. Hence
    the tests below cover all three read layers plus the empty case, not just
    the discordant one that surfaced it.

    The denial must now fire only on the actual absence of all four layers.

DEFECT B — window-size normalisation.
    MIN_ABSOLUTE_SUPPORT was a bare count, but the number of background
    reads in a window scales with the window's width. The same locus in the
    same BAM read "none" with 1 discordant pair at the default 500 bp and
    "weak" with 4 pairs at window_bp=1000 — the verdict moved because an
    argument moved, not because the evidence did.

    The threshold is now stated per MIN_ABSOLUTE_SUPPORT_WINDOW_BP and
    scaled to the window actually used, and both the window and the derived
    threshold are returned alongside the verdict.

Run: python3 stage1_igv_assistant/tests/test_subthreshold_observations.py
"""
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant.tools.bam_tools import (  # noqa: E402
    summarize_breakpoint_evidence,
    min_supporting_reads,
    MIN_ABSOLUTE_SUPPORT,
    MIN_ABSOLUTE_SUPPORT_WINDOW_BP,
)

FAILURES = []

BLANKET = "No discordant pairs, soft-clipping, split reads, or depth changes"

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


def _bam(path, disc_offsets=(), split_offsets=(), clip_offsets=()):
    """
    Uniform concordant background across the whole +/-2kb depth window (so the
    depth layer has nothing to say), plus signal reads at chosen offsets from
    the focus position 10,000.

    disc_offsets  reads whose mate is on chr8
    split_offsets reads carrying an SA tag pointing at chr8
    clip_offsets  reads with a 20bp soft clip, no SA tag
    """
    recs = [(p, f"n{p}", 0, p + 300, None, False)
            for p in range(7500, 12500, 4)]
    recs += [(10000 + o, f"d{i}", 1, 500000, None, False)
             for i, o in enumerate(disc_offsets)]
    recs += [(10000 + o, f"s{i}", 0, 10300, "chr8,500000,+,80M20S,60,0;", False)
             for i, o in enumerate(split_offsets)]
    recs += [(10000 + o, f"c{i}", 0, 10300, None, True)
             for i, o in enumerate(clip_offsets)]
    recs.sort(key=lambda x: x[0])
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        for pos, name, mate_ref, mate_pos, sa, clipped in recs:
            r = pysam.AlignedSegment(HEADER)
            r.query_name = name
            r.query_sequence = "A" * 100
            r.flag = 0x1 if (mate_ref == 1 or sa) else (0x1 | 0x2)
            r.reference_id = 0
            r.reference_start = pos
            r.mapping_quality = 60
            r.cigar = [(4, 20), (0, 80)] if (sa or clipped) else [(0, 100)]
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
    print("SUB-THRESHOLD OBSERVATION / WINDOW NORMALISATION SUITE")
    print("=" * 68)
    tmp = tempfile.mkdtemp()

    # ── DEFECT A ───────────────────────────────────────────────────────────
    print("\nA. a specific count and a blanket denial must never co-occur")
    for n in (1, 2):
        path = _bam(os.path.join(tmp, f"a_disc{n}.bam"),
                    disc_offsets=tuple(range(0, n * 50, 50)))
        res = summarize_breakpoint_evidence(path, "chr1", 10000)
        obs = res["supporting_observations"]
        joined = " ".join(obs)
        says_count = any("discordant pair" in o for o in obs)
        says_none = BLANKET in joined
        check(f"{n} discordant read(s): score is still 0 (strength 'none')",
              res["discordant_pair_score"] == 0.0
              and res["evidence_strength"] == "none",
              f"score={res['discordant_pair_score']} "
              f"strength={res['evidence_strength']!r}")
        check(f"{n} discordant read(s): the count is reported",
              says_count, f"observations={obs!r}")
        check(f"{n} discordant read(s): NOT also denied wholesale",
              not says_none, f"observations={obs!r}")
        check(f"{n} discordant read(s): named as sub-threshold, with the bar",
              "sub-threshold" in joined
              and f"{res['min_supporting_reads']}-read minimum" in joined,
              f"observations={obs!r}")
        # the invariant, stated once as the session stated it
        check(f"{n} discordant read(s): invariant — no count beside a denial",
              not (says_count and says_none),
              f"observations={obs!r}")

    print("\nA. split reads and soft clips get the same treatment")
    path = _bam(os.path.join(tmp, "a_split.bam"), split_offsets=(0, 40))
    res = summarize_breakpoint_evidence(path, "chr1", 10000)
    joined = " ".join(res["supporting_observations"])
    check("2 split reads: reported, not denied",
          "split read" in joined and BLANKET not in joined,
          f"observations={res['supporting_observations']!r}")
    check("2 split reads: named as sub-threshold",
          "sub-threshold" in joined and "split read(s), below the" in joined,
          f"observations={res['supporting_observations']!r}")

    path = _bam(os.path.join(tmp, "a_clip.bam"), clip_offsets=(0, 40))
    res = summarize_breakpoint_evidence(path, "chr1", 10000)
    joined = " ".join(res["supporting_observations"])
    check("2 scattered soft clips: reported, not denied",
          "soft-clipped" in joined and BLANKET not in joined,
          f"observations={res['supporting_observations']!r}")

    print("\nA. the denial survives where it is TRUE — a genuinely empty window")
    path = _bam(os.path.join(tmp, "a_empty.bam"))
    res = summarize_breakpoint_evidence(path, "chr1", 10000)
    obs = res["supporting_observations"]
    check("no support at all: strength is 'none'",
          res["evidence_strength"] == "none", f"got {res['evidence_strength']!r}")
    check("no support at all: the blanket denial IS emitted",
          BLANKET in " ".join(obs), f"observations={obs!r}")
    check("no support at all: nothing claimed as sub-threshold",
          "sub-threshold" not in " ".join(obs), f"observations={obs!r}")

    # ── DEFECT B ───────────────────────────────────────────────────────────
    print("\nB. the threshold scales with the window")
    check(f"reference width is exposed "
          f"(MIN_ABSOLUTE_SUPPORT_WINDOW_BP = {MIN_ABSOLUTE_SUPPORT_WINDOW_BP})",
          MIN_ABSOLUTE_SUPPORT_WINDOW_BP == 500,
          f"got {MIN_ABSOLUTE_SUPPORT_WINDOW_BP!r}")
    for window, expect in [(500, 3), (1000, 6), (1500, 9), (250, 2),
                           (200, 2), (100, 2), (50, 2)]:
        check(f"min_supporting_reads({window}) == {expect}",
              min_supporting_reads(window) == expect,
              f"got {min_supporting_reads(window)}")
    check("at the reference width it is exactly MIN_ABSOLUTE_SUPPORT",
          min_supporting_reads(MIN_ABSOLUTE_SUPPORT_WINDOW_BP) == MIN_ABSOLUTE_SUPPORT,
          f"got {min_supporting_reads(MIN_ABSOLUTE_SUPPORT_WINDOW_BP)}")
    check("it never falls below 2 — one read plus its duplicate",
          all(min_supporting_reads(w) >= 2 for w in (1, 10, 50, 100, 166)),
          "a window narrow enough would have let a single read score")

    print("\nB. widening the window must not flip the verdict on its own")
    # 2 discordant reads inside +/-500, 2 more between 500 and 1000. Under a
    # bare count of 3 this locus reads 'none' at 500 and 'weak' at 1000 —
    # exactly the P2 observation from the session.
    path = _bam(os.path.join(tmp, "b_scale.bam"),
                disc_offsets=(-250, 250, -750, 750))
    at500 = summarize_breakpoint_evidence(path, "chr1", 10000, window_bp=500)
    at1000 = summarize_breakpoint_evidence(path, "chr1", 10000, window_bp=1000)
    check("the fixture behaves as designed: 2 pairs at 500bp, 4 at 1000bp",
          at500["discordant_pairs"]["discordant_pairs"] == 2
          and at1000["discordant_pairs"]["discordant_pairs"] == 4,
          f"got {at500['discordant_pairs']['discordant_pairs']} and "
          f"{at1000['discordant_pairs']['discordant_pairs']}")
    check("500bp -> 'none' (2 below the 3-read bar)",
          at500["evidence_strength"] == "none",
          f"got {at500['evidence_strength']!r}")
    check("1000bp -> still 'none' (4 below the scaled 6-read bar)",
          at1000["evidence_strength"] == "none",
          f"got {at1000['evidence_strength']!r} "
          f"min={at1000['min_supporting_reads']}")
    check("the same locus does not change verdict with the window alone",
          at500["evidence_strength"] == at1000["evidence_strength"],
          f"{at500['evidence_strength']!r} vs {at1000['evidence_strength']!r}")

    print("\nB. a real signal still scores at a wide window")
    path = _bam(os.path.join(tmp, "b_real.bam"),
                disc_offsets=tuple(range(-800, 800, 100)))  # 16 pairs
    wide = summarize_breakpoint_evidence(path, "chr1", 10000, window_bp=1000)
    check("16 pairs at 1000bp clears the scaled bar of 6",
          wide["discordant_pair_score"] == 7.5
          and wide["evidence_strength"] == "weak",
          f"score={wide['discordant_pair_score']} "
          f"strength={wide['evidence_strength']!r} "
          f"n={wide['discordant_pairs']['discordant_pairs']}")

    print("\nB. the window is reported alongside the verdict")
    check("window_bp is returned",
          at1000["window_bp"] == 1000, f"got {at1000.get('window_bp')!r}")
    check("min_supporting_reads is returned, and matches the helper",
          at1000["min_supporting_reads"] == min_supporting_reads(1000) == 6,
          f"got {at1000.get('min_supporting_reads')!r}")
    check("the interpretation template states the window and the bar",
          "+/-1000bp" in at1000["interpretation_template"]
          and "6 supporting read(s)" in at1000["interpretation_template"],
          f"got {at1000['interpretation_template'][:400]!r}")
    check("and says the verdict is not reproducible without it",
          "not reproducible" in at1000["interpretation_template"],
          f"got {at1000['interpretation_template'][:400]!r}")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL SUB-THRESHOLD / WINDOW NORMALISATION TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
