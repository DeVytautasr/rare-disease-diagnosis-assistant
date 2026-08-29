#!/usr/bin/env python3
"""
Regression tests for the reciprocal-breakpoint verdict.

Named for the condition that exposed them. Real-data validation
(results/REAL_PATIENT_DATA_VALIDATION.md finding 3) found that
check_reciprocal_breakpoint could return two fields that contradict each
other, both of which look authoritative:

    verdict     : INSUFFICIENT EVIDENCE at both positions
    is_balanced : True

Every verdict branch required primary_disc >= 5, while is_balanced used
>= 3 on both sides. Anything in the 3-4 band on both sides fell through to
the catch-all. Reproduced on real data at two chr22 positions with 3
discordant pairs each.

A second defect in the same block: the catch-all message asserted "at both
positions" even when one side had signal — it fired on a case with 3
discordant pairs and 1 back-pointing read at the reciprocal position.

The verdict logic is now a pure function of three counts, which is what
lets this suite check every combination exhaustively instead of the handful
a BAM fixture happens to produce.

Run: python3 stage1_igv_assistant/tests/test_reciprocal_verdict.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant.tools.bam_tools import (  # noqa: E402
    _reciprocal_verdict,
    RECIPROCAL_STRONG,
    RECIPROCAL_PRESENT,
    RECIPROCAL_BALANCED,
    RECIPROCAL_BACKPOINT,
)

FAILURES = []


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def run_tests():
    print("=" * 68)
    print("RECIPROCAL VERDICT REGRESSION SUITE")
    print("=" * 68)

    # ── The invariant, over the whole input space ────────────────────────
    print("\nverdict and is_balanced can never contradict each other")
    contradictions = []
    checked = 0
    for p in range(0, 16):
        for r in range(0, 16):
            for b in range(0, 8):
                verdict, balanced = _reciprocal_verdict(p, r, b)
                checked += 1
                if balanced and verdict.startswith("INSUFFICIENT"):
                    contradictions.append((p, r, b, verdict))
    check(f"no is_balanced=True with an INSUFFICIENT verdict ({checked} combinations)",
          not contradictions,
          f"first offenders: {contradictions[:3]}")

    # ── The exact case reproduced on real data ───────────────────────────
    print("\nthe 3-4 band that produced the contradiction")
    for p, r in [(3, 3), (3, 4), (4, 3), (4, 4)]:
        verdict, balanced = _reciprocal_verdict(p, r, 0)
        check(f"primary={p} reciprocal={r}: balanced but not called INSUFFICIENT",
              balanced is True and not verdict.startswith("INSUFFICIENT"),
              f"verdict={verdict!r} is_balanced={balanced}")
        check(f"primary={p} reciprocal={r}: reported as POSSIBLE, not confirmed",
              verdict.startswith("RECIPROCAL POSSIBLE"),
              f"got {verdict!r}")

    # ── "at both positions" must mean both positions ─────────────────────
    print("\nthe catch-all message must not overclaim")
    v0, _ = _reciprocal_verdict(0, 0, 0)
    check("0/0 does say 'at both positions'",
          v0 == "INSUFFICIENT EVIDENCE at both positions", f"got {v0!r}")

    for p, r in [(5, 1), (2, 1), (1, 2), (4, 1), (1, 0), (0, 1)]:
        verdict, _ = _reciprocal_verdict(p, r, 0)
        if verdict.startswith("INSUFFICIENT"):
            check(f"primary={p} reciprocal={r}: does not claim 'at both positions'",
                  "at both positions" not in verdict, f"got {verdict!r}")

    check("one-sided signal names which side is short",
          "at the partner position" in _reciprocal_verdict(5, 0, 0)[0]
          or "NOT FOUND" in _reciprocal_verdict(5, 0, 0)[0],
          f"got {_reciprocal_verdict(5, 0, 0)[0]!r}")
    check("zero primary with partner signal names the primary",
          "at the primary position" in _reciprocal_verdict(0, 4, 0)[0],
          f"got {_reciprocal_verdict(0, 4, 0)[0]!r}")

    # ── The confident tiers still behave ─────────────────────────────────
    print("\nthe confident tiers are unchanged")
    v, bal = _reciprocal_verdict(5, 5, 3)
    check("strong both sides + back-pointing -> CONFIRMED",
          v.startswith("RECIPROCAL CONFIRMED") and bal is True, f"got {v!r}")
    v, _ = _reciprocal_verdict(5, 2, 0)
    check("strong primary + weak partner -> LIKELY",
          v.startswith("RECIPROCAL LIKELY"), f"got {v!r}")
    v, bal = _reciprocal_verdict(5, 0, 0)
    check("strong primary + silent partner -> NOT FOUND",
          v.startswith("RECIPROCAL NOT FOUND") and bal is False, f"got {v!r}")
    v, bal = _reciprocal_verdict(0, 0, 0)
    check("nothing anywhere -> INSUFFICIENT, not balanced",
          v.startswith("INSUFFICIENT") and bal is False, f"got {v!r}")

    # ── back_pointing must still gate CONFIRMED ──────────────────────────
    print("\nback_pointing still gates the top tier")
    v, _ = _reciprocal_verdict(9, 9, RECIPROCAL_BACKPOINT - 1)
    check("strong both sides but too few back-pointing -> not CONFIRMED",
          not v.startswith("RECIPROCAL CONFIRMED"), f"got {v!r}")

    # ── Thresholds are named, not buried ─────────────────────────────────
    print("\nthreshold provenance")
    check("thresholds exposed as named constants",
          (RECIPROCAL_STRONG, RECIPROCAL_PRESENT,
           RECIPROCAL_BALANCED, RECIPROCAL_BACKPOINT) == (5, 2, 3, 3),
          f"got {(RECIPROCAL_STRONG, RECIPROCAL_PRESENT, RECIPROCAL_BALANCED, RECIPROCAL_BACKPOINT)}")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL RECIPROCAL VERDICT TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
