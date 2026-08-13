"""
test_partner_distribution.py

Regression tests for _describe_partner_distribution and the observation
strings built from it.

These exist because summarize_breakpoint_evidence used to assert a
"predominant" partner chromosome that did not exist: it took
next(iter(mate_chromosomes)) -- the first-inserted dict key, not even the
maximum -- and unconditionally wrote "mates mapping predominantly to
<chrom>". A single discordant read therefore produced "mates mapping
predominantly to chr12" in the ADVERSARIAL benchmark case, whose prompt
falsely asserts a t(1;12) translocation. See
results/BENCHMARK_LOCAL_MODELS.md's correction section.

Run directly: python3 stage1_igv_assistant/tests/test_partner_distribution.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from stage1_igv_assistant.tools.bam_tools import (  # noqa: E402
    PARTNER_DOMINANCE_MIN_READS,
    PARTNER_DOMINANCE_MIN_SHARE,
    _describe_partner_distribution,
)

FAILURES = []


def check(label: str, condition: bool, detail: str = "") -> None:
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}{(' — ' + detail) if detail else ''}")
        FAILURES.append(label)


def test_single_read_never_predominant():
    """One read must never produce 'predominantly' — the original bug."""
    print("\nTEST 1: single discordant read")
    out = _describe_partner_distribution({"chr12": 1}, "mate")
    print(f"  -> {out!r}")
    check("no 'predominantly' for n=1", "predominantly" not in out, out)
    check("names the actual chromosome", "chr12" in out, out)
    check("says it is a single read", "single read" in out, out)


def test_seven_singletons_scattered():
    """Seven mates on seven chromosomes must read as scattered."""
    print("\nTEST 2: seven singletons across seven chromosomes")
    counts = {"chr4": 1, "chr14": 1, "chr12": 1, "chr8": 1,
              "chr11": 1, "chr19": 1, "chr5": 1}
    out = _describe_partner_distribution(counts, "mate")
    print(f"  -> {out!r}")
    check("no 'predominantly' for scattered mates", "predominantly" not in out, out)
    check("says 'scattered'", "scattered" in out, out)
    check("reports the chromosome count", "7 chromosomes" in out, out)
    check("says there is no dominant partner", "no dominant partner" in out, out)


def test_genuine_dominance_still_reported():
    """A real cluster must still be called predominant — the fix must not
    suppress true signal (e.g. the synthetic translocation fixture)."""
    print("\nTEST 3: genuine dominant partner (15/15 on chr8)")
    out = _describe_partner_distribution({"chr8": 15}, "mate")
    print(f"  -> {out!r}")
    check("says 'predominantly'", "predominantly" in out, out)
    check("names chr8", "chr8" in out, out)
    check("shows the supporting ratio", "15/15" in out, out)


def test_majority_threshold_boundary():
    """Dominance requires >= PARTNER_DOMINANCE_MIN_SHARE of mates."""
    print("\nTEST 4: majority threshold boundary")
    # 6/10 = 0.6 -> exactly at threshold, should qualify
    at = _describe_partner_distribution({"chr8": 6, "chr1": 4}, "mate")
    print(f"  6/10 -> {at!r}")
    check("60% share qualifies as predominant", "predominantly" in at, at)
    # 5/10 = 0.5 -> below threshold, must not qualify
    below = _describe_partner_distribution({"chr8": 5, "chr1": 5}, "mate")
    print(f"  5/10 -> {below!r}")
    check("50% share is not predominant", "predominantly" not in below, below)


def test_min_reads_threshold():
    """Dominance also requires >= PARTNER_DOMINANCE_MIN_READS total."""
    print("\nTEST 5: minimum-reads threshold")
    # 2/2 on one chromosome: 100% share but below the read-count floor
    out = _describe_partner_distribution({"chr8": 2}, "mate")
    print(f"  2/2 -> {out!r}")
    check(
        f"n=2 (< {PARTNER_DOMINANCE_MIN_READS}) is not predominant despite 100% share",
        "predominantly" not in out,
        out,
    )
    check("n=2 on one chromosome is not called 'scattered'",
          "scattered" not in out, out)
    check("no '1 chromosomes' grammatical artifact",
          "1 chromosomes" not in out, out)
    # 3/3 clears the floor
    out3 = _describe_partner_distribution({"chr8": 3}, "mate")
    print(f"  3/3 -> {out3!r}")
    check("n=3 with 100% share is predominant", "predominantly" in out3, out3)


def test_picks_max_not_insertion_order():
    """The specific original defect: first-inserted key was chosen, not the
    maximum. chr1 is inserted first but chr8 has the higher count."""
    print("\nTEST 6: chooses max by count, not dict insertion order")
    counts = {"chr1": 1, "chr8": 9}
    out = _describe_partner_distribution(counts, "mate")
    print(f"  -> {out!r}")
    check("names chr8 (the maximum)", "chr8" in out, out)
    check("does not name chr1 (merely first-inserted)", "chr1" not in out, out)


def test_deterministic_tie_break():
    """Ties must break deterministically, not by insertion order."""
    print("\nTEST 7: deterministic tie-breaking")
    a = _describe_partner_distribution({"chr1": 4, "chr8": 4, "chr2": 4}, "mate")
    b = _describe_partner_distribution({"chr8": 4, "chr2": 4, "chr1": 4}, "mate")
    print(f"  -> {a!r}")
    check("same counts in different insertion order give the same text", a == b,
          f"{a!r} != {b!r}")


def test_noun_is_parameterised():
    """Split reads use a different noun for the same shape of claim."""
    print("\nTEST 8: split-read noun")
    out = _describe_partner_distribution({"chr8": 8}, "supplementary alignment")
    print(f"  -> {out!r}")
    check("uses the supplied noun", "supplementary alignment" in out, out)
    check("pluralises it", "supplementary alignments" in out, out)


def test_empty_input():
    print("\nTEST 9: empty partner map")
    check("empty dict yields empty string",
          _describe_partner_distribution({}, "mate") == "")


def run_tests():
    print("=" * 66)
    print("_describe_partner_distribution regression tests")
    print(f"(thresholds: >= {PARTNER_DOMINANCE_MIN_READS} reads AND "
          f">= {PARTNER_DOMINANCE_MIN_SHARE:.0%} share)")
    print("=" * 66)

    test_single_read_never_predominant()
    test_seven_singletons_scattered()
    test_genuine_dominance_still_reported()
    test_majority_threshold_boundary()
    test_min_reads_threshold()
    test_picks_max_not_insertion_order()
    test_deterministic_tie_break()
    test_noun_is_parameterised()
    test_empty_input()

    print("\n" + "=" * 66)
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED:")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("All checks passed.")
    print("=" * 66)


if __name__ == "__main__":
    run_tests()
