#!/usr/bin/env python3
"""
Regression test: a parameter that appears to filter must filter.

From results/REAL_PATIENT_DATA_VALIDATION.md finding 7.

get_split_reads read only fields 0 and 1 of each SA record. Fields 2 (strand),
3 (CIGAR), 4 (mapQ) and 5 (NM) were parsed past and discarded. Two
consequences:

  * min_mapq filtered the PRIMARY alignment only. A caller asking for
    min_mapq=20 still had mapQ-0 partner segments counted at full weight.
    In one real window 144/619 and 224/316 SA records were mapQ 0 — 71% for
    one sample. Evidence that looked filtered was not.
  * Strand was dropped entirely, though 12.3% and 10.1% of real SA entries
    flip strand relative to the primary. Orientation is what distinguishes an
    inversion-type junction from a direct one.

Also fixed here: get_split_reads defaulted to min_mapq=0 while
count_discordant_pairs and count_soft_clipped_reads both defaulted to 20. A
bare call filtered nothing and silently disagreed with its neighbours;
summarize always passed 20 explicitly, which is why the asymmetry survived.

Run: python3 stage1_igv_assistant/tests/test_sa_mapq_strand.py
"""
import inspect
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant.tools.bam_tools import (  # noqa: E402
    get_split_reads,
    count_discordant_pairs,
    count_soft_clipped_reads,
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


def _bam(path, sa_tag, n=20, reverse=False):
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        for i in range(n):
            r = pysam.AlignedSegment(HEADER)
            r.query_name = f"r{i}"
            r.query_sequence = "A" * 100
            r.flag = 0x1 | (0x10 if reverse else 0)
            r.reference_id = 0
            r.reference_start = 10000 + i * 5
            r.mapping_quality = 60           # the PRIMARY is always high quality
            r.cigar = [(4, 20), (0, 80)]
            r.next_reference_id = 0
            r.next_reference_start = 10200
            r.query_qualities = pysam.qualitystring_to_array("I" * 100)
            r.set_tag("SA", sa_tag)
            bam.write(r)
    pysam.index(path)
    return path


def run_tests():
    print("=" * 68)
    print("SA MAPQ / STRAND REGRESSION SUITE")
    print("=" * 68)
    tmp = tempfile.mkdtemp()

    print("\nmin_mapq now reaches the SA record, not just the primary")
    for sa_mapq, expect_split in [(0, 0), (10, 0), (19, 0), (20, 20), (60, 20)]:
        path = _bam(os.path.join(tmp, f"q{sa_mapq}.bam"),
                    f"chr8,500000,+,80M20S,{sa_mapq},0;")
        r = get_split_reads(path, "chr1", 10050, min_mapq=20)
        check(f"SA mapQ {sa_mapq} with min_mapq=20 -> split_reads {expect_split}",
              r["split_reads"] == expect_split, f"got {r['split_reads']}")

    print("\ndropped SA records are counted, not silently discarded")
    low = get_split_reads(_bam(os.path.join(tmp, "low.bam"),
                               "chr8,500000,+,80M20S,0,0;"), "chr1", 10050, min_mapq=20)
    check("sa_entries_below_min_mapq reports the drops",
          low["sa_entries_below_min_mapq"] == 20,
          f"got {low['sa_entries_below_min_mapq']}")
    check("sa_entries_total still counts every raw record",
          low["sa_entries_total"] == 20, f"got {low['sa_entries_total']}")
    check("a read with no surviving partner is not split evidence",
          low["split_reads"] == 0 and low["partner_chromosomes"] == {},
          f"split={low['split_reads']} partners={low['partner_chromosomes']}")

    print("\nlowering min_mapq lets them back in (the parameter is real both ways)")
    permissive = get_split_reads(_bam(os.path.join(tmp, "perm.bam"),
                                      "chr8,500000,+,80M20S,0,0;"),
                                 "chr1", 10050, min_mapq=0)
    check("min_mapq=0 keeps mapQ-0 partners",
          permissive["split_reads"] == 20, f"got {permissive['split_reads']}")

    print("\nstrand is recorded rather than discarded")
    fwd = get_split_reads(_bam(os.path.join(tmp, "f.bam"),
                               "chr8,500000,+,80M20S,60,0;"), "chr1", 10050)
    rev = get_split_reads(_bam(os.path.join(tmp, "r.bam"),
                               "chr8,500000,-,80M20S,60,0;"), "chr1", 10050)
    check("same-strand SA segments count as concordant",
          fwd["partner_strand_concordant"] == 20 and fwd["partner_strand_flipped"] == 0,
          f"got {fwd['partner_strand_concordant']}/{fwd['partner_strand_flipped']}")
    check("opposite-strand SA segments count as flipped",
          rev["partner_strand_concordant"] == 0 and rev["partner_strand_flipped"] == 20,
          f"got {rev['partner_strand_concordant']}/{rev['partner_strand_flipped']}")

    print("\nstrand is relative to the PRIMARY, not to '+' absolutely")
    rev_primary = get_split_reads(
        _bam(os.path.join(tmp, "rp.bam"), "chr8,500000,-,80M20S,60,0;", reverse=True),
        "chr1", 10050)
    check("reverse primary with reverse SA is concordant, not flipped",
          rev_primary["partner_strand_concordant"] == 20,
          f"got {rev_primary['partner_strand_concordant']}/"
          f"{rev_primary['partner_strand_flipped']}")

    print("\nmalformed SA records do not crash the parse")
    for tag in ["chr8,500000;", "chr8,500000,+;", "chr8,500000,+,80M20S,notanint,0;"]:
        path = _bam(os.path.join(tmp, f"m{abs(hash(tag))%9999}.bam"), tag, n=5)
        try:
            r = get_split_reads(path, "chr1", 10050, min_mapq=20)
            check(f"malformed {tag!r} handled", isinstance(r, dict),
                  f"got {type(r).__name__}")
        except Exception as exc:
            check(f"malformed {tag!r} handled", False,
                  f"raised {type(exc).__name__}: {exc}")

    print("\nthe layer defaults now agree with each other")
    defaults = {
        "get_split_reads": inspect.signature(get_split_reads).parameters["min_mapq"].default,
        "count_discordant_pairs": inspect.signature(count_discordant_pairs).parameters["min_mapq"].default,
        "count_soft_clipped_reads": inspect.signature(count_soft_clipped_reads).parameters["min_mapq"].default,
    }
    check("all three layers default to the same min_mapq",
          len(set(defaults.values())) == 1, f"got {defaults}")
    check("that default is 20", defaults["get_split_reads"] == 20, f"got {defaults}")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL SA MAPQ / STRAND TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
