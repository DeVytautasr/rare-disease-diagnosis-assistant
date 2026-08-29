#!/usr/bin/env python3
"""
Regression test: partner_chromosomes counts reads, not SA entries.

From results/REAL_PATIENT_DATA_VALIDATION.md finding 9.

split_reads incremented once per READ; partner_chromosomes incremented once
per SA ENTRY. The two headline numbers sat next to each other in different
units, so sum(partner_chromosomes.values()) != split_reads (619 against 603
in one real window), and a read with two segments on the SAME contig counted
as two independent pieces of evidence for that contig. Multi-segment reads
are common: 3,807 and 3,282 reads with more than one SA entry across the two
samples, some with four. The docstring said "{chr_name: count}" without
saying count of what.

Partner counts are now deduplicated per read. One nuance, asserted below:
the sum equals split_reads only when no read names more than one DISTINCT
partner contig. A read genuinely pointing at two contigs contributes once to
each, so the sum exceeds split_reads — that is information, not
double-counting, and the docstring now says so. sa_entries_total preserves
the raw entry count.

Run: python3 stage1_igv_assistant/tests/test_split_partner_units.py
"""
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant.tools.bam_tools import get_split_reads  # noqa: E402

FAILURES = []

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 1000000},
           {"SN": "chr6", "LN": 900000},
           {"SN": "chr2", "LN": 800000}],
})


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def _bam(path, sa_tag, n=10):
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        for i in range(n):
            r = pysam.AlignedSegment(HEADER)
            r.query_name = f"r{i}"
            r.query_sequence = "A" * 100
            r.flag = 0x1
            r.reference_id = 0
            r.reference_start = 10000 + i * 10
            r.mapping_quality = 60
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
    print("SPLIT PARTNER UNIT REGRESSION SUITE")
    print("=" * 68)
    tmp = tempfile.mkdtemp()

    # ── The defect: repeated segments on ONE contig ──────────────────────
    print("\na read with several segments on one contig counts once for it")
    repeated = _bam(os.path.join(tmp, "rep.bam"),
                    "chr6,100,+,80M20S,60,0;chr6,200,+,80M20S,60,0;"
                    "chr6,300,+,80M20S,60,0;")
    r = get_split_reads(repeated, "chr1", 10050, min_mapq=20)
    check("split_reads counts reads", r["split_reads"] == 10,
          f"got {r['split_reads']}")
    check("chr6 counts 10 reads, not 30 entries",
          r["partner_chromosomes"].get("chr6") == 10,
          f"got {r['partner_chromosomes']}")
    check("sum(partner_chromosomes) == split_reads for single-partner reads",
          sum(r["partner_chromosomes"].values()) == r["split_reads"],
          f"sum={sum(r['partner_chromosomes'].values())} "
          f"split_reads={r['split_reads']}")
    check("sa_entries_total preserves the raw entry count",
          r["sa_entries_total"] == 30, f"got {r['sa_entries_total']}")

    # ── The nuance: genuinely multiple distinct partners ─────────────────
    print("\na read naming two DISTINCT contigs contributes once to each")
    multi = _bam(os.path.join(tmp, "multi.bam"),
                 "chr6,100,+,80M20S,60,0;chr6,200,+,80M20S,60,0;"
                 "chr2,300,+,80M20S,60,0;")
    m = get_split_reads(multi, "chr1", 10050, min_mapq=20)
    check("chr6 counts 10 (the repeat is deduped)",
          m["partner_chromosomes"].get("chr6") == 10,
          f"got {m['partner_chromosomes']}")
    check("chr2 counts 10",
          m["partner_chromosomes"].get("chr2") == 10,
          f"got {m['partner_chromosomes']}")
    check("sum exceeds split_reads, by design, for multi-partner reads",
          sum(m["partner_chromosomes"].values()) == 20 and m["split_reads"] == 10,
          f"sum={sum(m['partner_chromosomes'].values())} "
          f"split_reads={m['split_reads']}")
    check("sa_entries_total still counts every raw entry",
          m["sa_entries_total"] == 30, f"got {m['sa_entries_total']}")

    # ── The simple case must be unchanged ────────────────────────────────
    print("\nthe single-entry case is unaffected")
    single = _bam(os.path.join(tmp, "single.bam"), "chr6,100,+,80M20S,60,0;")
    s = get_split_reads(single, "chr1", 10050, min_mapq=20)
    check("one entry per read -> partner count equals split_reads",
          s["partner_chromosomes"].get("chr6") == s["split_reads"] == 10,
          f"got {s['partner_chromosomes']} split={s['split_reads']}")
    check("sa_entries_total equals split_reads here",
          s["sa_entries_total"] == 10, f"got {s['sa_entries_total']}")

    print("\nno reads means no counts, not a crash")
    empty = get_split_reads(single, "chr1", 500000, min_mapq=20)
    check("empty window returns zero entries",
          empty["sa_entries_total"] == 0, f"got {empty['sa_entries_total']}")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL SPLIT PARTNER UNIT TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
