#!/usr/bin/env python3
"""
Regression test: no per-read coordinate may be embedded in prose output.

From results/REAL_PATIENT_DATA_VALIDATION.md finding 12.

summarize_breakpoint_evidence's split-read observation interpolated a
concrete partner locus taken from a single read's SA tag:

    "3 split reads (2% of reads in window) with supplementary alignments
     mapping predominantly to chrN (3/3) (e.g. chrN:12345678)."

On patient data that trailing coordinate is a per-read position. It arrives
pre-embedded in the field a reporting assistant is most likely to copy
wholesale, and unlike a structured field it cannot be stripped without
rewriting the sentence.

The fix is a separation of concerns, not a removal of data:
example_partner_loci is still returned as a structured field for anyone who
wants it, and the caller decides whether it reaches a report.

Run: python3 stage1_igv_assistant/tests/test_prose_coordinates.py
"""
import os
import re
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant.tools.bam_tools import (  # noqa: E402
    get_split_reads,
    summarize_breakpoint_evidence,
)

FAILURES = []

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 1000000}, {"SN": "chr8", "LN": 900000}],
})
PARTNER_POS = 47000000 % 900000   # a concrete coordinate we can search for


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def _sa_bam(path):
    """Clipped reads on chr1 carrying SA tags that point at chr8."""
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        for i in range(30):
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
            r.set_tag("SA", f"chr8,{PARTNER_POS},+,80M20S,60,0;")
            bam.write(r)
    pysam.index(path)


def run_tests():
    print("=" * 68)
    print("PROSE COORDINATE REGRESSION SUITE")
    print("=" * 68)

    tmp = tempfile.mkdtemp()
    bam = os.path.join(tmp, "sa.bam")
    _sa_bam(bam)

    res = summarize_breakpoint_evidence(bam, "chr1", 10150)
    obs = " ".join(res["supporting_observations"])
    template = res.get("interpretation_template", "")
    prose = obs + " " + template

    print("\nthe split layer actually fired (else this proves nothing)")
    check("a split-read observation was emitted",
          "split read" in prose, f"observations={obs[:160]!r}")

    print("\nno per-read coordinate appears in any prose field")
    check("the partner position is absent from supporting_observations",
          str(PARTNER_POS) not in obs, f"got {obs[:220]!r}")
    check("the partner position is absent from interpretation_template",
          str(PARTNER_POS) not in template, f"got {template[:220]!r}")
    check("no '(e.g. ...)' example suffix remains",
          "e.g." not in prose, f"got {prose[:220]!r}")
    # Scoped to PARTNER contigs on purpose. Two coordinates legitimately
    # appear in prose and must not be flagged: the queried position, which is
    # the caller's own input echoed back, and the off-position dip position,
    # which is a bin-level aggregate. Neither is derived from an individual
    # read. What finding 12 concerns is a coordinate lifted from one read's
    # SA tag, and that always names a partner contig.
    partner_hits = re.findall(r"chr8:\d+", prose)
    check("no coordinate on a partner contig appears in prose",
          not partner_hits, f"matched {partner_hits[:3]}")
    check("the queried position is still stated (it is the caller's own input)",
          "chr1:10150" in prose, f"got {prose[:160]!r}")

    print("\nthe structured field is retained, not removed")
    split = get_split_reads(bam, "chr1", 10150, min_mapq=20)
    loci = split.get("example_partner_loci")
    check("example_partner_loci is still returned",
          isinstance(loci, list) and len(loci) > 0, f"got {loci!r}")
    check("it still carries the real partner locus",
          any(str(PARTNER_POS) in x for x in loci), f"got {loci!r}")
    check("a caller can drop it without touching the prose",
          "example_partner_loci" in split, f"keys={sorted(split)[:8]}")

    print("\nthe partner distribution is still described qualitatively")
    check("prose still names the partner contig",
          "chr8" in prose, f"got {prose[:200]!r}")

    print("\nthe retired interpolation is gone from the source")
    src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "tools", "bam_tools.py"), encoding="utf-8").read()
    live = [ln for ln in src.splitlines()
            if "example_suffix" in ln and not ln.strip().startswith("#")]
    check("no live code builds an example_suffix", not live, f"{live[:2]}")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL PROSE COORDINATE TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
