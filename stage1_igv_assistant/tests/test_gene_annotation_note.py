#!/usr/bin/env python3
"""
Regression test for the gene-annotation note that asserted a breakpoint.

From results/REAL_PATIENT_DATA_VALIDATION.md finding 13. It is prose, in a
field a reporting assistant is likely to quote verbatim, so a wrong word
travels further than a wrong number.

get_gene_at_locus returned, for ANY queried coordinate:
    "clinical_note": "Breakpoint directly disrupts N gene(s)."
It emitted that at arbitrary control positions where there was no breakpoint
and no evidence of one. A coordinate lookup establishes overlap; it cannot
establish that a breakpoint exists or that a gene is disrupted.

Ensembl is stubbed here: this suite tests wording, and must not fail because
a public API was slow.

Run: python3 stage1_igv_assistant/tests/test_gene_annotation_note.py
"""
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant.tools import bam_tools  # noqa: E402
from stage1_igv_assistant.tools.bam_tools import (  # noqa: E402
    get_gene_at_locus,
    summarize_breakpoint_evidence,
)

FAILURES = []


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


class _FakeResponse:
    def __init__(self, payload):
        self.status_code = 200
        self._payload = payload

    def json(self):
        return self._payload


class _StubEnsembl:
    """Replaces _requests.get so the wording test never touches the network."""

    def __init__(self, payload):
        self.payload = payload
        self._real = None

    def __enter__(self):
        self._real = bam_tools._requests.get
        bam_tools._requests.get = lambda *a, **k: _FakeResponse(self.payload)
        return self

    def __exit__(self, *exc):
        bam_tools._requests.get = self._real


GENE_HIT = [{
    "gene_id": "ENSG00000173218", "external_name": "VANGL1",
    "biotype": "protein_coding", "strand": 1,
    "start": 115641854, "end": 115698224,
}]

HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 248956422}],
})


def _depth_dip_bam(path):
    """High-MAPQ coverage with a gap at 10000 so the depth layer scores."""
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        n = 0
        for pos in range(8000, 12000, 25):
            if 9800 <= pos < 10200:
                continue
            r = pysam.AlignedSegment(HEADER)
            r.query_name = f"r{n}"
            r.query_sequence = "A" * 100
            r.flag = 0x1 | 0x2
            r.reference_id = 0
            r.reference_start = pos
            r.mapping_quality = 60
            r.cigar = [(0, 100)]
            r.next_reference_id = 0
            r.next_reference_start = pos + 200
            r.template_length = 300
            r.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(r)
            n += 1
    pysam.index(path)


def run_tests():
    print("=" * 68)
    print("GENE ANNOTATION NOTE REGRESSION SUITE")
    print("=" * 68)

    print("\ngene lookup states overlap, not disruption")
    with _StubEnsembl(GENE_HIT):
        hit = get_gene_at_locus("chr1", 115686862)
    with _StubEnsembl([]):
        miss = get_gene_at_locus("chr1", 120000000)

    for label, res in [("gene present", hit), ("intergenic", miss)]:
        check(f"{label}: clinical_note is gone",
              "clinical_note" not in res, f"keys={sorted(res)}")
        check(f"{label}: annotation_note is present",
              "annotation_note" in res, f"keys={sorted(res)}")
        note = res.get("annotation_note", "")
        check(f"{label}: does not assert a breakpoint",
              "breakpoint" not in note.lower(), f"got {note!r}")
        check(f"{label}: does not assert disruption",
              "disrupt" not in note.lower(), f"got {note!r}")

    check("gene present: names what actually overlaps",
          "VANGL1" in hit["annotation_note"] and "overlap" in hit["annotation_note"],
          f"got {hit['annotation_note']!r}")
    check("intergenic: says intergenic without claiming a breakpoint",
          "intergenic" in miss["annotation_note"].lower(),
          f"got {miss['annotation_note']!r}")
    check("gene_count/is_intergenic unchanged",
          hit["gene_count"] == 1 and miss["is_intergenic"] is True,
          f"hit={hit['gene_count']} miss={miss['is_intergenic']}")

    print("\nthe retired phrase is gone from the source")
    src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "tools", "bam_tools.py"), encoding="utf-8").read()
    disrupts = [ln for ln in src.splitlines()
                if "directly disrupts" in ln and not ln.strip().startswith("#")]
    check("no emitted string still says 'directly disrupts'",
          not disrupts, f"{disrupts[:2]}")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL GENE ANNOTATION NOTE TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
