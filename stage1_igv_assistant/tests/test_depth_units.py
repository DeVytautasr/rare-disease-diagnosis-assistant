#!/usr/bin/env python3
"""
Regression test for the stale depth unit label.

From results/REAL_PATIENT_DATA_VALIDATION.md finding 11. It sits in the
sentence a reporting assistant is most likely to paste verbatim.

The depth observation said "reads/window" for a quantity that
has been true per-base depth since the 2026-08-11 FIX-1, contradicting the
function's own docstring. The adjacent off-position branch printed the same
quantities with no unit at all.

Run: python3 stage1_igv_assistant/tests/test_depth_units.py
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
    print("DEPTH UNIT LABEL REGRESSION SUITE")
    print("=" * 68)

    print("\ndepth observations carry the right unit")
    tmp = tempfile.mkdtemp()
    bam = os.path.join(tmp, "dip.bam")
    _depth_dip_bam(bam)

    scored = summarize_breakpoint_evidence(bam, "chr1", 10000)
    obs = " ".join(scored["supporting_observations"]) \
        if isinstance(scored.get("supporting_observations"), list) \
        else str(scored.get("supporting_observations", ""))
    combined = obs + " " + scored.get("interpretation_template", "")

    check("depth layer actually fired (otherwise this test proves nothing)",
          "Read depth" in combined, f"observations={combined[:160]!r}")
    check("no stale 'reads/window' label anywhere in the output",
          "reads/window" not in combined, f"got {combined[:200]!r}")
    check("the depth quantity is labelled per-base",
          "per-base depth" in combined, f"got {combined[:200]!r}")

    off = summarize_breakpoint_evidence(bam, "chr1", 8600)
    off_obs = " ".join(off["supporting_observations"]) \
        if isinstance(off.get("supporting_observations"), list) \
        else str(off.get("supporting_observations", ""))
    if "off-position depth feature" in off_obs:
        check("off-position branch also carries a unit",
              "per-base depth" in off_obs, f"got {off_obs[:200]!r}")
        check("off-position branch has no 'reads/window'",
              "reads/window" not in off_obs, f"got {off_obs[:200]!r}")
    else:
        print("  (off-position branch not triggered here — scored branch covers the label)")

    print("\nthe retired label is gone from the source")
    src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "tools", "bam_tools.py"), encoding="utf-8").read()
    emitted = [ln for ln in src.splitlines()
               if "reads/window" in ln and not ln.strip().startswith("#")]
    check("no emitted string still says reads/window", not emitted, f"{emitted[:2]}")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL DEPTH UNIT TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
