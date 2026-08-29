#!/usr/bin/env python3
"""
Regression test for contig names invented by normalisation.

From results/REAL_PATIENT_DATA_VALIDATION.md finding 6.

The reference in use is hs38DH. Every contig class in it is chr-prefixed
except the 525 HLA-* contigs, which are not. The old normaliser prepended
"chr" unconditionally, so an SA tag pointing at HLA-A*01:01:01:01 was keyed
under chrHLA-A*01:01:01:01 -- a contig that does not exist in the header --
and emitted as if it were a real reference contig. In one MHC window, 16 of
44 partner keys named contigs absent from the header, carrying 102 SA
entries. No error was raised.

The load-bearing assertion here is the end-to-end one: every key in
partner_chromosomes and mate_chromosomes must exist in the BAM header. A
tool that reports a contig the reference does not contain is producing a
plausible-looking wrong answer, which is worse than a crash.

Run: python3 stage1_igv_assistant/tests/test_contig_naming.py
"""
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant.tools.bam_tools import (  # noqa: E402
    _canonical_chrom,
    _chrom_count,
    get_split_reads,
    count_discordant_pairs,
)

FAILURES = []

# Mirrors the real reference's shape: primary contigs are chr-prefixed,
# alt/random/Un are chr-prefixed, HLA contigs are NOT.
HLA_A = "HLA-A*01:01:01:01"
HLA_DRB1 = "HLA-DRB1*01:01:01"
HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [
        {"SN": "chr1", "LN": 248956422},
        {"SN": "chr6", "LN": 170805979},
        {"SN": "chr6_GL000256v2_alt", "LN": 4929269},
        {"SN": "chrUn_KI270302v1", "LN": 2274},
        {"SN": HLA_A, "LN": 3503},
        {"SN": HLA_DRB1, "LN": 11080},
    ],
})
REFS = {sq["SN"] for sq in HEADER.to_dict()["SQ"]}


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def _bam_with_hla_partners(path):
    """
    Reads on chr1 whose SA tags and mates point at contigs of every class,
    including the two HLA contigs.
    """
    partners = ["chr6", "chr6_GL000256v2_alt", "chrUn_KI270302v1", HLA_A, HLA_DRB1]
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        for i in range(40):
            partner = partners[i % len(partners)]
            r = pysam.AlignedSegment(HEADER)
            r.query_name = f"r{i}"
            r.query_sequence = "A" * 100
            r.flag = 0x1
            r.reference_id = 0
            r.reference_start = 10000 + i * 10
            r.mapping_quality = 60
            r.cigar = [(4, 20), (0, 80)]
            r.next_reference_id = HEADER.to_dict()["SQ"].index(
                next(sq for sq in HEADER.to_dict()["SQ"] if sq["SN"] == partner))
            r.next_reference_start = 500
            r.template_length = 0
            r.query_qualities = pysam.qualitystring_to_array("I" * 100)
            r.set_tag("SA", f"{partner},500,+,80M20S,60,0;")
            bam.write(r)
    pysam.index(path)


def run_tests():
    print("=" * 68)
    print("CONTIG NAMING REGRESSION SUITE")
    print("=" * 68)

    # ── The normaliser itself ────────────────────────────────────────────
    print("\nnames known to the header keep the header's spelling")
    for name in ["chr1", "chr6_GL000256v2_alt", "chrUn_KI270302v1", HLA_A, HLA_DRB1]:
        check(f"{name!r} unchanged", _canonical_chrom(name, REFS) == name,
              f"got {_canonical_chrom(name, REFS)!r}")

    print("\nchr-prefix equivalence still works (the property to preserve)")
    check("'1' resolves to 'chr1'", _canonical_chrom("1", REFS) == "chr1",
          f"got {_canonical_chrom('1', REFS)!r}")
    check("both spellings normalise identically",
          _canonical_chrom("1", REFS) == _canonical_chrom("chr1", REFS))

    print("\nHLA contigs are not given a fabricated chr prefix")
    for name in [HLA_A, HLA_DRB1]:
        out = _canonical_chrom(name, REFS)
        check(f"{name!r} does not become chr{name!r}",
              not out.startswith("chrHLA"), f"got {out!r}")
        check(f"{name!r} maps to a contig that exists", out in REFS, f"got {out!r}")

    print("\nthe fabricated form is repaired, not propagated")
    check("'chrHLA-A*01:01:01:01' resolves back to the real contig",
          _canonical_chrom("chr" + HLA_A, REFS) == HLA_A,
          f"got {_canonical_chrom('chr' + HLA_A, REFS)!r}")

    print("\nunknown names are returned unchanged, never invented")
    for name in ["chr99_nonexistent", "scaffold_xyz", "HLA-ZZZ*99:99"]:
        check(f"{name!r} returned as-is", _canonical_chrom(name, REFS) == name,
              f"got {_canonical_chrom(name, REFS)!r}")

    print("\nchr-prefix tolerant count lookup")
    check("'1' finds a 'chr1' key", _chrom_count({"chr1": 7}, "1") == 7)
    check("'chr1' finds a 'chr1' key", _chrom_count({"chr1": 7}, "chr1") == 7)
    check("absent contig counts zero", _chrom_count({"chr1": 7}, "chr9") == 0)

    # ── End to end: the assertion that actually matters ──────────────────
    print("\nEND TO END: no reported partner may be absent from the header")
    tmp = tempfile.mkdtemp()
    bam = os.path.join(tmp, "hla.bam")
    _bam_with_hla_partners(bam)

    split = get_split_reads(bam, "chr1", 10200, window_bp=500, min_mapq=20)
    partners = set(split["partner_chromosomes"])
    check("split_reads found SA partners (else this proves nothing)",
          len(partners) > 0, f"got {partners}")
    check("every partner_chromosomes key exists in the header",
          partners <= REFS, f"invented: {sorted(partners - REFS)}")
    check("no key was given a chrHLA prefix",
          not any(k.startswith("chrHLA") for k in partners),
          f"got {sorted(partners)}")
    check("both HLA contigs are reported under their real names",
          HLA_A in partners and HLA_DRB1 in partners, f"got {sorted(partners)}")

    disc = count_discordant_pairs(bam, "chr1", 10200, window_bp=500, min_mapq=20)
    mates = set(disc["mate_chromosomes"]) - {"unknown"}
    check("discordant found mates (else this proves nothing)",
          len(mates) > 0, f"got {mates}")
    check("every mate_chromosomes key exists in the header",
          mates <= REFS, f"invented: {sorted(mates - REFS)}")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL CONTIG NAMING TESTS PASSED")
    print("=" * 68)


if __name__ == "__main__":
    run_tests()
