#!/usr/bin/env python3
"""
Regression test: a supplementary alignment on a decoy contig is not a split-read
partner.

From results/decoy_partners_2026-09-27.json (Phase 13) and the fix registered in
results/decoy_fix_2026-09-27/registration.json.

get_split_reads kept every SA entry whose own mapQ passed, whatever its contig. An
entry naming a decoy (hs38DH's chrUn_*_decoy sequences, which absorb reads from
sequence missing from the assembly) made the read a split read and put the decoy
in partner_chromosomes. At the public NA12878 background breakend chr21:10,770,078,
266 of 269 split reads were partnered only with a decoy, and the summary rated the
position 72.5 "strong"; without them it is 55.0 "moderate".

Now a decoy entry is set aside: it is not a partner, a read left with no partner is
not a split read, and the decoy is reported in decoy_partners (reads per decoy
contig) with decoy_only_reads (reads dropped for having only decoy partners).

POSITIVE CONTROL: reads whose only SA entry names a decoy -- counted on the old
code, so this suite fails there. NEGATIVE CONTROL: an ordinary inter-chromosomal
entry (chr8) is still counted, including on a read that also names a decoy.

To run against another copy of bam_tools.py (e.g. the pre-fix file):
    BAM_TOOLS_FILE=/path/to/bam_tools.py python3 stage1_igv_assistant/tests/test_decoy_partners.py
Run: python3 stage1_igv_assistant/tests/test_decoy_partners.py
"""
import importlib.util
import os
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

if os.environ.get("BAM_TOOLS_FILE"):
    _spec = importlib.util.spec_from_file_location("bam_tools_under_test", os.environ["BAM_TOOLS_FILE"])
    bt = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(bt)
else:
    from stage1_igv_assistant.tools import bam_tools as bt  # noqa: E402

FAILURES = []
DECOY = "chrUn_JTFH01000001v1_decoy"
HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 1000000}, {"SN": "chr8", "LN": 900000}, {"SN": DECOY, "LN": 5000}],
})
LAYERS = ["discordant_pairs", "soft_clipped_reads", "split_reads", "read_depth"]


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def _bam(path, tags):
    """One primary read per SA tag in `tags` (None = no SA), MAPQ 60, at chr1:10,000+."""
    with pysam.AlignmentFile(path, "wb", header=HEADER) as bam:
        for i, tag in enumerate(tags):
            r = pysam.AlignedSegment(HEADER)
            r.query_name = f"r{i}"
            r.query_sequence = "A" * 100
            r.flag = 0x1
            r.reference_id = 0
            r.reference_start = 10000 + i * 5
            r.mapping_quality = 60
            r.cigar = [(4, 20), (0, 80)]
            r.next_reference_id = 0
            r.next_reference_start = 10200
            r.query_qualities = pysam.qualitystring_to_array("I" * 100)
            if tag:
                r.set_tag("SA", tag)
            bam.write(r)
    pysam.index(path)
    return path


def run_tests():
    print("=" * 68)
    print("DECOY PARTNER REGRESSION SUITE")
    print("=" * 68)
    tmp = tempfile.mkdtemp()
    decoy = f"{DECOY},1000,+,80M20S,60,0;"
    chr8 = "chr8,500000,+,80M20S,60,0;"

    print("\npositive control: a decoy-only SA entry is not a partner")
    r = bt.get_split_reads(_bam(os.path.join(tmp, "decoy.bam"), [decoy] * 20 + [None] * 20), "chr1", 10050)
    check("reads whose only partner is a decoy are not split reads", r["split_reads"] == 0,
          f"got split_reads={r['split_reads']}")
    check("the decoy is not in partner_chromosomes", DECOY not in r["partner_chromosomes"],
          f"got {r['partner_chromosomes']}")
    check("the decoy is reported in decoy_partners", r.get("decoy_partners") == {DECOY: 20},
          f"got {r.get('decoy_partners')}")
    check("decoy_only_reads counts the dropped reads", r.get("decoy_only_reads") == 20,
          f"got {r.get('decoy_only_reads')}")
    check("the fraction no longer counts them", r["split_read_fraction"] == 0.0,
          f"got {r['split_read_fraction']}")
    check("decoy entries are not partner strands",
          r["partner_strand_concordant"] + r["partner_strand_flipped"] == 0,
          f"got {r['partner_strand_concordant']}/{r['partner_strand_flipped']}")

    print("\npositive control through the summary: no split score, no decoy named")
    s = bt.summarize_breakpoint_evidence(os.path.join(tmp, "decoy.bam"), "chr1", 10050, applicable_layers=LAYERS)
    check("split_read_score is 0 when every SA partner is a decoy", s["split_read_score"] == 0.0,
          f"got {s['split_read_score']}")
    check("no observation names the decoy",
          not any(DECOY in o for o in s["supporting_observations"]),
          f"got {[o for o in s['supporting_observations'] if 'split' in o]}")

    print("\nnegative control: an ordinary inter-chromosomal entry is still counted")
    n = bt.get_split_reads(_bam(os.path.join(tmp, "interchrom.bam"), [chr8] * 20 + [None] * 20), "chr1", 10050)
    check("chr8 partners are counted", n["split_reads"] == 20, f"got {n['split_reads']}")
    check("chr8 is the partner", n["partner_chromosomes"] == {"chr8": 20}, f"got {n['partner_chromosomes']}")
    check("no decoy is reported", n.get("decoy_partners") == {} and n.get("decoy_only_reads") == 0,
          f"got {n.get('decoy_partners')} / {n.get('decoy_only_reads')}")

    print("\nnegative control: a read naming chr8 and a decoy is counted once, for chr8")
    m = bt.get_split_reads(_bam(os.path.join(tmp, "mixed.bam"), [chr8 + decoy] * 10 + [None] * 10), "chr1", 10050)
    check("mixed reads are still split reads", m["split_reads"] == 10, f"got {m['split_reads']}")
    check("only chr8 is a partner", m["partner_chromosomes"] == {"chr8": 10}, f"got {m['partner_chromosomes']}")
    check("the decoy is still reported", m.get("decoy_partners") == {DECOY: 10}, f"got {m.get('decoy_partners')}")
    check("mixed reads are not decoy-only", m.get("decoy_only_reads") == 0, f"got {m.get('decoy_only_reads')}")

    print("\nthe mapQ filter still comes first: a low-mapQ decoy entry is a mapQ drop, not a decoy")
    q = bt.get_split_reads(_bam(os.path.join(tmp, "lowq.bam"), [f"{DECOY},1000,+,80M20S,0,0;"] * 10), "chr1", 10050)
    check("sa_entries_below_min_mapq counts it", q["sa_entries_below_min_mapq"] == 10,
          f"got {q['sa_entries_below_min_mapq']}")
    check("and it is not reported as a decoy partner", q.get("decoy_partners") == {}, f"got {q.get('decoy_partners')}")

    print("\n" + "=" * 68)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S)")
        for f in FAILURES:
            print(f"  - {f}")
        sys.exit(1)
    print("ALL DECOY PARTNER TESTS PASSED")


if __name__ == "__main__":
    run_tests()
