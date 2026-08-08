"""
test_bam_tools.py
Tests for bam_tools.py using a synthetic BAM file.

We create a tiny BAM in memory with known reads, then verify
our tools return the expected numbers.
"""

import os
import pysam
import tempfile
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from stage1_igv_assistant.tools.bam_tools import (
    get_bam_stats_at_locus,
    count_discordant_pairs,
    count_soft_clipped_reads,
    get_split_reads,
    summarize_breakpoint_evidence,
)


def create_synthetic_bam(path: str):
    """
    Creates a minimal BAM with:
    - 20 normal concordant reads on chr1:1000-2000
    - 5 discordant reads at chr1:1500 whose mates are on chr8
    - 3 soft-clipped reads at chr1:1490
    """
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": "chr1", "LN": 248956422},
            {"SN": "chr8", "LN": 145138636},
        ]
    })

    with pysam.AlignmentFile(path, "wb", header=header) as bam:

        # 20 normal concordant reads on chr1
        for i in range(20):
            read = pysam.AlignedSegment(header)
            read.query_name = f"normal_read_{i}"
            read.query_sequence = "A" * 100
            read.flag = 0x1 | 0x2   # paired, proper pair
            read.reference_id = 0   # chr1
            read.reference_start = 1000 + i * 50
            read.mapping_quality = 60
            read.cigar = [(0, 100)]  # 100M
            read.next_reference_id = 0  # mate on chr1
            read.next_reference_start = 1100 + i * 50
            read.template_length = 200
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(read)

        # 5 discordant reads at chr1:1500 — mates on chr8
        for i in range(5):
            read = pysam.AlignedSegment(header)
            read.query_name = f"discordant_read_{i}"
            read.query_sequence = "T" * 100
            read.flag = 0x1         # paired, NOT proper pair
            read.reference_id = 0   # chr1
            read.reference_start = 1490 + i * 2
            read.mapping_quality = 55
            read.cigar = [(0, 100)]
            read.next_reference_id = 1  # mate on chr8
            read.next_reference_start = 47000000
            read.template_length = 0
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(read)

        # 3 soft-clipped reads at chr1:1490
        for i in range(3):
            read = pysam.AlignedSegment(header)
            read.query_name = f"clipped_read_{i}"
            read.query_sequence = "G" * 100
            read.flag = 0x1
            read.reference_id = 0
            read.reference_start = 1490
            read.mapping_quality = 50
            read.cigar = [(4, 15), (0, 85)]  # 15S 85M — 15 bases soft-clipped
            read.next_reference_id = 1
            read.next_reference_start = 47000000
            read.template_length = 0
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(read)

    # Sort and index
    sorted_path = path.replace(".bam", ".sorted.bam")
    pysam.sort("-o", sorted_path, path)
    pysam.index(sorted_path)
    return sorted_path


def create_translocation_bam(path: str):
    """
    Creates a synthetic BAM simulating a balanced translocation between
    chr1 and chr8, with SA tags on split reads — mimicking the output of a
    modern chimeric-alignment-aware aligner (unlike the 2018 HCC1143 BAM,
    which carries no SA tags anywhere in the file). Used to confirm the
    breakpoint evidence tools themselves work correctly end-to-end, given
    the alignment data they need.

    - 200 normal concordant read pairs on chr1:1,000,000-1,100,000, kept
      clear of the breakpoint window so they don't dilute the signal.
    - 15 discordant pairs at chr1:1,050,000 with mates on chr8:47,000,000
    - 8 split reads at chr1:1,050,000 carrying SA tags pointing to
      chr8:47,000,000 (ordinary discordant-pair mate info is NOT set for
      these — the SA tag is the signal being tested)
    - 5 soft-clipped reads at chr1:1,050,000 (ordinary clipping, no SA tag)
    """
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": "chr1", "LN": 248956422},
            {"SN": "chr8", "LN": 145138636},
        ]
    })

    with pysam.AlignmentFile(path, "wb", header=header) as bam:

        # 200 normal concordant read pairs, flanking the breakpoint window
        # on both sides so they never enter the inspection windows below
        for i in range(100):
            read = pysam.AlignedSegment(header)
            read.query_name = f"normal_read_L{i}"
            read.query_sequence = "A" * 100
            read.flag = 0x1 | 0x2
            read.reference_id = 0
            read.reference_start = 1000000 + i * 490
            read.mapping_quality = 60
            read.cigar = [(0, 100)]
            read.next_reference_id = 0
            read.next_reference_start = read.reference_start + 150
            read.template_length = 250
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(read)

        for i in range(100):
            read = pysam.AlignedSegment(header)
            read.query_name = f"normal_read_R{i}"
            read.query_sequence = "A" * 100
            read.flag = 0x1 | 0x2
            read.reference_id = 0
            read.reference_start = 1051000 + i * 490
            read.mapping_quality = 60
            read.cigar = [(0, 100)]
            read.next_reference_id = 0
            read.next_reference_start = read.reference_start + 150
            read.template_length = 250
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(read)

        # 15 discordant pairs at chr1:1,050,000 — mates cluster on chr8:47,000,000
        # Placed just past the 200bp inspection window (but inside the 500bp
        # one) so they don't dilute the soft-clip / split-read fractions.
        for i in range(15):
            read = pysam.AlignedSegment(header)
            read.query_name = f"discordant_read_{i}"
            read.query_sequence = "T" * 100
            read.flag = 0x1               # paired, NOT proper pair
            read.reference_id = 0
            read.reference_start = 1050250 + i * 5
            read.mapping_quality = 55
            read.cigar = [(0, 100)]
            read.next_reference_id = 1    # chr8
            read.next_reference_start = 47000000 + i * 20
            read.template_length = 0
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(read)

        # 8 split reads at chr1:1,050,000 — mate stays local on chr1 (so
        # they are NOT counted as discordant pairs), but each carries an SA
        # tag pointing to chr8:47,000,000, simulating a chimeric alignment.
        for i in range(8):
            read = pysam.AlignedSegment(header)
            read.query_name = f"split_read_{i}"
            read.query_sequence = "C" * 100
            read.flag = 0x1 | 0x2         # paired, proper pair (local mate)
            read.reference_id = 0
            read.reference_start = 1050000 + i * 2
            read.mapping_quality = 60
            read.cigar = [(0, 100)]
            read.next_reference_id = 0
            read.next_reference_start = read.reference_start + 150
            read.template_length = 250
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            read.set_tag("SA", f"chr8,{47000000 + i * 20},+,30M70S,60,0;", value_type="Z")
            bam.write(read)

        # 5 soft-clipped reads at chr1:1,050,000 — ordinary clipping, no SA tag
        for i in range(5):
            read = pysam.AlignedSegment(header)
            read.query_name = f"clipped_read_{i}"
            read.query_sequence = "G" * 100
            read.flag = 0x1
            read.reference_id = 0
            read.reference_start = 1049950 + i * 3
            read.mapping_quality = 50
            read.cigar = [(4, 15), (0, 85)]   # 15S 85M
            read.next_reference_id = 0
            read.next_reference_start = read.reference_start + 150
            read.template_length = 250
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            bam.write(read)

    # Sort and index
    sorted_path = path.replace(".bam", ".sorted.bam")
    pysam.sort("-o", sorted_path, path)
    pysam.index(sorted_path)
    return sorted_path


def run_tests():
    print("=" * 60)
    print("BAM TOOLS TEST SUITE")
    print("=" * 60)

    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
        tmp_path = f.name

    try:
        bam_path = create_synthetic_bam(tmp_path)
        print(f"\nSynthetic BAM created: {bam_path}\n")

        # ── Test 1: get_bam_stats_at_locus ──────────────────────────
        print("TEST 1: get_bam_stats_at_locus")
        stats = get_bam_stats_at_locus(bam_path, "chr1", 1000, 2000)
        print(f"  Total reads: {stats['total_reads']}")
        print(f"  Mean depth:  {stats['mean_depth']}")
        print(f"  Mean MAPQ:   {stats['mean_mapq']}")
        assert stats["total_reads"] >= 20, "Expected at least 20 reads"
        print("  PASSED ✓\n")

        # ── Test 2: count_discordant_pairs ───────────────────────────
        print("TEST 2: count_discordant_pairs")
        disc = count_discordant_pairs(
            bam_path, "chr1", position=1500, window_bp=200, min_mapq=0
        )
        print(f"  Total reads in window: {disc['total_reads_in_window']}")
        print(f"  Discordant pairs:      {disc['discordant_pairs']}")
        print(f"  Discordant fraction:   {disc['discordant_fraction']}")
        print(f"  Mate chromosomes:      {disc['mate_chromosomes']}")
        assert disc["discordant_pairs"] >= 5, "Expected 5 discordant pairs"
        assert "chr8" in disc["mate_chromosomes"], "Expected mates on chr8"
        print("  PASSED ✓\n")

        # ── Test 3: count_soft_clipped_reads ────────────────────────
        print("TEST 3: count_soft_clipped_reads")
        clips = count_soft_clipped_reads(
            bam_path, "chr1", position=1490, window_bp=100,
            min_clip_bases=10, min_mapq=0
        )
        print(f"  Total reads in window:  {clips['total_reads_in_window']}")
        print(f"  Soft-clipped reads:     {clips['soft_clipped_reads']}")
        print(f"  Clip fraction:          {clips['soft_clipped_fraction']}")
        print(f"  Consensus clip pos:     {clips['consensus_clip_position']}")
        assert clips["soft_clipped_reads"] >= 3, "Expected 3 clipped reads"
        print("  PASSED ✓\n")

    finally:
        os.unlink(tmp_path)
        if os.path.exists(bam_path):
            os.unlink(bam_path)
        if os.path.exists(bam_path + ".bai"):
            os.unlink(bam_path + ".bai")

    # ── TEST 6: full 5-tool pipeline on a synthetic translocation ──────────
    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
        tmp_path2 = f.name

    try:
        translocation_bam = create_translocation_bam(tmp_path2)
        print(f"Translocation BAM created: {translocation_bam}\n")

        print("TEST 6: full breakpoint evidence pipeline on synthetic translocation")
        chrom, pos = "chr1", 1050000

        disc = count_discordant_pairs(translocation_bam, chrom, pos)
        print(f"  Discordant pairs: {disc['discordant_pairs']} / {disc['total_reads_in_window']}"
              f"  mate_chromosomes={disc['mate_chromosomes']}")
        assert disc["discordant_pairs"] == 15, "Expected exactly 15 discordant pairs"
        assert set(disc["mate_chromosomes"].keys()) == {"chr8"}, "Expected all mates on chr8"

        clips = count_soft_clipped_reads(translocation_bam, chrom, pos)
        print(f"  Soft-clipped reads: {clips['soft_clipped_reads']} / {clips['total_reads_in_window']}")
        assert clips["soft_clipped_reads"] == 5, "Expected exactly 5 soft-clipped reads"

        split = get_split_reads(translocation_bam, chrom, pos)
        print(f"  Split reads: {split['split_reads']} / {split['total_reads_in_window']}"
              f"  partner_chromosomes={split['partner_chromosomes']}")
        assert split["split_reads"] == 8, "Expected exactly 8 split reads"
        assert set(split["partner_chromosomes"].keys()) == {"chr8"}, "Expected split-read partners on chr8"

        summary = summarize_breakpoint_evidence(translocation_bam, chrom, pos, label="synthetic_translocation")
        print(f"  Evidence score: {summary['evidence_score']}  strength: {summary['evidence_strength']}"
              f"  signal_layers: {summary['signal_layers']}")
        for obs in summary["supporting_observations"]:
            print(f"    - {obs}")
        assert summary["evidence_strength"] == "strong", "Expected STRONG evidence"
        assert summary["signal_layers"] == "3/3", "Expected all 3 evidence layers to show signal"

        print("  PASSED ✓\n")

    finally:
        os.unlink(tmp_path2)
        if os.path.exists(translocation_bam):
            os.unlink(translocation_bam)
        if os.path.exists(translocation_bam + ".bai"):
            os.unlink(translocation_bam + ".bai")

    print("=" * 60)
    print("ALL TESTS PASSED")
    print("=" * 60)


if __name__ == "__main__":
    run_tests()