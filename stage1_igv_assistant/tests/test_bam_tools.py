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

        print("=" * 60)
        print("ALL TESTS PASSED")
        print("=" * 60)

    finally:
        os.unlink(tmp_path)
        if os.path.exists(bam_path):
            os.unlink(bam_path)
        if os.path.exists(bam_path + ".bai"):
            os.unlink(bam_path + ".bai")


if __name__ == "__main__":
    run_tests()