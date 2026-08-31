"""
test_bam_tools.py
Tests for bam_tools.py using a synthetic BAM file.

We create a tiny BAM in memory with known reads, then verify
our tools return the expected numbers.
"""

import os
import pysam
import shutil
import tempfile
import time
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from stage1_igv_assistant.tools.bam_tools import (
    get_bam_stats_at_locus,
    count_discordant_pairs,
    count_soft_clipped_reads,
    get_split_reads,
    get_read_depth_profile,
    summarize_breakpoint_evidence,
    get_gene_at_locus,
    check_reciprocal_breakpoint,
    run_igv_screenshot,
    detect_applicable_layers,
    igv_evidence_panel,
    EVIDENCE_LAYER_NAMES,
    DEPTH_RATIO_DELETION_THRESHOLD,
    DEFAULT_PANEL_WINDOWS,
)
import stage1_igv_assistant.tools.bam_tools as bam_tools_module
from stage1_igv_assistant.case_object import BamCase


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
    - 15 reciprocal discordant pairs at chr8:47,000,000 with mates back on
      chr1:1,050,000 — a real balanced translocation shows discordant
      evidence on BOTH derivative chromosomes, not just one side
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

        # 15 reciprocal discordant pairs at chr8:47,000,000 — mates back on chr1
        for i in range(15):
            read = pysam.AlignedSegment(header)
            read.query_name = f"reciprocal_discordant_read_{i}"
            read.query_sequence = "T" * 100
            read.flag = 0x1               # paired, NOT proper pair
            read.reference_id = 1         # chr8
            read.reference_start = 47000000 + i * 5
            read.mapping_quality = 55
            read.cigar = [(0, 100)]
            read.next_reference_id = 0    # chr1
            read.next_reference_start = 1050000 + i * 20
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


def create_depth_dropout_bam(path: str):
    """
    Creates a synthetic BAM purpose-built for exercising get_read_depth_profile
    directly: dense, uniform tiled coverage on chr1:0-1100 and chr1:2000-3100,
    with a deliberate zero-read gap across chr1:1100-2000 in between. This is
    a pure depth-profile fixture — it makes no claim about representing any
    specific SV type (unlike create_translocation_bam's incidental read-free
    gaps, which TEST 6 already notes are a construction artifact, not a
    depth-profile test in their own right).
    """
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 248956422}],
    })
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        read_idx = 0
        for region_start in (0, 2000):
            for offset in range(0, 1000, 20):   # 50 tiled start positions
                for _ in range(5):               # 5x overlapping coverage each
                    read = pysam.AlignedSegment(header)
                    read.query_name = f"flank_read_{read_idx}"
                    read_idx += 1
                    read.query_sequence = "A" * 100
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = region_start + offset
                    read.mapping_quality = 60
                    read.cigar = [(0, 100)]
                    read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                    bam.write(read)

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
        print(f"  Low-MAPQ fraction: {stats['low_mapq_fraction']}")
        print(f"  Forward/reverse:   {stats['forward_reads']}/{stats['reverse_reads']}")
        assert stats["total_reads"] >= 20, "Expected at least 20 reads"
        # Exact values for this fixture (28 reads: 20 normal @ MAPQ 60,
        # 5 discordant @ MAPQ 55, 3 clipped @ MAPQ 50; none reverse-strand,
        # none below MAPQ 20) — previously only total_reads was asserted,
        # so a regression in mean_depth/mean_mapq/low_mapq_fraction/strand
        # counting would have passed silently.
        assert stats["total_reads"] == 28, f"Expected exactly 28 reads, got {stats['total_reads']}"
        assert stats["mean_depth"] == 2.71, f"Expected mean_depth=2.71, got {stats['mean_depth']}"
        assert stats["mean_mapq"] == 58.04, f"Expected mean_mapq=58.04, got {stats['mean_mapq']}"
        assert stats["low_mapq_fraction"] == 0.0, \
            f"Expected low_mapq_fraction=0.0 (no reads below MAPQ 20), got {stats['low_mapq_fraction']}"
        assert stats["forward_reads"] == 28, f"Expected all 28 reads forward, got {stats['forward_reads']}"
        assert stats["reverse_reads"] == 0, f"Expected 0 reverse reads, got {stats['reverse_reads']}"
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
        print(f"  Dominant clip side:     {clips['dominant_clip_side']}")
        assert clips["soft_clipped_reads"] >= 3, "Expected 3 clipped reads"
        # Exact values for this fixture: 3 reads, all left-clipped (15S85M),
        # all at the same reference_start (1490) — previously only the
        # raw count was asserted, so consensus_clip_position,
        # max_clips_at_position, and the left/right split could silently
        # regress (e.g. a left/right mixup, or a broken tie-break) without
        # failing this test.
        assert clips["consensus_clip_position"] == 1490, \
            f"Expected consensus_clip_position=1490, got {clips['consensus_clip_position']}"
        assert clips["max_clips_at_position"] == 3, \
            f"Expected max_clips_at_position=3 (all 3 reads pile up at 1490), got {clips['max_clips_at_position']}"
        assert clips["left_clip_reads"] == 3, f"Expected 3 left-clipped reads, got {clips['left_clip_reads']}"
        assert clips["right_clip_reads"] == 0, f"Expected 0 right-clipped reads, got {clips['right_clip_reads']}"
        assert clips["left_consensus_position"] == 1490, \
            f"Expected left_consensus_position=1490, got {clips['left_consensus_position']}"
        assert clips["left_max_clips_at_position"] == 3, \
            f"Expected left_max_clips_at_position=3, got {clips['left_max_clips_at_position']}"
        assert clips["right_consensus_position"] is None, \
            f"Expected right_consensus_position=None (no right clips), got {clips['right_consensus_position']}"
        assert clips["right_max_clips_at_position"] == 0, \
            f"Expected right_max_clips_at_position=0, got {clips['right_max_clips_at_position']}"
        assert clips["dominant_clip_side"] == "left", \
            f"Expected dominant_clip_side='left', got {clips['dominant_clip_side']}"
        print("  PASSED ✓\n")

    finally:
        os.unlink(tmp_path)
        if os.path.exists(bam_path):
            os.unlink(bam_path)
        if os.path.exists(bam_path + ".bai"):
            os.unlink(bam_path + ".bai")

    # get_read_depth_profile previously had zero direct test coverage — not
    # even imported in this file, only exercised indirectly (and unasserted
    # on its own fields) as a sub-call inside summarize_breakpoint_evidence.
    # TEST 4 and TEST 5 fill that gap directly, using a fixture with a real,
    # deliberate depth dropout rather than piggybacking on another test's BAM.
    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
        tmp_path_depth = f.name
    try:
        depth_bam = create_depth_dropout_bam(tmp_path_depth)
        print(f"Depth-dropout BAM created: {depth_bam}\n")

        # ── TEST 4: get_read_depth_profile on a uniform (no-dropout) region ──
        print("TEST 4: get_read_depth_profile (uniform region, no dropout)")
        flat_profile = get_read_depth_profile(depth_bam, "chr1", 0, 1000, window_size=200)
        assert "error" not in flat_profile, f"Unexpected error: {flat_profile}"
        assert len(flat_profile["windows"]) == 5, f"Expected 5 windows, got {len(flat_profile['windows'])}"
        assert flat_profile["summary"]["min_depth"] > 0, "Expected nonzero depth throughout this region"
        assert flat_profile["summary"]["depth_ratio_min_to_mean"] >= DEPTH_RATIO_DELETION_THRESHOLD, \
            f"Expected ratio >= {DEPTH_RATIO_DELETION_THRESHOLD} for a uniform region, " \
            f"got {flat_profile['summary']['depth_ratio_min_to_mean']}"
        assert flat_profile["summary"]["likely_deletion"] is False, \
            "Expected likely_deletion=False for a uniform-depth region"
        print(f"  min_depth={flat_profile['summary']['min_depth']} "
              f"mean_depth={flat_profile['summary']['mean_depth']} "
              f"ratio={flat_profile['summary']['depth_ratio_min_to_mean']} "
              f"likely_deletion={flat_profile['summary']['likely_deletion']}")
        print("  PASSED ✓\n")

        # ── TEST 5: get_read_depth_profile across the engineered dropout ──
        print("TEST 5: get_read_depth_profile (region spanning the depth dropout)")
        dropout_profile = get_read_depth_profile(depth_bam, "chr1", 0, 3000, window_size=200)
        assert "error" not in dropout_profile, f"Unexpected error: {dropout_profile}"
        assert len(dropout_profile["windows"]) == 15, f"Expected 15 windows, got {len(dropout_profile['windows'])}"
        assert dropout_profile["summary"]["min_depth"] == 0, \
            f"Expected min_depth=0 inside the engineered gap, got {dropout_profile['summary']['min_depth']}"
        # Every window fully inside chr1:1100-2000 should read exactly 0.
        zero_windows = [w for w in dropout_profile["windows"]
                         if w["window_start"] >= 1200 and w["window_end"] <= 2000]
        assert zero_windows and all(w["depth"] == 0 for w in zero_windows), \
            f"Expected all windows strictly inside the gap to have depth 0, got {zero_windows}"
        assert dropout_profile["summary"]["depth_ratio_min_to_mean"] < DEPTH_RATIO_DELETION_THRESHOLD, \
            f"Expected ratio < {DEPTH_RATIO_DELETION_THRESHOLD} across the dropout, " \
            f"got {dropout_profile['summary']['depth_ratio_min_to_mean']}"
        assert dropout_profile["summary"]["likely_deletion"] is True, \
            "Expected likely_deletion=True across the engineered dropout"
        print(f"  min_depth={dropout_profile['summary']['min_depth']} "
              f"mean_depth={dropout_profile['summary']['mean_depth']} "
              f"ratio={dropout_profile['summary']['depth_ratio_min_to_mean']} "
              f"likely_deletion={dropout_profile['summary']['likely_deletion']}")
        print("  PASSED ✓\n")
    finally:
        os.unlink(tmp_path_depth)
        if os.path.exists(depth_bam):
            os.unlink(depth_bam)
        if os.path.exists(depth_bam + ".bai"):
            os.unlink(depth_bam + ".bai")

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
        print(f"  Soft-clipped reads: {clips['soft_clipped_reads']} / {clips['total_reads_in_window']}"
              f"  max_clips_at_position={clips['max_clips_at_position']}")
        assert clips["soft_clipped_reads"] == 5, "Expected exactly 5 soft-clipped reads"
        # This fixture's 5 clipped reads each start 3bp apart (reference_start =
        # 1049950 + i*3) rather than piling up at one position like
        # create_synthetic_bam's clip fixture does — so max_clips_at_position
        # is 1, not 5. That matters below: it's below count_soft_clipped_reads'
        # own documented "< 3 = no real pileup" cutoff, which
        # summarize_breakpoint_evidence's soft_clip_score now scores on
        # directly (FIX 3) instead of soft_clipped_fraction.
        assert clips["max_clips_at_position"] == 1, \
            f"Expected max_clips_at_position=1 (reads staggered, no real pileup), got {clips['max_clips_at_position']}"

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
        # discordant_pairs and split_reads are the only layers that score here:
        # discordant_fraction and split_read_fraction both clear the top
        # (0.5/0.3) tier, so both score 25 full points — real, intentional
        # signal, unaffected by FIX 2/FIX 3.
        assert summary["discordant_pair_score"] == 25.0, \
            f"Expected discordant_pair_score=25.0, got {summary['discordant_pair_score']}"
        assert summary["split_read_score"] == 25.0, \
            f"Expected split_read_score=25.0, got {summary['split_read_score']}"
        # soft_clip_score is 0 (FIX 3): max_clips_at_position=1 (see above) is
        # below the pileup cutoff, regardless of soft_clipped_fraction being
        # high — this fixture never simulated a genuine same-position pileup,
        # it just happened to score full credit under the old fraction-only
        # scoring, which is exactly the failure mode FIX 3 closes.
        assert summary["soft_clip_score"] == 0.0, \
            f"Expected soft_clip_score=0.0 (no real pileup — max_clips_at_position=1), " \
            f"got {summary['soft_clip_score']}"
        # depth_score is 0 (FIX 2): the read-free gaps deliberately built into
        # this synthetic BAM (to keep discordant/clip/split fractions clean)
        # trip the depth-ratio threshold, but that dip sits well outside
        # dip_tolerance_bp of the queried position — an off-position depth
        # feature, not a signal AT this breakpoint. Confirm it's surfaced as
        # an observation instead of silently dropped.
        assert summary["depth_score"] == 0.0, \
            f"Expected depth_score=0.0 (off-position dip, gated by dip_is_at_focus), " \
            f"got {summary['depth_score']}"
        assert any("off-position depth feature" in obs for obs in summary["supporting_observations"]), \
            "Expected an off-position depth feature observation, got: " + repr(summary["supporting_observations"])
        assert summary["evidence_strength"] == "moderate", \
            f"Expected MODERATE evidence (25 discordant + 25 split only), got {summary['evidence_strength']!r}"
        assert summary["signal_layers"] == "2/4", \
            f"Expected 2/4 (discordant_pairs, split_reads only), got {summary['signal_layers']!r}"
        # Regression test for the bug in commit c093ed4: the 4 component
        # scores must always sum exactly to evidence_score. This broke in
        # production once (found by an LLM session, not by this test suite)
        # and there was no test that would have caught it.
        component_sum = (
            summary["discordant_pair_score"] + summary["soft_clip_score"]
            + summary["split_read_score"] + summary["depth_score"]
        )
        assert component_sum == summary["evidence_score"], (
            f"Component scores ({summary['discordant_pair_score']} + "
            f"{summary['soft_clip_score']} + {summary['split_read_score']} + "
            f"{summary['depth_score']} = {component_sum}) must sum exactly to "
            f"evidence_score ({summary['evidence_score']})"
        )
        print(f"  Component sum check: {summary['discordant_pair_score']} + "
              f"{summary['soft_clip_score']} + {summary['split_read_score']} + "
              f"{summary['depth_score']} == {summary['evidence_score']} — confirmed")

        print("  PASSED ✓\n")

        # ── TEST 8: check_reciprocal_breakpoint ─────────────────────────
        print("TEST 8: check_reciprocal_breakpoint")
        reciprocal = check_reciprocal_breakpoint(
            translocation_bam, "chr1", 1050000, "chr8", 47000000
        )
        print(f"  Verdict: {reciprocal['verdict']}")
        print(f"  is_balanced: {reciprocal['is_balanced']}")
        assert "RECIPROCAL" in reciprocal["verdict"], "Expected verdict to contain RECIPROCAL"
        assert reciprocal["is_balanced"] is True, "Expected is_balanced == True"
        print("  PASSED ✓\n")

    finally:
        os.unlink(tmp_path2)
        if os.path.exists(translocation_bam):
            os.unlink(translocation_bam)
        if os.path.exists(translocation_bam + ".bai"):
            os.unlink(translocation_bam + ".bai")

    # ── TEST 7: get_gene_at_locus (LIVE Ensembl) ─────────────────────────
    # This is the only assertion in the suite that depends on a third-party
    # service being up, and it used to fail the whole run when Ensembl was
    # not. Measured on 2026-08-31: 4 of 6 consecutive full-suite runs died
    # here on HTTP 500/503, while the same queries returned 200 from curl
    # seconds later — Ensembl throttles under repeated querying and expresses
    # it as a 5xx. The suite's exit code must not be a report on Ensembl's
    # weather.
    #
    # Two changes, both following conventions already used in this file and
    # in tests/test_gene_annotation_note.py:
    #   - unavailable or throttled now prints SKIPPED with the reason, the
    #     same way TEST 17 and TEST 19 already handle the GIAB BAM download
    #   - the retries now back off. The old loop made 4 calls with no delay,
    #     so all four landed inside one throttle window and failed together.
    #     Note this is NOT redundant with get_gene_at_locus's own backoff:
    #     that function retries 429 and connection exceptions, but on any
    #     other non-200 status it records the code and breaks immediately,
    #     so a 500 costs exactly one request and returns instantly.
    #
    # The logic this test covers is ALSO covered offline and unconditionally
    # by TEST 7b below, which stubs the HTTP layer. TEST 7 exists to confirm
    # the live contract still holds — the response shape, not the wording.
    print("TEST 7: get_gene_at_locus (live Ensembl — SKIPPED if unavailable)")

    def _live_ensembl(chromosome, position, attempts=4):
        """
        Live lookup with exponential backoff between attempts. Returns the
        result dict, or None if the service never responded — None means
        "could not test", never "test failed".
        """
        last_error = None
        for attempt in range(attempts):
            result = get_gene_at_locus(chromosome, position)
            if "error" not in result:
                return result
            last_error = result["error"]
            if attempt < attempts - 1:
                wait = 2 ** attempt          # 1s, 2s, 4s
                print(f"  (Ensembl attempt {attempt + 1} failed: {last_error}"
                      f" — retrying in {wait}s)")
                time.sleep(wait)
        print(f"  (Ensembl did not respond after {attempts} attempts:"
              f" {last_error})")
        return None

    gene_result = _live_ensembl("chr1", 115686862)
    intergenic_result = _live_ensembl("chr13", 15000000)

    if gene_result is None or intergenic_result is None:
        print("  SKIPPED — Ensembl REST API unavailable or throttled; "
              "the annotation logic is still covered offline by TEST 7b\n")
    else:
        print(f"  Gene(s) at chr1:115686862: {[g['gene_name'] for g in gene_result['genes']]}")
        assert gene_result["gene_count"] >= 1, "Expected at least 1 gene"
        assert gene_result["is_intergenic"] is False, "Expected is_intergenic == False"
        print(f"  chr13:15000000 is_intergenic: {intergenic_result['is_intergenic']}")
        assert intergenic_result["is_intergenic"] is True, "Expected is_intergenic == True"
        print("  PASSED ✓\n")

    # ── TEST 7b: get_gene_at_locus, offline ──────────────────────────────
    # The unconditional counterpart. Stubs the HTTP layer the same way
    # tests/test_gene_annotation_note.py and tests/test_depth_units.py do, so
    # the gene/intergenic branch is exercised on every run, on any machine,
    # with or without a network. If TEST 7 is skipped, this still ran.
    print("TEST 7b: get_gene_at_locus offline (stubbed Ensembl)")

    class _FakeResponse:
        def __init__(self, payload):
            self.status_code = 200
            self._payload = payload

        def json(self):
            return self._payload

    class _StubEnsembl:
        """Replaces _requests.get so this test never touches the network."""

        def __init__(self, payload):
            self.payload = payload
            self._real = None

        def __enter__(self):
            self._real = bam_tools_module._requests.get
            bam_tools_module._requests.get = lambda *a, **k: _FakeResponse(self.payload)
            return self

        def __exit__(self, *exc):
            bam_tools_module._requests.get = self._real

    _GENE_HIT = [{
        "gene_id": "ENSG00000173218", "external_name": "VANGL1",
        "biotype": "protein_coding", "strand": 1,
        "start": 115641854, "end": 115698224,
    }]

    with _StubEnsembl(_GENE_HIT):
        offline_hit = get_gene_at_locus("chr1", 115686862)
    with _StubEnsembl([]):
        offline_miss = get_gene_at_locus("chr13", 15000000)

    assert "error" not in offline_hit, f"stubbed lookup errored: {offline_hit}"
    assert offline_hit["gene_count"] == 1, f"got {offline_hit['gene_count']}"
    assert offline_hit["is_intergenic"] is False, "gene overlap must not read as intergenic"
    assert offline_hit["genes"][0]["gene_name"] == "VANGL1", \
        f"got {offline_hit['genes'][0]['gene_name']}"

    assert "error" not in offline_miss, f"stubbed lookup errored: {offline_miss}"
    assert offline_miss["gene_count"] == 0, f"got {offline_miss['gene_count']}"
    assert offline_miss["is_intergenic"] is True, "empty overlap must read as intergenic"

    # The two branches must not produce the same note, which is the failure
    # finding 13 was about.
    assert offline_hit["annotation_note"] != offline_miss["annotation_note"], \
        "gene-overlap and intergenic notes are identical"
    assert "overlap" in offline_hit["annotation_note"]
    assert "intergenic" in offline_miss["annotation_note"]
    print(f"  overlap note:    {offline_hit['annotation_note']}")
    print(f"  intergenic note: {offline_miss['annotation_note']}")
    print("  PASSED ✓\n")

    # ── TEST 9: BamCase object — save/load round-trip ─────────────────────
    print("TEST 9: BamCase object")
    case = BamCase.new("test_case_009", "test.bam")
    case.clinical.hpo_terms = ["HP:0001250", "HP:0000486"]
    case.clinical.suspected_sv_type = "balanced_translocation"
    case.clinical.cytogenetic_result = "t(1;8)(p36;q24)"
    bp1 = case.add_breakpoint("BP1", "chr1", 1050000)
    bp1.evidence_strength = "strong"
    bp1.evidence_score = 87.5
    bp1.signal_layers = "3/4"
    bp1.supporting_observations = [
        "12 discordant pairs (all mates on chr8)",
        "6 soft-clipped reads (consensus position 1050002)",
    ]
    case.add_breakpoint("BP2", "chr8", 47000000)
    assert len(case.breakpoints) == 2, "Expected 2 breakpoints"

    case_path = "/tmp/test_case_009.json"
    case.save(case_path)
    loaded = BamCase.load(case_path)
    assert loaded.case_id == "test_case_009", "Expected case_id to match after reload"
    assert len(loaded.breakpoints) == 2, "Expected 2 breakpoints after reload"

    # Clinical fields must round-trip, not just the identity fields.
    assert loaded.clinical.hpo_terms == ["HP:0001250", "HP:0000486"], \
        "Expected hpo_terms to survive save/load"
    assert loaded.clinical.suspected_sv_type == "balanced_translocation"
    assert loaded.clinical.cytogenetic_result == "t(1;8)(p36;q24)"

    # Breakpoint evidence fields (not just label/chromosome/position) must
    # also survive round-trip, including the supporting_observations list.
    loaded_bp1 = loaded.get_breakpoint("BP1")
    assert loaded_bp1 is not None, "Expected to retrieve BP1 by label after reload"
    assert loaded_bp1.evidence_strength == "strong"
    assert loaded_bp1.evidence_score == 87.5
    assert loaded_bp1.signal_layers == "3/4"
    assert loaded_bp1.supporting_observations == [
        "12 discordant pairs (all mates on chr8)",
        "6 soft-clipped reads (consensus position 1050002)",
    ], "Expected supporting_observations list to survive save/load exactly"
    os.unlink(case_path)

    assert loaded.sequencing.is_paired is True, "Default is_paired should be True"
    assert "discordant_pairs" in loaded.sequencing.applicable_evidence_layers, \
        "Expected discordant_pairs in applicable_evidence_layers for paired short-read data"
    print(f"  applicable_evidence_layers (default, paired short-read): "
          f"{loaded.sequencing.applicable_evidence_layers}")

    # Technology-dependent branch: unpaired long-read data should drop
    # discordant_pairs and pick up split_reads via has_sa_tags.
    loaded.sequencing.is_paired = False
    loaded.sequencing.technology = "pacbio_hifi"
    long_read_layers = loaded.sequencing.applicable_evidence_layers
    assert "discordant_pairs" not in long_read_layers, \
        "Expected discordant_pairs excluded for unpaired long-read data"
    assert "split_reads" in long_read_layers, \
        "Expected split_reads included for pacbio_hifi (has_sa_tags)"
    print(f"  applicable_evidence_layers (unpaired PacBio HiFi): {long_read_layers}")
    print("  PASSED ✓\n")

    # ── TEST 10: run_igv_screenshot ───────────────────────────────────────
    print("TEST 10: run_igv_screenshot")

    # Fake igv_path: must fail cleanly, not crash, regardless of whether
    # IGV happens to be installed on this machine.
    fake_result = run_igv_screenshot(
        bam_paths=["/tmp/does_not_matter.bam"],
        chromosome="chr1", start=1000, end=2000,
        output_path="/tmp/test_igv_screenshot_fake.png",
        igv_path="/nonexistent/path/igv.sh",
    )
    assert "error" in fake_result, "Expected an error dict for a nonexistent igv_path"
    assert "IGV not found" in fake_result["error"], \
        f"Expected 'IGV not found' in error, got: {fake_result['error']}"
    print(f"  Fake igv_path correctly reported: {fake_result['error']}")

    # Invalid color_by must be rejected immediately (well under a second),
    # not discovered only after IGV hangs for the full timeout. This must
    # hold regardless of whether IGV/java are actually installed, since the
    # validation happens before IGV is ever launched.
    t0 = time.time()
    bad_color_result = run_igv_screenshot(
        bam_paths=["/tmp/does_not_matter.bam"],
        chromosome="chr1", start=1000, end=2000,
        output_path="/tmp/test_igv_bad_color.png",
        igv_path="/nonexistent/path/igv.sh",
        color_by="MATE_CHROMOSOME",  # not a valid ColorOption in this IGV build
    )
    elapsed = time.time() - t0
    assert "error" in bad_color_result, "Expected an error dict for an invalid color_by"
    assert bad_color_result.get("error_type") == "invalid_parameters", \
        f"Expected error_type='invalid_parameters', got {bad_color_result.get('error_type')}"
    assert "UNEXPECTED_PAIR" in bad_color_result.get("valid_options", []), \
        "Expected the valid_options list to include UNEXPECTED_PAIR"
    assert elapsed < 5, f"color_by validation should fail fast, took {elapsed:.1f}s"
    print(f"  Invalid color_by 'MATE_CHROMOSOME' rejected in {elapsed:.2f}s "
          f"with valid_options listed")

    # Real IGV + java, if this machine has both — exercises an actual
    # screenshot run and asserts on success, not just batch-script content.
    # TEST 10 previously asserted only on the batch script's text, so it
    # kept reporting PASSED while every real run was silently failing
    # (color_by="MATE_CHROMOSOME" crashed IGV's batch executor). That is
    # worse than no test at all, so this now asserts on the real outcome.
    igv_candidates = [
        os.path.expanduser("~/IGV_2.17.4/igv.sh"),
        os.path.expanduser("~/igv/igv.sh"),
        "/opt/igv/igv.sh",
    ]
    real_igv_path = next((c for c in igv_candidates if os.path.exists(c)), None)
    java_available = shutil.which("java") is not None

    if real_igv_path is None or not java_available:
        reason = "IGV not installed" if real_igv_path is None else \
            "java not on PATH (IGV requires it — e.g. present in the 'rda' conda env)"
        print(f"  {reason} on this machine — skipping real-screenshot check.")
        print("  PASSED ✓ (error-path only)\n")
    else:
        with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
            tmp_path3 = f.name
        try:
            synthetic_bam = create_synthetic_bam(tmp_path3)
            snapshot_out = "/tmp/test_igv_screenshot_synthetic.png"
            real_result = run_igv_screenshot(
                bam_paths=[synthetic_bam],
                chromosome="chr1", start=1400, end=1600,
                output_path=snapshot_out,
                igv_path=real_igv_path,
                color_by="UNEXPECTED_PAIR",
                timeout_sec=120,
            )
            assert real_result.get("success") is True, \
                f"Expected success=True, got: {real_result}"
            assert os.path.exists(snapshot_out) and os.path.getsize(snapshot_out) > 0, \
                "Expected a non-empty PNG at the output path"
            batch = real_result["batch_script"]
            assert "goto" in batch, "Expected 'goto' command in batch script"
            assert "snapshot" in batch, "Expected 'snapshot' command in batch script"
            assert "UNEXPECTED_PAIR" in batch, "Expected coloring mode in batch script"
            print(f"  Real IGV run: success=True, "
                  f"file_size_bytes={real_result.get('file_size_bytes')}, "
                  f"shutdown_method={real_result.get('shutdown_method')}")
            print("  PASSED ✓\n")
        finally:
            os.unlink(tmp_path3)
            if os.path.exists(synthetic_bam):
                os.unlink(synthetic_bam)
            if os.path.exists(synthetic_bam + ".bai"):
                os.unlink(synthetic_bam + ".bai")
            if os.path.exists(snapshot_out):
                os.unlink(snapshot_out)

    # ── TEST 11: structured error handling for invalid inputs ────────────
    print("TEST 11: structured error handling for invalid inputs")

    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
        tmp_path4 = f.name
    try:
        err_bam = create_synthetic_bam(tmp_path4)  # header: chr1, chr8

        NONEXISTENT = "/nonexistent/path/does_not_exist.bam"

        def assert_error_dict(name, result, expected_error_type=None):
            assert isinstance(result, dict), f"{name}: expected a dict, got {type(result)}"
            assert "error" in result, f"{name}: expected an 'error' key, got keys {list(result.keys())}"
            if expected_error_type is not None:
                assert result.get("error_type") == expected_error_type, (
                    f"{name}: expected error_type={expected_error_type!r}, "
                    f"got {result.get('error_type')!r}"
                )

        # -- Nonexistent BAM path: every BAM-touching tool must return a
        #    clean error dict, never raise. --
        assert_error_dict(
            "get_bam_stats_at_locus/nonexistent",
            get_bam_stats_at_locus(NONEXISTENT, "chr1", 100, 200),
            "bam_access",
        )
        assert_error_dict(
            "count_discordant_pairs/nonexistent",
            count_discordant_pairs(NONEXISTENT, "chr1", 150),
            "bam_access",
        )
        assert_error_dict(
            "count_soft_clipped_reads/nonexistent",
            count_soft_clipped_reads(NONEXISTENT, "chr1", 150),
            "bam_access",
        )
        assert_error_dict(
            "get_split_reads/nonexistent",
            get_split_reads(NONEXISTENT, "chr1", 150),
            "bam_access",
        )
        assert_error_dict(
            "get_read_depth_profile/nonexistent",
            get_read_depth_profile(NONEXISTENT, "chr1", 100, 200),
            "bam_access",
        )
        assert_error_dict(
            "summarize_breakpoint_evidence/nonexistent",
            summarize_breakpoint_evidence(NONEXISTENT, "chr1", 150),
            "bam_access",
        )
        # check_reciprocal_breakpoint must propagate the failure, not mask
        # it behind a fabricated "INSUFFICIENT EVIDENCE" verdict.
        recip_err = check_reciprocal_breakpoint(NONEXISTENT, "chr1", 150, "chr8", 200)
        assert_error_dict("check_reciprocal_breakpoint/nonexistent", recip_err, "bam_access")
        assert "verdict" not in recip_err, \
            "A tool failure must not be reported alongside a verdict field"
        print("  Nonexistent BAM path: all tools return a clean error dict — PASSED")

        # -- Chromosome not in header --
        assert_error_dict(
            "get_bam_stats_at_locus/badchrom",
            get_bam_stats_at_locus(err_bam, "chrZZ", 100, 200),
            "invalid_region",
        )
        assert_error_dict(
            "count_discordant_pairs/badchrom",
            count_discordant_pairs(err_bam, "chrZZ", 150),
            "invalid_region",
        )
        assert_error_dict(
            "count_soft_clipped_reads/badchrom",
            count_soft_clipped_reads(err_bam, "chrZZ", 150),
            "invalid_region",
        )
        assert_error_dict(
            "get_split_reads/badchrom",
            get_split_reads(err_bam, "chrZZ", 150),
            "invalid_region",
        )
        assert_error_dict(
            "get_read_depth_profile/badchrom",
            get_read_depth_profile(err_bam, "chrZZ", 100, 200),
            "invalid_region",
        )
        assert_error_dict(
            "summarize_breakpoint_evidence/badchrom",
            summarize_breakpoint_evidence(err_bam, "chrZZ", 150),
            "invalid_region",
        )
        assert_error_dict(
            "check_reciprocal_breakpoint/badchrom",
            check_reciprocal_breakpoint(err_bam, "chrZZ", 150, "chr8", 200),
            "invalid_region",
        )
        print("  Chromosome not in header (chrZZ): all tools return a clean error dict — PASSED")

        # -- Contig-name normalisation: BAM header uses "chr1"/"chr8"; the
        #    unprefixed form must still resolve rather than error. --
        stats_noprefix = get_bam_stats_at_locus(err_bam, "1", 1000, 1100)
        assert "error" not in stats_noprefix, \
            f"Expected '1' to resolve to 'chr1' via alternate-form lookup, got {stats_noprefix}"
        print("  Contig normalisation ('1' resolves to 'chr1' in header) — PASSED")

        # -- Position beyond the end of the chromosome --
        assert_error_dict(
            "get_bam_stats_at_locus/beyond_end",
            get_bam_stats_at_locus(err_bam, "chr1", 999999999999, 999999999999 + 100),
            # Both of these were "invalid_region", which conflated two
            # different mistakes. A position past the contig end is now
            # out_of_range and names the contig length; an unknown contig
            # name is still invalid_region. See
            # REAL_PATIENT_DATA_VALIDATION.md finding 14.
            "out_of_range",
        )
        assert_error_dict(
            "count_discordant_pairs/beyond_end",
            count_discordant_pairs(err_bam, "chr1", 999999999999),
            "out_of_range",
        )
        print("  Position beyond end of chromosome: clean error dict — PASSED")

        # -- Negative position / negative start --
        assert_error_dict(
            "get_bam_stats_at_locus/negative",
            get_bam_stats_at_locus(err_bam, "chr1", -100, -50),
            "invalid_parameters",
        )
        print("  Negative start/end: clean error dict — PASSED")

        # -- start > end --
        assert_error_dict(
            "get_bam_stats_at_locus/start_gt_end",
            get_bam_stats_at_locus(err_bam, "chr1", 1000, 500),
            "invalid_parameters",
        )
        assert_error_dict(
            "get_read_depth_profile/start_gt_end",
            get_read_depth_profile(err_bam, "chr1", 1000, 500),
            "invalid_parameters",
        )
        print("  start > end: clean error dict — PASSED")

        # -- Empty bam_paths list for igv_screenshot: must fail fast, not
        #    launch IGV with zero load lines and wait out the full timeout.
        #    igv_path deliberately omitted (not pointed at a fake path) so
        #    this exercises the real "IGV is installed, bam_paths is just
        #    empty" scenario, not the separate "IGV not found" fast-path. --
        t0 = time.time()
        empty_bam_result = run_igv_screenshot(
            bam_paths=[], chromosome="chr1", start=100, end=200,
            output_path="/tmp/test_igv_empty_bams.png",
            timeout_sec=15,
        )
        empty_bam_elapsed = time.time() - t0
        assert isinstance(empty_bam_result, dict) and "error" in empty_bam_result, \
            "Expected a clean error dict for empty bam_paths"
        assert empty_bam_result.get("error_type") == "invalid_parameters", \
            f"Expected error_type='invalid_parameters', got {empty_bam_result.get('error_type')}"
        assert empty_bam_elapsed < 5, \
            f"Empty bam_paths should fail fast, took {empty_bam_elapsed:.1f}s " \
            f"(regression: this used to launch IGV with nothing loaded and " \
            f"wait out the full timeout)"
        print(f"  Empty bam_paths list for igv_screenshot: clean error dict "
              f"in {empty_bam_elapsed:.2f}s — PASSED")

        print("  PASSED ✓\n")
    finally:
        os.unlink(tmp_path4)
        if os.path.exists(err_bam):
            os.unlink(err_bam)
        if os.path.exists(err_bam + ".bai"):
            os.unlink(err_bam + ".bai")

    # ── TEST 12: detect_applicable_layers + applicable-layer normalisation ─
    print("TEST 12: detect_applicable_layers + applicable-layer normalisation")

    # -- Paired, no SA tags (create_synthetic_bam): discordant_pairs
    #    applicable, split_reads not — the same pattern found on the real
    #    HCC1143 2018-pipeline BAM in this repo. --
    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
        tmp_path5 = f.name
    try:
        paired_no_sa_bam = create_synthetic_bam(tmp_path5)

        detected = detect_applicable_layers(paired_no_sa_bam)
        assert "error" not in detected, f"Unexpected error: {detected}"
        assert detected["reads_sampled"] > 0, "Expected at least one sampled read"
        assert set(detected["applicable_layers"]) == {
            "discordant_pairs", "soft_clipped_reads", "read_depth"
        }, f"Expected 3 applicable layers (no SA tags in this fixture), got {detected['applicable_layers']}"
        print(f"  Paired/no-SA-tag BAM: applicable_layers={detected['applicable_layers']} — PASSED")

        # Unrecognised layer name must return a structured error, not raise.
        bad_layers_result = summarize_breakpoint_evidence(
            paired_no_sa_bam, "chr1", 1500, applicable_layers=["not_a_real_layer"]
        )
        assert "error" in bad_layers_result and bad_layers_result.get("error_type") == "invalid_parameters", \
            f"Expected invalid_parameters error, got {bad_layers_result}"
        print("  Unknown applicable_layers name rejected with structured error — PASSED")

        # Empty list must also be rejected (ambiguous: 0 applicable layers
        # would divide by zero if silently accepted).
        empty_layers_result = summarize_breakpoint_evidence(
            paired_no_sa_bam, "chr1", 1500, applicable_layers=[]
        )
        assert "error" in empty_layers_result and empty_layers_result.get("error_type") == "invalid_parameters", \
            f"Expected invalid_parameters error for empty applicable_layers, got {empty_layers_result}"
        print("  Empty applicable_layers list rejected with structured error — PASSED")
    finally:
        os.unlink(tmp_path5)
        if os.path.exists(paired_no_sa_bam):
            os.unlink(paired_no_sa_bam)
        if os.path.exists(paired_no_sa_bam + ".bai"):
            os.unlink(paired_no_sa_bam + ".bai")

    # -- Paired + SA tags (create_translocation_bam): all 4 layers
    #    applicable, and evidence_score should equal evidence_score_raw
    #    when applicable_layers covers all 4 (backward-compatible default). --
    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
        tmp_path6 = f.name
    try:
        all_layers_bam = create_translocation_bam(tmp_path6)

        detected_all = detect_applicable_layers(all_layers_bam)
        assert set(detected_all["applicable_layers"]) == set(EVIDENCE_LAYER_NAMES), \
            f"Expected all 4 layers applicable, got {detected_all['applicable_layers']}"

        default_summary = summarize_breakpoint_evidence(all_layers_bam, "chr1", 1050000)
        all4_summary = summarize_breakpoint_evidence(
            all_layers_bam, "chr1", 1050000, applicable_layers=detected_all["applicable_layers"]
        )
        assert default_summary["evidence_score"] == default_summary["evidence_score_raw"], \
            "Expected evidence_score == evidence_score_raw when applicable_layers defaults to all 4"
        assert all4_summary["evidence_score"] == default_summary["evidence_score"], \
            "Explicitly passing all 4 applicable layers should match the default"
        print(f"  All-4-applicable BAM: evidence_score == evidence_score_raw == "
              f"{default_summary['evidence_score']} — PASSED")

        # -- Restricting to a subset (e.g. as if split_reads/discordant_pairs
        #    were inapplicable) must raise evidence_score relative to
        #    evidence_score_raw when the applicable layers scored well,
        #    demonstrating the normalisation actually changes the number. --
        partial_summary = summarize_breakpoint_evidence(
            all_layers_bam, "chr1", 1050000,
            applicable_layers=["soft_clipped_reads", "read_depth"]
        )
        assert partial_summary["applicable_layers"] == ["soft_clipped_reads", "read_depth"]
        assert partial_summary["signal_layers"].endswith("/2"), \
            f"Expected signal_layers denominator of 2, got {partial_summary['signal_layers']}"
        expected_partial_score = round(
            (partial_summary["soft_clip_score"] + partial_summary["depth_score"]) * 100.0 / 50.0, 1
        )
        assert partial_summary["evidence_score"] == expected_partial_score, (
            f"Expected evidence_score={expected_partial_score} when normalised over "
            f"2 applicable layers, got {partial_summary['evidence_score']}"
        )
        assert partial_summary["evidence_score_raw"] == default_summary["evidence_score_raw"], \
            "evidence_score_raw must stay the same regardless of applicable_layers"
        print(f"  2-applicable-layer normalisation: evidence_score="
              f"{partial_summary['evidence_score']} vs evidence_score_raw="
              f"{partial_summary['evidence_score_raw']} — PASSED")
        print("  PASSED ✓\n")
    finally:
        os.unlink(tmp_path6)
        if os.path.exists(all_layers_bam):
            os.unlink(all_layers_bam)
        if os.path.exists(all_layers_bam + ".bai"):
            os.unlink(all_layers_bam + ".bai")

    # ── TEST 13: igv_evidence_panel ────────────────────────────────────────
    print("TEST 13: igv_evidence_panel")

    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
        tmp_path7 = f.name
    panel_dir = tempfile.mkdtemp(prefix="igv_evidence_panel_test_")
    try:
        panel_bam = create_synthetic_bam(tmp_path7)

        # Unknown applicable_layers value must return a structured error,
        # not raise, and must not attempt to launch IGV at all.
        bad_layers_panel = igv_evidence_panel(
            bam_paths=[panel_bam], chromosome="chr1", position=1500,
            output_dir=panel_dir, applicable_layers=["not_a_real_layer"],
        )
        assert "error" in bad_layers_panel and bad_layers_panel.get("error_type") == "invalid_parameters", \
            f"Expected invalid_parameters error, got {bad_layers_panel}"
        print("  Unknown applicable_layers value rejected with structured error — PASSED")

        # Unknown windows key must also return a structured error.
        bad_windows_panel = igv_evidence_panel(
            bam_paths=[panel_bam], chromosome="chr1", position=1500,
            output_dir=panel_dir, windows={"not_a_real_layer": 300},
        )
        assert "error" in bad_windows_panel and bad_windows_panel.get("error_type") == "invalid_parameters", \
            f"Expected invalid_parameters error for bad windows key, got {bad_windows_panel}"
        print("  Unknown windows key rejected with structured error — PASSED")

        # A missing IGV binary must produce a clean per-layer error, not a
        # crash, and must not block on layers that ARE inapplicable.
        fake_igv_panel = igv_evidence_panel(
            bam_paths=[panel_bam], chromosome="chr1", position=1500,
            output_dir=panel_dir, applicable_layers=["read_depth"],
            igv_path="/nonexistent/path/igv.sh",
        )
        assert "error" not in fake_igv_panel, f"Unexpected top-level error: {fake_igv_panel}"
        assert fake_igv_panel["panels"]["read_depth"].get("error") == "IGV not found", \
            f"Expected 'IGV not found' for the read_depth panel, got {fake_igv_panel['panels']['read_depth']}"
        for skipped_layer in ("discordant_pairs", "soft_clipped_reads", "split_reads"):
            assert fake_igv_panel["panels"][skipped_layer].get("skipped") is True, \
                f"Expected {skipped_layer} to be skipped (not in applicable_layers)"
        print("  Missing IGV binary: clean per-layer error, other layers correctly skipped — PASSED")

        # windows_used must reflect the caller-supplied start/end for
        # read_depth specifically, and the default half-window otherwise —
        # computed without needing to launch IGV.
        window_check_panel = igv_evidence_panel(
            bam_paths=[panel_bam], chromosome="chr1", position=1500,
            output_dir=panel_dir, start=100, end=3000,
            applicable_layers=["read_depth", "soft_clipped_reads"],
            windows={"soft_clipped_reads": 75},
            igv_path="/nonexistent/path/igv.sh",
        )
        wu = window_check_panel["windows_used"]
        assert wu["read_depth"] == {
            "start": 100, "end": 3000, "window_bp": None,
            "source": "caller-supplied start/end",
        }, f"Expected read_depth to use caller-supplied start/end, got {wu['read_depth']}"
        assert wu["soft_clipped_reads"] == {
            "start": 1425, "end": 1575, "window_bp": 75, "source": "override",
        }, f"Expected soft_clipped_reads to use the 75bp override, got {wu['soft_clipped_reads']}"
        assert wu["discordant_pairs"]["window_bp"] == DEFAULT_PANEL_WINDOWS["discordant_pairs"], \
            f"Expected discordant_pairs to use the default window, got {wu['discordant_pairs']}"
        print("  windows_used correctly reflects per-layer overrides and read_depth start/end — PASSED")

        # Real IGV, if this machine has it (and java) — exercises an actual
        # 4-panel run end-to-end on a BAM where all 4 layers are applicable.
        igv_candidates = [
            os.path.expanduser("~/IGV_2.17.4/igv.sh"),
            os.path.expanduser("~/igv/igv.sh"),
            "/opt/igv/igv.sh",
        ]
        real_igv_path = next((c for c in igv_candidates if os.path.exists(c)), None)
        java_available = shutil.which("java") is not None

        if real_igv_path is None or not java_available:
            reason = "IGV not installed" if real_igv_path is None else \
                "java not on PATH (IGV requires it — e.g. present in the 'rda' conda env)"
            print(f"  {reason} on this machine — skipping real 4-panel generation check.")
            print("  PASSED ✓ (error-path only)\n")
        else:
            with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
                tmp_path8 = f.name
            real_panel_dir = tempfile.mkdtemp(prefix="igv_evidence_panel_real_")
            try:
                real_panel_bam = create_translocation_bam(tmp_path8)
                real_panel = igv_evidence_panel(
                    bam_paths=[real_panel_bam], chromosome="chr1",
                    position=1050000, output_dir=real_panel_dir,
                    igv_path=real_igv_path, timeout_sec=120,
                )
                assert "error" not in real_panel, f"Unexpected top-level error: {real_panel}"
                assert set(real_panel["applicable_layers"]) == set(EVIDENCE_LAYER_NAMES), \
                    f"Expected all 4 layers applicable on this fixture, got {real_panel['applicable_layers']}"
                for layer in EVIDENCE_LAYER_NAMES:
                    panel = real_panel["panels"][layer]
                    assert panel.get("skipped") is not True, f"{layer} unexpectedly skipped: {panel}"
                    assert panel.get("success") is True, f"{layer} did not succeed: {panel}"
                    path = panel.get("screenshot_path")
                    assert path and os.path.exists(path) and os.path.getsize(path) > 0, \
                        f"{layer}: expected a non-empty PNG at {path}"
                    # Each layer must use its OWN window, not one shared region.
                    expected_half = DEFAULT_PANEL_WINDOWS[layer]
                    wu = real_panel["windows_used"][layer]
                    assert wu["window_bp"] == expected_half, \
                        f"{layer}: expected window_bp={expected_half}, got {wu}"
                    assert wu["end"] - wu["start"] == 2 * expected_half, \
                        f"{layer}: window span doesn't match window_bp: {wu}"
                    assert panel["window_bp"] == expected_half, \
                        f"{layer}: panel's own window_bp doesn't match windows_used: {panel.get('window_bp')}"
                # soft_clipped_reads must record how its sort side was chosen.
                clip_note = real_panel["panels"]["soft_clipped_reads"]["clip_side_determination"]
                assert clip_note and ("LEFT_CLIP" in clip_note or "RIGHT_CLIP" in clip_note or "defaulted" in clip_note), \
                    f"Expected a clip_side_determination note, got {clip_note!r}"
                print(f"  Real 4-layer panel: all layers succeeded, each with its own window — PASSED")
                print(f"  clip_side_determination: {clip_note}")
                print("  PASSED ✓\n")
            finally:
                os.unlink(tmp_path8)
                if os.path.exists(real_panel_bam):
                    os.unlink(real_panel_bam)
                if os.path.exists(real_panel_bam + ".bai"):
                    os.unlink(real_panel_bam + ".bai")
                shutil.rmtree(real_panel_dir, ignore_errors=True)
    finally:
        os.unlink(tmp_path7)
        if os.path.exists(panel_bam):
            os.unlink(panel_bam)
        if os.path.exists(panel_bam + ".bai"):
            os.unlink(panel_bam + ".bai")
        shutil.rmtree(panel_dir, ignore_errors=True)

    # ── TEST 14: summarize_breakpoint_evidence at a zero-read locus ──────
    # Regression test for the zero-coverage false-positive bug: a locus
    # with no reads at all must report NOT ASSESSABLE, never a fabricated
    # score (0 or otherwise) on any of the 4 components.
    print("TEST 14: summarize_breakpoint_evidence at a zero-read locus — must be NOT ASSESSABLE")
    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
        tmp_path9 = f.name
    try:
        zero_reads_bam = create_translocation_bam(tmp_path9)
        # chr1:5,000,000 is far outside every read group this fixture
        # places (all on chr1:1,000,000-1,051,320ish or chr8:47,000,000+),
        # so every layer's own window (±500bp default, ±2kb for depth)
        # is genuinely empty — not technologically inapplicable, just
        # nothing sequenced there.
        empty_locus_summary = summarize_breakpoint_evidence(zero_reads_bam, "chr1", 5000000)
        assert "error" not in empty_locus_summary, f"Unexpected error: {empty_locus_summary}"
        assert empty_locus_summary["evidence_strength"] == "NOT ASSESSABLE", \
            f"Expected NOT ASSESSABLE at a zero-read locus, got {empty_locus_summary['evidence_strength']!r}"
        assert empty_locus_summary["evidence_score"] is None, \
            f"Expected evidence_score=None (not 0) when nothing is assessable, got {empty_locus_summary['evidence_score']}"
        assert empty_locus_summary["signal_layers"] == "0/0", \
            f"Expected signal_layers='0/0', got {empty_locus_summary['signal_layers']!r}"
        assert set(empty_locus_summary["unassessable_layers"].keys()) == set(EVIDENCE_LAYER_NAMES), \
            f"Expected all 4 layers unassessable, got {empty_locus_summary['unassessable_layers']}"
        for component_field in ("discordant_pair_score", "soft_clip_score",
                                 "split_read_score", "depth_score"):
            assert empty_locus_summary[component_field] is None, \
                f"Expected {component_field}=None (never scored) at a zero-read locus, " \
                f"got {empty_locus_summary[component_field]}"
        # evidence_score_raw stays a plain number (unassessable layers
        # contribute 0), so it's still safely inspectable even though the
        # primary evidence_score field is None.
        assert empty_locus_summary["evidence_score_raw"] == 0.0, \
            f"Expected evidence_score_raw=0.0, got {empty_locus_summary['evidence_score_raw']}"
        print(f"  unassessable_layers: {empty_locus_summary['unassessable_layers']}")
        print(f"  evidence_score={empty_locus_summary['evidence_score']!r} "
              f"evidence_strength={empty_locus_summary['evidence_strength']!r} "
              f"evidence_score_raw={empty_locus_summary['evidence_score_raw']}")
        print("  PASSED ✓\n")
    finally:
        os.unlink(tmp_path9)
        if os.path.exists(zero_reads_bam):
            os.unlink(zero_reads_bam)
        if os.path.exists(zero_reads_bam + ".bai"):
            os.unlink(zero_reads_bam + ".bai")

    # ── TEST 15: check_reciprocal_breakpoint — naming-convention independence ─
    # Regression test for the contig-normalisation bug: mate_chromosomes'
    # keys come straight from the BAM header's own convention, so looking
    # them up with a caller-supplied chromosome string in the OPPOSITE
    # convention (bare digit vs chr-prefixed) must not silently zero out
    # back_pointing_to_primary and downgrade the verdict.
    print("TEST 15: check_reciprocal_breakpoint — contig-naming-convention independence")
    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
        tmp_path10 = f.name
    try:
        recip_bam = create_translocation_bam(tmp_path10)  # header uses "chr1"/"chr8"

        prefixed = check_reciprocal_breakpoint(recip_bam, "chr1", 1050000, "chr8", 47000000)
        bare = check_reciprocal_breakpoint(recip_bam, "1", 1050000, "8", 47000000)

        assert "error" not in prefixed, f"Unexpected error: {prefixed}"
        assert "error" not in bare, f"Unexpected error: {bare}"

        assert bare["verdict"] == prefixed["verdict"], (
            f"Expected identical verdict regardless of chr-prefix convention, got "
            f"bare={bare['verdict']!r} vs prefixed={prefixed['verdict']!r}"
        )
        assert bare["is_balanced"] == prefixed["is_balanced"], \
            "Expected is_balanced to match regardless of naming convention"
        assert bare["primary"]["discordant_pairs"] == prefixed["primary"]["discordant_pairs"] == 15, (
            f"Expected primary discordant_pairs=15 regardless of naming convention "
            f"(bare's local-mate split/clip reads must not be miscounted as discordant "
            f"just because '1' was compared against the header's 'chr1'), got "
            f"bare={bare['primary']['discordant_pairs']} vs prefixed={prefixed['primary']['discordant_pairs']}"
        )
        assert (
            bare["reciprocal"]["back_pointing_to_primary"]
            == prefixed["reciprocal"]["back_pointing_to_primary"]
            == 15
        ), (
            f"Expected back_pointing_to_primary=15 regardless of naming convention, got "
            f"bare={bare['reciprocal']['back_pointing_to_primary']} vs "
            f"prefixed={prefixed['reciprocal']['back_pointing_to_primary']}"
        )
        print(f"  Bare ('1'/'8') and chr-prefixed ('chr1'/'chr8') calls agree exactly: "
              f"verdict={bare['verdict']!r}, "
              f"primary_disc={bare['primary']['discordant_pairs']}, "
              f"back_pointing_to_primary={bare['reciprocal']['back_pointing_to_primary']}")
        print("  PASSED ✓\n")
    finally:
        os.unlink(tmp_path10)
        if os.path.exists(recip_bam):
            os.unlink(recip_bam)
        if os.path.exists(recip_bam + ".bai"):
            os.unlink(recip_bam + ".bai")

    # ── TEST 16: get_read_depth_profile bin-size invariance (FIX 1) ───────
    # Regression test for the read-count-vs-depth bug: before FIX 1,
    # get_read_depth_profile counted DISTINCT READS touching each bin
    # rather than true per-base depth, so mean_depth scaled with
    # window_size instead of being invariant to it (confirmed on real
    # HG002 data: 976.98 at window_size=100 vs 1354.25 at window_size=200
    # for the identical chr1:16,890,000 region — a ~39% difference for
    # nothing but a binning choice). True per-base depth must not depend
    # on how it's binned.
    print("TEST 16: get_read_depth_profile bin-size invariance (FIX 1)")
    with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as f:
        tmp_path11 = f.name
    try:
        invariance_bam = create_depth_dropout_bam(tmp_path11)
        # chr1:0-1000 is the uniform flank TEST 4 already confirmed has no
        # dropout — same region, three window_size values.
        means = {}
        for ws in (50, 100, 200):
            profile = get_read_depth_profile(invariance_bam, "chr1", 0, 1000, window_size=ws)
            assert "error" not in profile, f"Unexpected error at window_size={ws}: {profile}"
            means[ws] = profile["summary"]["mean_depth"]
        print(f"  mean_depth by window_size: {means}")
        baseline = means[100]
        for ws, m in means.items():
            pct_diff = abs(m - baseline) / baseline * 100
            assert pct_diff <= 5.0, (
                f"window_size={ws} mean_depth={m} differs from window_size=100's "
                f"{baseline} by {pct_diff:.1f}% (>5%) — depth must be bin-size invariant"
            )
        print(f"  All window sizes agree with window_size=100 within 5% — PASSED")
        print("  PASSED ✓\n")
    finally:
        os.unlink(tmp_path11)
        if os.path.exists(invariance_bam):
            os.unlink(invariance_bam)
        if os.path.exists(invariance_bam + ".bai"):
            os.unlink(invariance_bam + ".bai")

    # ── TEST 17/18: real HG002 BAM — FIX 1/2/3 verification ────────────────
    # These two tests hit a real, public, remote GIAB BAM over HTTPS (same
    # network-dependency precedent as TEST 7's Ensembl calls) at the exact
    # 3 coordinates from results/LLM_SESSION_3_BLIND.md, the blind session
    # that surfaced all 3 bugs these fixes address. Position A and B are
    # background/no-signal loci; position C is the real, GIAB-confirmed
    # HG002 deletion breakpoint at chr1:115,686,862 (VANGL1) also used
    # in TEST 7 and results/GIAB_PUBLIC_DATA_VALIDATION.md.
    REAL_HG002_BAM = (
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/"
        "AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/"
        "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam"
    )
    POSITION_A = ("chr1", 16890000)   # background — no real signal at this exact position
    POSITION_B = ("chr2", 96300000)   # background — cleanest of the 3, no nearby dip either
    POSITION_C = ("chr1", 115686862)  # real GIAB HG002 deletion breakpoint (VANGL1)

    print("TEST 17: get_read_depth_profile focus_position/dip_* on real HG002 BAM (FIX 1+2)")
    try:
        for label, (chrom, pos), expect_deletion, expect_dip_at_focus in (
            ("A", POSITION_A, False, False),
            ("C", POSITION_C, True, True),
        ):
            profile = get_read_depth_profile(
                REAL_HG002_BAM, chrom, pos - 3000, pos + 3000,
                window_size=100, focus_position=pos,
            )
            assert "error" not in profile, f"Unexpected error at position {label}: {profile}"
            s = profile["summary"]
            print(f"  Position {label} ({chrom}:{pos}): ratio={s['depth_ratio_min_to_mean']} "
                  f"likely_deletion={s['likely_deletion']} dip_position={s['dip_position']} "
                  f"dip_distance_from_focus={s['dip_distance_from_focus']} "
                  f"dip_is_at_focus={s['dip_is_at_focus']}")
            assert s["likely_deletion"] is expect_deletion, (
                f"Position {label}: expected likely_deletion={expect_deletion}, got {s['likely_deletion']}"
            )
            assert s["dip_is_at_focus"] is expect_dip_at_focus, (
                f"Position {label}: expected dip_is_at_focus={expect_dip_at_focus}, got {s['dip_is_at_focus']}"
            )
            if label == "A":
                # "Roughly 1500" per the FIX 2 spec; real measured value is
                # ~1400bp — asserted as a range, not pinned to the exact
                # integer, since the point being tested is "far enough to
                # correctly stay off-focus", not the precise offset.
                assert 1000 <= s["dip_distance_from_focus"] <= 2000, (
                    f"Position A: expected dip_distance_from_focus roughly 1500, "
                    f"got {s['dip_distance_from_focus']}"
                )

        # Bin-size invariance on real data too (window_size=100 vs 200,
        # ±2kb around position A, matching summarize_breakpoint_evidence's
        # own internal window) — before FIX 1: 976.98 vs 1354.25 (~39%
        # apart) for this same region.
        chrom, pos = POSITION_A
        p100 = get_read_depth_profile(REAL_HG002_BAM, chrom, pos - 2000, pos + 2000, window_size=100)
        p200 = get_read_depth_profile(REAL_HG002_BAM, chrom, pos - 2000, pos + 2000, window_size=200)
        m100, m200 = p100["summary"]["mean_depth"], p200["summary"]["mean_depth"]
        pct_diff = abs(m100 - m200) / m100 * 100
        print(f"  Real-data bin-size invariance at position A: window_size=100 mean_depth={m100}, "
              f"window_size=200 mean_depth={m200} ({pct_diff:.2f}% apart)")
        assert pct_diff <= 5.0, (
            f"Real-data mean_depth differs by {pct_diff:.1f}% between window sizes (>5%)"
        )
        # ...and consistent with bam_stats_at_locus's independently-computed
        # true per-base depth for a comparable span (500bp centered on A).
        stats_a = get_bam_stats_at_locus(REAL_HG002_BAM, chrom, pos - 250, pos + 250)
        assert "error" not in stats_a, f"Unexpected error: {stats_a}"
        pct_vs_stats = abs(m100 - stats_a["mean_depth"]) / stats_a["mean_depth"] * 100
        print(f"  vs bam_stats_at_locus mean_depth={stats_a['mean_depth']} ({pct_vs_stats:.2f}% apart)")
        assert pct_vs_stats <= 5.0, (
            f"read_depth_profile disagrees with bam_stats_at_locus by {pct_vs_stats:.1f}% (>5%)"
        )
        print("  PASSED ✓\n")
    except Exception as e:
        if "getaddrinfo" in str(e) or "Connection" in str(e) or "Network" in str(e):
            print(f"  SKIPPED — no network access to real GIAB BAM ({e})\n")
        else:
            raise

    print("TEST 18: soft-clip max_clips_at_position scoring on real HG002 BAM (FIX 3)")
    try:
        results = {}
        for label, (chrom, pos) in (("A", POSITION_A), ("B", POSITION_B), ("C", POSITION_C)):
            clips = count_soft_clipped_reads(REAL_HG002_BAM, chrom, pos)
            assert "error" not in clips, f"Unexpected error at position {label}: {clips}"
            results[label] = clips
            print(f"  Position {label}: soft_clipped_fraction={clips['soft_clipped_fraction']} "
                  f"max_clips_at_position={clips['max_clips_at_position']}")

        # The whole point of FIX 3: A and B have LOW max_clips_at_position
        # (no real pileup) despite nonzero fraction; C has a genuine pileup.
        assert results["A"]["max_clips_at_position"] < 3, \
            f"Expected position A max_clips_at_position < 3 (no pileup), got {results['A']['max_clips_at_position']}"
        assert results["B"]["max_clips_at_position"] < 3, \
            f"Expected position B max_clips_at_position < 3 (no pileup), got {results['B']['max_clips_at_position']}"
        assert results["C"]["max_clips_at_position"] >= 10, \
            f"Expected position C max_clips_at_position >= 10 (genuine pileup), got {results['C']['max_clips_at_position']}"

        # And the composite score reflects it: A/B score 0 on this layer, C scores full (25).
        applicable = detect_applicable_layers(REAL_HG002_BAM)["applicable_layers"]
        for label, (chrom, pos), expected_score in (
            ("A", POSITION_A, 0.0), ("B", POSITION_B, 0.0), ("C", POSITION_C, 25.0)
        ):
            summary = summarize_breakpoint_evidence(
                REAL_HG002_BAM, chrom, pos, label=f"position_{label}",
                applicable_layers=applicable,
            )
            assert "error" not in summary, f"Unexpected error at position {label}: {summary}"
            print(f"  Position {label}: soft_clip_score={summary['soft_clip_score']} "
                  f"depth_score={summary['depth_score']} evidence_score={summary['evidence_score']} "
                  f"strength={summary['evidence_strength']}")
            assert summary["soft_clip_score"] == expected_score, (
                f"Position {label}: expected soft_clip_score={expected_score}, "
                f"got {summary['soft_clip_score']}"
            )
        # C should also now score on depth (FIX 2's dip_is_at_focus gate
        # passes here — the real deletion's dip IS near the focus).
        summary_c = summarize_breakpoint_evidence(
            REAL_HG002_BAM, *POSITION_C, label="position_C", applicable_layers=applicable
        )
        assert summary_c["depth_score"] > 0, \
            f"Expected position C depth_score > 0 (real, on-position deletion signal), got {summary_c['depth_score']}"
        print("  PASSED ✓\n")
    except Exception as e:
        if "getaddrinfo" in str(e) or "Connection" in str(e) or "Network" in str(e):
            print(f"  SKIPPED — no network access to real GIAB BAM ({e})\n")
        else:
            raise

    print("=" * 60)
    print("ALL TESTS PASSED")
    print("=" * 60)


if __name__ == "__main__":
    run_tests()