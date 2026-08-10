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
    EVIDENCE_LAYER_NAMES,
    DEPTH_RATIO_DELETION_THRESHOLD,
)
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
        # signal_layers is 4/4 here rather than 3/4: the read-free gaps deliberately
        # built into this synthetic BAM (to keep discordant/clip/split fractions clean)
        # also trip the depth-drop threshold — an artifact of test construction, not a
        # real deletion signal. See TEST 6 output above: min depth is 0 reads/window.
        assert summary["signal_layers"] == "4/4", "Expected all 4 evidence layers to show signal"
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

    # ── TEST 7: get_gene_at_locus ────────────────────────────────────────
    print("TEST 7: get_gene_at_locus")
    gene_result = None
    for attempt in range(4):
        gene_result = get_gene_at_locus("chr1", 115686862)
        if "error" not in gene_result:
            break
        print(f"  (Ensembl attempt {attempt + 1} failed: {gene_result['error']} — retrying)")
    assert gene_result is not None and "error" not in gene_result, \
        "Ensembl REST API did not respond after retries"
    print(f"  Gene(s) at chr1:115686862: {[g['gene_name'] for g in gene_result['genes']]}")
    assert gene_result["gene_count"] >= 1, "Expected at least 1 gene"
    assert gene_result["is_intergenic"] is False, "Expected is_intergenic == False"

    intergenic_result = None
    for attempt in range(4):
        intergenic_result = get_gene_at_locus("chr13", 15000000)
        if "error" not in intergenic_result:
            break
        print(f"  (Ensembl attempt {attempt + 1} failed: {intergenic_result['error']} — retrying)")
    assert intergenic_result is not None and "error" not in intergenic_result, \
        "Ensembl REST API did not respond after retries"
    print(f"  chr13:15000000 is_intergenic: {intergenic_result['is_intergenic']}")
    assert intergenic_result["is_intergenic"] is True, "Expected is_intergenic == True"
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
            "invalid_region",
        )
        assert_error_dict(
            "count_discordant_pairs/beyond_end",
            count_discordant_pairs(err_bam, "chr1", 999999999999),
            "invalid_region",
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

        # -- Empty bam_paths list for igv_screenshot: must not raise --
        empty_bam_result = run_igv_screenshot(
            bam_paths=[], chromosome="chr1", start=100, end=200,
            output_path="/tmp/test_igv_empty_bams.png",
            igv_path="/nonexistent/path/igv.sh",
        )
        assert isinstance(empty_bam_result, dict) and "error" in empty_bam_result, \
            "Expected a clean error dict for empty bam_paths"
        print("  Empty bam_paths list for igv_screenshot: clean error dict — PASSED")

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

    print("=" * 60)
    print("ALL TESTS PASSED")
    print("=" * 60)


if __name__ == "__main__":
    run_tests()