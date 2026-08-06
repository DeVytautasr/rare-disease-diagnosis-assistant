"""
bam_tools.py
Core BAM inspection tools for the Stage 1 IGV breakpoint assistant.

Every function returns structured data only — no genomic interpretation.
The LLM reads these outputs and writes the report. It cannot add facts
that are not in the tool output.
"""

import pysam
from dataclasses import dataclass, asdict
from typing import Optional


# ── Data structures ──────────────────────────────────────────────────────────

@dataclass
class LocusStats:
    chromosome: str
    start: int
    end: int
    total_reads: int
    mean_depth: float
    mean_mapq: float
    low_mapq_fraction: float   # fraction of reads with MAPQ < 20
    forward_reads: int
    reverse_reads: int

@dataclass
class DiscordantPairResult:
    chromosome: str
    position: int
    window_bp: int
    total_reads_in_window: int
    discordant_pairs: int
    discordant_fraction: float
    mate_chromosomes: dict     # {chr_name: count}

@dataclass
class SoftClipResult:
    chromosome: str
    position: int
    window_bp: int
    total_reads_in_window: int
    soft_clipped_reads: int
    soft_clipped_fraction: float
    consensus_clip_position: Optional[int]   # position with most clipping
    max_clips_at_position: int


# ── Tool 1: Basic locus quality stats ────────────────────────────────────────

def get_bam_stats_at_locus(
    bam_path: str,
    chromosome: str,
    start: int,
    end: int
) -> dict:
    """
    Returns basic quality statistics for a genomic locus.
    Call this first before any other tool to check data quality.

    Args:
        bam_path:   Path to indexed BAM or CRAM file
        chromosome: Chromosome name (e.g. 'chr1' or '1')
        start:      Start position (0-based)
        end:        End position

    Returns:
        LocusStats as dict, or error dict if BAM cannot be read
    """
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": str(e), "bam_path": bam_path}

    total = 0
    mapq_sum = 0
    low_mapq = 0
    forward = 0
    reverse = 0
    depth_positions = {}

    for read in bam.fetch(chromosome, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        total += 1
        mapq_sum += read.mapping_quality
        if read.mapping_quality < 20:
            low_mapq += 1
        if read.is_forward:
            forward += 1
        else:
            reverse += 1
        # Track depth per position
        for pos in read.get_reference_positions():
            if start <= pos < end:
                depth_positions[pos] = depth_positions.get(pos, 0) + 1

    bam.close()

    region_length = end - start
    mean_depth = sum(depth_positions.values()) / region_length if region_length > 0 else 0
    mean_mapq = mapq_sum / total if total > 0 else 0
    low_mapq_fraction = low_mapq / total if total > 0 else 0

    result = LocusStats(
        chromosome=chromosome,
        start=start,
        end=end,
        total_reads=total,
        mean_depth=round(mean_depth, 2),
        mean_mapq=round(mean_mapq, 2),
        low_mapq_fraction=round(low_mapq_fraction, 3),
        forward_reads=forward,
        reverse_reads=reverse,
    )
    return asdict(result)


# ── Tool 2: Discordant read pairs ─────────────────────────────────────────────

def count_discordant_pairs(
    bam_path: str,
    chromosome: str,
    position: int,
    window_bp: int = 500,
    min_mapq: int = 20
) -> dict:
    """
    Counts read pairs where the mate maps to a different chromosome.
    This is the primary signal for inter-chromosomal translocations.

    A high fraction of discordant pairs at a locus strongly suggests
    a translocation breakpoint. The mate_chromosomes dict shows where
    the mates map — a cluster on one partner chromosome is most convincing.

    Args:
        bam_path:   Path to indexed BAM file
        chromosome: Chromosome of the candidate breakpoint
        position:   Centre of the inspection window
        window_bp:  Half-width of the inspection window (default 500bp)
        min_mapq:   Minimum mapping quality to include read (default 20)

    Returns:
        DiscordantPairResult as dict
    """
    start = max(0, position - window_bp)
    end = position + window_bp

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": str(e)}

    total = 0
    discordant = 0
    mate_chroms = {}

    for read in bam.fetch(chromosome, start, end):
        if read.is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.mate_is_unmapped:
            continue

        total += 1

        # Discordant = mate on a different chromosome
        if read.next_reference_name != chromosome:
            discordant += 1
            mate_chr = read.next_reference_name or "unknown"
            mate_chroms[mate_chr] = mate_chroms.get(mate_chr, 0) + 1

    bam.close()

    # Sort mate chromosomes by count descending
    mate_chroms_sorted = dict(
        sorted(mate_chroms.items(), key=lambda x: x[1], reverse=True)
    )

    result = DiscordantPairResult(
        chromosome=chromosome,
        position=position,
        window_bp=window_bp,
        total_reads_in_window=total,
        discordant_pairs=discordant,
        discordant_fraction=round(discordant / total, 3) if total > 0 else 0,
        mate_chromosomes=mate_chroms_sorted,
    )
    return asdict(result)


# ── Tool 3: Soft-clipped reads ────────────────────────────────────────────────

def count_soft_clipped_reads(
    bam_path: str,
    chromosome: str,
    position: int,
    window_bp: int = 200,
    min_clip_bases: int = 10,
    min_mapq: int = 20
) -> dict:
    """
    Counts reads with significant soft-clipping near a candidate breakpoint.
    Reads that span a junction are partially aligned and partially clipped.
    A pileup of clipped reads at the same position narrows the breakpoint
    to near-nucleotide resolution.

    Args:
        bam_path:        Path to indexed BAM file
        chromosome:      Chromosome of the candidate breakpoint
        position:        Centre of the inspection window
        window_bp:       Half-width of the window (default 200bp)
        min_clip_bases:  Minimum bases clipped to count (default 10)
        min_mapq:        Minimum mapping quality (default 20)

    Returns:
        SoftClipResult as dict
    """
    start = max(0, position - window_bp)
    end = position + window_bp

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": str(e)}

    total = 0
    clipped = 0
    clip_positions = {}

    for read in bam.fetch(chromosome, start, end):
        if read.is_unmapped or read.is_secondary:
            continue
        if read.mapping_quality < min_mapq:
            continue

        total += 1
        cigar = read.cigartuples
        if not cigar:
            continue

        # Check for soft-clip at read start (left clip)
        if cigar[0][0] == 4 and cigar[0][1] >= min_clip_bases:
            clipped += 1
            clip_pos = read.reference_start
            clip_positions[clip_pos] = clip_positions.get(clip_pos, 0) + 1

        # Check for soft-clip at read end (right clip)
        elif cigar[-1][0] == 4 and cigar[-1][1] >= min_clip_bases:
            clipped += 1
            clip_pos = read.reference_end
            clip_positions[clip_pos] = clip_positions.get(clip_pos, 0) + 1

    bam.close()

    consensus_pos = None
    max_clips = 0
    if clip_positions:
        consensus_pos = max(clip_positions, key=clip_positions.get)
        max_clips = clip_positions[consensus_pos]

    result = SoftClipResult(
        chromosome=chromosome,
        position=position,
        window_bp=window_bp,
        total_reads_in_window=total,
        soft_clipped_reads=clipped,
        soft_clipped_fraction=round(clipped / total, 3) if total > 0 else 0,
        consensus_clip_position=consensus_pos,
        max_clips_at_position=max_clips,
    )
    return asdict(result)