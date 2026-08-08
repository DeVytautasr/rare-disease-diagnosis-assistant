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

@dataclass
class SplitReadResult:
    chromosome: str
    position: int
    window_bp: int
    total_reads_in_window: int
    split_reads: int
    split_read_fraction: float
    partner_chromosomes: dict     # {chr_name: count}, parsed from SA tags
    example_partner_loci: list    # up to 5 "chrom:pos" strings for inspection

@dataclass
class BreakpointEvidenceSummary:
    label: str
    chromosome: str
    position: int
    evidence_score: float          # 0-100, normalized sum of the 3 components below
    evidence_strength: str         # "none" | "weak" | "moderate" | "strong"
    signal_layers: str             # e.g. "3/3" — how many of the 3 layers show any signal
    discordant_pair_score: float   # 0-50
    soft_clip_score: float         # 0-50
    split_read_score: float        # 0-50
    locus_stats: dict
    discordant_pairs: dict
    soft_clips: dict
    split_reads: dict
    supporting_observations: list


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
        # Unpaired reads (e.g. single-molecule long reads) have no mate at
        # all, so "mate_is_unmapped" defaults to False for them — without
        # this check they'd all be miscounted as discordant.
        if not read.is_paired:
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


# ── Tool 4: Split (chimeric) reads ─────────────────────────────────────────────

def get_split_reads(
    bam_path: str,
    chromosome: str,
    position: int,
    window_bp: int = 200,
    min_mapq: int = 0
) -> dict:
    """
    Counts split (chimeric) reads: primary alignments carrying an SA
    (supplementary alignment) tag, meaning part of the read aligns elsewhere
    in the genome. Unlike soft-clip counts, the SA tag gives the exact
    chromosome and position of the other segment, so split reads can pinpoint
    a breakpoint partner directly rather than just flagging that one exists.

    Args:
        bam_path:   Path to indexed BAM file
        chromosome: Chromosome of the candidate breakpoint
        position:   Centre of the inspection window
        window_bp:  Half-width of the inspection window (default 200bp)
        min_mapq:   Minimum mapping quality of the primary read (default 0)

    Returns:
        SplitReadResult as dict
    """
    start = max(0, position - window_bp)
    end = position + window_bp

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": str(e)}

    total = 0
    split = 0
    partner_chroms = {}
    example_loci = []

    for read in bam.fetch(chromosome, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue

        total += 1

        if not read.has_tag("SA"):
            continue

        split += 1
        # SA tag format: "rname,pos,strand,CIGAR,mapQ,NM;" (one or more entries)
        sa_entries = read.get_tag("SA").rstrip(";").split(";")
        first_rname, first_pos = None, None
        for entry in sa_entries:
            fields = entry.split(",")
            if len(fields) < 2:
                continue
            rname, pos_str = fields[0], fields[1]
            if first_rname is None:
                first_rname, first_pos = rname, pos_str
            partner_chroms[rname] = partner_chroms.get(rname, 0) + 1

        if first_rname is not None:
            locus = f"{first_rname}:{first_pos}"
            if locus not in example_loci and len(example_loci) < 5:
                example_loci.append(locus)

    bam.close()

    partner_chroms_sorted = dict(
        sorted(partner_chroms.items(), key=lambda x: x[1], reverse=True)
    )

    result = SplitReadResult(
        chromosome=chromosome,
        position=position,
        window_bp=window_bp,
        total_reads_in_window=total,
        split_reads=split,
        split_read_fraction=round(split / total, 3) if total > 0 else 0,
        partner_chromosomes=partner_chroms_sorted,
        example_partner_loci=example_loci,
    )
    return asdict(result)


# ── Tool 5: Combined breakpoint evidence summary ──────────────────────────────

def summarize_breakpoint_evidence(
    bam_path: str,
    chromosome: str,
    position: int,
    label: str = "",
    window_bp: int = 500,
    min_mapq: int = 20
) -> dict:
    """
    Combines discordant-pair, soft-clip, and split-read evidence into a single
    interpretable breakpoint evidence summary. Each of the 3 evidence layers is
    scored independently (0-50), summed, and normalized to a 0-100
    evidence_score, so the contribution of each layer stays visible rather
    than collapsed into a black-box number. signal_layers reports how many of
    the 3 layers show any signal at all (e.g. "3/3").

    This tool does not infer disease relevance — it only summarizes read-level
    structural-variant evidence at a candidate breakpoint.

    Args:
        bam_path:   Path to indexed BAM file
        chromosome: Chromosome of the candidate breakpoint
        position:   Candidate breakpoint position
        label:      Optional free-text label for this locus (e.g. a case/event name)
        window_bp:  Half-width of the inspection window (default 500bp)
        min_mapq:   Minimum mapping quality to include a read (default 20)

    Returns:
        BreakpointEvidenceSummary as dict, or error dict if underlying tools fail
    """
    stats = get_bam_stats_at_locus(
        bam_path, chromosome, max(0, position - window_bp), position + window_bp
    )
    disc = count_discordant_pairs(
        bam_path, chromosome, position, window_bp=window_bp, min_mapq=min_mapq
    )
    clips = count_soft_clipped_reads(
        bam_path, chromosome, position,
        window_bp=min(window_bp, 200), min_mapq=min_mapq
    )
    split = get_split_reads(
        bam_path, chromosome, position,
        window_bp=min(window_bp, 200), min_mapq=min_mapq
    )

    for result in (stats, disc, clips, split):
        if "error" in result:
            return {"error": result["error"], "bam_path": bam_path}

    observations = []

    if stats["total_reads"] == 0:
        observations.append("No reads found in the inspection window — evidence cannot be assessed.")

    # ── Discordant-pair component (0-50) ──
    disc_fraction = disc["discordant_fraction"]
    if disc_fraction >= 0.5:
        discordant_pair_score = 50.0
    elif disc_fraction >= 0.2:
        discordant_pair_score = 30.0
    elif disc_fraction > 0:
        discordant_pair_score = 15.0
    else:
        discordant_pair_score = 0.0

    if disc["discordant_pairs"] > 0:
        top_mate_chrom = next(iter(disc["mate_chromosomes"]))
        observations.append(
            f"{disc['discordant_pairs']} discordant pair(s) "
            f"({disc_fraction:.0%} of reads in window) with mates mapping "
            f"predominantly to {top_mate_chrom}."
        )

    # ── Soft-clip component (0-50) ──
    clip_fraction = clips["soft_clipped_fraction"]
    if clip_fraction >= 0.3:
        soft_clip_score = 50.0
    elif clip_fraction >= 0.1:
        soft_clip_score = 30.0
    elif clip_fraction > 0:
        soft_clip_score = 15.0
    else:
        soft_clip_score = 0.0

    if clips["soft_clipped_reads"] > 0:
        observations.append(
            f"{clips['soft_clipped_reads']} soft-clipped read(s) "
            f"({clip_fraction:.0%} of reads in window), "
            f"consensus clip position at {clips['consensus_clip_position']} "
            f"({clips['max_clips_at_position']} reads)."
        )

    # ── Split-read component (0-50) ──
    split_fraction = split["split_read_fraction"]
    if split_fraction >= 0.3:
        split_read_score = 50.0
    elif split_fraction >= 0.1:
        split_read_score = 30.0
    elif split_fraction > 0:
        split_read_score = 15.0
    else:
        split_read_score = 0.0

    if split["split_reads"] > 0:
        top_partner_chrom = next(iter(split["partner_chromosomes"]))
        observations.append(
            f"{split['split_reads']} split read(s) "
            f"({split_fraction:.0%} of reads in window) with supplementary "
            f"alignments mapping predominantly to {top_partner_chrom} "
            f"(e.g. {split['example_partner_loci'][:1]})."
        )

    # ── Combine: normalize 3 x (0-50) components onto a 0-100 scale ──
    raw_sum = discordant_pair_score + soft_clip_score + split_read_score
    evidence_score = round(raw_sum / 1.5, 1)
    signal_layers = f"{sum(1 for s in (discordant_pair_score, soft_clip_score, split_read_score) if s > 0)}/3"

    if evidence_score >= 70:
        evidence_strength = "strong"
    elif evidence_score >= 40:
        evidence_strength = "moderate"
    elif evidence_score > 0:
        evidence_strength = "weak"
    else:
        evidence_strength = "none"
        if stats["total_reads"] > 0:
            observations.append("No discordant pairs, soft-clipping, or split reads detected near this position.")

    result = BreakpointEvidenceSummary(
        label=label,
        chromosome=chromosome,
        position=position,
        evidence_score=evidence_score,
        evidence_strength=evidence_strength,
        signal_layers=signal_layers,
        discordant_pair_score=discordant_pair_score,
        soft_clip_score=soft_clip_score,
        split_read_score=split_read_score,
        locus_stats=stats,
        discordant_pairs=disc,
        soft_clips=clips,
        split_reads=split,
        supporting_observations=observations,
    )
    return asdict(result)