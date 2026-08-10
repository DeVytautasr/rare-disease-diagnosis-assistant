"""
bam_tools.py
Core BAM inspection tools for the Stage 1 IGV breakpoint assistant.

Every function returns structured data only — no genomic interpretation.
The LLM reads these outputs and writes the report. It cannot add facts
that are not in the tool output.
"""

import pysam
import requests as _requests
import time as _time
import subprocess as _subprocess
import os as _os
import tempfile as _tempfile
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
class ReadDepthProfile:
    chromosome: str
    start: int
    end: int
    window_size: int
    windows: list      # [{window_start, window_end, depth}, ...]
    summary: dict       # {min_depth, max_depth, mean_depth, depth_ratio_min_to_mean, likely_deletion}

@dataclass
class GeneLocus:
    chromosome: str
    position: int
    gene_count: int
    is_intergenic: bool
    genes: list

@dataclass
class BreakpointEvidenceSummary:
    label: str
    chromosome: str
    position: int
    evidence_score: float          # 0-100, direct sum of the 4 components below
    evidence_strength: str         # "none" | "weak" | "moderate" | "strong"
    signal_layers: str             # e.g. "3/4" — how many of the 4 layers show any signal
    discordant_pair_score: float   # 0-25
    soft_clip_score: float         # 0-25
    split_read_score: float        # 0-25
    depth_score: float             # 0-25
    locus_stats: dict
    discordant_pairs: dict
    soft_clips: dict
    split_reads: dict
    depth_profile: dict
    supporting_observations: list
    interpretation_template: str   # plain-language recap built only from the fields above


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


# ── Tool 5: Read depth profile ─────────────────────────────────────────────────

def get_read_depth_profile(
    bam_path: str,
    chromosome: str,
    start: int,
    end: int,
    window_size: int = 100
) -> dict:
    """
    Computes mean read depth in sliding windows across a region.
    Useful for detecting copy-number changes at SV breakpoints:
    - Deletions: depth drops inside the deleted region
    - Duplications: depth rises
    - Balanced events: depth stays flat (use other evidence layers)

    Depth here is approximated as the number of reads with actual aligned
    sequence in each window (not per-base pileup depth). Overlap is judged
    from each read's aligned reference positions rather than its
    reference_start/reference_end span, so a read carrying an internal CIGAR
    deletion is correctly NOT counted as covering the deleted region it
    spans but has no sequence in — this is exactly the case that matters
    for spotting a deletion inside a long read.

    Args:
        bam_path:    Path to indexed BAM file
        chromosome:  Chromosome of the region
        start:       Region start (0-based)
        end:         Region end
        window_size: Width of each sliding window in bp (default 100bp)

    Returns:
        ReadDepthProfile as dict, or error dict if BAM cannot be read
    """
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": str(e), "bam_path": bam_path}

    window_starts = list(range(start, end, window_size))
    counts = [0] * len(window_starts)

    for read in bam.fetch(chromosome, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        touched_windows = set()
        for pos in read.get_reference_positions():
            if start <= pos < end:
                touched_windows.add((pos - start) // window_size)
        for idx in touched_windows:
            counts[idx] += 1

    bam.close()

    windows = [
        {"window_start": w_start, "window_end": min(w_start + window_size, end), "depth": depth}
        for w_start, depth in zip(window_starts, counts)
    ]

    min_depth = min(counts) if counts else 0
    max_depth = max(counts) if counts else 0
    mean_depth = round(sum(counts) / len(counts), 2) if counts else 0.0
    depth_ratio_min_to_mean = round(min_depth / mean_depth, 3) if mean_depth > 0 else 0.0

    summary = {
        "min_depth": min_depth,
        "max_depth": max_depth,
        "mean_depth": mean_depth,
        "depth_ratio_min_to_mean": depth_ratio_min_to_mean,
        "likely_deletion": depth_ratio_min_to_mean < 0.6,
    }

    result = ReadDepthProfile(
        chromosome=chromosome,
        start=start,
        end=end,
        window_size=window_size,
        windows=windows,
        summary=summary,
    )
    return asdict(result)


# ── Tool 6: Combined breakpoint evidence summary ──────────────────────────────

def summarize_breakpoint_evidence(
    bam_path: str,
    chromosome: str,
    position: int,
    label: str = "",
    window_bp: int = 500,
    min_mapq: int = 20
) -> dict:
    """
    Combines discordant-pair, soft-clip, split-read, and read-depth evidence
    into a single interpretable breakpoint evidence summary. Each of the 4
    evidence layers is scored independently on a 0-25 scale (tiers: 0 / 7.5 /
    15 / 25 for discordant-pair, soft-clip, and split-read; 0 / 15 / 25 for
    depth, which has 3 tiers instead of 4). evidence_score is the direct,
    unweighted sum of the four component scores — discordant_pair_score +
    soft_clip_score + split_read_score + depth_score always equals
    evidence_score exactly, so the contribution of each layer stays visible
    rather than collapsed into a black-box number. signal_layers reports how
    many of the 4 layers show any signal at all (e.g. "3/4").

    This tool does not infer disease relevance — it only summarizes read-level
    structural-variant evidence at a candidate breakpoint.

    Depth profile uses a 4kb window (±2kb from position) to capture
    deletions larger than the short-read fragment size. Threshold 0.7
    calibrated against GIAB HG002 Illumina 300x validation.

    Args:
        bam_path:   Path to indexed BAM file
        chromosome: Chromosome of the candidate breakpoint
        position:   Candidate breakpoint position
        label:      Optional free-text label for this locus (e.g. a case/event name)
        window_bp:  Half-width of the inspection window for discordant-pair /
                    soft-clip / split-read evidence (default 500bp). Does NOT
                    affect the depth-profile window, which is fixed at ±2kb
                    (see above) regardless of this value.
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
    depth_profile = get_read_depth_profile(
        bam_path, chromosome, max(0, position - 2000), position + 2000,
        window_size=200
    )

    for result in (stats, disc, clips, split, depth_profile):
        if "error" in result:
            return {"error": result["error"], "bam_path": bam_path}

    observations = []

    if stats["total_reads"] == 0:
        observations.append("No reads found in the inspection window — evidence cannot be assessed.")

    # ── Discordant-pair component (0-25) ──
    disc_fraction = disc["discordant_fraction"]
    if disc_fraction >= 0.5:
        discordant_pair_score = 25.0
    elif disc_fraction >= 0.2:
        discordant_pair_score = 15.0
    elif disc_fraction > 0:
        discordant_pair_score = 7.5
    else:
        discordant_pair_score = 0.0

    if disc["discordant_pairs"] > 0:
        top_mate_chrom = next(iter(disc["mate_chromosomes"]))
        observations.append(
            f"{disc['discordant_pairs']} discordant pair(s) "
            f"({disc_fraction:.0%} of reads in window) with mates mapping "
            f"predominantly to {top_mate_chrom}."
        )

    # ── Soft-clip component (0-25) ──
    clip_fraction = clips["soft_clipped_fraction"]
    if clip_fraction >= 0.3:
        soft_clip_score = 25.0
    elif clip_fraction >= 0.1:
        soft_clip_score = 15.0
    elif clip_fraction > 0:
        soft_clip_score = 7.5
    else:
        soft_clip_score = 0.0

    if clips["soft_clipped_reads"] > 0:
        observations.append(
            f"{clips['soft_clipped_reads']} soft-clipped read(s) "
            f"({clip_fraction:.0%} of reads in window), "
            f"consensus clip position at {clips['consensus_clip_position']} "
            f"({clips['max_clips_at_position']} reads)."
        )

    # ── Split-read component (0-25) ──
    split_fraction = split["split_read_fraction"]
    if split_fraction >= 0.3:
        split_read_score = 25.0
    elif split_fraction >= 0.1:
        split_read_score = 15.0
    elif split_fraction > 0:
        split_read_score = 7.5
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

    # ── Read-depth component (0-25) ──
    # depth_ratio_min_to_mean < 0.7 flags a possible deletion (coverage drop);
    # the lower the ratio, the more pronounced the drop. 0.7 (raised from an
    # initial 0.6) is calibrated against the real GIAB HG002 deletion at
    # chr1:115,686,862, which measured 0.609 (PacBio HiFi) and 0.542
    # (Illumina 300x) using this tool's ±2kb/200bp-window depth profile —
    # see REAL_DATA_VALIDATION.md. Calibrated against a single confirmed
    # locus (replicated across 2 sequencing technologies, not 2 independent
    # loci) — treat as a starting point, not a validated general threshold.
    depth_ratio = depth_profile["summary"]["depth_ratio_min_to_mean"]
    if depth_ratio < 0.3:
        depth_score = 25.0
    elif depth_ratio < 0.7:
        depth_score = 15.0
    else:
        depth_score = 0.0

    if depth_score > 0:
        observations.append(
            f"Read depth drops to {depth_ratio:.0%} of the window mean "
            f"(min {depth_profile['summary']['min_depth']} vs mean "
            f"{depth_profile['summary']['mean_depth']} reads/window) — "
            f"consistent with a possible deletion."
        )

    # ── Combine: each component is already 0-25, so evidence_score is their
    # direct sum (0-100) — no separate normalization step. This means
    # discordant_pair_score + soft_clip_score + split_read_score + depth_score
    # always equals evidence_score exactly. ──
    component_scores = (discordant_pair_score, soft_clip_score, split_read_score, depth_score)
    evidence_score = round(sum(component_scores), 1)
    signal_layers = f"{sum(1 for s in component_scores if s > 0)}/4"

    if evidence_score >= 70:
        evidence_strength = "strong"
    elif evidence_score >= 40:
        evidence_strength = "moderate"
    elif evidence_score > 0:
        evidence_strength = "weak"
    else:
        evidence_strength = "none"
        if stats["total_reads"] > 0:
            observations.append("No discordant pairs, soft-clipping, split reads, or depth changes detected near this position.")

    # interpretation_template is built ONLY from values already computed above —
    # it adds no new facts, just restates them as one plain-language sentence.
    interpretation_template = (
        f"Breakpoint {label} at {chromosome}:{position}. "
        f"Evidence strength: {evidence_strength}. "
        f"Score: {evidence_score}/100 ({signal_layers} evidence layers showing signal). "
        f"Observations: {'; '.join(observations) if observations else 'none'}. "
        f"Technology note: discordant_pairs only valid for paired-end data; "
        f"split_reads only valid for modern-alignment BAMs with SA tags. "
        f"All values in this template come from tool outputs only."
    )

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
        depth_score=depth_score,
        locus_stats=stats,
        discordant_pairs=disc,
        soft_clips=clips,
        split_reads=split,
        depth_profile=depth_profile,
        supporting_observations=observations,
        interpretation_template=interpretation_template,
    )
    return asdict(result)


# ── Tool 7: Gene lookup at a locus ─────────────────────────────────────────────

def get_gene_at_locus(
    chromosome: str,
    position: int,
    genome_build: str = "GRCh38",
    max_retries: int = 3
) -> dict:
    """
    Query Ensembl REST API for genes overlapping a position.
    Tells the assistant whether a breakpoint disrupts a known gene,
    hits an intron, or falls in intergenic space — the key clinical question.

    Retries on HTTP 429 (rate-limited) with exponential backoff, and retries
    on network/timeout exceptions the same way. Ensembl's public REST API is
    observed to be intermittently slow (occasional 15s+ timeouts even on
    healthy requests) — this makes a single-attempt call unreliable in
    practice, so a transient failure is retried before giving up.

    Args:
        chromosome:   e.g. "chr1" or "1"
        position:     genomic position
        genome_build: "GRCh38" (default) or "GRCh37"
        max_retries:  attempts before giving up (default 3)

    Returns:
        dict with gene names, biotypes, and whether position is intergenic,
        or an error dict (with a "note" explaining gene annotation was skipped)
        if Ensembl could not be reached after max_retries attempts
    """
    chrom = chromosome.replace("chr", "")

    if genome_build == "GRCh37":
        server = "https://grch37.rest.ensembl.org"
    else:
        server = "https://rest.ensembl.org"

    url = f"{server}/overlap/region/human/{chrom}:{position}-{position}"
    headers = {"Content-Type": "application/json"}
    params = {"feature": "gene"}

    last_error = None
    for attempt in range(max_retries):
        try:
            r = _requests.get(url, headers=headers, params=params, timeout=15)
            if r.status_code == 200:
                genes = r.json()
                gene_list = [
                    {
                        "gene_id": g.get("gene_id", "unknown"),
                        "gene_name": g.get("external_name", "unknown"),
                        "biotype": g.get("biotype", "unknown"),
                        "strand": "+" if g.get("strand", 1) == 1 else "-",
                        "gene_start": g.get("start"),
                        "gene_end": g.get("end"),
                    }
                    for g in genes
                ]
                return {
                    "chromosome": chromosome,
                    "position": position,
                    "genome_build": genome_build,
                    "gene_count": len(gene_list),
                    "is_intergenic": len(gene_list) == 0,
                    "genes": gene_list,
                    "clinical_note": (
                        f"Breakpoint directly disrupts {len(gene_list)} gene(s)."
                        if gene_list
                        else "Breakpoint is intergenic — check nearby genes for positional effects."
                    )
                }
            elif r.status_code == 429:
                wait = 2 ** attempt
                _time.sleep(wait)
                last_error = f"HTTP 429 rate-limited (attempt {attempt + 1})"
            else:
                last_error = f"HTTP {r.status_code}"
                break
        except Exception as e:
            last_error = str(e)
            if attempt < max_retries - 1:
                _time.sleep(2 ** attempt)

    return {
        "error": last_error,
        "chromosome": chromosome,
        "position": position,
        "note": "Ensembl unavailable — gene annotation skipped"
    }


# ── Tool 8: Reciprocal breakpoint check ────────────────────────────────────────

def check_reciprocal_breakpoint(
    bam_path: str,
    primary_chromosome: str,
    primary_position: int,
    partner_chromosome: str,
    partner_position: int,
    window_bp: int = 500,
    min_mapq: int = 20
) -> dict:
    """
    For a suspected balanced translocation, verify the reciprocal breakpoint.

    If discordant pairs at primary_chromosome:primary_position point to
    partner_chromosome, this function checks partner_chromosome:partner_position
    and confirms whether discordant pairs there point BACK to primary_chromosome.

    True balanced translocation: both sides show reciprocal discordant signal.
    One-sided signal only: artifact or unbalanced event.

    Args:
        bam_path:           Path to indexed BAM
        primary_chromosome: First breakpoint chromosome (e.g. "chr1")
        primary_position:   First breakpoint position
        partner_chromosome: Partner chromosome (e.g. "chr8")
        partner_position:   Estimated partner breakpoint position
        window_bp:          Half-window for read counting
        min_mapq:           Minimum mapping quality

    Returns:
        dict with primary evidence, reciprocal evidence, and reciprocity verdict
    """
    # Check primary side
    primary = count_discordant_pairs(
        bam_path, primary_chromosome, primary_position, window_bp, min_mapq
    )

    # Check reciprocal side
    reciprocal = count_discordant_pairs(
        bam_path, partner_chromosome, partner_position, window_bp, min_mapq
    )

    # Reciprocity check: does partner side have discordant mates pointing back?
    primary_disc = primary.get("discordant_pairs", 0)
    reciprocal_disc = reciprocal.get("discordant_pairs", 0)

    # Check if reciprocal mates point back to primary chromosome
    reciprocal_mate_chroms = reciprocal.get("mate_chromosomes", {})
    back_pointing = reciprocal_mate_chroms.get(primary_chromosome, 0)

    if primary_disc >= 5 and reciprocal_disc >= 5 and back_pointing >= 3:
        verdict = "RECIPROCAL CONFIRMED — both breakpoints show concordant inter-chromosomal signal"
    elif primary_disc >= 5 and reciprocal_disc >= 2:
        verdict = "RECIPROCAL LIKELY — primary signal strong, partner signal present but weak"
    elif primary_disc >= 5 and reciprocal_disc == 0:
        verdict = "RECIPROCAL NOT FOUND — only one side shows signal; may be artifact or wrong partner coords"
    else:
        verdict = "INSUFFICIENT EVIDENCE at both positions"

    return {
        "primary": {
            "chromosome": primary_chromosome,
            "position": primary_position,
            "discordant_pairs": primary_disc,
            "mate_chromosomes": primary.get("mate_chromosomes", {}),
        },
        "reciprocal": {
            "chromosome": partner_chromosome,
            "position": partner_position,
            "discordant_pairs": reciprocal_disc,
            "back_pointing_to_primary": back_pointing,
            "mate_chromosomes": reciprocal_mate_chroms,
        },
        "verdict": verdict,
        "is_balanced": primary_disc >= 3 and reciprocal_disc >= 3,
    }


# ── Tool 9: IGV screenshot (visual evidence) ───────────────────────────────────

def run_igv_screenshot(
    bam_paths: list,
    chromosome: str,
    start: int,
    end: int,
    output_path: str,
    genome_build: str = "hg38",
    color_by: str = "MATE_CHROMOSOME",
    show_soft_clips: bool = True,
    max_coverage: int = None,
    squish: bool = True,
    igv_path: str = None,
    timeout_sec: int = 180
) -> dict:
    """
    Generate an IGV screenshot of a genomic region using headless batch mode.

    Produces the visual evidence a clinician would inspect manually:
    discordant pairs colored by mate chromosome, soft-clipped bases shown,
    and the region centred on the candidate breakpoint.

    Requires a usable display (real X11, or a Wayland/X compatibility layer
    such as WSLg) for IGV's Swing UI to render into — it is not invoked with
    an empty/overridden DISPLAY, since that crashes IGV's AWT event thread
    before it can take the snapshot. On a machine with no display at all,
    run this under `xvfb-run` (not handled internally, since that adds a
    dependency this tool doesn't otherwise need).

    Args:
        bam_paths:       List of BAM file paths or URLs to load as tracks
        chromosome:      Chromosome (e.g. "chr1")
        start:           Region start position
        end:             Region end position
        output_path:     Where to save the PNG
        genome_build:    IGV genome identifier ("hg38", "hg19")
        color_by:        IGV coloring mode. Options:
                         MATE_CHROMOSOME (translocations),
                         PAIR_ORIENTATION (inversions),
                         INSERT_SIZE (deletions/duplications),
                         NONE
        show_soft_clips: Display soft-clipped bases
        max_coverage:    If set, fixes the coverage track's max value (via
                         IGV's setDataRange) instead of autoscaling to the
                         tallest window in view. Set this slightly above the
                         observed max_depth from read_depth_profile so a
                         depth dip elsewhere in the region isn't flattened
                         by autoscaling to a taller peak.
        squish:          Compress read rows into a denser view when True
        igv_path:        Path to igv.sh (auto-detected if None)
        timeout_sec:     Max seconds to wait for IGV

    Returns:
        dict with screenshot_path, batch_script used, and success status
    """
    # Auto-detect IGV
    candidates = []
    if igv_path is None:
        candidates = [
            _os.path.expanduser("~/IGV_2.17.4/igv.sh"),
            _os.path.expanduser("~/igv/igv.sh"),
            "/opt/igv/igv.sh",
        ]
        for c in candidates:
            if _os.path.exists(c):
                igv_path = c
                break

    if igv_path is None or not _os.path.exists(igv_path):
        return {
            "error": "IGV not found",
            "searched": candidates if not igv_path else [igv_path],
            "note": "Install IGV or pass igv_path explicitly"
        }

    # Ensure output directory exists
    out_dir = _os.path.dirname(_os.path.abspath(output_path))
    _os.makedirs(out_dir, exist_ok=True)

    # Build IGV batch script
    lines = ["new", f"genome {genome_build}"]
    for bam in bam_paths:
        lines.append(f"load {bam}")
    lines.append("maxPanelHeight 600")
    if color_by and color_by != "NONE":
        lines.append(f"colorBy {color_by}")
    if show_soft_clips:
        lines.append("preference SAM.SHOW_SOFT_CLIPPED true")
    lines.append("preference SAM.SHOW_CENTER_LINE true")
    if max_coverage is not None:
        lines.append("preference SAM.MAX_VISIBLE_RANGE 1000")
    lines.append(f"goto {chromosome}:{start}-{end}")
    if max_coverage is not None:
        for bam in bam_paths:
            # The coverage sub-track's name is "<bam basename> Coverage" (with
            # a literal space), distinct from the alignment track's own name.
            # IGV's batch parser splits on unquoted spaces, so the track name
            # must be quoted as one token or the match silently fails and
            # autoscale stays on.
            track_name = f"{_os.path.basename(bam)} Coverage"
            lines.append(f'setDataRange 0,{max_coverage} "{track_name}"')
    lines.append("sort position")
    if squish:
        lines.append("squish")
    lines.append(f"snapshot {_os.path.abspath(output_path)}")
    lines.append("exit")

    batch_content = "\n".join(lines)

    # Write batch script to temp file
    with _tempfile.NamedTemporaryFile(mode="w", suffix=".bat",
                                       delete=False) as f:
        f.write(batch_content)
        batch_path = f.name

    try:
        # No env= override here: this inherits the parent process's DISPLAY.
        # Passing DISPLAY="" (rather than leaving it unset or omitting the
        # override entirely) crashes IGV's AWT EventDispatchThread before
        # the batch script runs, so no snapshot is ever produced — confirmed
        # directly against this IGV build, not assumed.
        result = _subprocess.run(
            [igv_path, "--batch", batch_path],
            capture_output=True,
            text=True,
            timeout=timeout_sec,
        )

        success = _os.path.exists(output_path)
        file_size = _os.path.getsize(output_path) if success else 0

        return {
            "success": success,
            "screenshot_path": _os.path.abspath(output_path) if success else None,
            "file_size_bytes": file_size,
            "region": f"{chromosome}:{start}-{end}",
            "color_by": color_by,
            "bam_tracks": len(bam_paths),
            "batch_script": batch_content,
            "igv_stdout": result.stdout[-500:] if result.stdout else "",
            "igv_stderr": result.stderr[-500:] if result.stderr else "",
        }
    except _subprocess.TimeoutExpired:
        return {"error": f"IGV timed out after {timeout_sec}s",
                "batch_script": batch_content}
    except Exception as e:
        return {"error": str(e), "batch_script": batch_content}
    finally:
        if _os.path.exists(batch_path):
            _os.unlink(batch_path)