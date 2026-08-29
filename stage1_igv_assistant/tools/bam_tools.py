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
import signal as _signal
import tempfile as _tempfile
import hashlib as _hashlib
import json as _json
from dataclasses import dataclass, asdict
from typing import Optional


# ── Threshold inventory and convention ──────────────────────────────────────
#
# "Threshold" means any numeric cutoff that changes what the assistant
# reports -- whether by altering a component score or by altering the prose a
# model reads and may quote. Strength bands are excluded: they only rename an
# already-computed score. Caller-overridable input filters are excluded but
# named below.
#
# Under that convention: 14 thresholds -- 11 scoring, 3 text-only -- of which
# 2 are empirically derived.
#
# SCORING (11), all in summarize_breakpoint_evidence:
#   discordant_pair_score   disc_fraction  >= 0.5 -> 25 | >= 0.2 -> 15 | > 0 -> 7.5
#   soft_clip_score         max_clips      >= 10  -> 25 | >= 3   -> 15
#   split_read_score        split_fraction >= 0.3 -> 25 | >= 0.1 -> 15 | > 0 -> 7.5
#   depth_score             depth_ratio    <  0.3 -> 25 | <  0.7 -> 15
#                           dip_tolerance_bp = 1000 zeroes a non-zero depth
#                           score when the dip is not localised to the focus
#
# TEXT-ONLY (3), which change no score but change what a model may quote:
#   PARTNER_DOMINANCE_MIN_SHARE = 0.6   } together gate whether a partner
#   PARTNER_DOMINANCE_MIN_READS = 3     } chromosome is called "predominant"
#   SOFT_CLIP_PILEUP_MIN_READS  = 3     gates "consensus clip position"
#                                       vs "no clip pileup"
#
# The text-only gates are counted deliberately. The predominance gate changed
# no score and caused a retraction: a model quoted the sentence it produced
# rather than the number behind it, and the published finding blamed the
# model. A convention that excluded these would have hidden the cutoff that
# did the most documented damage. See
# results/BENCHMARK_LOCAL_MODELS.md's correction notice.
#
# EMPIRICALLY DERIVED (2 of 14):
#   DEPTH_RATIO_DELETION_THRESHOLD = 0.7  one locus, two technologies
#   dip_tolerance_bp = 1000               two real loci, margin documented
#                                         on both sides (see
#                                         get_read_depth_profile's docstring)
# The remaining 12 are the author's judgement, documented as such.
#
# EXCLUDED but named: min_mapq = 20, the read-quality filter applied in every
# counting function. It is caller-overridable and not part of the scoring
# rubric, but it is the one judgement call that moves every fraction the
# scoring is built from -- change it and every tier above sees different
# input. low_mapq_fraction > 0.4 appears only as advisory text in a
# docstring; nothing in the code compares against it.

# ── Calibrated constants ────────────────────────────────────────────────────────
#
# Single source of truth for "does this depth_ratio_min_to_mean look like a
# deletion" — used identically by get_read_depth_profile's own
# `likely_deletion` flag and by summarize_breakpoint_evidence's depth_score.
# These used to be two independently hardcoded values (0.6 and 0.7) that
# drifted out of sync when only one was recalibrated, producing
# contradictory fields in the same summarize_breakpoint_evidence response.
#
# Calibrated against the real GIAB HG002 deletion at chr1:115,686,862,
# which measured 0.609 (PacBio HiFi) and 0.542 (Illumina 300x) using this
# tool's ±2kb/200bp-window depth profile — see REAL_DATA_VALIDATION.md.
# Calibrated against a single confirmed locus (replicated across 2
# sequencing technologies, not 2 independent loci) — treat as a starting
# point, not a validated general threshold.
#
# 2026-08-11 note: the 0.609/0.542 figures above were measured with
# get_read_depth_profile's pre-FIX-1 implementation, which computed a
# per-bin read count (inflated ~1.5-1.7x true depth, not bin-size
# invariant), not true per-base depth, and used the region's GLOBAL
# minimum rather than a focus-position-localized one (pre-FIX-2). Both
# bugs are now fixed. Re-measured at the same locus with the corrected
# tool (true depth, focus_position=115686862, dip_tolerance_bp=1000):
# 0.472 (window_size=100) / 0.502 (window_size=200) — still comfortably
# below 0.7, so THIS threshold's conclusion is unchanged; only the
# intermediate ratio values were ever wrong. The 0.609/0.542 numbers are
# left here as the historical calibration record, not rewritten — see
# REAL_DATA_VALIDATION.md's "Post-fix re-validation" addendum for the full
# before/after comparison.
DEPTH_RATIO_DELETION_THRESHOLD = 0.7

# Canonical evidence-layer names, shared between summarize_breakpoint_evidence's
# applicable_layers parameter, detect_applicable_layers' return value, and
# case_object.py's SequencingInfo.applicable_evidence_layers property (which
# uses this exact vocabulary already).
EVIDENCE_LAYER_NAMES = ("discordant_pairs", "soft_clipped_reads", "split_reads", "read_depth")

# A "predominant" partner chromosome may only be claimed when one actually
# dominates: at least PARTNER_DOMINANCE_MIN_READS partner records in total
# AND the top partner holding at least PARTNER_DOMINANCE_MIN_SHARE of them.
#
# This exists because the original implementation took
# next(iter(chrom_counts)) — the FIRST-INSERTED key, not even the maximum —
# and labelled it "predominantly" unconditionally, with no check that a
# dominant partner existed. Consequences found in committed benchmark data
# (see results/BENCHMARK_LOCAL_MODELS.md's correction section):
#   - ADVERSARIAL (prompt falsely asserts a t(1;12) translocation): a single
#     discordant read produced "mates mapping predominantly to chr12",
#     handing the model a sentence that reads as corroborating the false
#     premise.
#   - NEGATIVE (control locus, expected finding: no credible signal): 5
#     mates on 5 different chromosomes produced "predominantly to chr9".
# Both directly contradict this module's own documented semantics
# (discordant_pairs' docstring: "Mates scattered across many different
# chromosomes = background noise").
PARTNER_DOMINANCE_MIN_READS = 3
PARTNER_DOMINANCE_MIN_SHARE = 0.6

# Minimum reads sharing one clip position before the word "consensus" is
# used. Matches the first scoring tier in summarize_breakpoint_evidence
# (max_clips >= 3). Below it, describing a "consensus clip position" asserts
# agreement among reads that does not exist — the same defect class as the
# partner-dominance bug above, found in the same audit: the NEGATIVE control
# runs reported "consensus clip position at X (1 reads)" off a single read.
SOFT_CLIP_PILEUP_MIN_READS = 3


def _describe_partner_distribution(chrom_counts: dict, mate_noun: str) -> str:
    """
    Describe how partner chromosomes are distributed WITHOUT asserting a
    clustering pattern the counts don't show. Returns a clause to follow the
    read count, one of:

      "(mate on chr12) — single read, not a clustering signal"
      "with mates mapping predominantly to chr8 (12/15)"
      "with mates scattered across 7 chromosomes — no dominant partner"

    Args:
        chrom_counts: {chromosome: count}, e.g. discordant_pairs'
                      mate_chromosomes or split_reads' partner_chromosomes.
        mate_noun:    singular noun for one partner record ("mate",
                      "supplementary alignment"); plural adds "s".

    Ties are broken by chromosome name so the output is deterministic
    rather than dependent on dict insertion order — the specific flaw that
    made the original next(iter(...)) implementation pick an arbitrary
    chromosome and call it predominant.
    """
    if not chrom_counts:
        return ""
    total = sum(chrom_counts.values())
    top_chrom, top_count = max(chrom_counts.items(), key=lambda kv: (kv[1], kv[0]))

    n_chroms = len(chrom_counts)

    if total == 1:
        return f"({mate_noun} on {top_chrom}) — single read, not a clustering signal"
    if (total >= PARTNER_DOMINANCE_MIN_READS
            and top_count / total >= PARTNER_DOMINANCE_MIN_SHARE):
        return (f"with {mate_noun}s mapping predominantly to {top_chrom} "
                f"({top_count}/{total})")
    if n_chroms == 1:
        # All on one chromosome but under the read floor: not "scattered"
        # (there is nothing to scatter across), and not enough reads to call
        # a pattern either. Say exactly that.
        return (f"with all {total} {mate_noun}s on {top_chrom} — too few reads "
                f"to establish a clustering pattern")
    return (f"with {mate_noun}s scattered across {n_chroms} chromosomes "
            f"— no dominant partner")


# ── Data structures ──────────────────────────────────────────────────────────

@dataclass
class LocusStats:
    chromosome: str
    start: int                 # clamped to the contig, not the requested value
    end: int                   # clamped to the contig, not the requested value
    contig_length: int         # from the BAM header
    clamped_to_contig: bool    # True when the request overran the contig
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
    discordant_fraction: Optional[float]   # None when total_reads_in_window == 0 (undefined, not 0)
    mate_chromosomes: dict     # {chr_name: count}
    reads_below_min_mapq: int  # reads present but dropped by the MAPQ filter
    quality_limited: bool      # True when reads existed and ALL were filtered out;
                               # the layer is still assessed, and scores zero
    assessable: bool           # False ONLY when the window held no reads at all
    reason: Optional[str]      # explains why, when assessable is False

@dataclass
class SoftClipResult:
    chromosome: str
    position: int
    window_bp: int
    total_reads_in_window: int
    soft_clipped_reads: int
    soft_clipped_fraction: Optional[float]   # None when total_reads_in_window == 0 (undefined, not 0)
    consensus_clip_position: Optional[int]   # position with most clipping (left+right combined)
    max_clips_at_position: int
    # Left (5') and right (3') clips are tracked separately so callers can
    # tell which side actually dominates — the combined fields above alone
    # can't distinguish this, since a left-clip's reference_start and a
    # right-clip's reference_end can coincide at the same breakpoint.
    left_clip_reads: int
    right_clip_reads: int
    left_consensus_position: Optional[int]
    left_max_clips_at_position: int
    right_consensus_position: Optional[int]
    right_max_clips_at_position: int
    dominant_clip_side: Optional[str]   # "left" | "right" | "tied" | None (no clips at all)
    reads_below_min_mapq: int  # reads present but dropped by the MAPQ filter
    quality_limited: bool      # True when reads existed and ALL were filtered out;
                               # the layer is still assessed, and scores zero
    assessable: bool           # False ONLY when the window held no reads at all
    reason: Optional[str]      # explains why, when assessable is False

@dataclass
class SplitReadResult:
    chromosome: str
    position: int
    window_bp: int
    total_reads_in_window: int
    split_reads: int
    split_read_fraction: Optional[float]   # None when total_reads_in_window == 0 (undefined, not 0)
    partner_chromosomes: dict     # {chr_name: count}, parsed from SA tags
    example_partner_loci: list    # up to 5 "chrom:pos" strings for inspection
    reads_below_min_mapq: int  # reads present but dropped by the MAPQ filter
    quality_limited: bool      # True when reads existed and ALL were filtered out;
                               # the layer is still assessed, and scores zero
    assessable: bool           # False ONLY when the window held no reads at all
    reason: Optional[str]      # explains why, when assessable is False

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
    evidence_score: Optional[float]  # 0-100, normalised over applicable+assessable layers;
                                      # None when no layer could be scored at all (see evidence_strength)
    evidence_score_raw: float      # 0-100, direct sum of all 4 components regardless of applicability;
                                    # unassessable/inapplicable layers contribute 0, never None
    evidence_strength: str         # "none" | "weak" | "moderate" | "strong" | "NOT ASSESSABLE"
    signal_layers: str             # e.g. "2/3" — signal-showing / total applicable+assessable layers;
                                    # "0/0" when evidence_strength is "NOT ASSESSABLE"
    applicable_layers: list        # which of the 4 layers were counted in evidence_score's denominator
                                    # (caller-supplied applicable_layers, before assessability filtering)
    unassessable_layers: dict      # {layer_name: reason} for layers with zero reads in their window —
                                    # excluded from both evidence_score's numerator and denominator
    discordant_pair_score: Optional[float]   # 0-25, or None if discordant_pairs was unassessable
    soft_clip_score: Optional[float]         # 0-25, or None if soft_clipped_reads was unassessable
    split_read_score: Optional[float]        # 0-25, or None if split_reads was unassessable
    depth_score: Optional[float]             # 0-25, or None if read_depth was unassessable
    locus_stats: dict
    discordant_pairs: dict
    soft_clips: dict
    split_reads: dict
    depth_profile: dict
    supporting_observations: list
    interpretation_template: str   # plain-language recap built only from the fields above


# ── Error-handling helpers ─────────────────────────────────────────────────────

def _resolve_contig(bam, chromosome: str) -> Optional[str]:
    """
    Resolve a chromosome name against the BAM header, trying the
    chr-prefixed/unprefixed alternate form if the name as given isn't
    present. Returns the resolved name usable with bam.fetch(), or None
    if neither form is found in the header.
    """
    refs = bam.references
    if chromosome in refs:
        return chromosome
    alt = chromosome[3:] if chromosome.startswith("chr") else f"chr{chromosome}"
    if alt in refs:
        return alt
    return None


def _canonical_chrom(name: str, refs) -> str:
    """
    Canonical spelling of a contig name for use as a dict key or comparison,
    resolved against the BAM header.

    The previous helper (_normalize_chrom, since removed) prepended "chr"
    unconditionally. That is right for the 26
    primary contigs and wrong for the 525 HLA-* contigs in hs38DH, which carry
    no chr prefix: it turned HLA-A*01:01:01:01 into chrHLA-A*01:01:01:01, a
    name absent from the header, and emitted it as though it were a real
    reference contig. In one MHC window 16 of 44 partner keys named contigs
    that do not exist, with no error raised
    (REAL_PATIENT_DATA_VALIDATION.md finding 6).

    Here a name the header knows keeps the header's own spelling, and a name
    the header does not know is returned unchanged. Nothing is invented.
    Both operands of a comparison still normalise identically, which is the
    property that helper existed to provide.
    """
    if name in refs:
        return name
    alt = name[3:] if name.startswith("chr") else f"chr{name}"
    if alt in refs:
        return alt
    return name


def _chrom_count(counts: dict, name: str) -> int:
    """
    Look up a contig in a counts dict tolerant of chr-prefix spelling, since
    the dict is keyed by the header's spelling and the caller's string may
    use the other convention.
    """
    if name in counts:
        return counts[name]
    alt = name[3:] if name.startswith("chr") else f"chr{name}"
    return counts.get(alt, 0)


# Above this fraction of MAPQ<20 reads, a locus is reported as
# QUALITY-LIMITED instead of carrying a normalised evidence_score. 0.4 was
# already named as the advisory cutoff in this module's provenance notes but
# nothing compared against it; real-data validation made it operative. It is
# inherited judgement, not an independently validated figure.
LOW_MAPQ_QUALITY_GATE = 0.4


def _assess_window(passed: int, below_mapq: int):
    """
    Tri-state assessment of a counting window.

    Distinguishes a window with no reads at all -- structurally unassessable,
    and rightly excluded from the evidence denominator -- from a window whose
    reads were all removed by the MAPQ filter, which IS assessed and scores
    zero, and must stay in the denominator.

    Conflating the two made the evidence score rise as mapping quality fell:
    in an all-MAPQ-0 window the three filtered layers dropped out, the depth
    layer (which applies no MAPQ filter) was normalised over a denominator of
    1, and a 25/25 depth score became 100/100 "strong". See
    results/REAL_PATIENT_DATA_VALIDATION.md findings 1 and 2.

    Returns (assessable, reason, quality_limited).
    """
    if passed == 0 and below_mapq == 0:
        return False, "no reads in window", False
    if passed == 0:
        return True, None, True
    return True, None, False


def _contig_not_found_error(bam, chromosome: str) -> dict:
    return {
        "error": f"Chromosome '{chromosome}' not found in BAM header "
                 f"(also tried the alternate chr-prefix form).",
        "error_type": "invalid_region",
        "chromosome": chromosome,
        "contigs_in_header_sample": list(bam.references[:5]),
    }


def _validate_range(chromosome: str, start: int, end: int) -> Optional[dict]:
    """Returns an error dict for an invalid start/end pair, else None."""
    if start < 0 or end < 0:
        return {
            "error": f"start and end must be non-negative (got start={start}, end={end})",
            "error_type": "invalid_parameters",
            "chromosome": chromosome,
            "position": start,
        }
    if start > end:
        return {
            "error": f"start ({start}) must not be greater than end ({end})",
            "error_type": "invalid_parameters",
            "chromosome": chromosome,
            "position": start,
        }
    if start == end:
        # A zero-width window is a caller error, not a measurement. Left
        # unvalidated it silently produced "no reads in window" for every
        # counting layer, which dropped them from the evidence denominator and
        # let a clean locus score 60.0 "moderate" off the depth layer alone
        # (window_bp=0 -- REAL_PATIENT_DATA_VALIDATION.md finding 1).
        return {
            "error": f"zero-width region requested (start == end == {start}). "
                     f"A window with no width cannot be measured; if this came "
                     f"from window_bp=0, pass a positive window size.",
            "error_type": "invalid_parameters",
            "chromosome": chromosome,
            "position": start,
        }
    return None


def _clamp_or_error(bam, chromosome: str, resolved_chrom: str,
                    start: int, end: int):
    """
    Clamp a requested range to the contig's real bounds, or reject it as out
    of range.

    Two defects motivated this (REAL_PATIENT_DATA_VALIDATION.md findings 5
    and 14):

    * A window overrunning the contig end was neither clamped nor flagged,
      and the FULL requested width was reported. Asking for chrM:16000-17500
      on a 16,569 bp contig returned the same reads as 16000-16569 but a
      mean_depth of 1,011.06 instead of 2,665.36 -- the denominator included
      931 bases that do not exist. Depth profiles emitted whole bins past the
      contig end (20 of 20 in one case), and a window near the edge scored a
      maximum depth component off those empty bins.
    * A coordinate entirely past the end returned a clean zero result,
      indistinguishable from a genuine coverage gap. pysam does not raise for
      an out-of-range fetch, so nothing downstream noticed. A typo'd
      coordinate looked exactly like real data.

    Returns (start, end, clamp_info, error). Callers reassign start/end from
    this, so every span, denominator and echoed field is the clamped one.
    """
    length = bam.get_reference_length(resolved_chrom)
    if start >= length:
        return start, end, None, {
            "error": f"requested start {start} is at or past the end of "
                     f"{resolved_chrom} (contig length {length}). This is a "
                     f"coordinate error, not a coverage gap.",
            "error_type": "out_of_range",
            "chromosome": chromosome,
            "position": start,
            "contig_length": length,
        }
    new_start = max(0, start)
    new_end = min(end, length)
    return new_start, new_end, {
        "contig_length": length,
        "clamped_to_contig": (new_start != start or new_end != end),
        "requested_start": start,
        "requested_end": end,
    }, None


def _fetch_or_error(bam, chromosome: str, resolved_chrom: str, start: int, end: int):
    """
    Calls bam.fetch() and returns (iterator, None) on success or
    (None, error_dict) if the region is invalid (e.g. position beyond
    the end of the chromosome). pysam raises ValueError synchronously
    when fetch() is called with a bad region, before iteration starts.
    """
    try:
        return bam.fetch(resolved_chrom, start, end), None
    except ValueError as e:
        return None, {
            "error": str(e),
            "error_type": "invalid_region",
            "chromosome": chromosome,
            "position": start,
        }


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
        LocusStats as dict, or a structured error dict
        ({"error", "error_type", ...}) if the BAM cannot be read or the
        region is invalid. error_type is one of "bam_access",
        "invalid_region", "invalid_parameters".
    """
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": str(e), "error_type": "bam_access", "bam_path": bam_path}

    range_error = _validate_range(chromosome, start, end)
    if range_error:
        bam.close()
        return range_error

    resolved_chrom = _resolve_contig(bam, chromosome)
    if resolved_chrom is None:
        err = _contig_not_found_error(bam, chromosome)
        bam.close()
        return err

    start, end, clamp_info, range_error = _clamp_or_error(
        bam, chromosome, resolved_chrom, start, end
    )
    if range_error:
        bam.close()
        return range_error

    read_iter, fetch_error = _fetch_or_error(bam, chromosome, resolved_chrom, start, end)
    if fetch_error:
        bam.close()
        return fetch_error

    total = 0
    mapq_sum = 0
    low_mapq = 0
    forward = 0
    reverse = 0
    depth_positions = {}

    for read in read_iter:
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
        contig_length=clamp_info["contig_length"],
        clamped_to_contig=clamp_info["clamped_to_contig"],
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
        DiscordantPairResult as dict, or a structured error dict
        ({"error", "error_type", ...}) if the BAM cannot be read or the
        region is invalid. error_type is one of "bam_access",
        "invalid_region", "invalid_parameters".
    """
    start = max(0, position - window_bp)
    end = position + window_bp

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": str(e), "error_type": "bam_access", "bam_path": bam_path}

    range_error = _validate_range(chromosome, start, end)
    if range_error:
        bam.close()
        return range_error

    resolved_chrom = _resolve_contig(bam, chromosome)
    if resolved_chrom is None:
        err = _contig_not_found_error(bam, chromosome)
        bam.close()
        return err

    start, end, clamp_info, range_error = _clamp_or_error(
        bam, chromosome, resolved_chrom, start, end
    )
    if range_error:
        bam.close()
        return range_error

    read_iter, fetch_error = _fetch_or_error(bam, chromosome, resolved_chrom, start, end)
    if fetch_error:
        bam.close()
        return fetch_error

    total = 0
    discordant = 0
    mate_chroms = {}
    refs = set(bam.references)
    norm_chromosome = _canonical_chrom(chromosome, refs)

    below_mapq = 0
    for read in read_iter:
        if read.is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            below_mapq += 1
            continue
        # Unpaired reads (e.g. single-molecule long reads) have no mate at
        # all, so "mate_is_unmapped" defaults to False for them — without
        # this check they'd all be miscounted as discordant.
        if not read.is_paired:
            continue
        if read.mate_is_unmapped:
            continue

        total += 1

        # Discordant = mate on a different chromosome. Both sides are
        # normalised (chr-prefix-insensitive) before comparing: the mate's
        # name comes straight from the BAM header's own convention, while
        # `chromosome` is whatever convention the caller used — comparing
        # them raw silently downgrades every call where the two differ
        # (e.g. caller passes "1" against a "chr1"-header BAM, and every
        # same-chromosome mate gets miscounted as discordant).
        mate_name = read.next_reference_name
        mate_norm = _canonical_chrom(mate_name, refs) if mate_name is not None else None
        if mate_norm != norm_chromosome:
            discordant += 1
            mate_chr = mate_norm if mate_norm is not None else "unknown"
            mate_chroms[mate_chr] = mate_chroms.get(mate_chr, 0) + 1

    bam.close()

    # Sort mate chromosomes by count descending
    mate_chroms_sorted = dict(
        sorted(mate_chroms.items(), key=lambda x: x[1], reverse=True)
    )

    # total == 0 means the fraction is undefined, not 0 — a locus with no
    # reads at all is not the same claim as "reads present, none discordant".
    assessable, reason, quality_limited = _assess_window(total, below_mapq)

    result = DiscordantPairResult(
        chromosome=chromosome,
        position=position,
        window_bp=window_bp,
        total_reads_in_window=total,
        discordant_pairs=discordant,
        discordant_fraction=(round(discordant / total, 3) if total > 0
                             else (0.0 if quality_limited else None)),
        mate_chromosomes=mate_chroms_sorted,
        reads_below_min_mapq=below_mapq,
        quality_limited=quality_limited,
        assessable=assessable,
        reason=reason,
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

    consensus_clip_position/max_clips_at_position combine left-clips (5' end
    clipped, reads approaching a breakpoint from the right) and right-clips
    (3' end clipped, reads approaching from the left) into one count, since
    both can genuinely cluster at nearly the same position for a single
    clean breakpoint. left_/right_ prefixed fields report each side
    separately (position, max count, and total reads), and
    dominant_clip_side ("left" | "right" | "tied" | None) reports which
    side has more total clipped reads — e.g. to decide whether a
    visualization should sort by left-clip or right-clip length to make
    the pileup legible (see igv_evidence_panel).

    Args:
        bam_path:        Path to indexed BAM file
        chromosome:      Chromosome of the candidate breakpoint
        position:        Centre of the inspection window
        window_bp:       Half-width of the window (default 200bp)
        min_clip_bases:  Minimum bases clipped to count (default 10)
        min_mapq:        Minimum mapping quality (default 20)

    Returns:
        SoftClipResult as dict, or a structured error dict
        ({"error", "error_type", ...}) if the BAM cannot be read or the
        region is invalid. error_type is one of "bam_access",
        "invalid_region", "invalid_parameters".
    """
    start = max(0, position - window_bp)
    end = position + window_bp

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": str(e), "error_type": "bam_access", "bam_path": bam_path}

    range_error = _validate_range(chromosome, start, end)
    if range_error:
        bam.close()
        return range_error

    resolved_chrom = _resolve_contig(bam, chromosome)
    if resolved_chrom is None:
        err = _contig_not_found_error(bam, chromosome)
        bam.close()
        return err

    start, end, clamp_info, range_error = _clamp_or_error(
        bam, chromosome, resolved_chrom, start, end
    )
    if range_error:
        bam.close()
        return range_error

    read_iter, fetch_error = _fetch_or_error(bam, chromosome, resolved_chrom, start, end)
    if fetch_error:
        bam.close()
        return fetch_error

    total = 0
    clipped = 0
    clip_positions = {}
    left_clip_positions = {}
    right_clip_positions = {}

    below_mapq = 0
    for read in read_iter:
        if read.is_unmapped or read.is_secondary:
            continue
        if read.mapping_quality < min_mapq:
            below_mapq += 1
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
            left_clip_positions[clip_pos] = left_clip_positions.get(clip_pos, 0) + 1

        # Check for soft-clip at read end (right clip)
        elif cigar[-1][0] == 4 and cigar[-1][1] >= min_clip_bases:
            clipped += 1
            clip_pos = read.reference_end
            clip_positions[clip_pos] = clip_positions.get(clip_pos, 0) + 1
            right_clip_positions[clip_pos] = right_clip_positions.get(clip_pos, 0) + 1

    bam.close()

    def _consensus(positions):
        if not positions:
            return None, 0
        pos = max(positions, key=positions.get)
        return pos, positions[pos]

    consensus_pos, max_clips = _consensus(clip_positions)
    left_consensus_pos, left_max_clips = _consensus(left_clip_positions)
    right_consensus_pos, right_max_clips = _consensus(right_clip_positions)

    # total == 0 means the fraction is undefined, not 0.
    assessable, reason, quality_limited = _assess_window(total, below_mapq)

    # Which side dominates overall (not just at the combined consensus
    # position, since a left-clip cluster and a right-clip cluster from
    # reads approaching a breakpoint from opposite directions can coincide
    # at nearly the same position without either being individually
    # dominant there) — compares total clipped-read counts per side.
    left_total = sum(left_clip_positions.values())
    right_total = sum(right_clip_positions.values())
    if left_total == 0 and right_total == 0:
        dominant_side = None
    elif left_total > right_total:
        dominant_side = "left"
    elif right_total > left_total:
        dominant_side = "right"
    else:
        dominant_side = "tied"

    result = SoftClipResult(
        chromosome=chromosome,
        position=position,
        window_bp=window_bp,
        total_reads_in_window=total,
        soft_clipped_reads=clipped,
        soft_clipped_fraction=(round(clipped / total, 3) if total > 0
                               else (0.0 if quality_limited else None)),
        consensus_clip_position=consensus_pos,
        max_clips_at_position=max_clips,
        left_clip_reads=left_total,
        right_clip_reads=right_total,
        left_consensus_position=left_consensus_pos,
        left_max_clips_at_position=left_max_clips,
        right_consensus_position=right_consensus_pos,
        right_max_clips_at_position=right_max_clips,
        dominant_clip_side=dominant_side,
        reads_below_min_mapq=below_mapq,
        quality_limited=quality_limited,
        assessable=assessable,
        reason=reason,
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
        SplitReadResult as dict, or a structured error dict
        ({"error", "error_type", ...}) if the BAM cannot be read or the
        region is invalid. error_type is one of "bam_access",
        "invalid_region", "invalid_parameters".
    """
    start = max(0, position - window_bp)
    end = position + window_bp

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": str(e), "error_type": "bam_access", "bam_path": bam_path}

    range_error = _validate_range(chromosome, start, end)
    if range_error:
        bam.close()
        return range_error

    resolved_chrom = _resolve_contig(bam, chromosome)
    if resolved_chrom is None:
        err = _contig_not_found_error(bam, chromosome)
        bam.close()
        return err

    start, end, clamp_info, range_error = _clamp_or_error(
        bam, chromosome, resolved_chrom, start, end
    )
    if range_error:
        bam.close()
        return range_error

    read_iter, fetch_error = _fetch_or_error(bam, chromosome, resolved_chrom, start, end)
    if fetch_error:
        bam.close()
        return fetch_error

    total = 0
    split = 0
    partner_chroms = {}
    refs = set(bam.references)
    example_loci = []

    below_mapq = 0
    for read in read_iter:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            below_mapq += 1
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
            # Normalised so an aligner that writes "8" and one that writes
            # "chr8" for the same physical partner don't split into two
            # separate dict keys/counts.
            norm_rname = _canonical_chrom(rname, refs)
            partner_chroms[norm_rname] = partner_chroms.get(norm_rname, 0) + 1

        if first_rname is not None:
            # Display string keeps the SA tag's own original text (not
            # normalised) — this is for human inspection, not comparison.
            locus = f"{first_rname}:{first_pos}"
            if locus not in example_loci and len(example_loci) < 5:
                example_loci.append(locus)

    bam.close()

    partner_chroms_sorted = dict(
        sorted(partner_chroms.items(), key=lambda x: x[1], reverse=True)
    )

    # total == 0 means the fraction is undefined, not 0.
    assessable, reason, quality_limited = _assess_window(total, below_mapq)

    result = SplitReadResult(
        chromosome=chromosome,
        position=position,
        window_bp=window_bp,
        total_reads_in_window=total,
        split_reads=split,
        split_read_fraction=(round(split / total, 3) if total > 0
                             else (0.0 if quality_limited else None)),
        partner_chromosomes=partner_chroms_sorted,
        example_partner_loci=example_loci,
        reads_below_min_mapq=below_mapq,
        quality_limited=quality_limited,
        assessable=assessable,
        reason=reason,
    )
    return asdict(result)


# ── Tool 5: Read depth profile ─────────────────────────────────────────────────

def get_read_depth_profile(
    bam_path: str,
    chromosome: str,
    start: int,
    end: int,
    window_size: int = 100,
    focus_position: Optional[int] = None,
    dip_tolerance_bp: int = 1000,
) -> dict:
    """
    Computes true mean per-base read depth in sliding windows across a region.
    Useful for detecting copy-number changes at SV breakpoints:
    - Deletions: depth drops inside the deleted region
    - Duplications: depth rises
    - Balanced events: depth stays flat (use other evidence layers)

    Depth per window is the sum of aligned reference bases every read
    contributes within that window, divided by the window's actual width —
    the same true per-base coverage get_bam_stats_at_locus already computes
    (via the identical get_reference_positions() accumulation), just
    binned instead of averaged over one flat region. Overlap is judged
    from each read's aligned reference positions rather than its
    reference_start/reference_end span, so a read carrying an internal CIGAR
    deletion is correctly NOT counted as covering the deleted region it
    spans but has no sequence in — this is exactly the case that matters
    for spotting a deletion inside a long read.

    Before 2026-08-11 this instead counted the number of DISTINCT READS
    touching each bin (a read contributed a flat "+1" to every bin it had
    even one aligned base in, regardless of how much of it actually fell
    there). That is a read-count, not a depth: it scales with window_size
    (a read spanning a bin boundary gets double-counted more often in
    narrow bins) and runs systematically higher than true per-base
    coverage (confirmed on real data: ~1.5-1.7x bam_stats_at_locus's
    mean_depth for the same region, and not bin-size invariant — e.g.
    976.98 at window_size=100 vs 1354.25 at window_size=200 for the same
    span). See FIX 1 in git history / results/LLM_SESSION_3_BLIND.md for
    the worked example that surfaced this, and results docs predating this
    fix for numbers recorded under the old (read-count) implementation —
    those are annotated, not rewritten.

    summary["assessable"] is False when mean_depth is 0 (no reads touched
    any window in the queried region) — in that case
    depth_ratio_min_to_mean and likely_deletion are both None, not 0.0/True.
    A 0/0 ratio previously defaulted to 0.0, which reads as "depth dropped
    to 0% of the mean" — the strongest possible deletion signal — for a
    region that was never actually sequenced there. summary["reason"] is
    "no reads in queried region" when assessable is False, else None.

    focus_position (optional): when given, depth_ratio_min_to_mean (and so
    likely_deletion) is computed from the LOCAL depth around focus_position
    — the minimum across every bin within dip_tolerance_bp of it — against
    the window mean, instead of from the region's GLOBAL minimum. This
    matters because the global minimum of a wide scan can sit anywhere in
    the queried region, not necessarily anywhere near the position
    actually being investigated: a real but unrelated dip elsewhere in the
    window used to be indistinguishable from a dip AT the focus position
    (confirmed on real data: a dip ~1500bp from a queried position
    previously drove likely_deletion=True for a locus that was itself
    flat/at a local peak — see FIX 2 in results docs). The local search
    uses dip_tolerance_bp rather than a fixed one-bin neighborhood for the
    same reason dip_is_at_focus does (see that parameter's provenance note
    below) — a real transition's lowest sampled bin isn't always adjacent
    to where it starts. When given, the response also adds:
      dip_position:            window_start of the bin holding the
                                region's GLOBAL minimum depth (independent
                                of focus_position — always "where the
                                minimum actually is", for transparency)
      dip_distance_from_focus: abs(dip_position - focus_position)
      dip_is_at_focus:         True iff dip_distance_from_focus <=
                                dip_tolerance_bp

    dip_tolerance_bp (default 1000) — HEURISTIC, calibrated against exactly
    2 real loci, not validated beyond them: a real deletion breakpoint's
    depth transition is not always a clean step — it can decline gradually
    over several hundred bp past the breakpoint before reaching its lowest
    sampled bin, so requiring the global minimum to fall within a single
    bin (window_size, ~100-200bp) of focus_position is too strict and
    would wrongly zero out real signal. Calibrated against: chr1:16,890,000
    (a locus with NO real signal at the query position — the region's
    actual minimum is an unrelated fluctuation 1400bp away — correctly
    stays dip_is_at_focus=False at this tolerance) and chr1:115,686,862 (a
    real confirmed GIAB HG002 deletion breakpoint — the transition's
    steepest single-bin drop, -31%, is exactly at the query position, but
    its lowest sampled bin is 800-900bp downstream — correctly becomes
    dip_is_at_focus=True at this tolerance). 1000bp sits in the gap
    between those two real measurements with margin on both sides, but it
    is still only 2 data points: a real dip between roughly 1000-1400bp
    from a true breakpoint would currently be missed (dip_is_at_focus
    wrongly False), and an unrelated fluctuation within 1000bp of a flat
    locus would currently be miscounted as on-target (dip_is_at_focus
    wrongly True). Override via this parameter if a specific locus needs
    a different radius; revisit the default with more real loci if this
    becomes a bottleneck.

    All three dip_* fields are None when focus_position isn't given, or
    when the region isn't assessable (mirroring depth_ratio_min_to_mean/
    likely_deletion — an all-empty region has no meaningful dip location
    either). likely_deletion still compares against the same
    DEPTH_RATIO_DELETION_THRESHOLD; only its ratio's numerator changes.

    Args:
        bam_path:       Path to indexed BAM file
        chromosome:     Chromosome of the region
        start:          Region start (0-based)
        end:            Region end
        window_size:    Width of each sliding window in bp (default 100bp)
        focus_position: Optional genomic position to localize
                        depth_ratio_min_to_mean/likely_deletion/dip_* around,
                        instead of using the region's global minimum. Should
                        fall within [start, end) — outside that range it's
                        still accepted (clamped to the nearest scanned bin
                        for the local-ratio search) but dip_is_at_focus can
                        never be True, since no bin there was scanned.
        dip_tolerance_bp: Radius (bp) around focus_position used two ways:
                        (1) depth_ratio_min_to_mean's numerator is the
                        minimum depth among bins within this radius, and
                        (2) dip_is_at_focus is True iff the region's global
                        minimum also falls within this radius. Default
                        1000bp — see provenance note above. Ignored when
                        focus_position isn't given.

    Returns:
        ReadDepthProfile as dict, or a structured error dict
        ({"error", "error_type", ...}) if the BAM cannot be read or the
        region is invalid. error_type is one of "bam_access",
        "invalid_region", "invalid_parameters".
    """
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": str(e), "error_type": "bam_access", "bam_path": bam_path}

    range_error = _validate_range(chromosome, start, end)
    if range_error:
        bam.close()
        return range_error

    resolved_chrom = _resolve_contig(bam, chromosome)
    if resolved_chrom is None:
        err = _contig_not_found_error(bam, chromosome)
        bam.close()
        return err

    start, end, clamp_info, range_error = _clamp_or_error(
        bam, chromosome, resolved_chrom, start, end
    )
    if range_error:
        bam.close()
        return range_error

    read_iter, fetch_error = _fetch_or_error(bam, chromosome, resolved_chrom, start, end)
    if fetch_error:
        bam.close()
        return fetch_error

    window_starts = list(range(start, end, window_size))
    base_counts = [0] * len(window_starts)

    for read in read_iter:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        for pos in read.get_reference_positions():
            if start <= pos < end:
                base_counts[(pos - start) // window_size] += 1

    bam.close()

    # True per-base depth per bin: aligned bases summed into that bin,
    # divided by the bin's actual width — not window_size, since the last
    # bin in a region that doesn't divide evenly is narrower than the rest.
    windows = []
    depths = []
    for w_start, base_count in zip(window_starts, base_counts):
        w_end = min(w_start + window_size, end)
        bin_width = w_end - w_start
        depth = round(base_count / bin_width, 2) if bin_width > 0 else 0.0
        windows.append({"window_start": w_start, "window_end": w_end, "depth": depth})
        depths.append(depth)

    min_depth = min(depths) if depths else 0
    max_depth = max(depths) if depths else 0
    mean_depth = round(sum(depths) / len(depths), 2) if depths else 0.0

    # A region with zero reads (mean_depth == 0, whether from an empty
    # window list or windows that all read 0) has an UNDEFINED depth
    # ratio, not a ratio of 0. Previously depth_ratio_min_to_mean defaulted
    # to 0.0 here, which reads as 0/mean == 0% of mean depth — the single
    # most extreme possible "this looks like a deletion" value — and
    # likely_deletion (and downstream depth_score) took that at face
    # value, fabricating the strongest possible depth evidence from a
    # region with no data at all. assessable=False stops that: no ratio,
    # no likely_deletion verdict, and (in summarize_breakpoint_evidence)
    # no depth_score — a distinct "cannot be assessed" state rather than
    # a score.
    assessable = mean_depth > 0

    dip_position = None
    dip_distance_from_focus = None
    dip_is_at_focus = None

    if assessable:
        if focus_position is not None:
            # Same radius as dip_is_at_focus below (dip_tolerance_bp), not
            # just the single bin on either side of focus_position: a real
            # breakpoint's depth transition can decline gradually rather
            # than stepping cleanly, so its lowest sampled bin can sit
            # several hundred bp past the breakpoint itself (see
            # dip_tolerance_bp's provenance note above) — an immediate-
            # neighbor-only search would miss that real signal the same
            # way dip_is_at_focus would, for the same underlying reason,
            # so both use one consistent notion of "local".
            local_indices = [
                i for i, w_start in enumerate(window_starts)
                if abs(w_start - focus_position) <= dip_tolerance_bp
            ]
            if not local_indices:
                # focus_position farther than dip_tolerance_bp from every
                # scanned bin (e.g. outside [start, end) entirely) — fall
                # back to the single nearest bin so this never raises.
                nearest = (focus_position - start) // window_size
                local_indices = [max(0, min(nearest, len(depths) - 1))]
            local_min_depth = min(depths[i] for i in local_indices)
            depth_ratio_min_to_mean = round(local_min_depth / mean_depth, 3)
        else:
            depth_ratio_min_to_mean = round(min_depth / mean_depth, 3)
        likely_deletion = depth_ratio_min_to_mean < DEPTH_RATIO_DELETION_THRESHOLD

        min_idx = depths.index(min_depth)
        dip_position = windows[min_idx]["window_start"]
        if focus_position is not None:
            dip_distance_from_focus = abs(dip_position - focus_position)
            dip_is_at_focus = dip_distance_from_focus <= dip_tolerance_bp
    else:
        depth_ratio_min_to_mean = None
        likely_deletion = None

    summary = {
        "min_depth": min_depth,
        "max_depth": max_depth,
        "mean_depth": mean_depth,
        "depth_ratio_min_to_mean": depth_ratio_min_to_mean,
        "likely_deletion": likely_deletion,
        "dip_position": dip_position,
        "dip_distance_from_focus": dip_distance_from_focus,
        "dip_is_at_focus": dip_is_at_focus,
        "assessable": assessable,
        "reason": None if assessable else "no reads in queried region",
        # Clamped, not requested: bins past the contig end used to be emitted
        # as real zero-depth bins and pulled the mean down.
        "contig_length": clamp_info["contig_length"],
        "clamped_to_contig": clamp_info["clamped_to_contig"],
        # The depth layer deliberately applies NO MAPQ filter, while every
        # other counting layer filters at min_mapq=20. Reported rather than
        # changed: filtering here would silently shift every depth number and
        # invalidate what calibration the 0.7 deletion threshold has, which is
        # a separate open decision. Surfacing it lets a caller see that the
        # layers are not directly comparable.
        "min_mapq_applied": None,
        "mapq_note": ("depth counts ALL reads regardless of MAPQ; the "
                      "discordant, soft-clip and split layers filter at "
                      "min_mapq. In a low-MAPQ region this layer is the only "
                      "one still counting, so compare with care."),
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


# ── Composite scoring thresholds for summarize_breakpoint_evidence ────────────
#
# discordant-pair and split-read are scored 0/7.5/15/25 from a fraction
# cutoff; soft-clip is scored 0/15/25 (3 tiers, from 2026-08-11 — see below)
# from a max_clips_at_position cutoff; read-depth is scored 0/15/25 from a
# ratio cutoff, gated by a separate localization check. Honest accounting
# of where each number came from — "the numbers looked reasonable" is not
# the same as "the numbers were validated":
#
# DISCORDANT-PAIR (discordant_fraction cutoffs: 0.2, 0.5)              HEURISTIC
#   Chosen by judgement, not fit to data. No confirmed real balanced
#   translocation BAM was ever found to validate against (see
#   REAL_DATA_VALIDATION.md, "Bottom line") — every real BAM tested against
#   this tool was a deletion, not a translocation. The synthetic
#   translocation fixture in test_bam_tools.py hits the top tier by
#   construction (its discordant fraction was built to be clean), which
#   confirms the arithmetic works, not that 0.2/0.5 are correctly placed
#   for real data.
#
# SOFT-CLIP (max_clips_at_position cutoffs: 3, 10)                     HEURISTIC
#   Scored on max_clips_at_position — how many reads pile up at the SAME
#   position, i.e. count_soft_clipped_reads' own documented "< 3 = no real
#   pileup, treat as noise" line — not on soft_clipped_fraction (changed
#   2026-08-11). Fraction alone can't distinguish a genuine pileup from
#   scattered clipping: on real HG002 Illumina 300x data, three loci with
#   fraction 0.006 (max_clips=1), 0.008 (max_clips=2), and 0.027 WITH a
#   genuine 13-read pileup all sat under the old scheme's 0.1 floor and
#   scored identically — only the third was real signal by the tool's own
#   noise cutoff. max_clips-based tiers fix that specific failure, but 3
#   and 10 themselves are still chosen by judgement: the only real numbers
#   on record are HCC1143 (max_clips unrecorded, fraction 7.7%/22/287),
#   GIAB PacBio HiFi (max_clips=3, right at the new floor), and this same
#   GIAB Illumina 300x locus (max_clips=13, comfortably in the top tier) —
#   3 real data points, 2 of them the same underlying deletion on different
#   technologies, not an independently confirmed "pileup of size N should
#   score X" calibration.
#
# SPLIT-READ (split_read_fraction cutoffs: 0.1, 0.3)                   HEURISTIC
#   Same cutoff values as soft-clip used to have, same lack of calibration,
#   but weaker evidence even than that: no real split-read-positive locus
#   with known ground truth has been tested at all. Every real BAM checked
#   against split_reads so far (HCC1143, GIAB Illumina 300x/Novoalign)
#   happened to have zero SA tags anywhere in the file, so these fraction
#   cutoffs above 0 have never actually fired on real data — only on the
#   synthetic fixture, again by construction.
#
# READ-DEPTH (depth_ratio_min_to_mean cutoffs: 0.3, and
#             DEPTH_RATIO_DELETION_THRESHOLD=0.7)                       MIXED
#   The two cutoffs have different provenance:
#   - 0.7 (moderate/15pt tier boundary):                               EMPIRICAL
#     Calibrated against a real confirmed GIAB HG002 deletion
#     (chr1:115,686,862), replicated across 2 sequencing technologies —
#     PacBio HiFi and Illumina 300x both measured ratios below 0.7. The
#     exact ratio values on record from before 2026-08-11 (0.609, 0.542)
#     were computed with get_read_depth_profile's pre-FIX-1 read-count
#     implementation, not true per-base depth — see that function's
#     docstring and results/LLM_SESSION_3_BLIND.md. Re-measured after FIX
#     1+2 (true depth, focus-position-localized): 0.472 (window_size=100)
#     / 0.502 (window_size=200) at this same locus, still comfortably below
#     0.7. Still only 1 confirmed locus (not 2 independent loci), and not
#     checked against true-negative/normal-coverage regions for a
#     false-positive rate — see REAL_DATA_VALIDATION.md's "Calibration
#     update" section and this module's DEPTH_RATIO_DELETION_THRESHOLD
#     docstring for the full caveat.
#   - 0.3 (strong/25pt tier boundary):                                 HEURISTIC
#     Never actually confirmed by real data: the real ratios above land in
#     the 0.7 tier, not the 0.3 one — nothing has ever validated that 0.3
#     is where "strong" should start, because no confirmed-real locus has
#     ever scored there.
#   Since 2026-08-11 (FIX 2), a ratio below either cutoff is necessary but
#   not sufficient to score: get_read_depth_profile is now called with
#   focus_position=position, so depth_ratio_min_to_mean's numerator is
#   already restricted to bins within dip_tolerance_bp — but that alone
#   doesn't guarantee the region's actual lowest point (not just a nearby
#   one) is near the breakpoint. depth_score additionally requires
#   dip_is_at_focus (same dip_tolerance_bp radius, default 1000bp,
#   HEURISTIC — calibrated against exactly 2 real loci, see
#   get_read_depth_profile's docstring) before awarding points; otherwise
#   the dip is reported in supporting_observations as off-position, not
#   scored. This closes a real false-positive: a locus at a local peak
#   ~1500bp from an unrelated real dip elsewhere in the ±2kb window used to
#   score depth_score=15 for a deletion that wasn't there.
#
# The scale itself — 0-25 per layer, tiers of 0/7.5/15/25 (0/15/25 for
# soft-clip and depth) summing cleanly to a 0-100 composite — is also a
# HEURISTIC design choice (round numbers picked for readability), not
# derived from any statistical model, ROC analysis, or optimization
# against labeled data.
#
# Net: of the 9 distinct cutoff values used below (discordant-pair 0.2/0.5,
# soft-clip 3/10, split-read 0.1/0.3, depth 0.3/0.7, dip_tolerance_bp
# 1000), 1 (the depth layer's 0.7 threshold) is empirically grounded
# against a real confirmed-positive locus. The other 8 are unvalidated
# judgement calls. Treat evidence_score/evidence_strength as an
# interpretable, decomposed summary of what each tool found — not a
# calibrated probability of a true structural variant.

# ── Tool 6: Combined breakpoint evidence summary ──────────────────────────────

def summarize_breakpoint_evidence(
    bam_path: str,
    chromosome: str,
    position: int,
    label: str = "",
    window_bp: int = 500,
    min_mapq: int = 20,
    applicable_layers: list = None
) -> dict:
    """
    Combines discordant-pair, soft-clip, split-read, and read-depth evidence
    into a single interpretable breakpoint evidence summary. Each of the 4
    evidence layers is scored independently on a 0-25 scale: 0 / 7.5 / 15 /
    25 (4 tiers, from a fraction cutoff) for discordant-pair and split-read;
    0 / 15 / 25 (3 tiers) for soft-clip (from max_clips_at_position — how
    many reads pile up at the SAME position, not soft_clipped_fraction —
    see the provenance comment above this function) and for depth (from
    depth_ratio_min_to_mean, additionally gated on dip_is_at_focus — see
    below).

    evidence_score_raw is the direct, unweighted sum of all four component
    scores (0-100) regardless of whether a layer could possibly apply to this
    data — discordant_pair_score + soft_clip_score + split_read_score +
    depth_score always equals evidence_score_raw exactly.

    evidence_score is evidence_score_raw's sibling, normalised over only the
    layers listed in applicable_layers (all 4, by default). This matters
    because two of the four layers are structurally inapplicable to whole
    classes of real data — discordant_pairs is always 0 on unpaired
    long-read data, split_reads is always 0 on aligners that don't emit SA
    tags (e.g. Novoalign) — and scoring those as 0 while still dividing by
    100 systematically caps evidence_score_raw at 50-75 for those datasets
    regardless of how strong the applicable evidence actually is. Pass
    applicable_layers (e.g. from detect_applicable_layers) to exclude
    structurally-inapplicable layers from both the numerator and the
    denominator, so a real deletion on a Novoalign BAM (2 applicable layers:
    soft_clipped_reads, read_depth) can still reach "strong" on its own
    terms rather than being capped at "weak" by two layers that were never
    going to fire. evidence_strength and signal_layers are both derived from
    evidence_score (the normalised one), not evidence_score_raw.

    ASSESSABILITY (distinct from applicability): a layer can be structurally
    applicable to this technology but still have zero reads in its own
    window at this specific locus (e.g. a low-coverage region). That case
    used to default straight to a 0.0 fraction/ratio — which for the
    read-depth layer meant a 0/0 ratio silently read as "depth dropped to
    0% of the mean", awarding the single highest possible depth_score from
    a region with no data at all. Each layer's own tool now reports
    "assessable": false (with a "reason") when it has zero reads to work
    with, and this function:
      - never assigns that layer a component score (None, not 0.0 — "not
        scored" is a different claim from "scored 0")
      - excludes it from both evidence_score's numerator and denominator,
        exactly like a structurally-inapplicable layer
      - lists it in unassessable_layers (dict, {layer_name: reason})
      - still folds it into evidence_score_raw as a 0 contribution, since
        raw is deliberately the uncorrected, always-out-of-100 number —
        cite evidence_score for the corrected view
    If every layer this call would otherwise score is unassessable or
    inapplicable, evidence_score is None and evidence_strength is the
    distinct value "NOT ASSESSABLE" — never "none" (which means "we looked
    and found nothing", not "we couldn't look at all") and never a
    fabricated 0/100.

    This tool does not infer disease relevance — it only summarizes read-level
    structural-variant evidence at a candidate breakpoint.

    Depth profile uses a 4kb window (±2kb from position) to capture
    deletions larger than the short-read fragment size, and is queried with
    focus_position=position (see get_read_depth_profile) so
    depth_ratio_min_to_mean reflects depth near the breakpoint, not
    whatever the window's global minimum happens to be. Threshold 0.7
    (DEPTH_RATIO_DELETION_THRESHOLD) calibrated against GIAB HG002 Illumina
    300x validation.

    DEPTH LOCALIZATION: a ratio below threshold is necessary but not
    sufficient for depth_score > 0 — it also requires dip_is_at_focus
    (from the same get_read_depth_profile call) to be True. Without this,
    a real but unrelated dip elsewhere in the ±2kb window can pull the
    local ratio down even though the region's actual lowest point is
    nowhere near the queried position — this previously produced a
    concrete false positive (a locus at a local peak, ~1500bp from an
    unrelated real dip, scored depth_score=15 for a deletion that wasn't
    there). When the gate suppresses a would-be score, the dip is still
    surfaced in supporting_observations as an off-position depth feature,
    not silently dropped.

    THRESHOLD PROVENANCE — read before trusting evidence_strength as more
    than a decomposed summary (full detail in the comment block immediately
    above this function): of the 9 tier-cutoff values used across all 4
    layers, only ONE — the read-depth layer's 0.7 moderate-tier threshold —
    is empirically calibrated, and against a single confirmed real locus
    (GIAB HG002 deletion, replicated across 2 technologies, not 2
    independent loci). The other 8 — discordant-pair's 0.2/0.5, soft-clip's
    3/10 (on max_clips_at_position, not fraction — changed 2026-08-11),
    split-read's 0.1/0.3, the read-depth layer's own 0.3 strong-tier
    threshold, and dip_tolerance_bp's 1000bp localization radius — are
    HEURISTIC: chosen by judgement, never validated against real
    confirmed-positive or true-negative data beyond the 2-3 real loci noted
    in the comment block above this function. This is not a defect to work
    around silently — state it plainly in any report that cites
    evidence_strength: it reflects what fired, weighted by mostly-
    unvalidated cutoffs, not a calibrated probability.

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
        applicable_layers: Optional list of layer names to normalise
                    evidence_score over, drawn from EVIDENCE_LAYER_NAMES
                    ("discordant_pairs", "soft_clipped_reads", "split_reads",
                    "read_depth"). Defaults to all 4 (matching prior
                    behavior: evidence_score == evidence_score_raw). Get
                    this from detect_applicable_layers() rather than
                    guessing.

    Returns:
        BreakpointEvidenceSummary as dict, or error dict if underlying tools
        fail or applicable_layers contains an unrecognised name
        (error_type "invalid_parameters").
    """
    if applicable_layers is not None:
        unknown = [l for l in applicable_layers if l not in EVIDENCE_LAYER_NAMES]
        if unknown:
            return {
                "error": f"Unknown applicable_layers: {unknown}. "
                         f"Valid values: {list(EVIDENCE_LAYER_NAMES)}",
                "error_type": "invalid_parameters",
            }
        if not applicable_layers:
            return {
                "error": "applicable_layers must not be empty — pass None to use all 4 layers.",
                "error_type": "invalid_parameters",
            }
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
        window_size=200, focus_position=position
    )

    for result in (stats, disc, clips, split, depth_profile):
        if "error" in result:
            return {
                "error": result["error"],
                "error_type": result.get("error_type", "bam_access"),
                "bam_path": bam_path,
            }

    observations = []

    if stats["total_reads"] == 0:
        observations.append("No reads found in the inspection window — evidence cannot be assessed.")

    # ── Assessability: a layer with zero reads in its own window has an
    # undefined fraction/ratio, not a 0 one — see the "zero-coverage false
    # positive" note on get_read_depth_profile. Treated identically to a
    # structurally-inapplicable layer for scoring purposes: excluded from
    # both evidence_score's numerator and denominator, and never assigned
    # a component score (None, not 0.0) — a distinct "cannot be assessed"
    # state rather than a fabricated score. ──
    layer_assessable = {
        "discordant_pairs": disc["assessable"],
        "soft_clipped_reads": clips["assessable"],
        "split_reads": split["assessable"],
        "read_depth": depth_profile["summary"]["assessable"],
    }
    layer_reason = {
        "discordant_pairs": disc.get("reason"),
        "soft_clipped_reads": clips.get("reason"),
        "split_reads": split.get("reason"),
        "read_depth": depth_profile["summary"].get("reason"),
    }
    unassessable_layers = {
        layer: layer_reason[layer]
        for layer in EVIDENCE_LAYER_NAMES if not layer_assessable[layer]
    }
    for layer, reason in unassessable_layers.items():
        observations.append(f"{layer} could not be assessed: {reason}.")

    # A layer whose reads were ALL filtered out is assessed and scores zero --
    # it is not excluded. Say so explicitly, with the count, so the reason a
    # layer contributed nothing is never mistaken for an empty window.
    quality_limited_layers = {
        "discordant_pairs": disc,
        "soft_clipped_reads": clips,
        "split_reads": split,
    }
    for layer, payload in quality_limited_layers.items():
        if payload.get("quality_limited"):
            observations.append(
                f"{layer}: all {payload['reads_below_min_mapq']} reads in the "
                f"window were below MAPQ {min_mapq}; layer assessed as zero, "
                f"not excluded from the score."
            )

    # ── Discordant-pair component (0-25) ──
    if not layer_assessable["discordant_pairs"]:
        discordant_pair_score = None
    else:
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
        n_disc = disc["discordant_pairs"]
        partner_phrase = _describe_partner_distribution(disc["mate_chromosomes"], "mate")
        if n_disc == 1:
            # The fraction is omitted here deliberately: at n=1 it rounds to
            # "0% of reads in window", which reads as "no discordant pairs"
            # directly beside a sentence reporting one.
            observations.append(f"1 discordant pair {partner_phrase}.")
        else:
            observations.append(
                f"{n_disc} discordant pairs "
                f"({disc['discordant_fraction']:.0%} of reads in window) "
                f"{partner_phrase}."
            )

    # ── Soft-clip component (0-25) ──
    # Scored on max_clips_at_position — how many reads actually pile up at
    # the SAME position — not soft_clipped_fraction. fraction alone can't
    # tell a genuine pileup from scattered clipping: on real HG002 data, a
    # locus with fraction=0.006/max_clips=1, one with fraction=0.008/
    # max_clips=2, and one with fraction=0.027 AND a genuine 13-read pileup
    # all fell under the old fraction-only scheme's 0.1 tier and scored
    # identically, even though only the third was a real pileup by
    # count_soft_clipped_reads' own documented threshold (max_clips_at_
    # position < 3 = noise). Tiers below are HEURISTIC — chosen by
    # judgement (matching the depth layer's existing 0/15/25, 3-tier
    # shape), never validated against real confirmed-positive/negative
    # data; see the provenance comment above summarize_breakpoint_evidence.
    if not layer_assessable["soft_clipped_reads"]:
        soft_clip_score = None
    else:
        max_clips = clips["max_clips_at_position"]
        if max_clips >= 10:
            soft_clip_score = 25.0
        elif max_clips >= 3:
            soft_clip_score = 15.0
        else:
            soft_clip_score = 0.0

    if clips["soft_clipped_reads"] > 0:
        n_clips = clips["soft_clipped_reads"]
        max_clips_obs = clips["max_clips_at_position"]
        fraction_str = f"({clips['soft_clipped_fraction']:.0%} of reads in window)"
        if max_clips_obs >= SOFT_CLIP_PILEUP_MIN_READS:
            observations.append(
                f"{n_clips} soft-clipped read(s) {fraction_str}, "
                f"consensus clip position at {clips['consensus_clip_position']} "
                f"({max_clips_obs} reads)."
            )
        else:
            # "consensus" implies reads agreeing on a boundary; below the
            # pileup threshold they don't. Say so rather than dressing
            # scattered clipping up as a localized breakpoint signal.
            observations.append(
                f"{n_clips} soft-clipped read(s) {fraction_str}, but no clip pileup "
                f"— at most {max_clips_obs} read(s) share any single clip position, "
                f"so there is no consensus breakpoint."
            )

    # ── Split-read component (0-25) ──
    if not layer_assessable["split_reads"]:
        split_read_score = None
    else:
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
        n_split = split["split_reads"]
        partner_phrase = _describe_partner_distribution(
            split["partner_chromosomes"], "supplementary alignment"
        )
        examples = split.get("example_partner_loci") or []
        # Render the first example as a bare locus string; the previous
        # `[:1]` slice interpolated a Python list, so reports read
        # "(e.g. ['chr8:47000000'])" with brackets and quotes included.
        example_suffix = f" (e.g. {examples[0]})" if examples else ""
        if n_split == 1:
            observations.append(f"1 split read {partner_phrase}{example_suffix}.")
        else:
            observations.append(
                f"{n_split} split reads "
                f"({split['split_read_fraction']:.0%} of reads in window) "
                f"{partner_phrase}{example_suffix}."
            )

    # ── Read-depth component (0-25) ──
    # depth_ratio < DEPTH_RATIO_DELETION_THRESHOLD flags a possible deletion
    # (coverage drop) — see the module-level constant's docstring for
    # calibration provenance. depth_ratio here is already focus_position-
    # localized (get_read_depth_profile was called with focus_position=
    # position above), so it never disagrees with get_read_depth_profile's
    # own `likely_deletion` flag for this same call. 0.3 is a separate,
    # stricter cutoff for the "strong" scoring tier within this composite
    # only.
    #
    # A ratio below threshold is still gated on dip_is_at_focus before it's
    # allowed to score: depth_ratio's numerator is already restricted to
    # bins within dip_tolerance_bp of position (see get_read_depth_profile),
    # so a low ratio here means SOME nearby bin is low, but the region's
    # actual global minimum could still sit further out and dominate the
    # local mean/shape for reasons unrelated to this breakpoint. Requiring
    # dip_is_at_focus too means the region's single lowest point — not just
    # a nearby one — is within reach of the queried position before this
    # counts as evidence for THIS breakpoint. This is exactly the case that
    # produced a false positive before dip_is_at_focus existed: a locus
    # sitting at a local peak, ~1500bp from an unrelated real dip elsewhere
    # in the ±2kb window, previously scored depth_score=15 for a deletion
    # that wasn't there. When the gate suppresses a would-be score, the dip
    # is still surfaced in supporting_observations — off-position, not
    # absent.
    off_position_dip = False
    if not layer_assessable["read_depth"]:
        depth_score = None
    else:
        depth_ratio = depth_profile["summary"]["depth_ratio_min_to_mean"]
        dip_is_at_focus = depth_profile["summary"]["dip_is_at_focus"]
        if depth_ratio < 0.3:
            candidate_score = 25.0
        elif depth_ratio < DEPTH_RATIO_DELETION_THRESHOLD:
            candidate_score = 15.0
        else:
            candidate_score = 0.0

        if candidate_score > 0 and dip_is_at_focus is False:
            depth_score = 0.0
            off_position_dip = True
        else:
            depth_score = candidate_score

    if depth_score is not None and depth_score > 0:
        observations.append(
            f"Read depth drops to {depth_profile['summary']['depth_ratio_min_to_mean']:.0%} "
            f"of the window mean near the queried position (min {depth_profile['summary']['min_depth']} "
            f"vs mean {depth_profile['summary']['mean_depth']}, per-base depth) — "
            f"consistent with a possible deletion."
        )
    elif off_position_dip:
        observations.append(
            f"Read depth dips to {depth_profile['summary']['min_depth']} "
            f"(mean {depth_profile['summary']['mean_depth']}, per-base depth) at "
            f"{chromosome}:{depth_profile['summary']['dip_position']}, "
            f"{depth_profile['summary']['dip_distance_from_focus']}bp from the queried position — "
            f"an off-position depth feature, not counted as evidence for this breakpoint."
        )

    # ── Combine. evidence_score_raw is the direct sum over all 4 layers
    # (0-100), treating any unassessable/None component as a 0 contribution
    # — same convention it already used for structurally-inapplicable
    # layers, now extended to "no reads to assess" for the same reason:
    # raw is deliberately the uncorrected, always-out-of-100 number: The
    # individual component fields themselves stay None (never scored),
    # only the sum folds them in as 0. evidence_score normalises over only
    # the layers that are BOTH caller-applicable AND assessable, so it
    # isn't penalised by layers that were never going to fire, or that
    # happened to have no reads at this specific locus. ──
    scores_by_layer = {
        "discordant_pairs": discordant_pair_score,
        "soft_clipped_reads": soft_clip_score,
        "split_reads": split_read_score,
        "read_depth": depth_score,
    }
    all_component_scores = tuple(scores_by_layer[layer] for layer in EVIDENCE_LAYER_NAMES)
    evidence_score_raw = round(sum(s for s in all_component_scores if s is not None), 1)

    layers_used = list(applicable_layers) if applicable_layers is not None else list(EVIDENCE_LAYER_NAMES)
    scoreable_layers = [l for l in layers_used if layer_assessable[l]]

    if scoreable_layers:
        applicable_scores = [scores_by_layer[layer] for layer in scoreable_layers]
        applicable_max = len(scoreable_layers) * 25.0
        evidence_score = round(sum(applicable_scores) * (100.0 / applicable_max), 1)
        signal_layers = f"{sum(1 for s in applicable_scores if s > 0)}/{len(scoreable_layers)}"

        # ── Quality gate ──
        # Above LOW_MAPQ_QUALITY_GATE the underlying reads are mostly
        # multi-mapping, and a normalised score over them invites more
        # confidence than the data supports. Report the limitation instead of
        # a number. evidence_score_raw is still returned for reference.
        low_mapq_fraction = stats.get("low_mapq_fraction", 0) or 0
        if stats["total_reads"] > 0 and low_mapq_fraction > LOW_MAPQ_QUALITY_GATE:
            observations.append(
                f"Quality gate: {low_mapq_fraction:.1%} of reads at this locus "
                f"are below MAPQ 20 (threshold {LOW_MAPQ_QUALITY_GATE:.0%}). "
                f"Mapping here is mostly ambiguous, so no normalised evidence "
                f"score is reported."
            )
            evidence_score = None
            evidence_strength = "QUALITY-LIMITED"
        elif evidence_score >= 70:
            evidence_strength = "strong"
        elif evidence_score >= 40:
            evidence_strength = "moderate"
        elif evidence_score > 0:
            evidence_strength = "weak"
        else:
            evidence_strength = "none"
            if stats["total_reads"] > 0:
                observations.append("No discordant pairs, soft-clipping, split reads, or depth changes detected near this position.")
    else:
        # Every caller-applicable layer had zero reads to assess. A score
        # of 0/100 and "cannot be assessed" are different claims — this
        # must not collapse into evidence_strength == "none", which means
        # "we looked and found nothing", not "we couldn't look at all".
        evidence_score = None
        signal_layers = "0/0"
        evidence_strength = "NOT ASSESSABLE"

    # interpretation_template is built ONLY from values already computed above —
    # it adds no new facts, just restates them as one plain-language sentence.
    if evidence_strength == "QUALITY-LIMITED":
        score_line = (
            f"Score: withheld — {stats['low_mapq_fraction']:.1%} of reads here are "
            f"below MAPQ 20, above the {LOW_MAPQ_QUALITY_GATE:.0%} quality gate, so a "
            f"normalised score would overstate what the data supports (raw score over "
            f"all 4 layers was {evidence_score_raw}/100, shown for reference only). "
        )
    elif evidence_score is not None:
        score_line = (
            f"Score: {evidence_score}/100 normalised over {len(scoreable_layers)} applicable "
            f"layer(s) ({signal_layers} showing signal); raw score over all 4 layers "
            f"was {evidence_score_raw}/100. "
        )
    else:
        score_line = (
            f"Score: not assessable — no applicable layer had any reads to evaluate at this "
            f"locus (raw score over all 4 layers was {evidence_score_raw}/100, shown for "
            f"reference only, since it treats unassessable layers as 0). "
        )
    interpretation_template = (
        f"Breakpoint {label} at {chromosome}:{position}. "
        f"Evidence strength: {evidence_strength}. "
        f"{score_line}"
        f"Observations: {'; '.join(observations) if observations else 'none'}. "
        f"Technology note: discordant_pairs only valid for paired-end data; "
        f"split_reads only valid for modern-alignment BAMs with SA tags — use "
        f"detect_applicable_layers() to determine which apply to this BAM. "
        f"All values in this template come from tool outputs only."
    )

    result = BreakpointEvidenceSummary(
        label=label,
        chromosome=chromosome,
        position=position,
        evidence_score=evidence_score,
        evidence_score_raw=evidence_score_raw,
        evidence_strength=evidence_strength,
        signal_layers=signal_layers,
        applicable_layers=layers_used,
        unassessable_layers=unassessable_layers,
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


# ── Tool 10: Detect applicable evidence layers ─────────────────────────────────

def detect_applicable_layers(bam_path: str, sample_reads: int = 1000) -> dict:
    """
    Samples reads from a BAM to infer which of summarize_breakpoint_evidence's
    4 evidence layers can structurally produce signal in this data, without
    needing to already know the sequencing technology or aligner.

    soft_clipped_reads and read_depth are always applicable (any aligned BAM
    can show clipping and depth variation). discordant_pairs requires paired
    reads (false for single-molecule long reads: PacBio HiFi, ONT). split_reads
    requires at least one SA (supplementary alignment) tag anywhere in the
    sample (false for aligners that don't emit chimeric alignments, e.g. the
    2018-era Novoalign pipeline used for the HCC1143 validation BAM in this
    repo — see results/RESULTS_HCC1143.md).

    Call this once per BAM (not per locus) at the start of a session, before
    summarize_breakpoint_evidence, and pass its applicable_layers straight
    through so the composite score isn't penalised by layers that were never
    going to fire for this data.

    Args:
        bam_path:     Path to indexed BAM file
        sample_reads: Number of mapped, non-secondary/supplementary reads to
                      inspect before concluding (default 1000). Reads are
                      sampled in on-disk file order starting from the
                      beginning, not from a specific locus — this assumes a
                      technology/aligner is uniform across the whole BAM,
                      which holds for every BAM this project has validated
                      against.

    Returns:
        dict with applicable_layers (list, ready to pass to
        summarize_breakpoint_evidence), reads_sampled, and evidence (dict
        explaining each layer's determination), or a structured error dict
        ({"error", "error_type": "bam_access"}) if the BAM cannot be read.
    """
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return {"error": str(e), "error_type": "bam_access", "bam_path": bam_path}

    sampled = 0
    any_paired = False
    any_sa_tag = False
    try:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            sampled += 1
            if read.is_paired:
                any_paired = True
            if read.has_tag("SA"):
                any_sa_tag = True
            if sampled >= sample_reads:
                break
    finally:
        bam.close()

    if sampled == 0:
        # No mapped reads to sample at all — can't positively confirm or
        # rule out pairing/SA-tag support, so don't assert either way.
        # soft_clipped_reads/read_depth are kept applicable since they
        # require no reads to be structurally valid to attempt (they'd
        # just report zero signal), consistent with how those two tools
        # already behave on an empty window.
        return {
            "bam_path": bam_path,
            "reads_sampled": 0,
            "applicable_layers": ["read_depth", "soft_clipped_reads"],
            "evidence": {
                "soft_clipped_reads": "always applicable — any aligned BAM can show soft-clipping.",
                "read_depth": "always applicable — any aligned BAM can show a coverage drop.",
                "discordant_pairs": "inconclusive — no mapped reads found in the sample to check pairing.",
                "split_reads": "inconclusive — no mapped reads found in the sample to check for SA tags.",
            },
        }

    applicable_layers = ["soft_clipped_reads", "read_depth"]
    if any_paired:
        applicable_layers.append("discordant_pairs")
    if any_sa_tag:
        applicable_layers.append("split_reads")
    # Keep canonical ordering regardless of the order layers were appended above.
    applicable_layers = [l for l in EVIDENCE_LAYER_NAMES if l in applicable_layers]

    evidence = {
        "soft_clipped_reads": "always applicable — any aligned BAM can show soft-clipping.",
        "read_depth": "always applicable — any aligned BAM can show a coverage drop.",
        "discordant_pairs": (
            f"applicable — at least one paired read observed in {sampled} sampled reads."
            if any_paired else
            f"not applicable — no paired reads observed in {sampled} sampled reads "
            f"(consistent with unpaired long-read data, e.g. PacBio HiFi or ONT)."
        ),
        "split_reads": (
            f"applicable — at least one SA (supplementary alignment) tag observed "
            f"in {sampled} sampled reads."
            if any_sa_tag else
            f"not applicable — no SA tags observed in {sampled} sampled reads "
            f"(consistent with an aligner that doesn't emit chimeric alignments)."
        ),
    }

    return {
        "bam_path": bam_path,
        "reads_sampled": sampled,
        "applicable_layers": applicable_layers,
        "evidence": evidence,
    }


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
                    # Renamed from clinical_note, and reworded. A coordinate
                    # lookup cannot know a breakpoint exists, so it must not
                    # say one does: the old text asserted "Breakpoint directly
                    # disrupts N gene(s)" for ANY queried position, including
                    # arbitrary controls with no evidence of anything
                    # (REAL_PATIENT_DATA_VALIDATION.md finding 13). This states
                    # only what the Ensembl lookup established: overlap.
                    "annotation_note": (
                        f"{len(gene_list)} annotated gene(s) overlap this "
                        f"position: {', '.join(g['gene_name'] for g in gene_list)}."
                        if gene_list
                        else "No annotated gene overlaps this position "
                             "(intergenic). Nearby genes may still be relevant "
                             "to a positional effect."
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


# Discordant-pair thresholds for the reciprocal check. Named because verdict
# and is_balanced used to compute their own, and disagreed: every verdict
# branch required primary >= 5 while is_balanced used >= 3 on both sides, so
# anything in the 3-4 band returned "INSUFFICIENT EVIDENCE at both positions"
# together with is_balanced=True. Two callers reading the same object got
# opposite answers, and both fields look authoritative.
# (REAL_PATIENT_DATA_VALIDATION.md finding 3.)
RECIPROCAL_STRONG = 5      # a side counts as strong signal
RECIPROCAL_PRESENT = 2     # a partner side counts as present-but-weak
RECIPROCAL_BALANCED = 3    # both sides comparable -> is_balanced
RECIPROCAL_BACKPOINT = 3   # partner mates pointing back at the primary


def _reciprocal_verdict(primary_disc: int, reciprocal_disc: int,
                        back_pointing: int):
    """
    Single source of both the verdict and is_balanced, so they cannot
    contradict each other.

    Pure function of three counts — no BAM access — which is what lets the
    regression suite check every combination exhaustively rather than the
    handful a fixture happens to produce.

    Invariant, asserted in tests over the full grid: is_balanced is never True
    while the verdict says INSUFFICIENT. The old code had no branch for "both
    sides show weak but symmetric signal", so that whole region of the input
    space fell through to a catch-all that then contradicted is_balanced.

    Returns (verdict, is_balanced).
    """
    is_balanced = (primary_disc >= RECIPROCAL_BALANCED
                   and reciprocal_disc >= RECIPROCAL_BALANCED)

    if (primary_disc >= RECIPROCAL_STRONG
            and reciprocal_disc >= RECIPROCAL_STRONG
            and back_pointing >= RECIPROCAL_BACKPOINT):
        verdict = ("RECIPROCAL CONFIRMED — both breakpoints show concordant "
                   "inter-chromosomal signal")
    elif primary_disc >= RECIPROCAL_STRONG and reciprocal_disc >= RECIPROCAL_PRESENT:
        verdict = ("RECIPROCAL LIKELY — primary signal strong, partner signal "
                   "present but weak")
    elif is_balanced:
        # The branch the old code lacked. Both sides clear RECIPROCAL_BALANCED
        # but not RECIPROCAL_STRONG: real symmetry, too little of it to call.
        verdict = (f"RECIPROCAL POSSIBLE — both positions show weak but "
                   f"symmetric signal ({primary_disc} and {reciprocal_disc} "
                   f"discordant pairs); below the threshold for a confident call")
    elif primary_disc >= RECIPROCAL_STRONG and reciprocal_disc == 0:
        verdict = ("RECIPROCAL NOT FOUND — only one side shows signal; may be "
                   "artifact or wrong partner coords")
    elif primary_disc == 0 and reciprocal_disc == 0:
        verdict = "INSUFFICIENT EVIDENCE at both positions"
    elif primary_disc == 0:
        verdict = (f"INSUFFICIENT EVIDENCE at the primary position "
                   f"({reciprocal_disc} discordant pair(s) at the partner)")
    elif reciprocal_disc == 0:
        verdict = (f"INSUFFICIENT EVIDENCE at the partner position "
                   f"({primary_disc} discordant pair(s) at the primary)")
    else:
        verdict = (f"INSUFFICIENT EVIDENCE — {primary_disc} discordant pair(s) "
                   f"at the primary and {reciprocal_disc} at the partner, both "
                   f"below the threshold for a call")
    return verdict, is_balanced


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

    primary_chromosome/partner_chromosome may be given in either naming
    convention ("chr1" or "1") regardless of which one the BAM header
    itself uses — every chromosome comparison here (including matching a
    reciprocal read's mate back to primary_chromosome) is done on the
    normalised, chr-prefixed form of both sides, so the verdict does not
    depend on which convention the caller happened to use.

    Args:
        bam_path:           Path to indexed BAM
        primary_chromosome: First breakpoint chromosome (e.g. "chr1")
        primary_position:   First breakpoint position
        partner_chromosome: Partner chromosome (e.g. "chr8")
        partner_position:   Estimated partner breakpoint position
        window_bp:          Half-window for read counting
        min_mapq:           Minimum mapping quality

    Returns:
        dict with primary evidence, reciprocal evidence, and reciprocity
        verdict, or a structured error dict ({"error", "error_type",
        "side": "primary" | "reciprocal", ...}) if either side's
        underlying discordant-pair check failed. A tool failure (bad BAM
        path, invalid contig, etc.) is never reported as
        "INSUFFICIENT EVIDENCE" — those must stay visibly distinct.
    """
    # Check primary side
    primary = count_discordant_pairs(
        bam_path, primary_chromosome, primary_position, window_bp, min_mapq
    )
    if "error" in primary:
        return {
            "error": primary["error"],
            "error_type": primary.get("error_type", "bam_access"),
            "side": "primary",
            "chromosome": primary_chromosome,
            "position": primary_position,
        }

    # Check reciprocal side
    reciprocal = count_discordant_pairs(
        bam_path, partner_chromosome, partner_position, window_bp, min_mapq
    )
    if "error" in reciprocal:
        return {
            "error": reciprocal["error"],
            "error_type": reciprocal.get("error_type", "bam_access"),
            "side": "reciprocal",
            "chromosome": partner_chromosome,
            "position": partner_position,
        }

    # Reciprocity check: does partner side have discordant mates pointing back?
    primary_disc = primary.get("discordant_pairs", 0)
    reciprocal_disc = reciprocal.get("discordant_pairs", 0)

    # Check if reciprocal mates point back to primary chromosome.
    # mate_chromosomes' keys are already normalised (chr-prefixed) by
    # count_discordant_pairs, so primary_chromosome — whatever convention
    # the caller happened to use — must be normalised the same way before
    # this lookup, or a caller using the opposite convention to the BAM's
    # own header (e.g. "1" against a "chr1"-header BAM) silently gets
    # back_pointing=0 here even when the real count is nonzero, downgrading
    # the verdict below for no genomic reason.
    reciprocal_mate_chroms = reciprocal.get("mate_chromosomes", {})
    back_pointing = _chrom_count(reciprocal_mate_chroms, primary_chromosome)

    verdict, is_balanced = _reciprocal_verdict(
        primary_disc, reciprocal_disc, back_pointing
    )

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
        "is_balanced": is_balanced,
    }


def _signal_igv_process_group(proc, sig):
    """
    Send a signal to proc's whole process group, not just proc itself.

    igv.sh is a shell script that runs `java ...` as its last command
    without `exec`, so the shell stays alive as java's parent. Signaling
    proc.pid alone (as proc.terminate()/kill() do) only reaches that
    wrapper shell — the real IGV/java GUI process is a separate PID and
    is left running. proc is started with start_new_session=True so its
    pid is also its process group id, and killpg reaches both.
    """
    try:
        _os.killpg(_os.getpgid(proc.pid), sig)
    except ProcessLookupError:
        pass


# ── Tool 9: IGV screenshot (visual evidence) ───────────────────────────────────

# IGV's AlignmentTrack$ColorOption enum, as shipped in IGV desktop (verified
# against the installed IGV_2.17.4/lib/igv.jar via javap — this is NOT the
# same enum as AlignmentTrack$GroupOption, which is what "group by" batch
# commands use). MATE_CHROMOSOME is a GroupOption value, not a ColorOption
# value — it was never valid here, despite being a very natural name to
# reach for. UNEXPECTED_PAIR is the ColorOption that actually flags
# inter-chromosomal / anomalous-orientation / anomalous-insert-size pairs
# (labelled "insert size and pair orientation" in the IGV UI) and is the
# correct choice for translocation coloring.
VALID_COLOR_BY_OPTIONS = {
    "INSERT_SIZE", "READ_STRAND", "FIRST_OF_PAIR_STRAND", "PAIR_ORIENTATION",
    "READ_ORDER", "SAMPLE", "READ_GROUP", "LIBRARY", "MOVIE", "ZMW",
    "BISULFITE", "NOMESEQ", "TAG", "NONE", "UNEXPECTED_PAIR", "MAPPED_SIZE",
    "LINK_STRAND", "YC_TAG", "BASE_MODIFICATION", "BASE_MODIFICATION_2COLOR",
    "SMRT_SUBREAD_IPD", "SMRT_SUBREAD_PW", "SMRT_CCS_FWD_IPD",
    "SMRT_CCS_FWD_PW", "SMRT_CCS_REV_IPD", "SMRT_CCS_REV_PW",
}


def run_igv_screenshot(
    bam_paths: list,
    chromosome: str,
    start: int,
    end: int,
    output_path: str,
    genome_build: str = "hg38",
    color_by: str = "UNEXPECTED_PAIR",
    show_soft_clips: bool = True,
    max_coverage: int = None,
    coverage_height: int = 120,
    squish: bool = True,
    group_by: str = None,
    sort_by: str = "position",
    hide_alignment_tracks: bool = False,
    igv_path: str = None,
    timeout_sec: int = 180
) -> dict:
    """
    Generate an IGV screenshot of a genomic region using headless batch mode.

    Produces the visual evidence a clinician would inspect manually:
    discordant/anomalous pairs colored, soft-clipped bases shown, and the
    region centred on the candidate breakpoint.

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
        color_by:        IGV coloring mode, validated against this IGV
                         build's actual AlignmentTrack$ColorOption enum
                         before launching (an invalid value returns a
                         structured error immediately rather than hanging
                         IGV until timeout_sec). Recommended values:
                         UNEXPECTED_PAIR (translocations — flags
                         inter-chromosomal mates, and anomalous insert
                         size/orientation, in one coloring mode; this is
                         the default),
                         PAIR_ORIENTATION (inversions),
                         INSERT_SIZE (deletions/duplications),
                         NONE.
                         There is no "MATE_CHROMOSOME" ColorOption in this
                         IGV build — that name exists only as a *grouping*
                         option (AlignmentTrack$GroupOption), not a
                         coloring one. See VALID_COLOR_BY_OPTIONS in this
                         module for the full list this build supports.
        show_soft_clips: Display soft-clipped bases
        max_coverage:    If set, fixes the coverage track's max value (via
                         IGV's setDataRange) instead of autoscaling to the
                         tallest window in view. Set this slightly above the
                         observed max_depth from read_depth_profile so a
                         depth dip elsewhere in the region isn't flattened
                         by autoscaling to a taller peak.
        coverage_height: Coverage track height in pixels (IGV default is ~50px,
                         too short to render a depth dip legibly).
        squish:          Compress read rows into a denser view when True
        group_by:        IGV batch `group` argument, verified against this
                         build's AlignmentTrack$GroupOption menu (a
                         *different* enum from ColorOption — e.g. "TAG SA"
                         groups reads by their SA tag's value, putting every
                         split/chimeric read into its own labeled row group
                         showing the exact partner locus; not validated
                         against a fixed list the way color_by is, since
                         GroupOption also accepts free-form tag names).
                         None (default) skips the `group` command.
        sort_by:         IGV batch `sort` argument (default "position").
                         Verified values for this build include "base"
                         (SortOption.NUCLEOTIDE — sorts by the nucleotide at
                         the center line, useful for SNP/mismatch patterns,
                         less so for clip evidence specifically),
                         "insertSize", "strand", and others — see
                         org.broad.igv.sam.SortOption in igv.jar for the
                         full enum if a value here doesn't behave as
                         expected; unlike color_by this isn't validated
                         against a hardcoded list before launch.
        hide_alignment_tracks: If True, removes each bam_paths track's
                         alignment rows after loading (IGV batch `remove
                         <basename>`), leaving only its coverage sub-track
                         visible — confirmed empirically that IGV's
                         SAM.SHOW_ALIGNMENT_TRACK preference does NOT do
                         this despite the name (tested directly: reads stay
                         fully rendered), so `remove` is used instead.
        igv_path:        Path to igv.sh (auto-detected if None). Auto-detection
                         checks the IGV_PATH environment variable first, then
                         falls back to the hardcoded candidates below.
        timeout_sec:     Max seconds to wait for IGV

    Returns:
        dict with screenshot_path, batch_script used, and success status
    """
    # Validate inputs before doing anything else — both of these previously
    # only surfaced after IGV was launched and either crashed (bad
    # color_by) or sat with nothing loaded until timeout_sec (empty
    # bam_paths), confirmed directly: an empty bam_paths list produced a
    # batch script with zero `load` lines and burned the full timeout
    # before returning an error.
    if not bam_paths:
        return {
            "error": "bam_paths is empty — at least one BAM path or URL is required.",
            "error_type": "invalid_parameters",
        }
    if color_by and color_by not in VALID_COLOR_BY_OPTIONS:
        return {
            "error": f"'{color_by}' is not a valid color_by option for this IGV build.",
            "error_type": "invalid_parameters",
            "color_by": color_by,
            "valid_options": sorted(VALID_COLOR_BY_OPTIONS),
        }

    # Auto-detect IGV. IGV_PATH, if set, is checked before the hardcoded
    # fallback locations.
    candidates = []
    if igv_path is None:
        env_igv_path = _os.environ.get("IGV_PATH")
        candidates = ([env_igv_path] if env_igv_path else []) + [
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

    # Remove any pre-existing file at output_path before launching IGV.
    # Confirmed directly as a real bug, not a hypothetical: the completion
    # check below polls for output_path's size to be stable across two 1s
    # checks, then kills IGV. If a file already exists there (e.g. this
    # same locus was screenshotted in an earlier session), that condition
    # is trivially true from the first poll -- IGV gets SIGTERM'd before it
    # has done more than start the JVM, and the stale file's size is
    # reported back as if it were a fresh snapshot. On a remote-BAM-over-
    # HTTPS load this cuts a ~2-minute genuine render down to ~1 second
    # while still returning success=True. Deleting first means the file's
    # mere existence during polling is unambiguous proof of a fresh write.
    if _os.path.exists(output_path):
        _os.remove(output_path)

    # Build IGV batch script
    lines = ["new", f"genome {genome_build}"]
    for bam in bam_paths:
        lines.append(f"load {bam}")
    if hide_alignment_tracks:
        # `remove <track name>` deletes the alignment track but leaves its
        # "<basename> Coverage" sub-track untouched — confirmed empirically
        # this is the only way to get a coverage-only view (see docstring).
        for bam in bam_paths:
            lines.append(f"remove {_os.path.basename(bam)}")
    lines.append(f"preference SAM.COVERAGE_TRACK_HEIGHT {coverage_height}")
    lines.append("maxPanelHeight 600")
    if color_by:
        # Always issue this explicitly, including for "NONE", so the batch
        # script never silently depends on whatever colorBy state a
        # previous invocation happened to leave behind.
        #
        # IMPORTANT, confirmed empirically (not assumed): colorBy="NONE"
        # does NOT produce neutral/uncolored reads. Tested with a brand-new
        # BAM, an explicit track-name argument, and multiple colorBy values
        # (NONE, READ_STRAND) side by side — anomalous pairs (inter-
        # chromosomal mate, improper orientation) still render in their
        # characteristic red/orange/blue regardless of colorBy. This
        # appears to be baseline IGV behavior independent of the colorBy
        # selection, not a bug in this tool or evidence of a stale
        # preference — colorBy controls an *additional* coloring scheme on
        # top of that baseline, not a replacement for it. Anything relying
        # on a "clean, uncolored" screenshot of anomalous-pair data should
        # not assume color_by="NONE" achieves that.
        lines.append(f"colorBy {color_by}")
    if group_by:
        lines.append(f"group {group_by}")
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
    if sort_by:
        lines.append(f"sort {sort_by}")
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

    stdout_path = None
    stderr_path = None
    try:
        # No env= override here: this inherits the parent process's DISPLAY.
        # Passing DISPLAY="" (rather than leaving it unset or omitting the
        # override entirely) crashes IGV's AWT EventDispatchThread before
        # the batch script runs, so no snapshot is ever produced — confirmed
        # directly against this IGV build, not assumed.
        #
        # stdout/stderr go to temp files rather than PIPE: IGV logs steadily
        # while we poll below, and an unread PIPE fills its OS buffer and
        # blocks IGV's write — a second hang on top of the one this whole
        # poll loop exists to work around.
        stdout_f = _tempfile.NamedTemporaryFile(mode="w+", suffix=".igvout",
                                                  delete=False)
        stderr_f = _tempfile.NamedTemporaryFile(mode="w+", suffix=".igverr",
                                                  delete=False)
        stdout_path, stderr_path = stdout_f.name, stderr_f.name

        proc = _subprocess.Popen([igv_path, "--batch", batch_path],
                                  stdout=stdout_f, stderr=stderr_f,
                                  start_new_session=True)

        # IGV's `exit` batch command reliably writes the snapshot first but
        # has been observed to hang the AWT/JVM shutdown afterward, so we
        # don't wait for the process to exit on its own — we confirm the
        # PNG is fully written (two size checks, 1s apart, agreeing on a
        # non-zero size) and then terminate IGV ourselves.
        deadline = _time.time() + timeout_sec
        shutdown_method = None
        last_size = None

        while True:
            if proc.poll() is not None:
                shutdown_method = "clean_exit"
                break

            if _time.time() >= deadline:
                _signal_igv_process_group(proc, _signal.SIGTERM)
                try:
                    proc.wait(timeout=5)
                except _subprocess.TimeoutExpired:
                    _signal_igv_process_group(proc, _signal.SIGKILL)
                    try:
                        proc.wait(timeout=5)
                    except _subprocess.TimeoutExpired:
                        pass
                shutdown_method = "timeout"
                break

            if _os.path.exists(output_path):
                size = _os.path.getsize(output_path)
                if size > 0 and size == last_size:
                    _signal_igv_process_group(proc, _signal.SIGTERM)
                    try:
                        proc.wait(timeout=5)
                    except _subprocess.TimeoutExpired:
                        _signal_igv_process_group(proc, _signal.SIGKILL)
                        try:
                            proc.wait(timeout=5)
                        except _subprocess.TimeoutExpired:
                            pass
                    shutdown_method = "terminated_after_snapshot"
                    break
                last_size = size if size > 0 else None
                _time.sleep(1.0)
                continue

            _time.sleep(0.5)

        stdout_f.flush(); stderr_f.flush()
        stdout_f.seek(0); stderr_f.seek(0)
        stdout_text = stdout_f.read()
        stderr_text = stderr_f.read()
        stdout_f.close()
        stderr_f.close()

        success = _os.path.exists(output_path) and _os.path.getsize(output_path) > 0
        file_size = _os.path.getsize(output_path) if success else 0

        result = {
            "success": success,
            "screenshot_path": _os.path.abspath(output_path) if success else None,
            "file_size_bytes": file_size,
            "region": f"{chromosome}:{start}-{end}",
            "color_by": color_by,
            "bam_tracks": len(bam_paths),
            "batch_script": batch_content,
            "shutdown_method": shutdown_method,
            "igv_stdout": stdout_text[-500:] if stdout_text else "",
            "igv_stderr": stderr_text[-500:] if stderr_text else "",
            # The caller (an LLM agent) never receives the PNG's pixel data
            # through the MCP tool-call mechanism -- only this dict, as text.
            # Stated explicitly rather than left implicit, after a session
            # observed both claude-sonnet-5 and qwen2.5:7b write confident
            # descriptions of image contents neither had been shown (see
            # results/LLM_SESSION_4_VISUAL_*.md). evidence_panel's per-layer
            # panels embed this same dict, so this propagates there too.
            "image_content_available_to_caller": False,
            "note": (
                "This tool returns a file path only. The image itself has "
                "not been provided to you. You cannot describe its visual "
                "contents. Report that an image was generated and where, "
                "and state that visual interpretation requires a human "
                "reviewer or a vision-capable client."
            ),
        }
        if not success:
            result["error"] = (f"IGV timed out after {timeout_sec}s"
                                if shutdown_method == "timeout" else
                                "IGV exited without producing a screenshot")
        return result
    except Exception as e:
        return {"error": str(e), "batch_script": batch_content}
    finally:
        if _os.path.exists(batch_path):
            _os.unlink(batch_path)
        for p in (stdout_path, stderr_path):
            if p and _os.path.exists(p):
                _os.unlink(p)


# ── Image handles: keeping filesystem paths away from the model ─────────────
#
# The MCP tool layer (server.py) exposes NO path parameter and returns NO
# path: it assigns the output location itself and hands back an opaque
# handle. These helpers implement that.
#
# Why the tool signature had to change rather than just redacting the
# response: an earlier fix redacted screenshot_path from tool *results*, and
# qwen2.5:7b promptly cited "/tmp/igv_screenshot.png" in its report anyway --
# a path it had supplied itself as the output_path argument, which
# necessarily stays in the conversation history. You cannot redact away a
# path the model chose. Removing the parameter is the only way to make the
# path genuinely unavailable to it. See
# results/LLM_SESSION_4_VISUAL_qwen2.5-7b.md.
IMAGE_SESSION_DIR_ENV = "IGV_IMAGE_SESSION_DIR"
IMAGE_MANIFEST_NAME = "manifest.json"

# Fields that carry a filesystem path or an IGV process log (which also
# contains paths). Stripped from every model-visible screenshot result.
_PATH_BEARING_KEYS = ("screenshot_path", "batch_script", "igv_stdout", "igv_stderr")

IMAGE_HANDLE_NOTE = (
    "An image was generated but has NOT been provided to you: you have an "
    "opaque reference (image_ref) only -- no file path, and no pixel data. "
    "You cannot see this image and must not describe what it shows, "
    "contains, or looks like. Report that it was generated, which layer and "
    "region it covers, and what a human reviewer should check."
)


def image_handle(path: str) -> str:
    """Stable, non-path-shaped reference for one image file."""
    return "IMG_" + _hashlib.sha256(_os.path.abspath(path).encode("utf-8")).hexdigest()[:4]


def png_dimensions(path: str) -> Optional[str]:
    """"WxH" from the PNG IHDR chunk, or None. Describes the file, not its
    depicted contents, so it is safe to expose."""
    try:
        with open(path, "rb") as fh:
            header = fh.read(24)
        if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
            return None
        return (f"{int.from_bytes(header[16:20], 'big')}"
                f"x{int.from_bytes(header[20:24], 'big')}")
    except OSError:
        return None


def image_session_dir() -> str:
    """
    Directory this server process writes screenshots into. Honours
    IGV_IMAGE_SESSION_DIR (set by the benchmark harness so it knows where to
    read the handle manifest back from); otherwise falls back to a
    per-process directory under the repo's screenshots/sessions/.
    """
    configured = _os.environ.get(IMAGE_SESSION_DIR_ENV)
    if configured:
        base = configured
    else:
        base = _os.path.join(
            _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
            "screenshots", "sessions", f"session_{_os.getpid()}",
        )
    _os.makedirs(base, exist_ok=True)
    return base


def _record_handles(session_dir: str, mapping: dict) -> None:
    """Append to the session's handle->path manifest, so humans (and the
    harness) can resolve a handle back to a real file. The model never sees
    this file."""
    if not mapping:
        return
    manifest_path = _os.path.join(session_dir, IMAGE_MANIFEST_NAME)
    existing = {}
    if _os.path.exists(manifest_path):
        try:
            with open(manifest_path) as fh:
                existing = _json.load(fh)
        except (OSError, ValueError):
            existing = {}
    existing.update(mapping)
    with open(manifest_path, "w") as fh:
        _json.dump(existing, fh, indent=2)


def to_handle_result(result: dict, session_dir: str) -> dict:
    """
    Convert one run_igv_screenshot result into its model-visible form:
    every path-bearing field removed, an image_ref added, dimensions added.
    Records the handle->path mapping in the session manifest as a side
    effect. Results that produced no screenshot (errors, skips) pass through
    with path-bearing fields stripped but no handle.
    """
    if not isinstance(result, dict):
        return result
    out = {k: v for k, v in result.items() if k not in _PATH_BEARING_KEYS}
    path = result.get("screenshot_path")
    if isinstance(path, str) and path:
        handle = image_handle(path)
        _record_handles(session_dir, {handle: _os.path.abspath(path)})
        out["image_ref"] = handle
        dims = png_dimensions(path)
        if dims:
            out["image_dimensions"] = dims
    out["image_content_available_to_caller"] = False
    out["note"] = IMAGE_HANDLE_NOTE
    return out


# ── Tool 11: Per-layer visual evidence panel ────────────────────────────────

# Each entry: the run_igv_screenshot kwargs that isolate that one evidence
# layer visually, verified empirically against the installed IGV build
# (not assumed from option names — see the individual comments below and
# AUDIT_2026_08.md's follow-up on the color_by bug for why that verification
# step matters here specifically).
_PANEL_LAYER_SETTINGS = {
    "discordant_pairs": dict(
        color_by="UNEXPECTED_PAIR",   # confirmed valid ColorOption (Group 2)
        show_soft_clips=False,
    ),
    "split_reads": dict(
        # color_by is deliberately omitted (left as None, not "NONE") for
        # this layer: confirmed empirically that colorBy has NO effect at
        # all once `group TAG SA` is active — a screenshot with
        # color_by="NONE" and one with color_by="READ_STRAND" rendered
        # byte-for-byte the same colors (still pair-anomaly-style
        # red/orange/blue) under this grouping mode. This is a genuine IGV
        # behavior, not a bug in this tool — grouping-by-tag appears to
        # override colorBy rather than the reverse. Left unset rather than
        # forced to "NONE" so the batch script doesn't claim a setting that
        # has no effect.
        show_soft_clips=False,
        # Groups reads by their SA tag's literal value — confirmed
        # empirically this puts every SA-tagged (chimeric) read into its
        # own labeled row group showing the exact partner locus, visually
        # separating split-read evidence from ordinary pairs with no
        # equivalent color-by option (MATE_CHROMOSOME is a GroupOption,
        # not a ColorOption — see run_igv_screenshot's docstring).
        group_by="TAG SA",
        # squish (the default everywhere else) hides each group row's text
        # label — confirmed empirically the "chr8,47000000,+,..." partner-
        # locus label per group only renders in expanded mode. That label
        # is the single most informative part of this layer, so squish is
        # disabled here specifically even though it means fewer reads fit
        # in view.
        squish=False,
    ),
    "read_depth": dict(
        color_by="NONE",
        show_soft_clips=False,
        # SAM.SHOW_ALIGNMENT_TRACK=false does NOT hide the alignment rows
        # despite the name — confirmed empirically (reads stayed fully
        # rendered). `remove <track>` is the verified way to get a
        # coverage-only view; hide_alignment_tracks wraps that.
        hide_alignment_tracks=True,
        coverage_height=200,
    ),
    "soft_clipped_reads": dict(
        color_by="NONE",
        show_soft_clips=True,
        # sort_by is set dynamically per-call (LEFT_CLIP or RIGHT_CLIP,
        # chosen from count_soft_clipped_reads' dominant_clip_side) rather
        # than fixed here — see igv_evidence_panel. Confirmed empirically
        # against real GIAB data that sort=LEFT_CLIP/RIGHT_CLIP genuinely
        # reorders reads into a clean staircase pileup at the clip
        # boundary, unlike sort=base (SortOption.NUCLEOTIDE), which sorts
        # by the nucleotide at the center line — a SNP/mismatch view, not
        # a clip-length one, and was the original (less effective) choice
        # here before this was corrected.
    ),
}

# Per-layer window half-widths (position ± this many bp), used when the
# caller doesn't override via igv_evidence_panel's windows= parameter.
# discordant_pairs/split_reads need enough width to catch mate-pair/SA
# evidence a few hundred bp out; soft_clipped_reads is kept tight so the
# clip pileup itself isn't diluted by unrelated reads; read_depth is wide
# enough to show a multi-kb deletion/duplication span in context (and is
# overridden entirely by caller-supplied start/end when given, since a
# known SV span from a VCF is more accurate than a fixed guess).
DEFAULT_PANEL_WINDOWS = {
    "discordant_pairs": 1500,
    "split_reads": 1500,
    "soft_clipped_reads": 150,
    "read_depth": 3000,
}


def igv_evidence_panel(
    bam_paths: list,
    chromosome: str,
    position: int,
    output_dir: str,
    start: int = None,
    end: int = None,
    applicable_layers: list = None,
    windows: dict = None,
    genome_build: str = "hg38",
    igv_path: str = None,
    timeout_sec: int = 180,
) -> dict:
    """
    Generates one IGV screenshot per informative evidence layer, instead of
    a single image trying to show everything at once — each layer gets the
    coloring/grouping/sorting that actually isolates it visually (see
    _PANEL_LAYER_SETTINGS), so e.g. split-read evidence (invisible under
    plain pair coloring, since IGV colors by the primary alignment's actual
    mate, never by SA tags) gets its own SA-tag-grouped view instead of
    being silently absent from a single combined screenshot.

    Each layer also gets its OWN window around `position`, rather than one
    shared region for all four — confirmed directly (not assumed) that a
    single window can't serve every layer: a region wide enough to show a
    multi-kb deletion's depth dip renders individual soft-clip marks
    illegibly small, while a region tight enough for a clean clip pileup
    would crop a wide depth dip entirely. Defaults (see
    DEFAULT_PANEL_WINDOWS): discordant_pairs/split_reads ±1500bp,
    soft_clipped_reads ±150bp (tight, to keep the pileup uncluttered),
    read_depth uses the caller-supplied start/end if given (e.g. a known SV
    span from a VCF) or position±3000bp otherwise. Override any of these
    via windows={"soft_clipped_reads": 300, ...} (half-width in bp; keys
    from EVIDENCE_LAYER_NAMES). The exact region used for each layer is
    recorded in the returned "windows_used" dict (and echoed in each
    panel's own "region" field), so it can be cited directly rather than
    re-derived.

    The soft_clipped_reads layer's sort order is chosen dynamically, not
    fixed: calls count_soft_clipped_reads over that layer's own window and
    sorts by RIGHT_CLIP or LEFT_CLIP according to whichever side has more
    clipped reads (dominant_clip_side) — confirmed empirically against
    real data that this produces a clean staircase pileup at the clip
    boundary, unlike sorting by nucleotide/base. If the side can't be
    determined (tied, no clips found, or the lookup itself errors),
    defaults to LEFT_CLIP and records why in the panel's
    "clip_side_determination" field.

    If applicable_layers isn't given, calls detect_applicable_layers on
    bam_paths[0] first (assumes all tracks share one sequencing technology —
    true for every BAM this project has validated against) and uses its
    result. Layers found inapplicable are not screenshotted at all — each
    gets a {"skipped": True, "reason": ...} entry instead, using the same
    plain-language reason detect_applicable_layers already produces, so the
    caller (or the LLM) knows *why* a PNG is missing rather than just that
    it is.

    For the read_depth layer, automatically computes a fixed coverage-track
    scale from get_read_depth_profile's observed max_depth (+15% headroom)
    for that layer's own region, rather than requiring the caller to
    already know a good max_coverage value — this is the same reasoning
    run_igv_screenshot's own docstring already gives for setting
    max_coverage manually, done here automatically.

    Args:
        bam_paths:         List of BAM file paths or URLs (same technology
                           assumed across all of them — see above)
        chromosome:        Chromosome of the region
        position:          Candidate breakpoint — the center every layer's
                           window is built around (except read_depth when
                           start/end are both given)
        output_dir:        Directory to write "<layer>.png" files into
                           (created if it doesn't exist)
        start, end:        Optional explicit region, used ONLY by the
                           read_depth layer (e.g. a known deletion span
                           from a VCF); other layers always use
                           position ± their window regardless of these.
        applicable_layers: Optional list from EVIDENCE_LAYER_NAMES
                           ("discordant_pairs", "soft_clipped_reads",
                           "split_reads", "read_depth"). None (default)
                           auto-detects via detect_applicable_layers.
        windows:           Optional dict overriding DEFAULT_PANEL_WINDOWS'
                           half-widths per layer, e.g.
                           {"soft_clipped_reads": 300}. Keys must be from
                           EVIDENCE_LAYER_NAMES.
        genome_build:      IGV genome identifier, passed to every screenshot
        igv_path:          Path to igv.sh (auto-detected if None)
        timeout_sec:       Per-screenshot IGV timeout

    Returns:
        {
          "position": position, "chromosome": chromosome,
          "bam_paths": [...],
          "applicable_layers": [...],
          "applicable_layers_source": "detect_applicable_layers" | "caller-provided",
          "windows_used": {
            "<layer>": {"start": int, "end": int, "window_bp": int | None,
                        "source": "default" | "override" | "caller-supplied start/end"}
            for each of the 4 layers, regardless of skip status
          },
          "panels": {
            "<layer>": <run_igv_screenshot result dict, plus "window_bp"
                        and, for soft_clipped_reads, "clip_side_determination">
                       | {"skipped": True, "reason": "..."}  # if not applicable
            for each of the 4 layers
          }
        }
        or a structured error dict if detect_applicable_layers, an
        unrecognised applicable_layers value, or an unrecognised windows
        key fails first.
    """
    if windows is not None:
        unknown_window_keys = [l for l in windows if l not in EVIDENCE_LAYER_NAMES]
        if unknown_window_keys:
            return {
                "error": f"Unknown windows keys: {unknown_window_keys}. "
                         f"Valid values: {list(EVIDENCE_LAYER_NAMES)}",
                "error_type": "invalid_parameters",
            }

    if applicable_layers is None:
        detected = detect_applicable_layers(bam_paths[0])
        if "error" in detected:
            return detected
        applicable_layers = detected["applicable_layers"]
        applicable_layers_source = "detect_applicable_layers"
        layer_evidence = detected["evidence"]
    else:
        unknown = [l for l in applicable_layers if l not in EVIDENCE_LAYER_NAMES]
        if unknown:
            return {
                "error": f"Unknown applicable_layers: {unknown}. "
                         f"Valid values: {list(EVIDENCE_LAYER_NAMES)}",
                "error_type": "invalid_parameters",
            }
        applicable_layers_source = "caller-provided"
        layer_evidence = {}

    _os.makedirs(output_dir, exist_ok=True)

    # ── Compute each layer's own region ──
    windows_used = {}
    for layer in EVIDENCE_LAYER_NAMES:
        if layer == "read_depth" and start is not None and end is not None:
            windows_used[layer] = {
                "start": start, "end": end, "window_bp": None,
                "source": "caller-supplied start/end",
            }
            continue
        override = windows.get(layer) if windows else None
        half = override if override is not None else DEFAULT_PANEL_WINDOWS[layer]
        windows_used[layer] = {
            "start": max(0, position - half), "end": position + half,
            "window_bp": half,
            "source": "override" if override is not None else "default",
        }

    depth_max_coverage = None
    if "read_depth" in applicable_layers:
        dw = windows_used["read_depth"]
        depth_profile = get_read_depth_profile(bam_paths[0], chromosome, dw["start"], dw["end"])
        if "error" not in depth_profile:
            depth_max_coverage = int(depth_profile["summary"]["max_depth"] * 1.15) + 1
        # If this errors, fall through without a fixed scale rather than
        # failing the whole panel — the read_depth screenshot below will
        # just autoscale instead.

    soft_clip_sort_by = "LEFT_CLIP"
    clip_side_determination = None
    if "soft_clipped_reads" in applicable_layers:
        cw = windows_used["soft_clipped_reads"]
        clip_check = count_soft_clipped_reads(
            bam_paths[0], chromosome, position, window_bp=cw["window_bp"]
        )
        if "error" in clip_check:
            clip_side_determination = (
                f"Could not determine dominant clip side ({clip_check['error']}); "
                f"defaulted to LEFT_CLIP."
            )
        else:
            side = clip_check.get("dominant_clip_side")
            if side == "right":
                soft_clip_sort_by = "RIGHT_CLIP"
                clip_side_determination = (
                    f"Sorted by RIGHT_CLIP: {clip_check['right_clip_reads']} right-clipped "
                    f"vs {clip_check['left_clip_reads']} left-clipped reads observed "
                    f"in this window."
                )
            elif side == "left":
                soft_clip_sort_by = "LEFT_CLIP"
                clip_side_determination = (
                    f"Sorted by LEFT_CLIP: {clip_check['left_clip_reads']} left-clipped "
                    f"vs {clip_check['right_clip_reads']} right-clipped reads observed "
                    f"in this window."
                )
            else:
                clip_side_determination = (
                    f"Could not determine a dominant clip side "
                    f"(dominant_clip_side={side!r}: "
                    f"{clip_check['left_clip_reads']} left vs "
                    f"{clip_check['right_clip_reads']} right-clipped reads); "
                    f"defaulted to LEFT_CLIP."
                )

    panels = {}
    for layer in EVIDENCE_LAYER_NAMES:
        if layer not in applicable_layers:
            panels[layer] = {
                "skipped": True,
                "reason": layer_evidence.get(
                    layer, "excluded by caller-provided applicable_layers"
                ),
            }
            continue

        layer_kwargs = dict(_PANEL_LAYER_SETTINGS[layer])
        if layer == "read_depth" and depth_max_coverage is not None:
            layer_kwargs["max_coverage"] = depth_max_coverage
        if layer == "soft_clipped_reads":
            layer_kwargs["sort_by"] = soft_clip_sort_by

        lw = windows_used[layer]
        output_path = _os.path.join(output_dir, f"{layer}.png")
        result = run_igv_screenshot(
            bam_paths, chromosome, lw["start"], lw["end"], output_path,
            genome_build=genome_build, igv_path=igv_path, timeout_sec=timeout_sec,
            **layer_kwargs,
        )
        result["window_bp"] = lw["window_bp"]
        if layer == "soft_clipped_reads":
            result["clip_side_determination"] = clip_side_determination
        panels[layer] = result

    return {
        "position": position,
        "chromosome": chromosome,
        "bam_paths": bam_paths,
        "applicable_layers": applicable_layers,
        "applicable_layers_source": applicable_layers_source,
        "windows_used": windows_used,
        "panels": panels,
    }