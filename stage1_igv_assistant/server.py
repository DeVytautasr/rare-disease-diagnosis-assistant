"""
server.py
FastMCP server exposing 6 BAM inspection tools to an LLM.
Anti-hallucination design: the LLM reads only tool output,
never adds genomic facts from its own training data.

Run: python stage1_igv_assistant/server.py
"""

from fastmcp import FastMCP
from stage1_igv_assistant.tools.bam_tools import (
    get_bam_stats_at_locus,
    count_discordant_pairs,
    count_soft_clipped_reads,
    get_split_reads,
    get_read_depth_profile,
    summarize_breakpoint_evidence,
    get_gene_at_locus,
    check_reciprocal_breakpoint,
)

mcp = FastMCP(
    name="IGV Breakpoint Assistant",
    instructions="""
You are a structural variant breakpoint inspection assistant.

RULES:
1. Call bam_stats_at_locus FIRST to check data quality.
2. Call tools in order: stats → discordant_pairs → soft_clips → split_reads → depth_profile → summarize.
3. Your final report must cite ONLY values the tools returned in this session.
4. Do NOT use prior knowledge about genes, cell lines, or variants.
5. If all evidence types return 0, state that clearly. Do not invent signal.
6. State the sequencing technology at the start — it determines which evidence layers apply.
7. For balanced translocations: flat depth is EXPECTED. Do not interpret it as negative evidence.
""",
)

@mcp.tool()
def bam_stats_at_locus(bam_path: str, chromosome: str, start: int, end: int) -> dict:
    """
    Basic quality check for a genomic locus. ALWAYS call this first.
    Returns depth, mean MAPQ, low-MAPQ fraction, strand balance.
    Low-MAPQ fraction > 0.4 means repetitive region — interpret all other signals with caution.
    """
    return get_bam_stats_at_locus(bam_path, chromosome, start, end)

@mcp.tool()
def discordant_pairs(bam_path: str, chromosome: str, position: int,
                     window_bp: int = 500, min_mapq: int = 20) -> dict:
    """
    Count reads whose mates map to a different chromosome.
    PRIMARY signal for balanced translocations in short-read paired-end data.
    NOT applicable to unpaired long reads (PacBio HiFi, Oxford Nanopore).
    Real translocation = many discordant pairs clustering on ONE partner chromosome.
    Mates scattered across many different chromosomes = background noise.
    """
    return count_discordant_pairs(bam_path, chromosome, position, window_bp, min_mapq)

@mcp.tool()
def soft_clipped_reads(bam_path: str, chromosome: str, position: int,
                       window_bp: int = 200, min_clip_bases: int = 10,
                       min_mapq: int = 20) -> dict:
    """
    Count reads with soft-clipped overhangs near a breakpoint.
    A pileup at the SAME position (consensus_clip_position) narrows the breakpoint precisely.
    max_clips_at_position < 3 = no real pileup, treat as noise.
    """
    return count_soft_clipped_reads(bam_path, chromosome, position, window_bp, min_clip_bases, min_mapq)

@mcp.tool()
def split_reads(bam_path: str, chromosome: str, position: int,
                window_bp: int = 300, min_mapq: int = 20) -> dict:
    """
    Find reads with supplementary SA tags spanning a breakpoint junction.
    Most direct evidence of a structural variant. Works best with long reads
    and modern BWA-MEM alignments. If the whole BAM has zero SA tags
    (2018-era pipelines), this tool cannot contribute regardless of locus.
    Partner positions in SA tags reveal the other side of the breakpoint.
    """
    return get_split_reads(bam_path, chromosome, position, window_bp, min_mapq)

@mcp.tool()
def read_depth_profile(bam_path: str, chromosome: str, start: int,
                       end: int, window_size: int = 100) -> dict:
    """
    Sliding-window read depth across a region.
    Deletions: ~50% depth drop inside deleted region (heterozygous).
    Duplications: depth rises.
    Balanced translocations: depth stays FLAT — this is expected and correct.
    depth_ratio_min_to_mean < 0.6 suggests deletion.
    """
    return get_read_depth_profile(bam_path, chromosome, start, end, window_size)

@mcp.tool()
def breakpoint_evidence_summary(bam_path: str, chromosome: str,
                                 position: int, label: str = "candidate_BP") -> dict:
    """
    Integrates all 4 evidence layers into one structured report. Call this LAST.
    Returns: evidence_strength, evidence_score (0-100), signal_layers (N/4),
    supporting_observations, and interpretation_template.
    The interpretation_template contains ONLY facts from tool outputs.
    """
    return summarize_breakpoint_evidence(bam_path, chromosome, position, label)

@mcp.tool()
def gene_at_locus(chromosome: str, position: int, genome_build: str = "GRCh38") -> dict:
    """
    Look up which gene (if any) is at a breakpoint position using Ensembl.
    This answers the key clinical question: does the breakpoint disrupt a gene?
    Call this after finding strong discordant/split-read evidence.
    Returns gene name, biotype (protein_coding / lncRNA / etc), strand,
    and whether the position is intergenic.
    Requires internet access to query the Ensembl REST API.
    """
    return get_gene_at_locus(chromosome, position, genome_build)

@mcp.tool()
def reciprocal_breakpoint(bam_path: str, primary_chromosome: str, primary_position: int,
                          partner_chromosome: str, partner_position: int,
                          window_bp: int = 500, min_mapq: int = 20) -> dict:
    """
    Verify both sides of a suspected balanced translocation.
    After finding discordant pairs at one breakpoint pointing to a partner chromosome,
    call this to check whether the partner location shows reciprocal signal back.
    True balanced translocation: both sides show inter-chromosomal discordant pairs.
    One-sided signal: likely artifact or wrong partner coordinates.
    Provide estimated partner position from the mate_chromosomes output of discordant_pairs.
    """
    return check_reciprocal_breakpoint(
        bam_path, primary_chromosome, primary_position,
        partner_chromosome, partner_position, window_bp, min_mapq
    )

if __name__ == "__main__":
    mcp.run()
