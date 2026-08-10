"""
server.py
FastMCP server exposing 11 BAM/breakpoint inspection tools to an LLM.
Anti-hallucination design: the LLM reads only tool output,
never adds genomic facts from its own training data.

Run from repo root: python -m stage1_igv_assistant.server
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
    run_igv_screenshot,
    detect_applicable_layers,
    igv_evidence_panel,
)

mcp = FastMCP(
    name="IGV Breakpoint Assistant",
    instructions="""
You are a structural variant breakpoint inspection assistant.

RULES:
1. Call applicable_layers FIRST, once per BAM, to learn which of the
   4 evidence layers can structurally produce signal for this data (pairing,
   SA-tag support). Pass its applicable_layers field straight through to
   every breakpoint_evidence_summary call at that BAM's loci.
2. Call bam_stats_at_locus next, per locus, to check data quality.
3. Call tools in order: stats → discordant_pairs → soft_clips → split_reads → depth_profile → summarize.
4. Your final report must cite ONLY values the tools returned in this session.
5. Do NOT use prior knowledge about genes, cell lines, or variants.
6. If all evidence types return 0, state that clearly. Do not invent signal.
7. State the sequencing technology at the start — it determines which evidence layers apply.
8. For balanced translocations: flat depth is EXPECTED. Do not interpret it as negative evidence.
9. breakpoint_evidence_summary returns both evidence_score (normalised over
   applicable_layers) and evidence_score_raw (always over all 4 layers).
   Report evidence_score/evidence_strength as the primary finding — it's the
   one that isn't artificially capped by layers this data could never
   produce — but evidence_score_raw is available if asked for the
   unnormalised view.
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
                window_bp: int = 200, min_mapq: int = 0) -> dict:
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
    depth_ratio_min_to_mean < 0.7 suggests deletion (summary.likely_deletion
    reflects this same threshold — see bam_tools.DEPTH_RATIO_DELETION_THRESHOLD,
    shared with breakpoint_evidence_summary's depth_score).
    """
    return get_read_depth_profile(bam_path, chromosome, start, end, window_size)

@mcp.tool()
def breakpoint_evidence_summary(bam_path: str, chromosome: str, position: int,
                                 label: str = "", applicable_layers: list = None) -> dict:
    """
    Integrates evidence layers into one structured report. Call this LAST.

    Pass applicable_layers (from detect_applicable_layers, called once per
    BAM) so the composite score isn't penalised by layers that can't
    structurally produce signal for this data — e.g. discordant_pairs on
    unpaired long reads, or split_reads on an aligner with no SA-tag
    support. Without it, evidence_score is capped well below "strong" for
    any technology missing a layer, even when the applicable evidence is
    overwhelming.

    Returns exactly these fields:
      label, chromosome, position — echoed back from the call
      evidence_score (float, 0-100, normalised over applicable_layers only —
        this is the primary score; report this one)
      evidence_score_raw (float, 0-100, direct sum over all 4 layers
        regardless of applicability — always ≤ evidence_score)
      evidence_strength ("none"|"weak"|"moderate"|"strong", derived from
        evidence_score, not evidence_score_raw)
      applicable_layers (list[str] — which layers evidence_score was
        normalised over; defaults to all 4 if not passed)
      signal_layers (str, "N/M" — M is len(applicable_layers), N is how many
        of those showed any signal)
      discordant_pair_score, soft_clip_score, split_read_score, depth_score (each 0-25,
        the decomposed per-layer scores — these always sum exactly to evidence_score_raw,
        and to evidence_score only when all 4 layers are applicable)
      locus_stats, discordant_pairs, soft_clips, split_reads, depth_profile
        (the full raw dict returned by each underlying tool, for inspection)
      supporting_observations (list[str] — plain-language notes on what fired)
      interpretation_template (str — one sentence restating the above fields;
        adds no new facts, only recombines values already in this same dict)
    """
    return summarize_breakpoint_evidence(bam_path, chromosome, position, label,
                                         applicable_layers=applicable_layers)

@mcp.tool()
def applicable_layers(bam_path: str, sample_reads: int = 1000) -> dict:
    """
    Samples reads from a BAM to determine which of the 4 evidence layers can
    structurally produce signal for this data — without needing to already
    know the sequencing technology or aligner. Call this FIRST, once per
    BAM, before breakpoint_evidence_summary.

    discordant_pairs requires paired reads (false for unpaired long reads:
    PacBio HiFi, ONT). split_reads requires at least one SA (supplementary
    alignment) tag anywhere in the sample (false for aligners that don't
    emit chimeric alignments). soft_clipped_reads and read_depth are always
    applicable.

    Returns applicable_layers (list[str], pass this straight through to
    breakpoint_evidence_summary's applicable_layers parameter),
    reads_sampled (int), and evidence (dict explaining each layer's
    determination in plain language).
    """
    return detect_applicable_layers(bam_path, sample_reads)

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

@mcp.tool()
def igv_screenshot(bam_paths: list, chromosome: str, start: int, end: int,
                   output_path: str, genome_build: str = "hg38",
                   color_by: str = "UNEXPECTED_PAIR",
                   max_coverage: int = None,
                   coverage_height: int = 120) -> dict:
    """
    Generate a visual IGV screenshot of a breakpoint region.
    Call this AFTER gathering evidence, to produce visual confirmation
    a clinician can review.

    Choose color_by based on the suspected event:
    - UNEXPECTED_PAIR for translocations (shows inter-chromosomal pairs,
      plus anomalous insert size/orientation — this is the default)
    - PAIR_ORIENTATION for inversions
    - INSERT_SIZE for deletions and duplications

    color_by is validated against this IGV build's actual coloring options
    before launch; an invalid value returns a structured error listing the
    valid options immediately rather than hanging until timeout. Note there
    is no "MATE_CHROMOSOME" coloring option in this IGV build — that name
    only exists as a "group by" option, not a "color by" one.

    Use a window of roughly 2-5kb around the breakpoint for readable output.
    Set max_coverage slightly above the observed max_depth from
    read_depth_profile so the coverage track is not clipped — clipping
    hides deletions, since IGV otherwise autoscales the coverage track to
    the tallest window in view, flattening a real depth dip elsewhere in
    the region.
    coverage_height sets the coverage track's pixel height (IGV default
    is ~50px, too short to render a depth dip legibly).
    Returns the path to the PNG image and the exact IGV batch script used,
    so the visualization is fully reproducible.
    """
    return run_igv_screenshot(bam_paths, chromosome, start, end,
                              output_path, genome_build, color_by,
                              max_coverage=max_coverage,
                              coverage_height=coverage_height)

@mcp.tool()
def evidence_panel(bam_paths: list, chromosome: str, start: int, end: int,
                   output_dir: str, applicable_layers: list = None) -> dict:
    """
    Generates one screenshot PER evidence layer, instead of a single image
    trying to show everything at once — each layer gets the IGV settings
    that actually isolate it visually:
      - discordant_pairs: UNEXPECTED_PAIR coloring (inter-chromosomal /
        anomalous pairs highlighted red)
      - split_reads: grouped by SA tag value — every chimeric read gets its
        own labeled row showing the exact partner locus text (e.g.
        "chr8,47000000,+,..."). This is the ONLY way to see split-read
        evidence visually: ordinary pair coloring never reflects SA tags,
        since IGV colors by the primary alignment's actual mate, not by
        tag content.
      - read_depth: alignment rows removed entirely, coverage track only,
        tall, at a fixed scale auto-computed from this exact region's
        observed max depth (+15% headroom) — not autoscaled, so a real dip
        isn't flattened by a taller peak elsewhere in view.
      - soft_clipped_reads: soft-clipped bases shown, sorted by base.
        Caveat found during testing: at a wide region (the kind read_depth
        needs to show a whole deletion span) individual soft-clip marks
        are not legible — this layer specifically benefits from a much
        narrower window (a few hundred bp) than the others. Consider
        calling this tool twice with different windows if both a depth
        overview and legible clip detail are needed.

    If applicable_layers isn't given, calls applicable_layers (the other
    tool) on bam_paths[0] first and uses its result — call that tool
    yourself first if you want to inspect the reasoning before generating
    screenshots. A layer found inapplicable (e.g. split_reads on a BAM with
    no SA tags) is not screenshotted — its entry in the result is
    {"skipped": true, "reason": "..."} instead of a PNG path, with the same
    plain-language reason applicable_layers itself would give.

    Returns:
      region, bam_paths, applicable_layers, applicable_layers_source, and
      panels: a dict with one entry per evidence layer, each either a
      run_igv_screenshot-shaped result (screenshot_path, success, etc.) or
      a {"skipped": true, "reason": ...} entry.
    """
    return igv_evidence_panel(bam_paths, chromosome, start, end, output_dir,
                              applicable_layers=applicable_layers)

if __name__ == "__main__":
    mcp.run()
