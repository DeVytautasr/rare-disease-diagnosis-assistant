"""
server.py
FastMCP server exposing 11 BAM/breakpoint inspection tools to an LLM.
Anti-hallucination design: the LLM reads only tool output,
never adds genomic facts from its own training data.

Run from repo root: python -m stage1_igv_assistant.server
"""

import os

from fastmcp import FastMCP
from stage1_igv_assistant.tools.bam_tools import (
    image_session_dir,
    to_handle_result,
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
10. You have not seen any image produced by the screenshot tools. Never
    describe what an image shows, contains, or looks like. State that it
    was generated, give its path, and say what a reviewer should check.
    Describing image contents you have not been shown is fabrication,
    equivalent to citing a number no tool returned.
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
    "assessable": false (with a "reason") and discordant_fraction: null mean
    zero reads fell in this window — a different claim from "0 discordant
    pairs found among reads that were there".
    """
    return count_discordant_pairs(bam_path, chromosome, position, window_bp, min_mapq)

@mcp.tool()
def soft_clipped_reads(bam_path: str, chromosome: str, position: int,
                       window_bp: int = 200, min_clip_bases: int = 10,
                       min_mapq: int = 20) -> dict:
    """
    Count reads with soft-clipped overhangs near a breakpoint.
    A pileup at the SAME position (consensus_clip_position) narrows the breakpoint precisely.
    max_clips_at_position < 3 = no real pileup, treat as noise. This is also
    exactly what breakpoint_evidence_summary's soft_clip_score now tiers on
    (3-9 = partial, >=10 = full) instead of soft_clipped_fraction, since
    fraction alone can't tell a genuine pileup from scattered clipping.
    "assessable": false (with a "reason") and soft_clipped_fraction: null mean
    zero reads fell in this window — a different claim from "0 clipped reads
    found among reads that were there".
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
    "assessable": false (with a "reason") and split_read_fraction: null mean
    zero reads fell in this window — a different claim from "0 split reads
    found among reads that were there".
    """
    return get_split_reads(bam_path, chromosome, position, window_bp, min_mapq)

@mcp.tool()
def read_depth_profile(bam_path: str, chromosome: str, start: int,
                       end: int, window_size: int = 100,
                       focus_position: int = None,
                       dip_tolerance_bp: int = 1000) -> dict:
    """
    Sliding-window TRUE per-base read depth across a region (each window's
    "depth" is aligned bases summed into that bin, divided by bin width —
    not a read count; confirmed bin-size invariant to within ~1-2% on real
    data, unlike the pre-2026-08-11 implementation).
    Deletions: ~50% depth drop inside deleted region (heterozygous).
    Duplications: depth rises.
    Balanced translocations: depth stays FLAT — this is expected and correct.
    depth_ratio_min_to_mean < 0.7 suggests deletion (summary.likely_deletion
    reflects this same threshold — see bam_tools.DEPTH_RATIO_DELETION_THRESHOLD,
    shared with breakpoint_evidence_summary's depth_score).
    If summary.assessable is false, the region had zero reads at all —
    depth_ratio_min_to_mean and likely_deletion are both null, not 0.0/true.
    A 0/0 ratio is undefined, not "depth dropped to 0% of the mean".

    Pass focus_position (a candidate breakpoint within [start, end)) to
    localize depth_ratio_min_to_mean/likely_deletion to the depth AROUND
    that position instead of the region's global minimum, which can sit
    anywhere in a wide scan — including nowhere near the position actually
    being investigated. When set, the response also adds dip_position
    (where the region's true minimum actually is), dip_distance_from_focus,
    and dip_is_at_focus (true when the minimum is within dip_tolerance_bp,
    default 1000bp — HEURISTIC, see the function's own docstring for
    calibration provenance) — use these to tell "real signal at this
    position" apart from "an unrelated dip elsewhere dragged the ratio
    down", which breakpoint_evidence_summary's depth_score now requires
    before awarding points.
    """
    return get_read_depth_profile(bam_path, chromosome, start, end, window_size,
                                  focus_position=focus_position,
                                  dip_tolerance_bp=dip_tolerance_bp)

@mcp.tool()
def breakpoint_evidence_summary(bam_path: str, chromosome: str, position: int,
                                 label: str = "", applicable_layers: list = None,
                                 window_bp: int = 500, min_mapq: int = 20) -> dict:
    """
    Integrates evidence layers into one structured report. Call this LAST.

    THRESHOLD PROVENANCE — state this plainly in any report citing
    evidence_strength, don't present it as calibrated: of the 9 tier
    cutoffs behind the 4 component scores, only ONE (the read-depth
    layer's 0.7 moderate-tier threshold) is empirically calibrated, and
    against a single confirmed real locus. The other 8 (discordant-pair
    0.2/0.5, soft-clip 3/10 — on max_clips_at_position, not fraction —,
    split-read 0.1/0.3, read-depth's own 0.3 strong-tier threshold, and
    the depth layer's dip_tolerance_bp=1000bp localization radius) are
    heuristic judgement calls, never validated against real data beyond a
    handful of loci. evidence_strength is an interpretable decomposition
    of what fired, not a calibrated probability.

    DEPTH SCORING is now two-part: depth_ratio_min_to_mean (localized to
    focus_position=position, not the window's global minimum — see
    read_depth_profile) must cross the tier threshold, AND dip_is_at_focus
    must be True — the region's actual lowest point, not just a nearby
    one, has to be near the breakpoint. A dip that fails the second check
    is reported in supporting_observations as off-position, not scored.
    This closes a real false positive: a locus at a local peak ~1500bp
    from an unrelated real dip elsewhere in the window used to score full
    depth points for a deletion that wasn't there.

    Pass applicable_layers (from detect_applicable_layers, called once per
    BAM) so the composite score isn't penalised by layers that can't
    structurally produce signal for this data — e.g. discordant_pairs on
    unpaired long reads, or split_reads on an aligner with no SA-tag
    support. Without it, evidence_score is capped well below "strong" for
    any technology missing a layer, even when the applicable evidence is
    overwhelming.

    ASSESSABILITY: a layer with zero reads in its own window at this locus
    (distinct from being structurally inapplicable to the technology) is
    never assigned a component score — None, not 0.0, since "not scored"
    and "scored 0" are different claims — and is excluded from
    evidence_score's numerator/denominator the same way an inapplicable
    layer is. It still contributes 0 to evidence_score_raw (which stays a
    plain number, never None). If every layer this call would otherwise
    score turns out unassessable or inapplicable, evidence_score is None
    and evidence_strength is the distinct value "NOT ASSESSABLE" — never
    "none" (which means real evidence was checked and found absent) and
    never a bare 0/100.

    Returns exactly these fields:
      label, chromosome, position — echoed back from the call
      evidence_score (float 0-100, or None if NOT ASSESSABLE — normalised
        over applicable AND assessable layers only; this is the primary
        score when not None; report this one)
      evidence_score_raw (float, 0-100, direct sum over all 4 layers
        regardless of applicability/assessability, treating any excluded
        layer as a 0 contribution — never None, but see ASSESSABILITY
        above before treating a low raw score as informative on its own)
      evidence_strength ("none"|"weak"|"moderate"|"strong"|"NOT ASSESSABLE",
        derived from evidence_score, not evidence_score_raw)
      applicable_layers (list[str] — which layers were candidates for
        evidence_score's denominator, before assessability filtering;
        defaults to all 4 if not passed)
      unassessable_layers (dict[str, str] — {layer_name: reason} for any
        layer excluded for having zero reads in its window; empty dict
        when every applicable layer had reads to assess)
      signal_layers (str, "N/M" — M is the count of applicable-AND-assessable
        layers, N is how many of those showed any signal; "0/0" when
        evidence_strength is "NOT ASSESSABLE")
      discordant_pair_score, soft_clip_score, split_read_score, depth_score
        (each 0-25, or None for a layer in unassessable_layers — the
        non-None ones always sum exactly to evidence_score_raw's
        corresponding contribution, and to evidence_score only when all 4
        layers are both applicable and assessable)
      locus_stats, discordant_pairs, soft_clips, split_reads, depth_profile
        (the full raw dict returned by each underlying tool, for inspection —
        each of the latter 4 carries its own "assessable"/"reason" fields)
      supporting_observations (list[str] — plain-language notes on what fired,
        including a note for each unassessable layer)
      interpretation_template (str — one sentence restating the above fields;
        adds no new facts, only recombines values already in this same dict)

    window_bp (default 500) sets the half-width for the discordant-pair /
    soft-clip / split-read windows — narrow it (e.g. 150-200) for a tight,
    precisely-localized breakpoint, widen it if the exact position is
    uncertain. Does NOT affect the depth-profile window, which is always a
    fixed ±2kb regardless of this value. min_mapq (default 20) sets the
    minimum mapping quality across all layers — lower it for regions with
    known poor mappability (check low_mapq_fraction from bam_stats_at_locus
    first) rather than silently accepting an all-zero result there.
    """
    return summarize_breakpoint_evidence(bam_path, chromosome, position, label,
                                         applicable_layers=applicable_layers,
                                         window_bp=window_bp, min_mapq=min_mapq)

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
    Requires internet access to query the Ensembl REST API. Retries on
    HTTP 429 (rate-limited) and transient network errors before giving up.

    Returns exactly these fields on success:
      chromosome, position, genome_build — echoed back from the call
      gene_count (int — number of genes overlapping this exact position)
      is_intergenic (bool — True iff gene_count == 0)
      genes (list[dict] — one entry per overlapping gene, each with
        gene_id, gene_name, biotype (protein_coding / lncRNA / etc),
        strand ("+"|"-"), gene_start, gene_end)
      annotation_note (str — plain-language recap of what the lookup
        established: which annotated genes overlap this position, or
        that it is intergenic. States overlap only — a coordinate
        lookup cannot establish that a breakpoint exists or that any
        gene is disrupted)

    On failure (Ensembl unreachable after retries), returns instead:
      error (str), chromosome, position,
      note ("Ensembl unavailable — gene annotation skipped")
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
                   genome_build: str = "hg38",
                   color_by: str = "UNEXPECTED_PAIR",
                   max_coverage: int = None,
                   coverage_height: int = 120) -> dict:
    """
    Generate a visual IGV screenshot of a breakpoint region.
    Call this AFTER gathering evidence, to produce visual confirmation
    a clinician can review.

    You do NOT choose where the image is written and you do NOT receive a
    file path. The server assigns the location and returns an opaque
    image_ref (e.g. "IMG_a3f9") plus metadata you can legitimately reason
    about: region, coloring mode, pixel dimensions, and success/failure.
    The image itself is never provided to you, so you cannot and must not
    describe what it shows — report that it was generated, which region it
    covers, and what a reviewer should check.

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
    """
    session_dir = image_session_dir()
    output_path = os.path.join(
        session_dir, f"{chromosome}_{start}_{end}_{color_by}.png"
    )
    result = run_igv_screenshot(bam_paths, chromosome, start, end,
                                output_path, genome_build, color_by,
                                max_coverage=max_coverage,
                                coverage_height=coverage_height)
    return to_handle_result(result, session_dir)

@mcp.tool()
def evidence_panel(bam_paths: list, chromosome: str, position: int,
                   start: int = None, end: int = None,
                   applicable_layers: list = None, windows: dict = None) -> dict:
    """
    Generates one screenshot PER evidence layer, instead of a single image
    trying to show everything at once.

    You do NOT choose where images are written and you do NOT receive file
    paths. The server assigns the location; each panel entry carries an
    opaque image_ref (e.g. "IMG_a3f9") plus region, coloring mode, pixel
    dimensions, and success/failure. The images themselves are never
    provided to you, so you cannot and must not describe what they show —
    report which layers were generated, the region each covers, and what a
    reviewer should check. Each layer gets BOTH the IGV
    settings that isolate it visually AND its own window around `position`
    — a single shared window can't serve every layer (confirmed directly:
    a region wide enough for a depth dip renders soft-clip marks
    illegibly small, and vice versa):
      - discordant_pairs (±1500bp default): UNEXPECTED_PAIR coloring
        (inter-chromosomal / anomalous pairs highlighted red)
      - split_reads (±1500bp default): grouped by SA tag value — every
        chimeric read gets its own labeled row showing the exact partner
        locus text (e.g. "chr8,47000000,+,..."). This is the ONLY way to
        see split-read evidence visually: ordinary pair coloring never
        reflects SA tags, since IGV colors by the primary alignment's
        actual mate, not by tag content.
      - read_depth (caller-supplied start/end if given, else ±3000bp):
        alignment rows removed entirely, coverage track only, tall, at a
        fixed scale auto-computed from that region's observed max depth
        (+15% headroom) — not autoscaled, so a real dip isn't flattened
        by a taller peak elsewhere in view.
      - soft_clipped_reads (±150bp default, tight): soft-clipped bases
        shown, sorted by RIGHT_CLIP or LEFT_CLIP — whichever side actually
        has more clipped reads at this locus (from
        count_soft_clipped_reads' dominant_clip_side), not a fixed choice.
        Confirmed empirically this produces a clean staircase pileup at
        the clip boundary. If the dominant side can't be determined,
        defaults to LEFT_CLIP and says so in the panel's
        "clip_side_determination" field.

    Override any layer's half-window via windows={"soft_clipped_reads": 300}
    (bp; keys from applicable_layers' vocabulary). The exact region used
    per layer is in the returned "windows_used" dict — cite that rather
    than assuming the defaults were used, since overrides and the
    read_depth start/end special-case both change it.

    If applicable_layers isn't given, calls applicable_layers (the other
    tool) on bam_paths[0] first and uses its result — call that tool
    yourself first if you want to inspect the reasoning before generating
    screenshots. A layer found inapplicable (e.g. split_reads on a BAM with
    no SA tags) is not screenshotted — its entry in the result is
    {"skipped": true, "reason": "..."} instead of a PNG path, with the same
    plain-language reason applicable_layers itself would give.

    Returns:
      position, chromosome, bam_paths, applicable_layers,
      applicable_layers_source, windows_used (per-layer region actually
      used, for all 4 layers regardless of skip status), and panels: a
      dict with one entry per evidence layer, each either a screenshot
      result (image_ref, image_dimensions, success, region, color_by,
      window_bp and — for soft_clipped_reads — clip_side_determination)
      or a {"skipped": true, "reason": ...} entry.
    """
    session_dir = image_session_dir()
    result = igv_evidence_panel(bam_paths, chromosome, position, session_dir,
                                start=start, end=end,
                                applicable_layers=applicable_layers,
                                windows=windows)
    panels = result.get("panels")
    if isinstance(panels, dict):
        result["panels"] = {
            layer: (entry if isinstance(entry, dict) and entry.get("skipped")
                    else to_handle_result(entry, session_dir))
            for layer, entry in panels.items()
        }
    return result

if __name__ == "__main__":
    mcp.run()
