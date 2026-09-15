"""
candidate_server.py
FastMCP server exposing the candidate-set bridge: 4 tools that read an
offline SV caller's VCF/BCF and hand normalised junctions to the assistant.

Deliberately a SEPARATE server from server.py. The eleven evidence tools and
their contract are unchanged, and tests/test_server.py still asserts exactly
those eleven. A client configures both servers.

Anti-hallucination design, unchanged from server.py: the assistant may state
a genomic fact only if a tool returned it in this session. These tools are the
only route by which a candidate coordinate can enter the conversation.

Run from repo root: python -m stage1_igv_assistant.candidate_server
"""

from fastmcp import FastMCP
from stage1_igv_assistant.tools.vcf_tools import (
    load_candidate_set as _load_candidate_set,
    list_candidates as _list_candidates,
    get_candidate as _get_candidate,
    compare_candidate_sets as _compare_candidate_sets,
    DEDUP_TOLERANCE_BP,
    RECURRENCE_TOLERANCE_BP,
)

mcp = FastMCP(
    name="SV Candidate Set Bridge",
    instructions="""
You read candidate structural variants produced by an offline SV caller.

RULES:
1. Call load_candidate_set FIRST for each candidate file. It returns a set_id
   and aggregate counts only -- no coordinates.
2. Obtain candidates ONLY from list_candidates or get_candidate. Never state a
   candidate coordinate that a tool did not return in this session, and never
   pass a coordinate you composed yourself to the evidence tools.
3. Every list is filtered unless you filtered nothing. Each response reports
   total_matching AND total_in_set. When you quote a count, quote the filters
   with it -- they are in filters_applied, each with a provenance label.
4. provenance "author judgement" means the threshold is the author's choice,
   not a recommendation from the caller or the reference. Say so when you
   quote a number that depended on one.
5. Recurrence between two samples is an ARTIFACT INDICATOR, not a finding.
   Do not interpret it as shared biology.
6. These tools return no sample name and no file path. You do not know which
   patient a set came from and must not speculate.
7. Do not interpret the biology of any candidate. Report machinery and counts.
""",
)


@mcp.tool()
def load_candidate_set(path: str, label: str,
                       dedup_tolerance_bp: int = DEDUP_TOLERANCE_BP) -> dict:
    """
    Register an SV caller's VCF/BCF for this session and report what is in it.

    Detects the breakend convention from the file's header AND records:
    single-record CHR2/POS2 (DELLY) or MATEID-paired records (manta, VCF spec).
    Fails loudly if BND records use neither, rather than guessing.

    Returns set_id, detected convention, total records, counts by SVTYPE and
    by FILTER, and deduplication counts. Returns NO coordinates, NO sample
    name and NO file path. Use the set_id in every later call.

    dedup_tolerance_bp merges junctions whose both breakends agree within it;
    it is AUTHOR JUDGEMENT and is echoed back in thresholds_applied.
    """
    return _load_candidate_set(path, label, dedup_tolerance_bp)


@mcp.tool()
def list_candidates(set_id: str, svtype: str = None, filter_pass: bool = False,
                    min_pe: int = None, min_sr: int = None,
                    primary_only: bool = False, mask_path: str = None,
                    limit: int = 50, offset: int = 0) -> dict:
    """
    Candidates from a registered set, with every applied threshold echoed back.

    Each entry in filters_applied carries the threshold's value, its
    provenance ('tool-defined', 'reference-defined' or 'author judgement'),
    the cumulative survivors after that step, AND how many that filter would
    remove from the unfiltered set on its own. Both are given because a step
    can look inert in a cumulative chain purely because an earlier step
    already removed the records it would have caught.

    The response always reports total_matching AND total_in_set. Do not
    present a filtered list as the complete set.

    Each candidate carries an opaque candidate_id, stable across calls, so you
    can refer back to one without restating its coordinates.
    """
    return _list_candidates(set_id, svtype, filter_pass, min_pe, min_sr,
                            primary_only, mask_path, limit, offset)


@mcp.tool()
def get_candidate(set_id: str, candidate_id: str) -> dict:
    """
    Full detail for one junction, including breakend_1 and breakend_2 shaped
    as ready-to-use arguments for the Stage 1 evidence tools.

    Pass those straight through. The coordinates originate from the candidate
    file; you must not retype or adjust them.
    """
    return _get_candidate(set_id, candidate_id)


@mcp.tool()
def compare_candidate_sets(set_a: str, set_b: str,
                           tolerance_bp: int = RECURRENCE_TOLERANCE_BP) -> dict:
    """
    Junction-level recurrence between two registered sets.

    Two uses: recurrence between two samples (an artifact indicator, because
    the same junction in two unrelated genomes is more likely a mapping
    artifact or common variant than a finding), and concordance between two
    callers on one sample once a second caller is added.

    Returns matched and unmatched counts for both sets plus the matched
    candidate_id pairs, and echoes tolerance_bp with its provenance
    (AUTHOR JUDGEMENT). Reports counts; does not interpret them.
    """
    return _compare_candidate_sets(set_a, set_b, tolerance_bp)


if __name__ == "__main__":
    mcp.run()
