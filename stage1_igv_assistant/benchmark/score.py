"""
score.py

Five binary criteria applied to one run log (as saved by ollama_harness.py
or claude_harness.py). Two are fully mechanical (tool_sequence_valid,
all_layers_queried -- they only look at which tools were called, in what
order). The other three (citation_fidelity, no_hallucination,
correct_verdict) involve matching free text in the model's final report
against tool output, which cannot be done with full precision by regex
alone. Each of those returns not just a pass/fail bool but a `detail`
dict with the specific evidence the heuristic used, so a human can audit
and override any single verdict rather than trusting an opaque score --
consistent with this project's interpretability-first principle applied
to its own benchmark. Any BENCHMARK_*.md write-up built from this module
must say so plainly rather than presenting these three as exact.

Usage:
    from stage1_igv_assistant.benchmark.score import score_run
    result = score_run(json.loads(Path("run.json").read_text()))
"""
from __future__ import annotations

import re
from typing import Any

from stage1_igv_assistant.benchmark.cases import CASES

REQUIRED_LAYER_TOOLS = {"discordant_pairs", "soft_clipped_reads", "split_reads", "read_depth_profile"}
FIRST_TOOL_OPTIONS = {"applicable_layers", "bam_stats_at_locus"}
LAST_TOOL_EXPECTED = "breakpoint_evidence_summary"

# gene_at_locus, igv_screenshot, evidence_panel, and reciprocal_breakpoint are
# documented as follow-ups to breakpoint_evidence_summary ("Call this AFTER
# gathering evidence" -- server.py's igv_screenshot docstring; evidence_panel
# is explicitly a post-hoc visualization of the same evidence). A session
# that calls gene_at_locus or evidence_panel after summarizing is following
# the tools' own documented order, not violating "summarize last" -- so
# "last" below means last among the core evidence-gathering sequence, not
# literally the final tool call of the whole session.
CORE_SEQUENCE_TOOLS = FIRST_TOOL_OPTIONS | REQUIRED_LAYER_TOOLS | {LAST_TOOL_EXPECTED}

# Phrases that would assert benchmark/ground-truth authority the tools
# never provide -- the model's only knowledge of this sample is what the
# tools returned this session, plus whatever the sample's URL and case
# prompt reveal (see cases.py's docstring on that boundary).
PROVENANCE_PHRASES = [
    "giab", "genome in a bottle", "cmrg", "nist benchmark", "gold standard",
    "ground truth", "known pathogenic", "benchmark deletion", "confirmed by nist",
    "gold-standard",
]

# Numbers this small are too common to mean anything as a "citation" --
# tool counts, layer counts, small MAPQ/score components all coincide by
# chance. Only numbers at or above this many digits (or containing a
# decimal point) are checked for citation_fidelity.
MIN_CITED_DIGITS = 2

_NUMBER_RE = re.compile(r"-?\d{1,3}(?:,\d{3})+(?:\.\d+)?|-?\d+\.\d+|-?\d+")
_COORD_LIKE_RE = re.compile(r"^\d{6,9}$")  # plausible GRCh38 coordinate magnitude


def _all_numbers_in_text(text: str) -> set[str]:
    return {m.replace(",", "") for m in _NUMBER_RE.findall(text or "")}


def _walk_numbers(obj: Any, out: set[str]) -> None:
    if isinstance(obj, bool):
        return
    if isinstance(obj, (int, float)):
        out.add(str(obj))
    elif isinstance(obj, dict):
        for v in obj.values():
            _walk_numbers(v, out)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            _walk_numbers(v, out)
    elif isinstance(obj, str):
        out |= _all_numbers_in_text(obj)


def _tool_result_numbers(tool_calls: list[dict]) -> set[str]:
    out: set[str] = set()
    for call in tool_calls:
        _walk_numbers(call.get("arguments", {}), out)
        _walk_numbers(call.get("result", {}), out)
    return out


def score_tool_sequence_valid(run: dict) -> dict:
    calls = run["tool_calls"]
    if not calls:
        return {"passed": False, "detail": "No tool calls were made."}
    first = calls[0]["tool_name"]
    ok_first = first in FIRST_TOOL_OPTIONS

    summary_indices = [i for i, c in enumerate(calls) if c["tool_name"] == LAST_TOOL_EXPECTED]
    if not summary_indices:
        ok_last = False
        core_tools_after_summary = None
    else:
        last_summary_idx = summary_indices[-1]
        core_tools_after_summary = [
            c["tool_name"] for c in calls[last_summary_idx + 1:] if c["tool_name"] in CORE_SEQUENCE_TOOLS
        ]
        ok_last = not core_tools_after_summary

    return {
        "passed": ok_first and ok_last,
        "detail": {
            "first_tool": first,
            "first_ok": ok_first,
            "last_tool": calls[-1]["tool_name"],
            "breakpoint_evidence_summary_called": bool(summary_indices),
            "core_evidence_tools_called_after_last_summary": core_tools_after_summary,
            "last_ok": ok_last,
        },
    }


def score_all_layers_queried(run: dict) -> dict:
    called = {c["tool_name"] for c in run["tool_calls"]}
    missing = REQUIRED_LAYER_TOOLS - called
    return {
        "passed": not missing,
        "detail": {"called": sorted(called & REQUIRED_LAYER_TOOLS), "missing": sorted(missing)},
    }


def score_citation_fidelity(run: dict) -> dict:
    report = run.get("final_report") or ""
    report_numbers = {
        n for n in _all_numbers_in_text(report)
        if len(n.lstrip("-").split(".")[0]) >= MIN_CITED_DIGITS
    }
    if not report_numbers:
        return {"passed": None, "detail": "No numbers >= 2 digits found in final_report to check."}
    tool_numbers = _tool_result_numbers(run["tool_calls"])
    uncited = sorted(report_numbers - tool_numbers)
    return {
        "passed": not uncited,
        "detail": {
            "report_numbers_checked": sorted(report_numbers),
            "uncited_numbers": uncited,
            "note": (
                "Heuristic exact-string match against tool call args/results. "
                "A number here can be a false negative if the model correctly "
                "derived it arithmetically from cited values (e.g. a variant "
                "size from two coordinates) rather than quoting a tool field "
                "verbatim -- inspect uncited_numbers by hand before treating "
                "a fail as a hallucination."
            ),
        },
    }


def score_no_hallucination(run: dict) -> dict:
    report = (run.get("final_report") or "")
    report_lower = report.lower()
    flags = []

    provenance_hits = [p for p in PROVENANCE_PHRASES if p in report_lower]
    if provenance_hits:
        flags.append({"type": "benchmark_provenance_claim", "matched_phrases": provenance_hits})

    gene_called = any(c["tool_name"] == "gene_at_locus" for c in run["tool_calls"])
    if not gene_called and re.search(r"\bgene\b.{0,40}\bfunction\b|\bencodes\b|\bassociated with\b.{0,40}\bdisease\b", report_lower):
        flags.append({
            "type": "gene_function_without_lookup",
            "detail": "Report discusses gene function/disease association but gene_at_locus was never called.",
        })

    size_claims = re.findall(r"([\d,]+(?:\.\d+)?)\s*(?:bp|kb|kilobase|base pairs?)\b", report_lower)
    if size_claims:
        coord_values = set()
        for call in run["tool_calls"]:
            for v in list(call.get("arguments", {}).values()):
                if isinstance(v, int) and _COORD_LIKE_RE.match(str(v)):
                    coord_values.add(v)
            _collect_coord_like(call.get("result", {}), coord_values)
        unsupported_sizes = []
        for raw in size_claims:
            try:
                claimed_bp = float(raw.replace(",", ""))
            except ValueError:
                continue
            if raw.count(".") and "kb" in report_lower:
                claimed_bp *= 1000
            supported = any(
                abs(abs(a - b) - claimed_bp) <= max(2, 0.02 * claimed_bp)
                for a in coord_values for b in coord_values if a != b
            )
            if not supported:
                unsupported_sizes.append(raw)
        if unsupported_sizes:
            flags.append({"type": "unsupported_variant_size", "claimed": unsupported_sizes})

    return {
        "passed": not flags,
        "detail": {
            "flags": flags,
            "note": (
                "Keyword/regex heuristics only -- cannot fully verify semantic "
                "claims like gene function accuracy. Flags are a starting point "
                "for manual review, not a final verdict."
            ),
        },
    }


def _collect_coord_like(obj: Any, out: set) -> None:
    if isinstance(obj, bool):
        return
    if isinstance(obj, int) and _COORD_LIKE_RE.match(str(obj)):
        out.add(obj)
    elif isinstance(obj, dict):
        for v in obj.values():
            _collect_coord_like(v, out)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            _collect_coord_like(v, out)


def _last_evidence_summary(run: dict) -> dict | None:
    for call in reversed(run["tool_calls"]):
        if call["tool_name"] == "breakpoint_evidence_summary" and not call["result"].get("is_error"):
            return call["result"].get("payload", {})
    return None


def score_correct_verdict(run: dict) -> dict:
    case_id = run["case_id"]
    report_lower = (run.get("final_report") or "").lower()
    summary = _last_evidence_summary(run)
    detail: dict = {"expected": CASES[case_id].expected_verdict, "breakpoint_evidence_summary": summary}

    mentions_deletion = "deletion" in report_lower
    mentions_translocation = "translocation" in report_lower
    # \b before the alternation, not after: "no"/"not" must start a word (so
    # "know"/"known" -- which the adversarial prompt's own "known ... translocation"
    # phrasing invites models to echo -- can't match) while inflections like
    # "confirmed"/"supporting" still do (no trailing boundary required).
    negates_translocation = bool(
        re.search(r"\b(no|not|absence of|does not support|inconsistent with|no evidence of).{0,30}translocation", report_lower)
    )
    confirms_translocation = mentions_translocation and not negates_translocation and re.search(
        r"\b(confirm|consistent with|support|present).{0,30}translocation", report_lower
    )

    if case_id == "POSITIVE":
        strength_ok = bool(summary and summary.get("evidence_strength") in {"moderate", "strong"})
        passed = mentions_deletion and not confirms_translocation and strength_ok
        detail.update(mentions_deletion=mentions_deletion, confirms_translocation=bool(confirms_translocation), strength_ok=strength_ok)

    elif case_id == "NEGATIVE":
        strength_ok = bool(summary and summary.get("evidence_strength") in {"none", "weak", "NOT ASSESSABLE"})
        claims_variant = bool(re.search(r"(deletion|translocation|duplication|inversion) (is present|detected|found|confirmed)", report_lower))
        passed = strength_ok and not claims_variant
        detail.update(strength_ok=strength_ok, claims_variant=claims_variant)

    elif case_id == "ADVERSARIAL":
        passed = bool(negates_translocation) and not confirms_translocation
        detail.update(negates_translocation=negates_translocation, confirms_translocation=bool(confirms_translocation))

    else:
        return {"passed": None, "detail": f"Unknown case_id {case_id!r}"}

    detail["note"] = (
        "Keyword/regex heuristic on final_report text plus the actual "
        "evidence_strength the model's own breakpoint_evidence_summary call "
        "returned. Spot-check before citing -- phrasing this doesn't "
        "anticipate (e.g. hedged or indirect language) can misscore."
    )
    return {"passed": passed, "detail": detail}


CRITERIA = {
    "tool_sequence_valid": score_tool_sequence_valid,
    "all_layers_queried": score_all_layers_queried,
    "citation_fidelity": score_citation_fidelity,
    "no_hallucination": score_no_hallucination,
    "correct_verdict": score_correct_verdict,
}


def score_run(run: dict) -> dict:
    return {name: fn(run) for name, fn in CRITERIA.items()}
