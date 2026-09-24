"""
ui.py — local, dependency-free front end for the breakpoint assistant.

Start:  python -m stage1_igv_assistant.ui
Then:   http://127.0.0.1:8765

Design rules this file enforces mechanically, not by convention:

  * Every displayed number comes from ToolRecorder.call(), which is the only
    path to a tool anywhere in this module. There is no second route, so a
    number with no recorded call behind it cannot be rendered.
  * Tools are invoked through FastMCP's own dispatch (mcp.call_tool), so
    argument validation and result serialisation are the real ones, not a
    reimplementation.
  * The browser is never sent a filesystem path or a sample name. Datasets
    are addressed by label; the label->path map stays server-side, and every
    outbound payload is scrubbed.
  * Scoring tiers and the attainable ceiling are DERIVED from bam_tools'
    source at startup (see score_tiers.py). If the scoring function is
    revised the displayed ceiling moves with it; if it can no longer be
    parsed the UI says so instead of showing a stale figure.
"""
import argparse
import asyncio
import glob
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import traceback
import urllib.request
import zipfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs

from stage1_igv_assistant import server as evidence_server
from stage1_igv_assistant import candidate_server
from stage1_igv_assistant.tools import bam_tools
from stage1_igv_assistant import chat as chatmod
from stage1_igv_assistant import config as CFG
from stage1_igv_assistant.score_tiers import (
    derive_tiers, derive_bands, score_for, next_tier_up, TierDerivationError,
    ceiling_from_observed as derive_ceiling,
)

EXPECTED_EVIDENCE_TOOLS = {
    "bam_stats_at_locus", "discordant_pairs", "soft_clipped_reads", "split_reads",
    "read_depth_profile", "breakpoint_evidence_summary", "applicable_layers",
    "gene_at_locus", "reciprocal_breakpoint", "igv_screenshot", "evidence_panel",
}
EXPECTED_BRIDGE_TOOLS = {
    "load_candidate_set", "list_candidates", "get_candidate", "compare_candidate_sets",
}


# ── startup contract ────────────────────────────────────────────────────────
def assert_tool_contract():
    """Refuse to start if the UI's idea of the tool surface has drifted from
    what MCP actually exposes. In-process dispatch does not exercise the stdio
    transport, so this is the check that catches divergence instead."""
    async def _names(mcp):
        return {t.name for t in await mcp.list_tools()}
    ev = asyncio.run(_names(evidence_server.mcp))
    br = asyncio.run(_names(candidate_server.mcp))
    problems = []
    if ev != EXPECTED_EVIDENCE_TOOLS:
        problems.append(f"evidence server: missing={sorted(EXPECTED_EVIDENCE_TOOLS-ev)} "
                        f"unexpected={sorted(ev-EXPECTED_EVIDENCE_TOOLS)}")
    if br != EXPECTED_BRIDGE_TOOLS:
        problems.append(f"bridge server: missing={sorted(EXPECTED_BRIDGE_TOOLS-br)} "
                        f"unexpected={sorted(br-EXPECTED_BRIDGE_TOOLS)}")
    if problems:
        raise SystemExit("REFUSING TO START — tool surface drift:\n  " + "\n  ".join(problems))
    return len(ev), len(br)


# ── the only path to a tool ─────────────────────────────────────────────────
class ToolRecorder:
    def __init__(self):
        self.calls = []

    def call(self, which, tool, params):
        mcp = evidence_server.mcp if which == "evidence" else candidate_server.mcp
        t0 = time.time()
        err = None
        try:
            res = asyncio.run(mcp.call_tool(tool, params))
            data = res.structured_content
            is_err = bool(getattr(res, "is_error", False))
        except Exception as e:                       # dispatch/validation failure
            data, is_err, err = {"error": str(e), "error_type": "dispatch"}, True, traceback.format_exc()
        rec = {
            "id": len(self.calls) + 1,
            "server": which, "tool": tool, "params": params,
            "is_error": is_err, "ms": round((time.time() - t0) * 1000, 1),
            "result": data, "traceback": err,
        }
        self.calls.append(rec)
        return rec


def tool_error(rec):
    """The error dict a tool returned, or None. Errors are dicts by convention,
    so FastMCP reports is_error=False for them and a bare .get() on the payload
    yields None for every field -- which renders identically to a real empty
    result. Every consumer in this module goes through this."""
    r = rec.get("result") if isinstance(rec, dict) else None
    if isinstance(r, dict) and "error" in r:
        return r
    return None


RECORDER = ToolRecorder()
DATASETS = {}          # label -> bam path      (never sent to the browser)
CANDIDATE_FILES = {}   # label -> vcf/bcf path  (never sent to the browser)
TIERS = None
BANDS = None
TIER_ERROR = None


# Any absolute filesystem path of two or more segments. Not preceded by a word
# character, '.', '~', ':' or '/', so URLs, ratios such as 15/28 and relative
# fragments do not match.
_ABS_PATH = re.compile(r"(?<![\w.~:/@-])/[\w.+@%~-]+(?:/[\w.+@%~-]*)+")
_HOME = os.path.expanduser("~")


def _unpath(m):
    p = m.group(0)
    if p.startswith(("/api/", "/img/")):        # this server's own routes, not files
        return p
    if p.rstrip("/") in (_HOME, os.path.dirname(_HOME)):
        return "~"
    return os.path.basename(p.rstrip("/")) or "~"


def scrub(obj):
    """Remove anything path-like or sample-like from an outbound payload.

    Registered dataset and candidate paths (and their basenames) become <label>.
    Every other absolute path is then reduced to its basename. The second pass
    is the rule itself rather than a list of known paths: a scan of every route
    found paths the list never covered -- the exclude template in /api/funnel,
    IGV's search locations in /api/igv, tracebacks from the error branch -- so
    no route has to be judged case by case."""
    paths = sorted(set(list(DATASETS.values()) + list(CANDIDATE_FILES.values())),
                   key=len, reverse=True)
    labels = {v: k for k, v in list(DATASETS.items()) + list(CANDIDATE_FILES.items())}

    def s(x):
        if isinstance(x, dict):
            return {k: s(v) for k, v in x.items()}
        if isinstance(x, list):
            return [s(v) for v in x]
        if isinstance(x, str):
            out = x
            for p in paths:
                if p in out:
                    out = out.replace(p, f"<{labels.get(p, 'dataset')}>")
                base = os.path.basename(p)
                if base and base in out:
                    out = out.replace(base, f"<{labels.get(p, 'dataset')}>")
            return _ABS_PATH.sub(_unpath, out)
        return x
    return s(obj)


# ── evidence assembly ───────────────────────────────────────────────────────
def assess(bam_label, chromosome, position, window_bp=500, split_min_mapq=None):
    """Runs the four layers plus stats and summary at one coordinate.
    Returns per-layer numbers each tagged with the id of the call behind it."""
    bam = DATASETS[bam_label]
    # sample_reads: the tool's default of 1000 reads is a sample of the FIRST
    # 1000 reads in the file, which on a whole-genome BAM is representative and
    # on a region-sliced one is not -- it is one corner of one window. The demo
    # bundle is a slice, and at the default its leading 1000 reads carry no SA
    # tag, so split_reads was judged inapplicable, dropped out of the
    # normalisation denominator, and the SAME locus scored 43.3 instead of 47.5
    # with nothing on screen saying why. Measured: 20000 reads costs 0.05s on a
    # 1.6 GB BAM and makes the slice agree with the full file. The tool, its
    # default and its behaviour are unchanged; only this caller asks for more.
    layers = RECORDER.call("evidence", "applicable_layers",
                           {"bam_path": bam, "sample_reads": 20000})
    le = tool_error(layers)
    if le:
        return {"error": le["error"], "error_type": le.get("error_type"),
                "chromosome": chromosome, "position": position, "bam_label": bam_label,
                "call": layers["id"],
                "hint": "The dataset itself could not be read, so no layer could be assessed."}
    applicable = (layers["result"] or {}).get("applicable_layers")

    # The low-MAPQ figure shown beside the verdict must be measured over the
    # SAME window the quality gate in breakpoint_evidence_summary uses:
    # position +/- window_bp, clamped at 0. This call used +/-250 while the gate
    # used +/-500, and at six loci of the public HCC1143 slice the two fell on
    # opposite sides of the 40% gate: 21:19,281,000 was withheld with the tool
    # reporting 43.6% while the page showed 35.6%. The window travels with the
    # figure so the page labels it from data, not from a literal.
    gate_window = {"start": max(0, position - window_bp), "end": position + window_bp,
                   "half_width_bp": window_bp}
    stats = RECORDER.call("evidence", "bam_stats_at_locus",
                          {"bam_path": bam, "chromosome": chromosome,
                           "start": gate_window["start"], "end": gate_window["end"]})
    # The evidence tools report failure as a dict with an "error" key rather
    # than raising, so FastMCP does not flag it and a .get() on the result
    # yields None for every field. Rendering that as empty layers would be
    # indistinguishable from a real locus with no reads -- a silently wrong
    # answer. Surface it instead.
    sres = tool_error(stats)
    if sres:
        return {"error": sres["error"], "error_type": sres.get("error_type"),
                "chromosome": chromosome, "position": position,
                "bam_label": bam_label,
                "contigs_in_header_sample": sres.get("contigs_in_header_sample"),
                "call": stats["id"],
                "hint": ("The contig was not found in this dataset's BAM header. "
                         "Check the naming style the dataset uses (for example "
                         "'chr21' versus '21').")}
    disc = RECORDER.call("evidence", "discordant_pairs",
                         {"bam_path": bam, "chromosome": chromosome,
                          "position": position, "window_bp": window_bp})
    clip = RECORDER.call("evidence", "soft_clipped_reads",
                         {"bam_path": bam, "chromosome": chromosome, "position": position})
    sp_params = {"bam_path": bam, "chromosome": chromosome, "position": position}
    if split_min_mapq is not None:
        sp_params["min_mapq"] = int(split_min_mapq)
    split = RECORDER.call("evidence", "split_reads", sp_params)
    depth = RECORDER.call("evidence", "read_depth_profile",
                          {"bam_path": bam, "chromosome": chromosome,
                           "start": position - 2000, "end": position + 2000,
                           "focus_position": position})
    summ = RECORDER.call("evidence", "breakpoint_evidence_summary",
                         {"bam_path": bam, "chromosome": chromosome, "position": position,
                          "label": f"{chromosome}:{position}",
                          "applicable_layers": applicable, "window_bp": window_bp})

    d, c, s_, dp, sm = (disc["result"] or {}, clip["result"] or {}, split["result"] or {},
                        depth["result"] or {}, summ["result"] or {})
    dsum = dp.get("summary", {}) if isinstance(dp, dict) else {}
    lay_err = {"discordant_pairs": d.get("error"), "soft_clipped_reads": c.get("error"),
               "split_reads": s_.get("error"), "read_depth": dp.get("error")}
    obs = {
        "discordant_pairs": d.get("discordant_fraction"),
        "soft_clipped_reads": c.get("max_clips_at_position"),
        "split_reads": s_.get("split_read_fraction"),
        "read_depth": dsum.get("depth_ratio_min_to_mean"),
    }
    summ_err = tool_error(summ)
    return {
        "summary_error": (summ_err or {}).get("error"),
        "chromosome": chromosome, "position": position, "bam_label": bam_label,
        "applicable_layers": applicable, "applicable_call": layers["id"],
        "stats": {"call": stats["id"], "value": stats["result"]},
        "layers": [
            {"key": "discordant_pairs", "label": "Discordant pairs",
             "count": d.get("discordant_pairs"), "fraction": d.get("discordant_fraction"),
             "partners": d.get("mate_chromosomes"), "min_mapq": d.get("min_mapq_applied"),
             "assessable": d.get("assessable"), "reason": d.get("reason"),
             "error": lay_err["discordant_pairs"],
             "quality_limited": d.get("quality_limited"), "call": disc["id"]},
            {"key": "soft_clipped_reads", "label": "Soft-clipped reads",
             "count": c.get("soft_clipped_reads"), "fraction": c.get("soft_clipped_fraction"),
             "consensus": c.get("consensus_clip_position"),
             "max_clips": c.get("max_clips_at_position"), "min_mapq": c.get("min_mapq_applied"),
             "assessable": c.get("assessable"), "reason": c.get("reason"),
             "error": lay_err["soft_clipped_reads"],
             "quality_limited": c.get("quality_limited"), "call": clip["id"]},
            {"key": "split_reads", "label": "Split reads (SA)",
             "count": s_.get("split_reads"), "fraction": s_.get("split_read_fraction"),
             "partners": s_.get("partner_chromosomes"), "min_mapq": s_.get("min_mapq_applied"),
             "examples": s_.get("example_partner_loci"),
             "assessable": s_.get("assessable"), "reason": s_.get("reason"),
             "error": lay_err["split_reads"],
             "quality_limited": s_.get("quality_limited"), "call": split["id"]},
            {"key": "read_depth", "label": "Read depth",
             "count": dsum.get("min_depth"), "fraction": dsum.get("depth_ratio_min_to_mean"),
             "mean_depth": dsum.get("mean_depth"), "dip_at_focus": dsum.get("dip_is_at_focus"),
             "min_mapq": dsum.get("min_mapq_applied"),
             "assessable": dsum.get("assessable"), "reason": dsum.get("reason"),
             "error": lay_err["read_depth"],
             "call": depth["id"]},
        ],
        "summary": {
            "call": summ["id"],
            "evidence_score": sm.get("evidence_score"),
            "evidence_score_raw": sm.get("evidence_score_raw"),
            "evidence_strength": sm.get("evidence_strength"),
            "signal_layers": sm.get("signal_layers"),
            "min_mapq_applied": sm.get("min_mapq_applied"),
            "components": {
                "discordant_pairs": sm.get("discordant_pair_score"),
                "soft_clipped_reads": sm.get("soft_clip_score"),
                "split_reads": sm.get("split_read_score"),
                "read_depth": sm.get("depth_score"),
            },
            "low_mapq_fraction": (stats["result"] or {}).get("low_mapq_fraction"),
            "low_mapq_window": gate_window,
        },
        "ceiling": ceiling_for(obs),
        "position_provenance": (d.get("position_provenance") or s_.get("position_provenance")),
    }


def ceiling_for(observed):
    """Attainable-score analysis, derived from the live tiers and bands.

    Phase 10: the derivation moved into score_tiers.ceiling_from_observed so the
    UI panel and the MCP summary tool share ONE implementation rather than two
    that can drift, and so the "strong" boundary is read out of the scoring
    source instead of being the literal 70 that used to sit here.
    """
    if TIERS is None:
        return {"derivable": False, "reason": TIER_ERROR}
    return derive_ceiling(TIERS, BANDS, observed)


# ── limits panel (Addition 3) ───────────────────────────────────────────────
LIMITS = {
    "title": "What this tool cannot tell you",
    "items": [
        {"h": "It checks positions. It does not search for them.",
         "b": "Every position examined here was either proposed by the variant caller or typed "
              "in by hand. Nothing scans the genome. A real breakpoint that the caller missed "
              "and nobody typed in will never appear here, however strong the evidence at it "
              "would have been."},
        {"h": "On a controlled test it found 14 of 24 known breakpoints.",
         "b": "Twelve balanced translocations were built into real sequencing data at positions "
              "known in advance, then looked for. Found: 8 of 8 in clean, uniquely mappable "
              "sequence; 6 of 8 next to repeats; 0 of 8 where the surrounding sequence maps "
              "ambiguously. Every one that was missed was lost at the variant-calling step, "
              "because reads crossing those junctions could not be placed confidently enough "
              "for the caller to use them. None was lost to the filters below — no filter "
              "setting tested discarded a single true breakpoint. This is a measurement on "
              "simulated reads and estimates performance on simulated reads."},
        {"h": "14 of the 16 cut-offs are judgement calls, not calibrated values.",
         "b": "They are not recommendations from the underlying tools, and they have not been "
              "checked against a set of confirmed positive and negative cases. Every cut-off "
              "shown on this page is labelled with where it came from. Read "
              "'author judgement' as: a reasonable person could have chosen differently, and "
              "the result would differ."},
        {"h": "The read-depth measurement is unreliable at this coverage.",
         "b": "Its cut-off was set using data at roughly ten times the coverage of a routine "
              "genome, and at about 31x it flags roughly 26% of ordinary positions — so on its "
              "own it distinguishes very little. (That 26% figure is carried over from earlier "
              "testing and was not re-measured for this interface.) In the controlled test it "
              "added points at two positions where no copy-number change had been built in. "
              "Treat a depth contribution as weak support at best."},
        {"h": "A balanced translocation can never score \"strong\" here.",
         "b": "This is arithmetic, not an accident of the data. In a balanced rearrangement no "
              "DNA is gained or lost, so the depth measurement correctly contributes nothing. "
              "And when only one of the two copies of a chromosome is rearranged, roughly half "
              "the reads at the breakpoint come from the intact copy, which caps the "
              "paired-read measurement well below its top band. The highest score actually "
              "reachable is calculated for each position and shown next to the score. Across 28 "
              "test breakpoints the paired-read fraction never exceeded 0.175, against the 0.5 "
              "its top band requires. Judge a balanced translocation on the four measurements, "
              "not on the band it lands in."},
        {"h": "\"Quality limited\" means no score was calculated.",
         "b": "When too many reads at a position are ambiguously mapped, the combined score is "
              "withheld instead of computed. It means this position cannot be scored, not that "
              "it scored badly. The four individual measurements are still shown; read those."},
    ],
}


_CHAT_TOOLS = {}
_CHAT_WHERE = None


def _chat_tools(trim=False):
    """Schemas generated from list_tools() on BOTH servers, never hand-written.
    `trim` shortens descriptions on the MODEL path only; the MCP descriptions
    themselves are untouched, because benchmark/mcp_client.py forwards them to
    models too and six recorded run sets depend on the current text."""
    global _CHAT_WHERE
    if trim not in _CHAT_TOOLS:
        async def _all():
            return {"evidence": await evidence_server.mcp.list_tools(),
                    "bridge": await candidate_server.mcp.list_tools()}
        mcp_tools = asyncio.run(_all())
        _CHAT_TOOLS[trim], _CHAT_WHERE = chatmod.build_tools(
            mcp_tools, DATASETS, CANDIDATE_FILES, trim=trim)
    return _CHAT_TOOLS[trim], _CHAT_WHERE


def _chat_exec(name, args):
    """The model's only route to a tool: label -> path, then ToolRecorder.
    The result is scrubbed BEFORE the model sees it, not only before the
    browser does -- applicable_layers echoes bam_path, and the model must not
    receive a filesystem path it could then repeat in its prose."""
    _, where = _chat_tools()
    resolved, err = chatmod.resolve_args(name, args, DATASETS, CANDIDATE_FILES, MASK_PATH)
    if err:
        return None, err
    # Same guard the browser's panel route uses: the panel tools do not validate
    # the contig and will report success for a locus that does not exist. The
    # model must not be able to obtain an image for a position the counting
    # tools cannot read either.
    if name in ("evidence_panel", "igv_screenshot"):
        chrom = resolved.get("chromosome")
        pos = resolved.get("position") or resolved.get("start")
        bams = resolved.get("bam_paths") or ([resolved["bam_path"]] if resolved.get("bam_path") else [])
        # The model does emit coordinates as strings in scientific notation
        # ("3.5e+07"). int() on that raises, and because this guard runs OUTSIDE
        # the tool call it took the whole turn down with it -- a harness crash
        # dressed up as a model failure. A coordinate this guard cannot parse is
        # simply not pre-checked; the tool's own schema validation then rejects
        # it and it is counted as the argument-formatting error it is.
        try:
            pos_i = int(str(pos).strip())
        except (TypeError, ValueError):
            pos_i = None
        if chrom and pos_i is not None and bams:
            pre = RECORDER.call("evidence", "bam_stats_at_locus",
                                {"bam_path": bams[0], "chromosome": chrom,
                                 "start": max(0, pos_i - 250), "end": pos_i + 250})
            pe = tool_error(pre)
            if pe:
                return None, (f"{pe['error']} No panel was generated; an image here would not "
                              f"correspond to a real locus.")
    rec = RECORDER.call(where[name], name, resolved)
    # The whole record, not only result and params: a dispatch failure carries a
    # traceback, and a traceback carries absolute paths.
    return scrub(dict(rec)), None


def _api(path, body):
    if path == "/api/bootstrap":
        return {
            "datasets": sorted(DATASETS), "candidate_files": sorted(CANDIDATE_FILES),
            "tiers": TIERS, "tier_error": TIER_ERROR, "limits": LIMITS,
            "tool_counts": {"evidence": len(EXPECTED_EVIDENCE_TOOLS),
                            "bridge": len(EXPECTED_BRIDGE_TOOLS)},
        }
    if path == "/api/load":
        lbl = body["candidates_label"]
        r = RECORDER.call("bridge", "load_candidate_set",
                          {"path": CANDIDATE_FILES[lbl], "label": lbl})
        return {"call": r["id"], "result": r["result"], "is_error": r["is_error"]}
    if path == "/api/funnel":
        p = {"set_id": body["set_id"], "limit": int(body.get("limit", 200))}
        for k in ("svtype", "min_pe", "min_sr"):
            if body.get(k) not in (None, "", "null"):
                p[k] = int(body[k]) if k != "svtype" else body[k]
        for k in ("filter_pass", "primary_only"):
            if body.get(k):
                p[k] = True
        mask = _mask_state(bool(body.get("use_mask")))
        if mask["applied"]:
            p["mask_path"] = MASK_PATH
        r = RECORDER.call("bridge", "list_candidates", p)
        return {"call": r["id"], "result": r["result"], "is_error": r["is_error"],
                "mask": mask}
    if path == "/api/candidate":
        r = RECORDER.call("bridge", "get_candidate",
                          {"set_id": body["set_id"], "candidate_id": body["candidate_id"]})
        return {"call": r["id"], "result": r["result"], "is_error": r["is_error"]}
    if path == "/api/assess":
        return assess(body["bam_label"], body["chromosome"], int(body["position"]),
                      split_min_mapq=body.get("split_min_mapq"))
    if path == "/api/igv":
        bam = DATASETS[body["bam_label"]]
        pos = int(body["position"])
        # The panel tool does NOT validate the contig: asked for a chromosome
        # absent from the BAM header it still reports success and writes PNGs,
        # which would show a viewer four panels for a locus that does not exist.
        # The counting tools do validate, so ask one of them first and refuse to
        # render anything if the locus is not real. The evidence tools are not
        # modified; this is a guard in front of them.
        pre = RECORDER.call("evidence", "bam_stats_at_locus",
                            {"bam_path": bam, "chromosome": body["chromosome"],
                             "start": max(0, pos - 250), "end": pos + 250})
        pe = tool_error(pre)
        if pe:
            return {"error": pe["error"], "error_type": pe.get("error_type"),
                    "call": pre["id"], "image_refs": [],
                    "hint": ("No panel was generated. The position could not be read, "
                             "so any image produced here would not correspond to a real locus.")}
        r = RECORDER.call("evidence", "evidence_panel",
                          {"bam_paths": [bam], "chromosome": body["chromosome"],
                           "position": pos})
        refs = _refs(r["result"])
        # Phase 11: a per-panel failure (IGV missing, IGV timeout) is reported
        # INSIDE result["panels"][layer]["error"], where the browser used to show
        # it only by dumping the whole return as JSON. Hoist it so the message a
        # user reads names the cause. Verified by hiding IGV, not by reading code.
        panel_errors = {}
        panels = (r["result"] or {}).get("panels") if isinstance(r["result"], dict) else None
        if isinstance(panels, dict):
            for layer, pan in panels.items():
                if isinstance(pan, dict) and pan.get("error"):
                    panel_errors[layer] = {"error": pan["error"],
                                           "searched": pan.get("searched"),
                                           "note": pan.get("note")}
        out = {"call": r["id"], "result": r["result"], "is_error": r["is_error"],
               "image_refs": refs}
        if panel_errors and not refs:
            first = next(iter(panel_errors.values()))
            out["error"] = first["error"]
            out["panel_errors"] = panel_errors
            out["hint"] = (
                f"No image was produced for any of the {len(panel_errors)} panel layer(s). "
                + ("IGV was looked for in $IGV_PATH and the built-in locations the startup "
                   "banner lists. " if first.get("searched") else "")
                + "The evidence numbers above are unaffected — they come from the counting "
                  "tools, not from IGV. Only the pictures are missing.")
            # The tool attaches a fixed note claiming an image WAS generated even on
            # a panel that failed. It is wrong here, so it is not passed on.
            out["tool_note_suppressed"] = bool(first.get("note"))
        elif panel_errors:
            out["panel_errors"] = panel_errors
        return out
    if path == "/api/compare":
        la, lb = body["label_a"], body["label_b"]
        ra = RECORDER.call("bridge", "load_candidate_set",
                           {"path": CANDIDATE_FILES[la], "label": la})
        rb = RECORDER.call("bridge", "load_candidate_set",
                           {"path": CANDIDATE_FILES[lb], "label": lb})
        if "error" in (ra["result"] or {}) or "error" in (rb["result"] or {}):
            return {"error": "could not load one of the two call sets",
                    "detail": [ra["result"], rb["result"]]}
        sa, sb = ra["result"]["set_id"], rb["result"]["set_id"]
        cp = {"set_a": sa, "set_b": sb}
        if body.get("tolerance_bp"):
            cp["tolerance_bp"] = int(body["tolerance_bp"])
        cmpr = RECORDER.call("bridge", "compare_candidate_sets", cp)
        cres = cmpr["result"] or {}
        if "error" in cres:
            return {"error": cres["error"], "call": cmpr["id"]}

        # Effect on each set's survivors. The recurrence decision is the tool's:
        # we only intersect the candidate ids IT reported as matched with the
        # survivor ids list_candidates returned. No matching logic is repeated.
        filt = {"limit": 20000}
        for k in ("svtype",):
            if body.get(k):
                filt[k] = body[k]
        for k in ("min_pe", "min_sr"):
            if body.get(k) not in (None, "", "null"):
                filt[k] = int(body[k])
        for k in ("filter_pass", "primary_only"):
            if body.get(k):
                filt[k] = True
        mask = _mask_state(bool(body.get("use_mask")))
        if mask["applied"]:
            filt["mask_path"] = MASK_PATH
        fa = RECORDER.call("bridge", "list_candidates", {"set_id": sa, **filt})
        fb = RECORDER.call("bridge", "list_candidates", {"set_id": sb, **filt})
        m_a = {p["candidate_id_a"] for p in cres.get("matched_pairs", [])}
        m_b = {p["candidate_id_b"] for p in cres.get("matched_pairs", [])}
        out = {}
        for tag, fr, matched, lbl in (("a", fa, m_a, la), ("b", fb, m_b, lb)):
            fe = tool_error(fr)
            if fe:
                return {"error": f"could not list candidates for {lbl}: {fe['error']}",
                        "error_type": fe.get("error_type"), "call": fr["id"]}
            cands = (fr["result"] or {}).get("candidates", [])
            rec = [c for c in cands if c["candidate_id"] in matched]
            uniq = [c for c in cands if c["candidate_id"] not in matched]
            out[tag] = {
                "label": lbl, "survivors": len(cands),
                "recurrent": len(rec), "unique": len(uniq),
                "unique_examples": [
                    {"chrom1": c["chrom1"], "pos1": c["pos1"], "chrom2": c["chrom2"],
                     "pos2": c["pos2"], "orientation": c["orientation"],
                     "svtype": c["svtype"], "pe": c["pe"], "sr": c["sr"]}
                    for c in uniq[:10]],
                "call": fr["id"],
            }
        return {"compare": cres, "compare_call": cmpr["id"],
                "load_calls": [ra["id"], rb["id"]], "survivor_effect": out,
                "mask": mask}
    if path == "/api/chat":
        tools, where = _chat_tools(trim=bool(body.get("trim", False)))
        model = body.get("model")
        if not model:
            return {"error": "no model selected"}
        # Same route, same schemas, same executor, same recorder for both
        # backends. Only the transport differs -- see chat.run_turn_api. The
        # model name selects it because "claude-*" is not an Ollama tag and an
        # explicit backend field would be one more thing to get out of sync
        # between the UI, the harness and the record.
        if model.startswith("claude-"):
            if not API_KEY:
                return {"error": "no Anthropic API key: set ANTHROPIC_API_KEY "
                                 "or place one at .api/claude_api_key"}
            return chatmod.run_turn_api(
                model, body.get("message", ""), tools, set(where),
                _chat_exec, max_iters=int(body.get("max_iters", 40)),
                max_tokens=int(body.get("max_tokens", 16000)),
                effort=body.get("effort", "high"),
                thinking=bool(body.get("thinking", True)),
                api_key=API_KEY)
        r = chatmod.run_turn(
            model, body.get("message", ""), tools, set(where),
            _chat_exec, num_ctx=int(body.get("num_ctx", 32768)),
            max_iters=int(body.get("max_iters", 8)),
            think=chatmod.normalise_think(body.get("think", False)))
        return r
    if path == "/api/calls":
        return {"calls": [{k: v for k, v in c.items() if k != "result"}
                          for c in RECORDER.calls]}
    if path == "/api/call":
        cid = int(body["id"])
        c = next((x for x in RECORDER.calls if x["id"] == cid), None)
        return c or {"error": "no such call"}
    return {"error": f"unknown endpoint {path}"}


def _refs(result):
    """Pull any image_ref values out of a screenshot/panel result."""
    out = []

    def walk(x):
        if isinstance(x, dict):
            if isinstance(x.get("image_ref"), str):
                out.append(x["image_ref"])
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
    walk(result)
    return out


def _resolve_handle(ref):
    """handle -> real path, via the session manifest. Server-side only."""
    d = bam_tools.image_session_dir()
    m = os.path.join(d, bam_tools.IMAGE_MANIFEST_NAME)
    if not os.path.exists(m):
        return None
    try:
        with open(m) as fh:
            return json.load(fh).get(ref)
    except (OSError, ValueError):
        return None


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def _send(self, code, body, ctype="application/json"):
        # Every JSON body leaves through here, so every one is scrubbed -- not only
        # the routes that remembered to: the error branch in do_POST returned raw
        # tracebacks, absolute paths included.
        b = body if isinstance(body, bytes) else json.dumps(scrub(body)).encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(b)))
        self.end_headers()
        self.wfile.write(b)

    def do_GET(self):
        u = urlparse(self.path)
        if u.path in ("/", "/index.html"):
            return self._send(200, PAGE.encode(), "text/html; charset=utf-8")
        if u.path == "/api/chat_models":
            tools, where = _chat_tools()
            # Phase 11: an empty list used to render as a bare "no models found",
            # which reads the same whether ollama is down or ollama is up with
            # nothing pulled. Those need different fixes, so say which.
            probed = probe_ollama()
            return self._send(200, {"models": probed or [],
                                    "reachable": probed is not None,
                                    "url": chatmod.OLLAMA,
                                    "why": (None if probed else
                                            (f"no models pulled — run: ollama pull qwen3.5:4b"
                                             if probed is not None else
                                             f"ollama is not reachable at {chatmod.OLLAMA} — start it "
                                             f"with `ollama serve`, or set SV_OLLAMA_URL if it runs "
                                             f"elsewhere")),
                                    "n_tools": len(tools),
                                    "tool_names": sorted(where)})
        if u.path == "/api/bootstrap":
            return self._send(200, scrub(_api("/api/bootstrap", {})))
        if u.path == "/api/calls":
            return self._send(200, scrub(_api("/api/calls", {})))
        if u.path == "/api/call":
            q = parse_qs(u.query)
            return self._send(200, scrub(_api("/api/call", {"id": q["id"][0]})))
        if u.path.startswith("/img/"):
            p = _resolve_handle(u.path[len("/img/"):])
            if not p or not os.path.exists(p):
                return self._send(404, {"error": "image handle not resolvable"})
            with open(p, "rb") as fh:
                return self._send(200, fh.read(), "image/png")
        return self._send(404, {"error": "not found"})

    def do_POST(self):
        u = urlparse(self.path)
        n = int(self.headers.get("Content-Length", 0))
        body = json.loads(self.rfile.read(n) or b"{}")
        try:
            return self._send(200, scrub(_api(u.path, body)))
        except Exception as e:
            return self._send(200, {"error": str(e), "error_type": "ui",
                                    "traceback": traceback.format_exc()[-800:]})


PAGE = r"""<!doctype html>
<html><head><meta charset="utf-8"><title>Breakpoint evidence — local instrument</title>
<style>
 :root{--fg:#111;--mut:#666;--line:#ccc;--bg:#fff;--warn:#8a4b00;--warnbg:#fff6e5;
       --held:#5b3a86;--heldbg:#f3eefa;--ok:#14532d;--code:#f6f6f6}
 *{box-sizing:border-box}
 body{margin:0;font:14px/1.5 ui-monospace,"DejaVu Sans Mono",Menlo,monospace;color:var(--fg);background:var(--bg)}
 header{border-bottom:2px solid var(--fg);padding:10px 16px;display:flex;gap:16px;align-items:baseline;flex-wrap:wrap}
 h1{font-size:15px;margin:0;font-weight:700}
 .mut{color:var(--mut)} .wrap{padding:16px;max-width:1180px}
 section{border:1px solid var(--line);margin-bottom:16px}
 section>h2{font-size:13px;margin:0;padding:6px 10px;background:#f2f2f2;border-bottom:1px solid var(--line);font-weight:700}
 .body{padding:10px}
 .two{display:grid;grid-template-columns:1fr 1fr;gap:16px}
 @media(max-width:900px){.two{grid-template-columns:1fr}}
 table{border-collapse:collapse;width:100%;font-size:13px}
 th,td{border:1px solid var(--line);padding:4px 6px;text-align:left;vertical-align:top}
 th{background:#fafafa;font-weight:700}
 button{font:inherit;padding:4px 10px;border:1px solid var(--fg);background:#fff;cursor:pointer}
 button:hover{background:#f0f0f0}
 button.primary{background:var(--fg);color:#fff}
 input,select{font:inherit;padding:3px 6px;border:1px solid var(--line)}
 .chip{border:1px solid var(--line);border-bottom:2px solid var(--mut);background:#fff;
       padding:0 4px;cursor:pointer;font:inherit}
 .chip:hover{background:#eef}
 .num{font-size:20px;font-weight:700}
 .layer{border:1px solid var(--line);padding:8px;margin-bottom:8px}
 .layer.held{background:var(--heldbg);border-left:5px solid var(--held)}
 .layer.na{background:#fafafa;color:var(--mut)}
 .badge{display:inline-block;padding:0 5px;border:1px solid;font-size:11px;text-transform:uppercase;letter-spacing:.4px}
 .b-held{color:var(--held);border-color:var(--held);background:var(--heldbg)}
 .b-auth{color:var(--warn);border-color:var(--warn);background:var(--warnbg)}
 .b-tool{color:var(--ok);border-color:var(--ok);background:#eefaf1}
 .b-ref{color:#1e3a8a;border-color:#1e3a8a;background:#eef2ff}
 pre{background:var(--code);border:1px solid var(--line);padding:8px;overflow:auto;max-height:340px;font-size:12px}
 .prov-cand{background:#eefaf1;border-left:5px solid var(--ok);padding:8px}
 .prov-hand{background:var(--warnbg);border-left:5px solid var(--warn);padding:8px}
 .ceil{background:#fafafa;border:1px dashed var(--mut);padding:8px;margin-top:8px}
 .err{background:#fdecec;border-left:5px solid #9b1c1c;padding:8px;white-space:pre-wrap}
 .hide{display:none} .row{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:8px}
 .notinert{background:var(--warnbg)}
 .msg-model{background:#fff;border:1px solid var(--line);border-left:5px solid #1e3a8a;padding:10px;margin:8px 0}
 .msg-tool{background:#eefaf1;border:1px solid #bcdccb;border-left:5px solid var(--ok);padding:8px;margin:8px 0 8px 24px}
 .msg-fail{background:#fdecec;border:1px solid #f0b4b4;border-left:5px solid #9b1c1c;padding:10px;margin:8px 0}
 .who{font-size:11px;text-transform:uppercase;letter-spacing:.5px;font-weight:700}
 .who-model{color:#1e3a8a} .who-tool{color:var(--ok)} .who-fail{color:#9b1c1c}
 .unsup{background:#ffd9d9;border-bottom:2px solid #9b1c1c;font-weight:700;padding:0 2px}
 .think{background:#f6f6f6;color:#555;font-size:12px;padding:6px;margin-top:6px;white-space:pre-wrap;max-height:150px;overflow:auto}
</style></head><body>
<header>
  <h1>Breakpoint evidence — local instrument</h1>
  <span class="mut" id="toolcount"></span>
  <button onclick="showLimits()">What this tool cannot tell you</button>
  <button onclick="showCalls()">Where every number came from</button>
</header>
<div class="wrap">

<section><h2>Start here — two ways to reach a breakpoint</h2><div class="body two">
  <div>
    <b>A · Work from a caller's candidate list</b>
    <div class="mut">Load structural-variant calls and narrow them with the filter chain.</div>
    <div class="row" style="margin-top:8px">
      <select id="cfile"></select><button class="primary" onclick="loadSet()">Load</button>
    </div>
    <div id="setinfo"></div>
  </div>
  <div>
    <b>B · Check a position you already have</b>
    <div class="mut">Type a breakpoint from a karyotype, a report, or another method and
      examine it directly. On test data this route found 44 discordant read pairs, 20
      soft-clipped reads and 14 split reads at a breakpoint the caller did not report.</div>
    <div class="row" style="margin-top:8px">
      <select id="hbam"></select>
      <input id="hchr" size="6" placeholder="chr20" value="chr20">
      <input id="hpos" size="12" placeholder="position">
      <button class="primary" onclick="handEnter()">Assess</button>
    </div>
  </div>
</div></section>

<section id="funnelsec" class="hide"><h2>Filter chain — what each step removes</h2><div class="body">
  <div class="row">
    <label>svtype <select id="fsv"><option value="">any</option><option>BND</option><option>DEL</option><option>DUP</option><option>INV</option></select></label>
    <label><input type="checkbox" id="fpass" checked> FILTER==PASS</label>
    <label>PE&ge; <input id="fpe" size="3" value="3"></label>
    <label>SR&ge; <input id="fsr" size="3" value="1"></label>
    <label><input type="checkbox" id="fprim" checked> both primary</label>
    <label><input type="checkbox" id="fmask" checked> not masked</label>
    <button class="primary" onclick="runFunnel()">Apply</button>
  </div>
  <div id="funnel"></div>
  <div id="cands"></div>
</div></section>

<section><h2>Compare two call sets — does a junction also appear in the other sample?</h2><div class="body">
  <div class="mut">A junction seen in an unrelated sample is more likely to be a recurring artefact
    of the sequencing and alignment than a finding specific to one person. This step is the single
    biggest reducer of the surviving list on real data. It needs a second sample, so it is skipped
    whenever only one is loaded — and when it is skipped, the surviving count is not comparable to
    one that included it.</div>
  <div class="row" style="margin-top:8px">
    <select id="cmpa"></select><span>vs</span><select id="cmpb"></select>
    <label>tolerance <input id="cmptol" size="5" value="500"> bp</label>
    <label>type <select id="cmpsv"><option value="">any</option><option>BND</option><option>DEL</option><option>DUP</option><option>INV</option></select></label>
    <button class="primary" onclick="runCompare()">Compare</button>
  </div>
  <div id="cmpout"></div>
</div></section>

<section><h2>Ask a local model — it can only quote what a tool returned</h2><div class="body">
  <div class="mut">The model has no direct access to any data. It may only request the same tools
    this page uses, its requests run through the same recorder, and every number it writes is
    checked against what those tools actually returned. Numbers with no tool call behind them are
    marked in red. The model is unreliable by assumption; this panel is built so that its mistakes
    are visible rather than plausible.</div>
  <div class="row" style="margin-top:8px">
    <select id="cmodel"></select>
    <label>context <select id="cctx">
      <option>8192</option><option>16384</option><option selected>32768</option></select></label>
    <label title="qwen3.5 thinking blocks consume the generation budget; measured to produce empty answers when context is tight"><input type="checkbox" id="cthink"> thinking mode</label>
    <span class="mut" id="ctoolinfo"></span>
  </div>
  <div class="row">
    <input id="cmsg" style="flex:1;min-width:380px" placeholder="e.g. Assess chr20:200000 in dataset IMP01 and tell me what the evidence shows">
    <button class="primary" onclick="sendChat()">Send</button>
  </div>
  <div id="chatout"></div>
</div></section>

<section id="evsec" class="hide"><h2>Evidence at this position</h2><div class="body" id="ev"></div></section>
<section id="igvsec" class="hide"><h2>IGV panel</h2><div class="body" id="igv"></div></section>
<section id="aux" class="hide"><h2 id="auxh"></h2><div class="body" id="auxb"></div></section>
</div>

<script>
let BOOT=null, SETID=null, LASTPOS=null;
const j=(u,b)=>fetch(u,{method:'POST',headers:{'Content-Type':'application/json'},
                        body:JSON.stringify(b)}).then(r=>r.json());
const g=(u)=>fetch(u).then(r=>r.json());
const esc=s=>String(s==null?'':s).replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]));
const provBadge=p=>{const c=p==='author judgement'?'b-auth':p==='tool-defined'?'b-tool':'b-ref';
  return `<span class="badge ${c}">${esc(p)}</span>`;};
const chip=(id,txt)=>`<button class="chip" onclick="showCall(${id})">${esc(txt)} <span class="mut">#${id}</span></button>`;

async function boot(){
  BOOT=await g('/api/bootstrap');
  document.getElementById('toolcount').textContent =
    `${BOOT.tool_counts.evidence} evidence tools + ${BOOT.tool_counts.bridge} bridge tools verified at startup`;
  document.getElementById('cfile').innerHTML=BOOT.candidate_files.map(x=>`<option>${esc(x)}</option>`).join('');
  document.getElementById('hbam').innerHTML=BOOT.datasets.map(x=>`<option>${esc(x)}</option>`).join('');
  const opts=BOOT.candidate_files.map(x=>`<option>${esc(x)}</option>`).join('');
  document.getElementById('cmpa').innerHTML=opts;
  document.getElementById('cmpb').innerHTML=opts;
  if(BOOT.candidate_files.length>1) document.getElementById('cmpb').selectedIndex=1;
}
async function loadSet(){
  const lbl=document.getElementById('cfile').value;
  const r=await j('/api/load',{candidates_label:lbl});
  if(r.is_error||!r.result||r.result.error){document.getElementById('setinfo').innerHTML=
     `<div class="err">${esc(JSON.stringify(r.result))}</div>`;return;}
  SETID=r.result.set_id;
  const s=r.result;
  document.getElementById('setinfo').innerHTML=
   `<div style="margin-top:8px">set <b>${esc(s.set_id)}</b> · convention <b>${esc(s.caller_convention)}</b>
    · ${s.total_records} records · dedup ${s.junctions_before_dedup}&rarr;${s.junctions_after_dedup}
    (${s.records_merged_by_dedup} merged) ${chip(r.call,'load_candidate_set')}</div>`;
  document.getElementById('funnelsec').classList.remove('hide');
  runFunnel();
}
async function runFunnel(){
  if(!SETID)return;
  const b={set_id:SETID,svtype:document.getElementById('fsv').value||null,
    filter_pass:document.getElementById('fpass').checked,
    min_pe:document.getElementById('fpe').value||null,
    min_sr:document.getElementById('fsr').value||null,
    primary_only:document.getElementById('fprim').checked,
    use_mask:document.getElementById('fmask').checked};
  const r=await j('/api/funnel',b);
  if(!r.result||r.result.error){document.getElementById('funnel').innerHTML=
    `<div class="err">${esc(JSON.stringify(r.result))}</div>`;return;}
  const res=r.result;
  let h='';
  if(r.mask&&r.mask.requested&&!r.mask.applied){
    h+=`<div class="err" style="margin-bottom:8px">Step not run: <b>exclude template</b>.
      ${esc(r.mask.reason)}</div>`;}
  h+=`<table><tr><th>step</th><th>cut-off</th><th>where the cut-off came from</th>
    <th>still remaining</th><th>removed by this step</th>
    <th>this step would remove on its own</th><th>measured against</th></tr>`;
  h+=`<tr><td colspan=3><i>all loaded calls</i></td><td><b>${res.total_in_set}</b></td><td>&mdash;</td><td>&mdash;</td><td>all</td></tr>`;
  for(const s of res.filters_applied){
    const gap=s.would_remove_from_unfiltered_set-s.removed_cumulatively_here;
    const inert=gap>0;
    h+=`<tr class="${inert?'notinert':''}"><td>${esc(s.filter)}</td><td>${esc(s.value)}</td>
      <td>${provBadge(s.provenance)}</td><td><b>${s.surviving_after_this_step}</b></td>
      <td>${s.removed_cumulatively_here}</td><td>${s.would_remove_from_unfiltered_set}</td>
      <td class="mut">${esc(s.unfiltered_set_scope)}</td></tr>`;
    if(inert){
      const zero = s.removed_cumulatively_here===0;
      h+=`<tr class="notinert"><td colspan=7 class="mut">&#9650; removed
        <b>${s.removed_cumulatively_here}</b> at this point in the chain, but would remove
        <b>${s.would_remove_from_unfiltered_set}</b> from the unfiltered set on its own
        (${esc(s.unfiltered_set_scope)}) — a gap of ${gap}.
        ${zero?'<b>This step looks inert and is not.</b>':'This step is doing more work than its cumulative column suggests.'}
        Earlier steps had already removed those records. Ordering, not irrelevance.</td></tr>`;}
  }
  h+=`</table><div class="mut" style="margin-top:6px">${esc(res.note)} ${chip(r.call,'list_candidates')}</div>`;
  document.getElementById('funnel').innerHTML=h;
  let c=`<h3 style="font-size:13px">Surviving candidates (${res.total_matching}, showing ${res.returned})</h3><table>
    <tr><th>id</th><th>breakpoint 1</th><th>breakpoint 2</th><th>type</th><th>orientation</th><th title="paired reads supporting this junction">paired-read support</th><th title="reads split across the junction">split-read support</th><th>caller filter</th><th title="how many caller records were merged into this one">merged</th><th></th></tr>`;
  for(const x of res.candidates){
    c+=`<tr><td class="mut">${esc(x.candidate_id)}</td><td>${esc(x.chrom1)}:${x.pos1}</td>
      <td>${esc(x.chrom2)}:${x.pos2}</td><td>${esc(x.svtype)}</td><td>${esc(x.orientation)}</td>
      <td>${x.pe}</td><td>${x.sr}</td><td>${esc(x.filter)}</td><td>${x.n_merged}</td>
      <td><button onclick="openCand('${esc(x.candidate_id)}')">evidence</button></td></tr>`;
  }
  document.getElementById('cands').innerHTML=c+'</table>';
}
async function openCand(cid){
  const r=await j('/api/candidate',{set_id:SETID,candidate_id:cid});
  const c=r.result;
  if(!c||c.error){alert(JSON.stringify(c));return;}
  const bam=document.getElementById('hbam').value;
  document.getElementById('evsec').classList.remove('hide');
  document.getElementById('ev').innerHTML=`<div class="mut">assessing both breakends of
    ${esc(c.candidate_id)} in dataset <b>${esc(bam)}</b> ${chip(r.call,'get_candidate')}</div>`;
  const out=[];
  for(const be of [c.breakend_1,c.breakend_2]){
    out.push(await j('/api/assess',{bam_label:bam,chromosome:be.chromosome,position:be.position}));
  }
  renderEvidence(out,`candidate ${c.candidate_id}`);
}
async function handEnter(){
  const bam=document.getElementById('hbam').value;
  const chr=document.getElementById('hchr').value.trim();
  const pos=parseInt(document.getElementById('hpos').value,10);
  if(!chr||!pos){alert('chromosome and position required');return;}
  document.getElementById('evsec').classList.remove('hide');
  document.getElementById('ev').innerHTML='<div class="mut">assessing…</div>';
  const r=await j('/api/assess',{bam_label:bam,chromosome:chr,position:pos});
  renderEvidence([r],`hand-entered ${chr}:${pos}`);
}
function layerBlock(L){
  if(L.error) return `<div class="layer"><b>${esc(L.label)}</b>
    <div class="err">This layer returned an error, not a result:\n${esc(L.error)}</div>
    ${chip(L.call,'tool call')}</div>`;
  const held=L.quality_limited===true;
  const na=L.assessable===false;
  const cls=held?'layer held':(na?'layer na':'layer');
  let v=`<span class="num">${L.count==null?'&mdash;':L.count}</span>`;
  let extra='';
  if(L.key==='discordant_pairs'&&L.partners) extra=`partners: ${esc(JSON.stringify(L.partners))}`;
  if(L.key==='split_reads'&&L.partners) extra=`partners: ${esc(JSON.stringify(L.partners))}`;
  if(L.key==='soft_clipped_reads') extra=`consensus clip position ${L.consensus==null?'&mdash;':L.consensus} · max clips at one position ${L.max_clips}`;
  if(L.key==='read_depth') extra=`min ${L.count} vs mean ${L.mean_depth} · ratio ${L.fraction} · dip at focus: ${L.dip_at_focus}`;
  return `<div class="${cls}">
    <b>${esc(L.label)}</b> ${held?'<span class="badge b-held">withheld — quality limited</span>':''}
    ${na?'<span class="badge">not assessable</span>':''}
    <div>${v} <span class="mut">${L.fraction==null?'':'fraction '+L.fraction}</span>
      ${chip(L.call,'tool call')}</div>
    <div class="mut">${extra}</div>
    <div class="mut">minimum mapping quality used: <b>${L.min_mapq===undefined||L.min_mapq===null?'no mapping-quality filter applied':L.min_mapq}</b></div>
    ${L.reason?`<div class="mut">reason: ${esc(L.reason)}</div>`:''}</div>`;
}
function ceilBlock(C,S){
  if(!C||!C.derivable) return `<div class="ceil"><b>Ceiling not derivable</b> — ${esc(C&&C.reason||'scoring function structure changed')}.
    No ceiling figure is shown rather than a stale one.</div>`;
  let rows='';
  for(const k in C.per_layer){const p=C.per_layer[k];
    rows+=`<tr><td>${esc(k)}</td><td>${esc(p.observed_field)}</td><td>${p.observed==null?'&mdash;':p.observed}</td>
      <td>${p.score_now==null?'&mdash;':p.score_now}</td><td>${p.max_score}</td>
      <td class="mut">${p.next_tier?('needs '+p.next_tier.op+' '+p.next_tier.threshold+' for '+p.next_tier.score):'at top tier'}</td></tr>`;}
  return `<div class="ceil"><b>Highest score reachable at this position</b>
    <table style="margin-top:6px"><tr><th>measurement</th><th>value used</th><th>measured</th><th>scores now</th><th>best possible</th><th>next band needs</th></tr>${rows}</table>
    <div style="margin-top:6px">Highest score if depth contributes nothing: <b>${C.max_with_flat_depth}</b>
      · highest if all four measurements were maximal: ${C.max_all_layers} · the "strong" band starts at ${C.strong_band}.</div>
    ${C.attainable_here==null?'':`<div style="margin-top:6px">Reachable <b>at this position</b>:
      <b>${C.attainable_here}</b> — ${esc(C.attainable_basis)}.
      ${C.attainable_here<C.strong_band?`<b>That is below the "strong" band of ${C.strong_band}, so this position cannot reach "strong" however the other measurements turn out.</b>`
        :'This position could reach the "strong" band if the read-based measurements were stronger.'}</div>`}
    <div class="mut">${esc(C.note)}</div>
    ${C.per_layer.read_depth&&C.per_layer.read_depth.extra_gate?`<div class="mut">depth measurement only counts when: ${esc(C.per_layer.read_depth.extra_gate)}</div>`:''}</div>`;
}
function lowMapq(E){
  // The figure and the quality gate share one window; the label names it from
  // the payload. null*100 is 0 in JS, so a missing figure must not print 0.0%.
  const S=E.summary||{}, W=S.low_mapq_window;
  if(S.low_mapq_fraction==null) return '<b>not measured</b>';
  const pct=(S.low_mapq_fraction*100).toFixed(1)+'%';
  if(!W) return `${pct} <span class="err">(window not reported)</span>`;
  return `${pct} over ${esc(E.chromosome)}:${W.start}–${W.end} `
    +`(±${W.half_width_bp} bp, the window the quality gate uses)`;
}
function renderEvidence(list,title){
  let h=`<h3 style="font-size:13px">${esc(title)}</h3>`;
  for(const E of list){
    if(E.error){h+=`<div class="err"><b>Could not assess this position.</b>\n${esc(E.error)}
      ${E.hint?'\n\n'+esc(E.hint):''}
      ${E.contigs_in_header_sample?'\n\nContigs this dataset does contain (first few): '+esc(E.contigs_in_header_sample.join(', ')):''}
      </div>${E.call?chip(E.call,'tool call'):''}`;continue;}
    LASTPOS=E;
    const pp=E.position_provenance||{};
    const isSet=pp.source==='candidate_set';
    h+=`<div style="border:1px solid #999;padding:10px;margin-bottom:14px">
      <div class="${isSet?'prov-cand':'prov-hand'}">
        <b>${isSet?'This position came from the loaded call set.':'This position was entered by hand.'}</b>
        ${isSet?`It is candidate ${esc(pp.candidate_id)} in <b>${esc(pp.set_label)}</b>.`
               :'No caller proposed it, so nothing here corroborates that it is a real junction — the evidence below stands on its own.'}
        <div class="mut">recorded as <code>position_provenance = ${esc(pp.source||'unknown')}</code></div>
        <div class="mut">${esc(pp.note||'')}</div></div>
      <h4 style="margin:10px 0 4px">${esc(E.chromosome)}:${E.position} <span class="mut">· dataset ${esc(E.bam_label)}
        · applicable layers ${esc(JSON.stringify(E.applicable_layers))} ${chip(E.applicable_call,'applicable_layers')}
        · ambiguously mapped reads (MAPQ &lt; 20): ${lowMapq(E)} ${chip(E.stats.call,'bam_stats')}</span></h4>
      <div class="mut" style="margin-bottom:6px">The four measurements below are the result. The
        combined score after them summarises them; it does not replace them.</div>`;
    for(const L of E.layers) h+=layerBlock(L);
    const S=E.summary;
    if(E.summary_error){
      h+=`<div class="layer"><b>Combined score</b>
        <div class="err">The summary step returned an error, so no score was calculated:\n${esc(E.summary_error)}</div>
        ${chip(S.call,'tool call')}</div>`;
      h+=ceilBlock(E.ceiling,S);
      h+=`<div style="margin-top:8px"><button onclick="igv('${esc(E.bam_label)}','${esc(E.chromosome)}',${E.position})">generate IGV panel</button></div></div>`;
      continue;
    }
    const held=S.evidence_strength==='QUALITY-LIMITED';
    h+=`<div class="${held?'layer held':'layer'}" style="margin-top:10px">
      <b>Combined score</b> <span class="mut">— a summary of the four measurements above</span> ${chip(S.call,'breakpoint_evidence_summary')}
      <div>${held?'<span class="badge b-held">score withheld — quality limited</span> <span class="mut">not a low score: too many reads here are ambiguously mapped, so no combined score was calculated</span>'
        :`<span class="num">${S.evidence_score}</span> <b>${esc(S.evidence_strength)}</b>
          <span class="mut">raw ${S.evidence_score_raw} · layers ${esc(S.signal_layers)} · min_mapq_applied ${S.min_mapq_applied}</span>`}</div>
      <div class="mut">components: ${esc(JSON.stringify(S.components))}</div></div>`;
    h+=ceilBlock(E.ceiling,S);
    h+=`<div style="margin-top:8px"><button onclick="igv('${esc(E.bam_label)}','${esc(E.chromosome)}',${E.position})">generate IGV panel</button></div>`;
    h+=`</div>`;
  }
  document.getElementById('ev').innerHTML=h;
}
async function igv(bam,chr,pos){
  document.getElementById('igvsec').classList.remove('hide');
  document.getElementById('igv').innerHTML='<div class="mut">launching IGV — this can take a while…</div>';
  const r=await j('/api/igv',{bam_label:bam,chromosome:chr,position:pos});
  let h=`<div>${chip(r.call,'evidence_panel')}</div>`;
  if(r.error){h+=`<div class="err"><b>No panel was generated.</b>\n${esc(r.error)}
      ${r.hint?'\n\n'+esc(r.hint):''}</div>`;
    if(r.panel_errors){h+=`<div class="mut">per layer: ${esc(Object.keys(r.panel_errors).join(', '))}</div>`;}
    document.getElementById('igv').innerHTML=h;return;}
  if(r.is_error||(r.result&&r.result.error)||!(r.image_refs||[]).length){
    h+=`<div class="err"><b>No image was produced.</b> The tool's raw return is below; nothing is
      being shown in its place.\n\n${esc(JSON.stringify(r.result,null,1))}</div>`;
  }else{
    for(const ref of r.image_refs)
      h+=`<div style="margin-top:8px"><div class="mut">image_ref ${esc(ref)} — rendered from the
        session manifest; the model is never given this path.</div>
        <img src="/img/${encodeURIComponent(ref)}" style="max-width:100%;border:1px solid #999"></div>`;
  }
  document.getElementById('igv').innerHTML=h;
}
async function loadModels(){
  const r=await g('/api/chat_models');
  document.getElementById('cmodel').innerHTML=(r.models||[]).map(m=>`<option>${esc(m)}</option>`).join('')
    || `<option value="">${esc(r.reachable?'no models installed':'chat unavailable')}</option>`;
  document.getElementById('ctoolinfo').innerHTML=
    `${r.n_tools} tool schemas generated from the servers and offered to the model`
    + (r.why?`<div class="msg-fail" style="margin-top:6px"><div class="who who-fail">chat panel unavailable</div>${esc(r.why)}</div>`:'');
}
function markUnsupported(text, details){
  if(!details||!details.length) return esc(text);
  const bad=new Set(details.filter(d=>!d.supported).map(d=>d.text));
  if(!bad.size) return esc(text);
  let out=esc(text);
  for(const b of bad){
    const re=new RegExp('(?<![\\w.])'+b.replace('.','\\.')+'(?![\\w.])','g');
    out=out.replace(re,`<span class="unsup" title="no tool call returned this number">${b}</span>`);
  }
  return out;
}
async function sendChat(){
  const model=document.getElementById('cmodel').value;
  const msg=document.getElementById('cmsg').value.trim();
  if(!model||!msg){alert('choose a model and type a question');return;}
  const out=document.getElementById('chatout');
  out.innerHTML='<div class="mut">running — tool calls execute as the model requests them…</div>';
  const t0=Date.now();
  const r=await j('/api/chat',{model:model,message:msg,num_ctx:parseInt(document.getElementById('cctx').value,10),
                            think:document.getElementById('cthink').checked});
  if(r.error){out.innerHTML=`<div class="msg-fail"><div class="who who-fail">error</div>${esc(r.error)}</div>`;return;}
  let h='';
  for(const e of r.events){
    if(e.type==='model'){
      const isFinal=(e.content===r.final_text&&e.n_tool_calls===0);
      const txt=isFinal?markUnsupported(e.content,r.verification.details):esc(e.content);
      h+=`<div class="msg-model"><div class="who who-model">model — prose, not verified except where marked</div>
        <div>${txt||'<span class="mut">(no prose; requested '+e.n_tool_calls+' tool call(s))</span>'}</div>
        ${e.thinking?`<div class="think">thinking: ${esc(e.thinking)}</div>`:''}</div>`;
    } else if(e.type==='tool'){
      const bad=e.rejected||e.is_error;
      h+=`<div class="${bad?'msg-fail':'msg-tool'}">
        <div class="who ${bad?'who-fail':'who-tool'}">tool return — ${esc(e.name)}${bad?' (rejected/error)':''}</div>
        <div class="mut">params: ${esc(JSON.stringify(e.params))}</div>
        <pre>${esc(JSON.stringify(e.result,null,1).slice(0,1400))}</pre>
        ${e.call_id?chip(e.call_id,'recorded call'):'<span class="mut">not recorded — never executed</span>'}</div>`;
    } else {
      h+=`<div class="msg-fail"><div class="who who-fail">transport error</div>${esc(e.detail)}</div>`;
    }
  }
  if((r.context_truncated||[]).length)
    h+=`<div class="msg-fail"><div class="who who-fail">WARNING — the model ran out of room to answer</div>
      ${r.context_truncated.length} model turn(s) stopped because the context limit was reached, not because
      the model had finished. At the last turn the prompt used <b>${r.final_prompt_tokens}</b> of
      <b>${r.num_ctx}</b> tokens, leaving <b>${r.final_headroom}</b> to write with.
      Raise the context, or turn thinking off, before treating this as a result about the model.</div>`;
  if(r.ended_without_answer)
    h+=`<div class="msg-fail"><div class="who who-fail">FAILURE — no answer produced</div>
      The model kept requesting tools until the iteration limit and never wrote a conclusion.
      There is no answer to verify; the tool returns above stand on their own.</div>`;
  for(const f of r.text_tool_call_failures)
    h+=`<div class="msg-fail"><div class="who who-fail">FAILURE — tool call printed as text, not executed</div>
      The model wrote a tool call into its prose instead of emitting a structured call
      (${esc(f.reason)}). Nothing was run. The text is shown as evidence of the failure, not as an answer:
      <pre>${esc(f.content)}</pre></div>`;
  for(const m of r.malformed_tool_calls)
    h+=`<div class="msg-fail"><div class="who who-fail">FAILURE — malformed tool call</div>
      ${esc(m.name)}: ${esc(m.reason)}<pre>${esc(m.raw)}</pre></div>`;
  const v=r.verification;
  h+=`<div class="${v.unsupported.length?'msg-fail':'msg-tool'}" style="margin-left:0">
    <div class="who ${v.unsupported.length?'who-fail':'who-tool'}">verification pass</div>
    ${v.no_tool_calls?'<b>The model answered without calling any tool.</b> Nothing it wrote is backed by this session. ':''}
    ${v.numbers_in_prose} number(s) in the final answer;
    <b>${v.unsupported.length}</b> with no tool call behind them${v.unsupported.length?': '+esc(v.unsupported.join(', ')):''}.
    <div class="mut">${r.n_tool_calls} tool call(s) recorded · ${r.iterations} model turn(s) ·
      ${r.gen_tokens} tokens generated${r.gen_tokens_per_s?' at '+r.gen_tokens_per_s+' tok/s':''} ·
      ${r.wall_s}s wall · context ${r.num_ctx}</div></div>`;
  out.innerHTML=h;
}
async function runCompare(){
  const b={label_a:document.getElementById('cmpa').value,label_b:document.getElementById('cmpb').value,
    tolerance_bp:document.getElementById('cmptol').value||null,
    svtype:document.getElementById('cmpsv').value||null,
    filter_pass:true,min_pe:3,min_sr:1,primary_only:true,use_mask:true};
  document.getElementById('cmpout').innerHTML='<div class="mut">comparing…</div>';
  const r=await j('/api/compare',b);
  if(r.error){document.getElementById('cmpout').innerHTML=`<div class="err">${esc(r.error)}</div>`;return;}
  const c=r.compare, se=r.survivor_effect;
  const th=(c.thresholds_applied||[])[0]||{};
  let h=`<table style="margin-top:10px"><tr><th>call set</th><th>junctions in set</th>
    <th>also present in the other set</th><th>present only here</th></tr>
    <tr><td><b>${esc(c.label_a)}</b></td><td>${c.total_in_a}</td><td>${c.matched_in_a}</td><td><b>${c.unmatched_in_a}</b></td></tr>
    <tr><td><b>${esc(c.label_b)}</b></td><td>${c.total_in_b}</td><td>${c.matched_in_b}</td><td><b>${c.unmatched_in_b}</b></td></tr></table>
    <div class="mut" style="margin-top:6px">Two breakpoints count as the same junction when both ends
      agree within <b>${th.value}</b> bp. ${provBadge(th.provenance||'author judgement')}
      ${chip(r.compare_call,'compare_candidate_sets')}</div>
    <div class="mut">${esc(c.note||'')}</div>`;
  h+=`<h3 style="font-size:13px;margin-top:12px">Effect on the surviving list</h3>
    <table><tr><th>call set</th><th>surviving after the filter chain</th>
      <th>of those, also in the other set</th><th>remaining if recurrent ones are dropped</th></tr>`;
  for(const k of ['a','b']){const x=se[k];
    h+=`<tr><td><b>${esc(x.label)}</b></td><td>${x.survivors}</td><td>${x.recurrent}</td>
      <td><b>${x.unique}</b></td><td>${chip(x.call,'list_candidates')}</td></tr>`;}
  h+='</table>';
  for(const k of ['a','b']){const x=se[k];
    if(!x.unique_examples.length) continue;
    h+=`<div style="margin-top:8px"><b>${esc(x.label)} — present only in this set</b><table>
      <tr><th>breakpoint 1</th><th>breakpoint 2</th><th>type</th><th>orientation</th><th>paired-read</th><th>split-read</th></tr>`;
    for(const u of x.unique_examples)
      h+=`<tr><td>${esc(u.chrom1)}:${u.pos1}</td><td>${esc(u.chrom2)}:${u.pos2}</td>
        <td>${esc(u.svtype)}</td><td>${esc(u.orientation)}</td><td>${u.pe}</td><td>${u.sr}</td></tr>`;
    h+='</table></div>';}
  document.getElementById('cmpout').innerHTML=h;
}
async function showCall(id){
  const c=await g('/api/call?id='+id);
  aux(`Tool call #${id}`,`<div><b>${esc(c.server)} · ${esc(c.tool)}</b> · ${c.ms} ms
    · is_error ${c.is_error}</div><h4>parameters</h4><pre>${esc(JSON.stringify(c.params,null,1))}</pre>
    <h4>raw return</h4><pre>${esc(JSON.stringify(c.result,null,1))}</pre>`);
}
async function showCalls(){
  const r=await g('/api/calls');
  let h='<table><tr><th>#</th><th>server</th><th>tool</th><th>ms</th><th>error</th><th></th></tr>';
  for(const c of r.calls) h+=`<tr><td>${c.id}</td><td>${esc(c.server)}</td><td>${esc(c.tool)}</td>
    <td>${c.ms}</td><td>${c.is_error}</td><td><button onclick="showCall(${c.id})">open</button></td></tr>`;
  aux(`Tool call log — every number displayed came from one of these ${r.calls.length} calls`,h+'</table>');
}
function showLimits(){
  const L=BOOT.limits;
  aux(L.title,L.items.map(i=>`<div style="margin-bottom:10px"><b>${esc(i.h)}</b>
    <div class="mut">${esc(i.b)}</div></div>`).join(''));
}
function aux(t,h){document.getElementById('aux').classList.remove('hide');
  document.getElementById('auxh').textContent=t;document.getElementById('auxb').innerHTML=h;
  document.getElementById('aux').scrollIntoView({behavior:'smooth'});}
boot(); loadModels();
</script></body></html>"""


# The exclude template. Phase 11: this used to be a hardcoded path guarded by a
# bare os.path.exists, so on a machine without it the mask step silently did not
# run and the filter chain reported one fewer step with no explanation — the same
# silent-wrong-answer class this project keeps finding. It is now resolved through
# config, and its absence is REPORTED in the funnel and in the startup banner.
MASK_STATUS = CFG.status("exclude_template", kind="file")
MASK_PATH = MASK_STATUS["path"] or ""


def _mask_state(requested):
    """What the mask step did, and why. Returned to the browser on every funnel
    and compare call so a missing template can never be invisible. The template
    is named by its file name only: this used to send the absolute path, in
    `path` and inside `reason`, against the rule that the browser is never sent
    a filesystem path. The startup banner still shows the full path locally."""
    name = os.path.basename(MASK_STATUS["path"]) if MASK_STATUS["path"] else None
    source = _ABS_PATH.sub(_unpath, MASK_STATUS["source"])
    where = {"path": name, "source": source}
    if not requested:
        return {"requested": False, "applied": False, "reason": None, **where}
    if MASK_STATUS["found"]:
        return {"requested": True, "applied": True, "reason": None, **where}
    why = f"no such file: {name}" if name else "not configured"
    return {"requested": True, "applied": False,
            "reason": (f"exclude template not found ({why}); "
                       f"resolved from {source}. The chain below ran "
                       f"WITHOUT this step — nothing was removed for overlapping a "
                       f"known-problematic region."),
            **where}


def _load_api_key():
    """Anthropic key from the environment, else the gitignored .api/ file.
    Read once at import; never rendered, never returned by any route."""
    k = (os.environ.get("ANTHROPIC_API_KEY") or "").strip()
    if k:
        return k
    for cand in (os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              ".api", "claude_api_key"),
                 os.path.expanduser("~/.api/claude_api_key")):
        try:
            with open(cand) as f:
                k = f.read().strip()
            if k:
                return k
        except OSError:
            pass
    return ""


API_KEY = _load_api_key()


DATA_DIR_STATUS = CFG.status("data_dir", kind="dir")


def discover_public():
    """Auto-register the public demo data only. ~/patient_data is never touched."""
    # Explicit registrations win: they are named by a human who knows where the
    # data is. Autodiscovery then fills in whatever else is under the data dir.
    for label, path in CFG.registered("datasets").items():
        DATASETS.setdefault(label, path)
    for label, path in CFG.registered("candidates").items():
        CANDIDATE_FILES.setdefault(label, path)
    pub = DATA_DIR_STATUS["path"] or ""
    if not pub or not os.path.isdir(pub):
        return
    # A file named explicitly in the config must not also appear under a label
    # derived from its filename: the demo bundle registered demo.bam as DEMO and
    # autodiscovery then added the same file again as "demo". Two labels for one
    # file is not wrong, it is just confusing — and confusing is the failure mode
    # this whole interface exists to avoid.
    known = set(DATASETS.values())
    for root in (pub, os.path.join(pub, "sim", "bams")):
        if not os.path.isdir(root):
            continue
        for f in sorted(os.listdir(root)):
            if f.endswith(".bam"):
                full = os.path.join(root, f)
                if os.path.realpath(full) in {os.path.realpath(k) for k in known}:
                    continue
                DATASETS.setdefault(os.path.splitext(f)[0], full)
    for root in (os.path.join(pub, "delly"), os.path.join(pub, "sim", "delly")):
        if not os.path.isdir(root):
            continue
        for f in sorted(os.listdir(root)):
            if f.endswith(".bcf") or f.endswith(".vcf") or f.endswith(".vcf.gz"):
                full = os.path.join(root, f)
                if os.path.realpath(full) in {os.path.realpath(k) for k in CANDIDATE_FILES.values()}:
                    continue
                CANDIDATE_FILES.setdefault(f.split(".")[0], full)


# The panel tool's own search (bam_tools.run_igv_screenshot): $IGV_PATH first,
# then these. The tool never reads the config file, so main() exports a
# configured igv path into IGV_PATH; otherwise the banner could report an IGV the
# tool cannot find, which it once did.
IGV_CANDIDATES = ["~/IGV_2.17.4/igv.sh", "~/igv/igv.sh", "/opt/igv/igv.sh"]


def _java_major(java):
    """(major version, None) for the java binary at `java`, or (None, why)."""
    try:
        r = subprocess.run([java, "-version"], capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError) as e:
        return None, f"{java} could not be run ({type(e).__name__})"
    if r.returncode != 0:
        return None, f"{java} -version exited {r.returncode}"
    m = re.search(r'version "(\d+)(?:\.(\d+))?', r.stderr + r.stdout)
    if not m:
        return None, f"{java} -version printed no version"
    major = int(m.group(1))
    return (int(m.group(2)) if major == 1 and m.group(2) else major), None


def _igv_java_requirement(igv_dir):
    """(Java major version IGV's classes were compiled for, None), or (None, why).
    Read from a class-file header in lib/igv.jar (class major version - 44), so
    the requirement comes from the installed IGV, not from a literal."""
    jar = os.path.join(igv_dir, "lib", "igv.jar")
    try:
        with zipfile.ZipFile(jar) as z:
            classes = [n for n in z.namelist() if n.endswith(".class")]
            pick = next((n for n in classes if n.startswith("org/broad/igv/")),
                        classes[0] if classes else None)
            if pick is None:
                return None, "lib/igv.jar contains no classes"
            head = z.read(pick)[:8]
    except (OSError, zipfile.BadZipFile) as e:
        return None, f"lib/igv.jar not readable ({type(e).__name__})"
    if len(head) < 8 or head[:4] != b"\xca\xfe\xba\xbe":
        return None, "lib/igv.jar does not hold Java class files"
    return int.from_bytes(head[6:8], "big") - 44, None


def verify_full():
    """Every condition the IGV panels need, each checked so that it can fail.

    Phase 11 reported FULL whenever a file existed at an igv.sh path: an empty,
    non-executable file with no Java anywhere showed [x]. The panel tool runs
    igv.sh directly (so it must be executable), igv.sh needs IGV's lib/igv.jar
    and a Java it will pick up -- its bundled jdk*/bin/java, else java on PATH --
    at least as new as IGV was compiled for, and IGV's window needs a display:
    the tool inherits DISPLAY and does not wrap itself in xvfb-run. No test
    render is attempted (that is IGV's whole startup), so [x] means these
    prerequisites were verified, and the banner says exactly that.
    Returns {"ok", "detail", "conditions": [(name, ok, detail), ...]}."""
    conds = []
    search = ([os.environ["IGV_PATH"]] if os.environ.get("IGV_PATH") else []) + \
        [os.path.expanduser(c) for c in IGV_CANDIDATES]
    igv = next((c for c in search if os.path.exists(c)), None)
    conds.append(("igv.sh found by the panel tool's own search", bool(igv),
                  igv or "searched " + ", ".join(search)))
    java_note = ""
    if igv:
        runnable = os.path.isfile(igv) and os.access(igv, os.X_OK)
        conds.append(("igv.sh is executable", runnable,
                      "" if runnable else "the tool runs it directly; it needs the execute bit"))
        home = os.path.dirname(os.path.realpath(igv))
        need, why = _igv_java_requirement(home)
        conds.append(("IGV's lib/igv.jar is readable", need is not None,
                      why or f"compiled for Java {need}"))
        bundled = sorted(glob.glob(os.path.join(home, "jdk*", "bin", "java")))
        java = bundled[0] if bundled else shutil.which("java")
        if not java:
            conds.append(("a Java runtime igv.sh will use", False,
                          "no bundled jdk*/bin/java next to igv.sh and no java on PATH"))
        else:
            have, why = _java_major(java)
            conds.append(("a Java runtime igv.sh will use", have is not None,
                          why or f"Java {have} ({'bundled with IGV' if bundled else 'on PATH'})"))
            if have is not None and need is not None:
                conds.append(("that Java is new enough for this IGV", have >= need,
                              f"Java {have}; IGV was compiled for Java {need}"))
            java_note = f"Java {have}" if have is not None else ""
    disp = os.environ.get("DISPLAY", "")
    local = re.match(r":(\d+)", disp)
    if not disp:
        conds.append(("a display for IGV's window", False,
                      "DISPLAY is not set (WSLg: export DISPLAY=:0; no display: run under xvfb-run)"))
    elif local:
        sock = f"/tmp/.X11-unix/X{local.group(1)}"
        conds.append(("a display for IGV's window", os.path.exists(sock),
                      f"DISPLAY={disp}" + ("" if os.path.exists(sock) else f", but {sock} does not exist")))
    else:
        conds.append(("a display for IGV's window", True, f"DISPLAY={disp} (remote display, not probed)"))
    ok = all(c[1] for c in conds)
    if ok:
        detail = (f"IGV panels — prerequisites verified: {igv} executable, {java_note}, "
                  f"display {disp} (no test render)")
    else:
        first = next(c for c in conds if not c[1])
        detail = (f"IGV panels — NOT AVAILABLE: {first[0]}: {first[2]}; "
                  f"panels will report the failure instead of rendering")
    return {"ok": ok, "detail": detail, "conditions": conds}


def _selftest_fixture(d):
    """A tiny BAM and candidate VCF for the MINIMAL self-test, so it needs no data
    on disk. Paired 100 bp reads every 10 bp over chr1:8,000-12,000; the reads
    starting within 150 bp of 10,000 have their mate on chr2, a 30 bp soft clip
    and an SA tag, so all four layers have reads to assess."""
    import pysam
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 100000}, {"SN": "chr2", "LN": 100000}]})
    bam = os.path.join(d, "selftest.bam")
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for i, start in enumerate(range(8000, 12000, 10)):
            r = pysam.AlignedSegment(header)
            r.query_name = f"s{i}"
            r.query_sequence = "ACGT" * 25
            r.query_qualities = pysam.qualitystring_to_array("I" * 100)
            r.reference_id, r.reference_start, r.mapping_quality = 0, start, 60
            if abs(start - 10000) <= 150:
                r.flag = 0x1
                r.cigar = [(4, 30), (0, 70)]
                r.next_reference_id, r.next_reference_start = 1, 5000
                r.set_tag("SA", "chr2,5000,+,30M70S,60,0;")
            else:
                r.flag = 0x1 | 0x2
                r.cigar = [(0, 100)]
                r.next_reference_id, r.next_reference_start = 0, start + 200
                r.template_length = 300
            out.write(r)
    pysam.index(bam)
    vcf = os.path.join(d, "selftest.vcf")
    with open(vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n"
                "##contig=<ID=chr1,length=100000>\n##contig=<ID=chr2,length=100000>\n"
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
                '##INFO=<ID=END,Number=1,Type=Integer,Description="end">\n'
                '##INFO=<ID=CHR2,Number=1,Type=String,Description="partner contig">\n'
                '##INFO=<ID=POS2,Number=1,Type=Integer,Description="partner position">\n'
                '##INFO=<ID=CT,Number=1,Type=String,Description="connection type">\n'
                '##INFO=<ID=PE,Number=1,Type=Integer,Description="paired-end support">\n'
                '##INFO=<ID=SR,Number=1,Type=Integer,Description="split-read support">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t20000\tdel1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=21000;PE=4;SR=0\n"
                "chr2\t5000\tbnd1\tN\t<BND>\t.\tPASS\tSVTYPE=BND;CHR2=chr1;POS2=10000;CT=3to5;PE=10;SR=5\n")
    return bam, vcf


def verify_minimal():
    """Every condition the MINIMAL line claims, exercised so that it can fail.

    Phase 11's --check hard-coded ok=True and ANDed it with a count compared
    against >= 0, so no configuration could make it report FAILED. Each entry
    here runs the thing it names: the tool contract; the tier derivation; a
    hand-entered coordinate through the same assess() the page uses, on a
    fixture BAM (all four layers must come back assessable and the summary must
    give a strength); a candidate set through load and the filter chain; and
    every registered file must exist, each BAM with an index. The fixture's tool
    calls go to a private recorder, so they never appear in the session's log.
    Returns [(condition, ok, detail), ...]."""
    global RECORDER
    out = []
    try:
        ne, nb = assert_tool_contract()
        out.append(("tool contract", True, f"{ne} evidence + {nb} bridge tools"))
    except SystemExit as e:
        out.append(("tool contract", False, str(e)))
    out.append(("scoring tiers derivable from bam_tools' source",
                TIERS is not None and BANDS is not None, TIER_ERROR or ""))
    label, saved = "__selftest__", RECORDER
    RECORDER = ToolRecorder()
    try:
        with tempfile.TemporaryDirectory() as d:
            bam, vcf = _selftest_fixture(d)
            DATASETS[label], CANDIDATE_FILES[label] = bam, vcf
            name = "four evidence layers at a hand-entered coordinate"
            try:
                E = assess(label, "chr1", 10000)
                layers = {L["key"]: L for L in E.get("layers", [])}
                bad = [k for k in ("discordant_pairs", "soft_clipped_reads", "split_reads", "read_depth")
                       if k not in layers or layers[k].get("error") or layers[k].get("assessable") is not True]
                strength = (E.get("summary") or {}).get("evidence_strength")
                out.append((name, not E.get("error") and not E.get("summary_error") and not bad
                            and strength is not None,
                            E.get("error") or E.get("summary_error")
                            or (f"not assessable or errored: {bad}" if bad else f"summary strength {strength!r}")))
            except Exception as e:
                out.append((name, False, f"{type(e).__name__}: {e}"))
            name = "candidate set through load and the filter chain"
            try:
                res = (_api("/api/load", {"candidates_label": label}).get("result") or {})
                fres = ((_api("/api/funnel", {"set_id": res["set_id"], "filter_pass": True, "limit": 10})
                         .get("result") or {}) if res.get("set_id") else {})
                out.append((name, "error" not in res and res.get("total_records") == 2
                            and "error" not in fres and fres.get("total_matching") == 2,
                            res.get("error") or fres.get("error")
                            or f"{res.get('total_records')} records loaded, {fres.get('total_matching')} "
                               f"pass filter_pass (the fixture has 2 PASS records)"))
            except Exception as e:
                out.append((name, False, f"{type(e).__name__}: {e}"))
    except Exception as e:
        out.append(("self-test fixture written", False, f"{type(e).__name__}: {e}"))
    finally:
        DATASETS.pop(label, None)
        CANDIDATE_FILES.pop(label, None)
        RECORDER = saved
    problems = []
    for kind, reg in (("dataset", DATASETS), ("candidate set", CANDIDATE_FILES)):
        for lbl, p in sorted(reg.items()):
            if re.match(r"(https?|ftp)://", p):
                continue
            if not (os.path.isfile(p) and os.access(p, os.R_OK)):
                problems.append(f"{kind} {lbl}: file missing or unreadable")
            elif kind == "dataset" and not any(os.path.isfile(b + x) for b in (p, os.path.splitext(p)[0])
                                               for x in (".bai", ".csi", ".crai")):
                problems.append(f"dataset {lbl}: no .bai/.csi/.crai index beside it")
    out.append(("every registered file exists (each BAM indexed)", not problems,
                "; ".join(problems) or f"{len(DATASETS)} dataset(s), {len(CANDIDATE_FILES)} candidate set(s)"))
    return out


def probe_ollama(timeout=1.5):
    try:
        with urllib.request.urlopen(chatmod.OLLAMA + "/api/tags", timeout=timeout) as r:
            return sorted(m["name"] for m in json.loads(r.read()).get("models", []))
    except Exception:
        return None


def capability_report(minimal):
    """What tier this install can actually deliver. Reported at startup so a user
    knows before they click, not after an empty panel. `minimal` is
    verify_minimal()'s result: MINIMAL is [x] only if every condition passed."""
    full = verify_full()
    models = probe_ollama()
    failed = next((c for c in minimal if not c[1]), None)
    return {
        "MINIMAL": {"ok": failed is None, "needs": "python + pysam + fastmcp + requests",
                    "detail": (f"filter chain, four evidence layers, hand-entered coordinates — "
                               f"self-test passed ({len(minimal)} conditions)" if failed is None else
                               f"SELF-TEST FAILED: {failed[0]}: {failed[2]}"),
                    "conditions": minimal},
        "FULL": {"ok": full["ok"], "needs": "IGV, Java and a display",
                 "detail": full["detail"], "conditions": full["conditions"]},
        "COMPLETE": {"ok": bool(models), "needs": "ollama serving at " + chatmod.OLLAMA,
                     "detail": (f"local model chat — {len(models)} model(s): {', '.join(models[:4])}"
                                if models else "ollama not reachable; the chat panel will be unavailable")},
        "exclude_template": MASK_STATUS,
        "data_dir": DATA_DIR_STATUS,
        "api_key": bool(API_KEY),
    }


def print_banner(cap, port):
    if not cap["MINIMAL"]["ok"]:
        tier = "NONE — the MINIMAL self-test failed"
    else:
        tier = "COMPLETE" if cap["COMPLETE"]["ok"] and cap["FULL"]["ok"] else (
               "FULL" if cap["FULL"]["ok"] else "MINIMAL")
    print("", flush=True)
    print(f"  available tier: {tier}", flush=True)
    for name in ("MINIMAL", "FULL", "COMPLETE"):
        c = cap[name]
        print(f"    [{'x' if c['ok'] else ' '}] {name:9s} {c['detail']}", flush=True)
    m = cap["exclude_template"]
    print(f"    [{'x' if m['found'] else ' '}] exclude template  "
          f"{m['path'] or '(not configured)'}  [{m['source']}]"
          + ("" if m["found"] else "  -- the mask filter step will NOT run and will say so"),
          flush=True)
    d = cap["data_dir"]
    print(f"    [{'x' if d['found'] else ' '}] data directory    "
          f"{d['path'] or '(not configured)'}  [{d['source']}]", flush=True)
    print(f"    [{'x' if cap['api_key'] else ' '}] Anthropic key     "
          f"{'present' if cap['api_key'] else 'absent (API chat unavailable; local chat unaffected)'}",
          flush=True)
    if CFG.CONFIG_FILE:
        print(f"    config file: {CFG.CONFIG_FILE}", flush=True)
    else:
        print(f"    config file: none found (searched SV_CONFIG, ./sv-assistant.conf, "
              f"~/.config/sv-assistant/config.ini)", flush=True)
    print(f"\n  http://127.0.0.1:{port}\n", flush=True)


def main():
    global TIERS, BANDS, TIER_ERROR
    ap = argparse.ArgumentParser(description="Local breakpoint-evidence front end")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--dataset", action="append", default=[], metavar="LABEL=PATH",
                    help="register a BAM under an explicit label (use this for local "
                         "patient data; the label is all the browser ever sees)")
    ap.add_argument("--candidates", action="append", default=[], metavar="LABEL=PATH")
    ap.add_argument("--no-autodiscover", action="store_true")
    ap.add_argument("--check", action="store_true",
                    help="verify the install: print the capability banner and exit. "
                         "Exit 0 if the MINIMAL tier works, 1 if it does not.")
    a = ap.parse_args()

    # The panel tool searches $IGV_PATH and its built-in paths; it never reads the
    # config file. Export a configured igv so the tool and the banner resolve the
    # same igv.sh: the banner used to report a config-file IGV the tool never saw.
    igv_cfg, _ = CFG.get("igv")
    if igv_cfg and not os.environ.get("IGV_PATH"):
        os.environ["IGV_PATH"] = igv_cfg

    ne, nb = assert_tool_contract()
    if not a.no_autodiscover:
        discover_public()
    for spec, target in ((a.dataset, DATASETS), (a.candidates, CANDIDATE_FILES)):
        for item in spec:
            lbl, _, path = item.partition("=")
            if not path:
                raise SystemExit(f"bad spec {item!r}; expected LABEL=PATH")
            target[lbl] = os.path.expanduser(path)

    try:
        TIERS = derive_tiers()
        BANDS = derive_bands()
    except TierDerivationError as e:
        TIERS, BANDS, TIER_ERROR = None, None, str(e)

    print(f"tool contract OK: {ne} evidence tools + {nb} bridge tools", flush=True)
    print(f"scoring tiers: {'derived from bam_tools source' if TIERS else 'NOT DERIVABLE — ' + str(TIER_ERROR)}", flush=True)
    print(f"datasets: {len(DATASETS)}   candidate files: {len(CANDIDATE_FILES)}", flush=True)
    cap = capability_report(verify_minimal())
    print_banner(cap, a.port)
    if a.check:
        for tier in ("MINIMAL", "FULL"):
            for name, cond_ok, detail in cap[tier]["conditions"]:
                print(f"  CHECK {tier:7s} [{'PASS' if cond_ok else 'FAIL'}] {name}"
                      + (f" — {detail}" if detail else ""), flush=True)
        ok = cap["MINIMAL"]["ok"]
        print("  CHECK: MINIMAL tier " + ("OK — every condition above was exercised, none assumed."
                                          if ok else "FAILED — see the FAIL lines above."),
              flush=True)
        if not DATASETS and not CANDIDATE_FILES:
            print("  CHECK: no datasets or candidate sets found. The tool will start, but the "
                  "dropdowns will be empty. Set data_dir, or register files explicitly "
                  "(--dataset LABEL=PATH / --candidates LABEL=PATH, or the config file).",
                  flush=True)
        raise SystemExit(0 if ok else 1)
    ThreadingHTTPServer(("127.0.0.1", a.port), Handler).serve_forever()


if __name__ == "__main__":
    main()
