"""
save_blind_run.py

Converts one instruction-blinded Claude Code session (run via the Agent
tool, not this harness's own MCP loop) into the same RunLog JSON shape
ollama_harness.py / claude_harness.py produce, so score.py can score all
three arms identically.

This arm's tool-call trace is self-reported by the subagent (asked to
transcribe its own tool calls as a JSON block at the end of its report) --
there is no independent capture of the real MCP wire traffic the way the
other two harnesses have via mcp_client.py. Treat the trace as a
plausibility-checked transcript, not a ground-truth log: the self-reported
"chr1:115686862" values in the pilot run matched the independently-verified
figures already committed in results/EVIDENCE_PANEL_VALIDATION.md, which is
reassuring but not a per-run guarantee.

Usage: called once per completed Agent-tool session; not a CLI entry point.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from stage1_igv_assistant.benchmark.mcp_client import RunLog, ToolCallRecord

RUNS_DIR = Path(__file__).resolve().parent / "runs" / "blind_arm"


def save_blind_run(
    model: str,
    case_id: str,
    run_index: int,
    wall_clock_seconds: float,
    tool_calls_raw: list[dict],
    final_report: str,
    agent_id: str,
) -> Path:
    tool_calls = []
    for entry in tool_calls_raw:
        result = entry.get("result", {})
        if not (isinstance(result, dict) and "is_error" in result and "payload" in result):
            result = {"is_error": False, "payload": result}
        tool_calls.append(
            ToolCallRecord(
                turn=entry.get("turn", 0),
                tool_name=entry["tool_name"],
                arguments=entry.get("arguments", {}),
                result=result,
                wall_clock_seconds=0.0,
                malformed=False,
                malformed_reason=None,
            ).to_dict()
        )

    log = RunLog(
        model=model,
        backend="anthropic-instruction-blind",
        case_id=case_id,
        run_index=run_index,
        timestamp=datetime.now(timezone.utc).isoformat(),
        system_prompt=(
            "(not captured verbatim for this arm -- the subagent was instructed "
            "to use ONLY the igv-breakpoint-assistant MCP tools, not read source "
            "files, and not rely on prior knowledge; see cases.py's prompt template)"
        ),
        messages=[{"role": "agent_id", "content": agent_id}],
        tool_calls=tool_calls,
        final_report=final_report,
        hit_max_turns=False,
        wall_clock_seconds=round(wall_clock_seconds, 3),
        error=None,
    )
    RUNS_DIR.mkdir(parents=True, exist_ok=True)
    safe_model = model.replace(":", "-").replace("/", "-")
    path = RUNS_DIR / f"{safe_model}-blind__{case_id}__run{run_index}.json"
    log.save(path)
    return path
