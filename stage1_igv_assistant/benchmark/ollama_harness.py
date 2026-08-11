"""
ollama_harness.py

Runs a local Ollama model as the tool-calling agent against the IGV
breakpoint assistant MCP server (stage1_igv_assistant/server.py), over
the same stdio transport a real MCP client uses. The MCP session
handling, tool schema conversion, and run-log format are shared with
claude_harness.py via mcp_client.py -- only the model-specific chat call
differs between the two backends.

Usage:
    python -m stage1_igv_assistant.benchmark.ollama_harness \\
        --model qwen2.5:7b --case POSITIVE --runs 3 \\
        --output stage1_igv_assistant/benchmark/runs
"""
from __future__ import annotations

import argparse
import asyncio
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import ollama

from stage1_igv_assistant.benchmark.cases import CASES
from stage1_igv_assistant.benchmark.mcp_client import (
    MAX_TURNS,
    RunLog,
    ToolCallRecord,
    call_mcp_tool,
    mcp_schema_to_ollama_tool,
    mcp_session,
)

# Default Ollama context (4096) is too small for this workload -- 11 tool
# schemas plus growing conversation history over up to MAX_TURNS turns
# overflows it and silently truncates. Verified empirically on this 8GB-VRAM
# GPU: qwen2.5:7b at num_ctx=16384 still loads 100% GPU-resident (~7.3GB used,
# ~0.9GB headroom) -- no CPU offload, no slowdown.
NUM_CTX = 16384


async def run_once(model: str, case_id: str, run_index: int) -> RunLog:
    case = CASES[case_id]
    start = time.monotonic()
    error: str | None = None

    async with mcp_session() as (session, init_result):
        system_prompt = init_result.instructions or ""
        tools_result = await session.list_tools()
        ollama_tools = [mcp_schema_to_ollama_tool(t) for t in tools_result.tools]

        messages: list[dict] = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": case.prompt},
        ]
        tool_call_log: list[ToolCallRecord] = []
        hit_max_turns = True
        final_report: str | None = None

        for turn in range(1, MAX_TURNS + 1):
            try:
                response = ollama.chat(
                    model=model,
                    messages=messages,
                    tools=ollama_tools,
                    options={"num_ctx": NUM_CTX},
                )
            except Exception as exc:  # noqa: BLE001 -- record and stop, don't crash the sweep
                error = f"ollama.chat failed on turn {turn}: {exc}"
                break

            msg = response.message
            assistant_entry: dict = {"role": "assistant", "content": msg.content or ""}
            if msg.tool_calls:
                assistant_entry["tool_calls"] = [
                    {
                        "function": {
                            "name": tc.function.name,
                            "arguments": dict(tc.function.arguments),
                        }
                    }
                    for tc in msg.tool_calls
                ]
            messages.append(assistant_entry)

            if not msg.tool_calls:
                final_report = msg.content or ""
                hit_max_turns = False
                break

            for tc in msg.tool_calls:
                name = tc.function.name
                args = dict(tc.function.arguments)
                call_start = time.monotonic()
                malformed = False
                malformed_reason = None
                try:
                    result = await call_mcp_tool(session, name, args)
                except Exception as exc:  # noqa: BLE001 -- a malformed call is a result, not a crash
                    malformed = True
                    malformed_reason = str(exc)
                    result = {"is_error": True, "payload": {"error": str(exc)}}
                call_elapsed = time.monotonic() - call_start
                tool_call_log.append(
                    ToolCallRecord(
                        turn=turn,
                        tool_name=name,
                        arguments=args,
                        result=result,
                        wall_clock_seconds=round(call_elapsed, 3),
                        malformed=malformed,
                        malformed_reason=malformed_reason,
                    )
                )
                messages.append(
                    {
                        "role": "tool",
                        "tool_name": name,
                        "content": json.dumps(result["payload"], default=str),
                    }
                )

        elapsed = time.monotonic() - start
        return RunLog(
            model=model,
            backend="ollama",
            case_id=case_id,
            run_index=run_index,
            timestamp=datetime.now(timezone.utc).isoformat(),
            system_prompt=system_prompt,
            messages=messages,
            tool_calls=[t.to_dict() for t in tool_call_log],
            final_report=final_report,
            hit_max_turns=hit_max_turns,
            wall_clock_seconds=round(elapsed, 3),
            error=error,
        )


async def run_sweep(model: str, case_ids: list[str], runs: int, output_dir: Path) -> list[Path]:
    saved: list[Path] = []
    for case_id in case_ids:
        for run_index in range(1, runs + 1):
            print(f"[{model}] {case_id} run {run_index}/{runs} ...", flush=True)
            log = await run_once(model, case_id, run_index)
            safe_model = model.replace(":", "-").replace("/", "-")
            path = output_dir / f"{safe_model}__{case_id}__run{run_index}.json"
            log.save(path)
            saved.append(path)
            status = "MAX_TURNS" if log.hit_max_turns else "done"
            n_tools = len(log.tool_calls)
            print(
                f"    -> {status} in {log.wall_clock_seconds:.1f}s, "
                f"{n_tools} tool calls, saved to {path}",
                flush=True,
            )
    return saved


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, help="Ollama model tag, e.g. qwen2.5:7b")
    parser.add_argument(
        "--case",
        default="ALL",
        help="Case id (POSITIVE, NEGATIVE, ADVERSARIAL) or ALL (default)",
    )
    parser.add_argument("--runs", type=int, default=3, help="Runs per case (default 3)")
    parser.add_argument(
        "--output",
        default="stage1_igv_assistant/benchmark/runs",
        help="Directory to write per-run JSON logs to",
    )
    args = parser.parse_args()

    case_ids = list(CASES.keys()) if args.case == "ALL" else [args.case]
    for cid in case_ids:
        if cid not in CASES:
            raise SystemExit(f"Unknown case '{cid}'. Known cases: {list(CASES.keys())}")

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    asyncio.run(run_sweep(args.model, case_ids, args.runs, output_dir))


if __name__ == "__main__":
    main()
