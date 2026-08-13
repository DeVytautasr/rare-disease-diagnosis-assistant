"""
claude_harness.py

Runs a Claude model as the tool-calling agent against the IGV breakpoint
assistant MCP server, via the Anthropic Messages API directly rather than
through the Claude Code CLI -- this is the "same code" Claude baseline:
the MCP session handling, tool schema conversion, and run-log format are
shared with ollama_harness.py via mcp_client.py, so the only thing that
differs between the local-model runs and this baseline is which chat API
is called. Contrast with the instruction-blinded Claude Code sessions
(results/BENCHMARK_CLAUDE_BASELINE.md's other arm), which reuse this same
MCP server but through the actual Claude Code client rather than this
harness.

Requires ANTHROPIC_API_KEY, read from the environment or a repo-root
.env file (already .gitignore'd) via python-dotenv. Calls are billed to
that key's Anthropic account -- separate from any Claude Code
subscription.

Usage:
    python -m stage1_igv_assistant.benchmark.claude_harness \\
        --model claude-sonnet-5 --case POSITIVE --runs 3 \\
        --output stage1_igv_assistant/benchmark/runs
"""
from __future__ import annotations

import argparse
import asyncio
import os
import time
from datetime import datetime, timezone
from pathlib import Path

from dotenv import load_dotenv

from stage1_igv_assistant.benchmark.cases import CASES
from stage1_igv_assistant.benchmark.mcp_client import (
    MAX_TURNS,
    REPO_ROOT,
    RunLog,
    ToolCallRecord,
    call_mcp_tool,
    mcp_schema_to_anthropic_tool,
    mcp_session,
    read_image_manifest,
    redact_image_payload,
)

load_dotenv(REPO_ROOT / ".env")

MAX_TOKENS = 8192

# claude-sonnet-5 runs adaptive thinking by default when `thinking` is
# omitted (unlike Opus 4.8/4.7, where omitting it meant no thinking) --
# thinking tokens bill as output at the same rate as the response text.
# This task is a bounded tool-call sequence + evidence-cited report, not
# open-ended reasoning, so thinking is turned off to keep the API-metered
# arm's cost predictable; see results/BENCHMARK_CLAUDE_BASELINE.md for the
# rationale and the actual measured cost this produced.
THINKING: dict = {"type": "disabled"}

# Sonnet 5 pricing as of 2026-08 (introductory, through 2026-08-31):
# $2.00 / MTok input, $10.00 / MTok output. Update if pricing changes.
PRICE_PER_MTOK_INPUT_USD = 2.00
PRICE_PER_MTOK_OUTPUT_USD = 10.00


async def run_once(client, model: str, case_id: str, run_index: int, max_turns: int = MAX_TURNS) -> RunLog:
    case = CASES[case_id]
    start = time.monotonic()
    error: str | None = None

    async with mcp_session() as (session, init_result, image_session_dir):
        system_prompt = init_result.instructions or ""
        tools_result = await session.list_tools()
        anthropic_tools = [mcp_schema_to_anthropic_tool(t) for t in tools_result.tools]

        messages: list[dict] = [{"role": "user", "content": case.prompt}]
        tool_call_log: list[ToolCallRecord] = []
        hit_max_turns = True
        final_report: str | None = None
        input_tokens_total = 0
        output_tokens_total = 0
        image_handles: dict = {}

        for turn in range(1, max_turns + 1):
            try:
                response = client.messages.create(
                    model=model,
                    max_tokens=MAX_TOKENS,
                    system=system_prompt,
                    messages=messages,
                    tools=anthropic_tools,
                    thinking=THINKING,
                )
            except Exception as exc:  # noqa: BLE001 -- record and stop, don't crash the sweep
                error = f"messages.create failed on turn {turn}: {exc}"
                break

            input_tokens_total += response.usage.input_tokens
            output_tokens_total += response.usage.output_tokens

            assistant_blocks = [block.model_dump() for block in response.content]
            messages.append({"role": "assistant", "content": assistant_blocks})

            tool_use_blocks = [b for b in response.content if b.type == "tool_use"]
            if response.stop_reason != "tool_use" or not tool_use_blocks:
                text_parts = [b.text for b in response.content if b.type == "text"]
                final_report = "\n".join(text_parts)
                hit_max_turns = False
                break

            tool_result_blocks = []
            for block in tool_use_blocks:
                name = block.name
                args = dict(block.input)
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
                # The server already returns handles rather than paths (FIX
                # C), so this is now a no-op safety net rather than the
                # primary control -- kept so a future tool that forgets to
                # redact still can't leak a path into the transcript.
                visible_payload, new_handles = redact_image_payload(result["payload"])
                image_handles.update(new_handles)
                tool_result_blocks.append(
                    {
                        "type": "tool_result",
                        "tool_use_id": block.id,
                        "content": str(visible_payload),
                        "is_error": bool(result["is_error"]),
                    }
                )
            messages.append({"role": "user", "content": tool_result_blocks})

        # Server-assigned handles: the only route from an image_ref the
        # model cited back to a real file a human can open.
        image_handles.update(read_image_manifest(image_session_dir))

        elapsed = time.monotonic() - start
        estimated_cost_usd = (
            input_tokens_total * PRICE_PER_MTOK_INPUT_USD
            + output_tokens_total * PRICE_PER_MTOK_OUTPUT_USD
        ) / 1_000_000
        return RunLog(
            model=model,
            backend="anthropic",
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
            usage={
                "input_tokens": input_tokens_total,
                "output_tokens": output_tokens_total,
                "estimated_cost_usd": round(estimated_cost_usd, 4),
            },
            image_handles=image_handles or None,
        )


async def run_sweep(model: str, case_ids: list[str], runs: int, output_dir: Path) -> list[Path]:
    import anthropic

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        raise SystemExit(
            "ANTHROPIC_API_KEY not set. Add it to a repo-root .env file "
            "(ANTHROPIC_API_KEY=sk-ant-...) or export it before running."
        )
    client = anthropic.Anthropic(api_key=api_key)

    saved: list[Path] = []
    running_cost_usd = 0.0
    for case_id in case_ids:
        for run_index in range(1, runs + 1):
            print(f"[{model}] {case_id} run {run_index}/{runs} ...", flush=True)
            log = await run_once(client, model, case_id, run_index)
            safe_model = model.replace(":", "-").replace("/", "-")
            path = output_dir / f"{safe_model}__{case_id}__run{run_index}.json"
            log.save(path)
            saved.append(path)
            status = "MAX_TURNS" if log.hit_max_turns else "done"
            n_tools = len(log.tool_calls)
            run_cost = (log.usage or {}).get("estimated_cost_usd", 0.0)
            running_cost_usd += run_cost
            print(
                f"    -> {status} in {log.wall_clock_seconds:.1f}s, "
                f"{n_tools} tool calls, ${run_cost:.4f} this run "
                f"(${running_cost_usd:.4f} running total), saved to {path}",
                flush=True,
            )
    return saved


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", default="claude-sonnet-5", help="Anthropic model id")
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
