"""
mcp_client.py

Shared MCP mechanics for both benchmark harnesses (ollama_harness.py,
claude_harness.py) so the session handling, tool schema conversion, and
run-log format are identical regardless of which model backend is under
test -- only the model-specific chat call differs between the two.
"""
from __future__ import annotations

import dataclasses
import json
import os
import sys
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any, Optional

from mcp import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client
from mcp.types import Tool

REPO_ROOT = Path(__file__).resolve().parents[2]
MAX_TURNS = 20


def build_server_env() -> dict:
    """
    Env for the MCP server subprocess. Must include the active Python's
    bin/ directory in PATH or java (needed by igv_screenshot/
    evidence_panel) is unreachable through MCP even though it works from
    an activated shell -- the same PATH issue documented in
    stage1_igv_assistant/README.md's Known limitations and hit again in
    results/EVIDENCE_PANEL_VALIDATION.md's "Note on this session's
    environment".
    """
    env = dict(os.environ)
    conda_bin = str(Path(sys.executable).parent)
    env["PATH"] = conda_bin + os.pathsep + env.get("PATH", "")
    return env


@asynccontextmanager
async def mcp_session():
    """Launch stage1_igv_assistant/server.py over stdio and yield (session, initialize_result)."""
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "stage1_igv_assistant.server"],
        cwd=str(REPO_ROOT),
        env=build_server_env(),
    )
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            init_result = await session.initialize()
            yield session, init_result


def mcp_schema_to_ollama_tool(tool: Tool) -> dict:
    """OpenAI-style function-tool dict, the format ollama-python expects."""
    schema = tool.inputSchema or {"type": "object", "properties": {}}
    return {
        "type": "function",
        "function": {
            "name": tool.name,
            "description": tool.description or "",
            "parameters": {
                "type": schema.get("type", "object"),
                "properties": schema.get("properties", {}),
                "required": schema.get("required", []),
            },
        },
    }


def mcp_schema_to_anthropic_tool(tool: Tool) -> dict:
    """Anthropic Messages-API tool dict."""
    schema = tool.inputSchema or {"type": "object", "properties": {}}
    return {
        "name": tool.name,
        "description": tool.description or "",
        "input_schema": schema,
    }


async def call_mcp_tool(session: ClientSession, name: str, arguments: dict) -> dict:
    """
    Execute one tool call against the live MCP session and return a
    JSON-serialisable {"is_error": bool, "payload": ...} dict, using
    structuredContent when the server provided it (it always does here,
    since every stage1 tool returns a plain dict) and falling back to
    parsing the text content otherwise.
    """
    result = await session.call_tool(name, arguments)
    if result.structuredContent is not None:
        payload = result.structuredContent
    else:
        text = "".join(getattr(block, "text", "") for block in result.content)
        try:
            payload = json.loads(text)
        except json.JSONDecodeError:
            payload = {"raw_text": text}
    return {"is_error": bool(result.isError), "payload": payload}


@dataclasses.dataclass
class ToolCallRecord:
    turn: int
    tool_name: str
    arguments: dict
    result: dict
    wall_clock_seconds: float
    malformed: bool = False
    malformed_reason: Optional[str] = None

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)


@dataclasses.dataclass
class RunLog:
    model: str
    backend: str            # "ollama" | "anthropic"
    case_id: str
    run_index: int
    timestamp: str
    system_prompt: str
    messages: list           # full chat history, backend-native message shape
    tool_calls: list         # list[dict] (ToolCallRecord.to_dict()), in call order
    final_report: Optional[str]
    hit_max_turns: bool
    wall_clock_seconds: float
    error: Optional[str] = None

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2, default=str))


def load_run_log(path: Path) -> dict:
    return json.loads(path.read_text())
