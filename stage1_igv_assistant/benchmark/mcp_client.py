"""
mcp_client.py

Shared MCP mechanics for both benchmark harnesses (ollama_harness.py,
claude_harness.py) so the session handling, tool schema conversion, and
run-log format are identical regardless of which model backend is under
test -- only the model-specific chat call differs between the two.
"""
from __future__ import annotations

import copy
import dataclasses
import hashlib
import json
import os
import sys
import tempfile
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any, Optional

from mcp import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client
from mcp.types import Tool

REPO_ROOT = Path(__file__).resolve().parents[2]
MAX_TURNS = 20


def read_image_manifest(image_session_dir: str) -> dict:
    """
    {image_ref: real path} for images the server generated this run.

    The server assigns output paths itself and returns only handles (see
    bam_tools.image_session_dir / to_handle_result), so this manifest is the
    only way back to the files. The model never sees it.
    """
    manifest = Path(image_session_dir) / "manifest.json"
    if not manifest.exists():
        return {}
    try:
        return json.loads(manifest.read_text())
    except (OSError, ValueError):
        return {}


def build_server_env(image_session_dir: Optional[str] = None) -> dict:
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
    if image_session_dir:
        # Tells the server where to write screenshots, so this harness knows
        # where to read the handle->path manifest back from without the
        # model ever seeing a path.
        env["IGV_IMAGE_SESSION_DIR"] = image_session_dir
    return env


@asynccontextmanager
async def mcp_session(image_session_dir: Optional[str] = None):
    """
    Launch stage1_igv_assistant/server.py over stdio and yield
    (session, initialize_result, image_session_dir).

    image_session_dir is where the server writes screenshots; a fresh
    temporary directory is created when none is given, so each run's
    handle manifest is isolated.
    """
    if image_session_dir is None:
        image_session_dir = tempfile.mkdtemp(prefix="igv_images_")
    os.makedirs(image_session_dir, exist_ok=True)
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "stage1_igv_assistant.server"],
        cwd=str(REPO_ROOT),
        env=build_server_env(image_session_dir),
    )
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            init_result = await session.initialize()
            yield session, init_result, image_session_dir


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


# ── Image-path redaction (harness-level constraint) ────────────────────────
#
# The screenshot tools return a filesystem path. The MCP tool-call mechanism
# carries only text, so the model never receives the image's pixels -- yet in
# testing both claude-sonnet-5 and qwen2.5:7b wrote confident descriptions of
# image contents they had not been shown (results/LLM_SESSION_4_VISUAL_*.md).
#
# Adding an advisory field and a system-prompt rule fixed this for
# claude-sonnet-5 and did NOT fix it for qwen2.5:7b, which kept writing
# declarative claims ("Read Depth Profile Image: Demonstrate a significant
# depth drop...") with the advisory text sitting in its context. An
# instruction the model can simply not follow is not a constraint.
#
# So the path is removed before the tool result enters the message history
# and replaced with an opaque handle. The model cannot describe pixels it has
# no path to, no image of, and nothing path-shaped to pattern-match on. It
# retains everything it can legitimately reason about: the region, layer,
# coloring mode, dimensions, and whether generation succeeded.
#
# IMPORTANT SCOPE NOTE: this is a HARNESS-level constraint, not a model-level
# one. It works because this harness feeds the model text only. A
# vision-capable client that genuinely passes image content would need a
# different approach entirely -- there, the model can see the image, and the
# correct requirement is that its description match the pixels, not that it
# abstain. Nothing here generalises to that case.
_IMAGE_PATH_KEY = "screenshot_path"

# Dropped from the model-visible copy only. batch_script embeds the output
# path verbatim in its `snapshot <path>` line (so redacting screenshot_path
# alone would leak it straight back); igv_stdout/igv_stderr are IGV process
# logs that also carry filesystem paths and have no analytical use. All three
# are preserved in full in the run log.
_MODEL_HIDDEN_KEYS = ("batch_script", "igv_stdout", "igv_stderr")

_REDACTED_IMAGE_NOTE = (
    "An image was generated but has NOT been provided to you: you have "
    "received an opaque reference (image_ref) only -- no file path, and no "
    "pixel data. You cannot see this image and must not describe what it "
    "shows, contains, or looks like. Report that it was generated, which "
    "layer and region it covers, and what a human reviewer should check. "
    "The image_ref maps to a real file via the run log's image_handles."
)


def _image_handle(path: str) -> str:
    """Stable, non-path-shaped reference for one image file."""
    return "IMG_" + hashlib.sha256(path.encode("utf-8")).hexdigest()[:4]


def _png_dimensions(path: str) -> Optional[str]:
    """
    "WxH" read from the PNG IHDR chunk, or None. Dimensions are legitimate
    metadata (they describe the file, not its depicted contents), so they
    stay visible to the model.
    """
    try:
        with open(path, "rb") as fh:
            header = fh.read(24)
        if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
            return None
        width = int.from_bytes(header[16:20], "big")
        height = int.from_bytes(header[20:24], "big")
        return f"{width}x{height}"
    except OSError:
        return None


def redact_image_payload(payload: Any) -> tuple[Any, dict]:
    """
    Return (model_visible_payload, {handle: real_path}).

    Walks the payload and rewrites every dict carrying a screenshot_path --
    which covers igv_screenshot's top-level result and each of
    evidence_panel's per-layer panel entries, without either tool needing to
    be special-cased. Non-image payloads pass through untouched.

    The caller keeps the ORIGINAL payload in the run log and sends only the
    returned copy to the model, so images stay retrievable by a human while
    being unreachable by the agent.
    """
    handles: dict[str, str] = {}

    def _walk(node: Any) -> Any:
        if isinstance(node, dict):
            if isinstance(node.get(_IMAGE_PATH_KEY), str):
                path = node[_IMAGE_PATH_KEY]
                handle = _image_handle(path)
                handles[handle] = path
                out = {
                    k: _walk(v)
                    for k, v in node.items()
                    if k != _IMAGE_PATH_KEY and k not in _MODEL_HIDDEN_KEYS
                }
                out["image_ref"] = handle
                dims = _png_dimensions(path)
                if dims:
                    out["image_dimensions"] = dims
                out["image_content_available_to_caller"] = False
                out["note"] = _REDACTED_IMAGE_NOTE
                return out
            return {k: _walk(v) for k, v in node.items()}
        if isinstance(node, list):
            return [_walk(v) for v in node]
        return node

    return _walk(copy.deepcopy(payload)), handles


@dataclasses.dataclass
class ToolCallRecord:
    turn: int
    tool_name: str
    arguments: dict
    result: dict
    wall_clock_seconds: float
    malformed: bool = False
    malformed_reason: Optional[str] = None
    via_text_fallback: bool = False

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
    usage: Optional[dict] = None  # {"input_tokens", "output_tokens", "estimated_cost_usd"} -- Anthropic backend only
    # {image_ref: real filesystem path} for every screenshot generated this
    # run. The model saw only the refs (see redact_image_payload); this is
    # how a human gets back to the actual files.
    image_handles: Optional[dict] = None
    # Which tool/scorer fixes were in place when this run executed. Provenance
    # lives here rather than only in the directory name, so a log stays
    # self-describing if it is ever moved -- see benchmark/runs/README.md.
    # Backfilled for existing logs; set it explicitly on new sweeps.
    stage: Optional[str] = None

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2, default=str))


def load_run_log(path: Path) -> dict:
    return json.loads(path.read_text())
