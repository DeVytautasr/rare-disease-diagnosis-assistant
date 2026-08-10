"""
test_server.py
Confirms the FastMCP server actually starts and exposes the tools it
claims to. Previously nothing exercised server.py at all — README's tool
count and the @mcp.tool() decorator count could drift from what a real
MCP client would see, with no test to catch it.
"""

import asyncio
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from stage1_igv_assistant.server import mcp

EXPECTED_TOOL_NAMES = {
    "applicable_layers",
    "bam_stats_at_locus",
    "breakpoint_evidence_summary",
    "discordant_pairs",
    "evidence_panel",
    "gene_at_locus",
    "igv_screenshot",
    "read_depth_profile",
    "reciprocal_breakpoint",
    "soft_clipped_reads",
    "split_reads",
}


def run_tests():
    print("=" * 60)
    print("SERVER STARTUP TEST SUITE")
    print("=" * 60)

    print("TEST 1: MCP server starts and exposes the expected tools")

    async def _list_tools():
        return await mcp.list_tools()

    tools = asyncio.run(_list_tools())
    tool_names = {t.name for t in tools}

    print(f"  Tools exposed ({len(tool_names)}): {sorted(tool_names)}")

    assert len(tools) == len(EXPECTED_TOOL_NAMES), (
        f"Expected exactly {len(EXPECTED_TOOL_NAMES)} tools, got {len(tools)}: "
        f"{sorted(tool_names)}"
    )
    assert tool_names == EXPECTED_TOOL_NAMES, (
        f"Tool set mismatch.\n  Missing: {EXPECTED_TOOL_NAMES - tool_names}\n"
        f"  Unexpected: {tool_names - EXPECTED_TOOL_NAMES}"
    )
    # Every tool must have a non-empty docstring — this server's whole design
    # depends on the LLM reading tool descriptions rather than guessing.
    undocumented = [t.name for t in tools if not (t.description or "").strip()]
    assert not undocumented, f"Tools with no docstring/description: {undocumented}"

    print("  PASSED ✓\n")

    print("=" * 60)
    print("ALL TESTS PASSED")
    print("=" * 60)


if __name__ == "__main__":
    run_tests()
