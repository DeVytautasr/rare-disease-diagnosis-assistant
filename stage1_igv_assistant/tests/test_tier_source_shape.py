#!/usr/bin/env python3
"""
Regression test: the scoring tiers must derive from the scoring function's
source whatever the layout of its signature, and a source that genuinely
cannot be parsed must still report "not derivable" rather than crash.

Named for the condition that exposed it. score_tiers dedented the source with
inspect.cleandoc, a docstring tool that measures the margin from the SECOND
line on. The real signature spans several lines and closes with ") -> dict:"
at column 0, which pinned that margin to zero. Written on one line, the body's
own indentation became the margin and was stripped, parsing raised
IndentationError, and the interface reported the ceiling as not derivable for
a source that was perfectly valid Python.

Every case feeds a variant of the REAL source through derive_tiers() and
derive_bands() by replacing inspect.getsource for the duration of the case:

  CONDITION   the one-line variant has its whole signature on its first line, is
              valid Python, and is broken by cleandoc -- so the test cannot pass
              vacuously if the real source's shape ever changes
  ONE-LINE    derives tiers and bands identical to the real source's
  INDENTED    the whole function indented four spaces, as getsource returns a
              method, derives the same
  BROKEN      a source with an unmatched ")" raises TierDerivationError
              ("could not parse")
  NO SOURCE   getsource raising OSError raises TierDerivationError
  NO LADDER   a valid function without the tier chains raises TierDerivationError
              for both tiers and bands, not an unhandled exception

Run: python3 stage1_igv_assistant/tests/test_tier_source_shape.py
"""
import ast
import inspect
import os
import re
import sys
import textwrap
from unittest import mock

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from stage1_igv_assistant import score_tiers  # noqa: E402
from stage1_igv_assistant.score_tiers import TierDerivationError  # noqa: E402
from stage1_igv_assistant.tools import bam_tools  # noqa: E402

FAILURES = []


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}" + (f" -- {detail}" if detail else ""))
        FAILURES.append(label)


def derive_with(source=None, raises=None):
    """(tiers, bands) or (exception, exception) when getsource returns `source`."""
    def fake(obj):
        if raises is not None:
            raise raises
        return source
    out = []
    with mock.patch.object(score_tiers.inspect, "getsource", fake):
        for fn in (score_tiers.derive_tiers, score_tiers.derive_bands):
            try:
                out.append(fn())
            except Exception as e:          # the test classifies; it must not crash
                out.append(e)
    return tuple(out)


def one_line(src):
    m = re.search(r"def summarize_breakpoint_evidence\((.*?)\)\s*->\s*dict:", src, re.S)
    params = re.sub(r"\s*\n\s*", " ", m.group(1)).strip().rstrip(",")
    return src[:m.start()] + f"def summarize_breakpoint_evidence({params}) -> dict:" + src[m.end():]


def main():
    real_src = inspect.getsource(bam_tools.summarize_breakpoint_evidence)
    real = derive_with(real_src)
    print("\nREAL")
    check("the unmodified source derives tiers and bands",
          isinstance(real[0], dict) and isinstance(real[1], list), repr(real))

    print("\nCONDITION")
    ol = one_line(real_src)
    first = ol.splitlines()[0]
    check("the variant's first line holds the whole signature",
          first.startswith("def summarize_breakpoint_evidence(") and first.endswith("-> dict:"))
    try:
        ast.parse(ol)
        valid = True
    except SyntaxError:
        valid = False
    check("the variant is valid Python", valid)
    try:
        ast.parse(inspect.cleandoc(ol))
        broken_by_cleandoc = False
    except SyntaxError:                      # IndentationError is a SyntaxError
        broken_by_cleandoc = True
    check("cleandoc breaks the variant (the defect this test is named for)", broken_by_cleandoc)

    print("\nONE-LINE")
    got = derive_with(ol)
    check("tiers identical to the real source's", got[0] == real[0], repr(got[0])[:200])
    check("bands identical to the real source's", got[1] == real[1], repr(got[1])[:200])

    print("\nINDENTED")
    got = derive_with(textwrap.indent(real_src, "    "))
    check("tiers identical to the real source's", got[0] == real[0], repr(got[0])[:200])
    check("bands identical to the real source's", got[1] == real[1], repr(got[1])[:200])

    print("\nBROKEN")
    # An unmatched ")" as the body's first line. (Truncating the source is not a
    # reliable breakage: cut at a statement boundary, it still parses -- this
    # test's own condition check caught exactly that in its first version.)
    head, sep, body = real_src.partition(") -> dict:\n")
    cut = head + sep + "    )\n" + body
    try:
        ast.parse(textwrap.dedent(cut))
        cut_parses = True
    except SyntaxError:
        cut_parses = False
    check("condition: the damaged source really does not parse", bool(sep) and not cut_parses)
    got = derive_with(cut)
    for name, r in zip(("tiers", "bands"), got):
        check(f"{name}: TierDerivationError saying it could not parse",
              isinstance(r, TierDerivationError) and "could not parse" in str(r), repr(r)[:200])

    print("\nNO SOURCE")
    got = derive_with(raises=OSError("could not get source code"))
    for name, r in zip(("tiers", "bands"), got):
        check(f"{name}: TierDerivationError", isinstance(r, TierDerivationError), repr(r)[:200])

    print("\nNO LADDER")
    got = derive_with("def summarize_breakpoint_evidence(bam_path):\n    return {}\n")
    for name, r in zip(("tiers", "bands"), got):
        check(f"{name}: TierDerivationError, not a crash", isinstance(r, TierDerivationError),
              repr(r)[:200])

    print("\n" + "=" * 60)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("ALL TIER SOURCE-SHAPE TESTS PASSED")


if __name__ == "__main__":
    main()
