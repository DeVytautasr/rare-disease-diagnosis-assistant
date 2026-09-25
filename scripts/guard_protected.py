#!/usr/bin/env python3
"""AST-diff guard: the evidence tools, their thresholds and the scoring must be
unchanged, and any script about to run must be exactly what is committed.

    guard_protected.py [--base COMMIT] [--also PATH ...]
    guard_protected.py --self-test

PROTECTED files are compared three ways: the working tree against HEAD (bytes
and AST), and HEAD against --base (bytes; default 9ae73bc, the last commit on
the public remote before the 2026-09 recovery work). Any difference fails --
a comment-only edit is still an edit to a file nobody may edit.

--also PATH adds a file that must be identical to its committed version (bytes
and, for Python, AST) before it is run: what runs is then what is recorded.

Exit 0 all identical, 1 a difference, 2 could not check (e.g. untracked file).
"""
import argparse
import ast
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROTECTED = [
    "stage1_igv_assistant/tools/bam_tools.py",      # the eleven evidence tools, thresholds, scoring
    "stage1_igv_assistant/server.py",               # their MCP wrappers
    "stage1_igv_assistant/tools/vcf_tools.py",      # candidate-set parsing, dedup and recurrence
    "stage1_igv_assistant/candidate_server.py",     # the four candidate-set tools
]
DEFAULT_BASE = "9ae73bc"


def git_bytes(rev, path):
    p = subprocess.run(["git", "-C", REPO, "show", f"{rev}:{path}"], capture_output=True)
    return p.stdout if p.returncode == 0 else None


def ast_dump(src, path):
    try:
        return ast.dump(ast.parse(src, filename=path), include_attributes=False)
    except SyntaxError as e:
        return f"<SyntaxError {e.msg}>"


def compare(path, work, head, base=None):
    """[problems] for one file."""
    probs = []
    if head is None:
        return [f"{path}: not in HEAD (untracked or new) -- cannot vouch for it"]
    if work != head:
        what = "bytes differ"
        if path.endswith(".py"):
            same_ast = ast_dump(work, path) == ast_dump(head, path)
            what += ", AST identical (formatting/comments)" if same_ast else ", AST DIFFERS"
        probs.append(f"{path}: working tree vs HEAD: {what}")
    if base is not None and head != base:
        probs.append(f"{path}: HEAD vs base: bytes differ")
    return probs


def run(base, also):
    probs, checked = [], 0
    for path in PROTECTED:
        with open(os.path.join(REPO, path), "rb") as f:
            work = f.read()
        b = git_bytes(base, path)
        if b is None:
            print(f"cannot read {path} at {base}")
            return 2
        probs += compare(path, work, git_bytes("HEAD", path), b)
        checked += 1
    for path in also:
        rel = os.path.relpath(os.path.abspath(path), REPO)
        with open(os.path.join(REPO, rel), "rb") as f:
            work = f.read()
        head = git_bytes("HEAD", rel)
        if head is None:
            print(f"{rel}: not committed -- commit it before running it")
            return 2
        probs += compare(rel, work, head)
        checked += 1
    for p in probs:
        print(f"  DIFFERS  {p}")
    print(f"guard: {checked} file(s) checked against HEAD"
          f" ({len(PROTECTED)} protected, also against {base}): "
          + ("ALL IDENTICAL" if not probs else f"{len(probs)} DIFFERENCE(S)"))
    return 1 if probs else 0


def self_test():
    """Positive and negative controls on in-memory copies; nothing is written."""
    path = PROTECTED[0]
    head = git_bytes("HEAD", path)
    src = head.decode()
    ok = True

    def expect(name, probs, want_fire):
        nonlocal ok
        fired = bool(probs)
        good = fired == want_fire
        ok &= good
        print(f"  {'PASS' if good else 'FAIL'}  {name}: {'fired' if fired else 'silent'}"
              + (f" ({probs[0].split(': ', 1)[1]})" if probs else ""))

    expect("unchanged file", compare(path, head, head, head), False)
    # The assignment itself, at the start of a line: the same text also occurs in
    # a comment, and mutating a comment is (correctly) not an AST change.
    marker = re.compile(r"^LOW_MAPQ_QUALITY_GATE = 0\.4$", re.M)
    if len(marker.findall(src)) != 1:
        print(f"  FAIL  self-test needs exactly one assignment LOW_MAPQ_QUALITY_GATE = 0.4 in {path}")
        return 1
    changed = marker.sub("LOW_MAPQ_QUALITY_GATE = 0.5", src).encode()
    expect("a threshold changed 0.4 -> 0.5", compare(path, changed, head), True)
    expect("... and it is reported as an AST difference",
           [p for p in compare(path, changed, head) if "AST DIFFERS" in p], True)
    comment = (src + "\n# a comment\n").encode()
    expect("a comment appended (still an edit)", compare(path, comment, head), True)
    expect("HEAD differs from base", compare(path, head, head, changed), True)
    expect("untracked file", compare("new.py", b"x = 1\n", None), True)
    print("self-test " + ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=DEFAULT_BASE)
    ap.add_argument("--also", nargs="*", default=[])
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    sys.exit(self_test() if a.self_test else run(a.base, a.also))
