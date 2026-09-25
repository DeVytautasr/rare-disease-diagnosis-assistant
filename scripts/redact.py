#!/usr/bin/env python3
"""Mask patient identifiers in anything that is about to be displayed.

    some_command 2>&1 | scripts/redact.sh
    scripts/redact.sh -- some_command args ...   (stdout+stderr filtered, exit code kept)

Reads the private term list ~/patient_data/.redact_terms ("term<TAB>label" per
line; --terms FILE overrides). Every term is matched case-insensitively, a
separator in a term (. _ - or space) matches any other separator, longer terms
win over shorter ones, and each match is replaced by [label], e.g.
[SAMPLE_A:SM]. Standard library only; streams line by line, so it can sit
behind tail -f.

It FAILS CLOSED. If the term list is missing, unreadable, empty, or readable by
anyone but its owner, it writes nothing from its input and exits 2. A filter
that passes text through unchanged when it cannot load its terms looks exactly
like a working one, which is worse than having no filter.
"""
import argparse
import os
import re
import stat
import sys

DEFAULT_TERMS = os.path.expanduser("~/patient_data/.redact_terms")
SEPARATORS = "-_. "


class TermListError(Exception):
    pass


def load_terms(path):
    st = os.stat(path)
    if stat.S_IMODE(st.st_mode) & 0o077:
        raise TermListError(f"the term list is readable beyond its owner "
                              f"(mode {oct(stat.S_IMODE(st.st_mode))}); expected 0o600")
    terms = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            term, _, label = line.partition("\t")
            if term.strip():
                terms.append((term, label.strip() or "REDACTED"))
    if not terms:
        raise TermListError("the term list contains no terms")
    return terms


def compile_terms(terms):
    terms = sorted(terms, key=lambda t: -len(t[0]))
    alts = []
    for term, _ in terms:
        alts.append("(" + "".join("[-_.\\s]" if ch in SEPARATORS else re.escape(ch)
                                  for ch in term) + ")")
    pat = re.compile("|".join(alts), re.IGNORECASE)
    labels = [label for _, label in terms]
    return pat, (lambda m: f"[{labels[m.lastindex - 1]}]")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--terms", default=DEFAULT_TERMS)
    a = ap.parse_args()
    try:
        pat, repl = compile_terms(load_terms(a.terms))
    except Exception as e:
        # The reason never quotes the path: the default path is generic, but a
        # --terms path could be anything.
        why = os.strerror(e.errno) if isinstance(e, OSError) and e.errno else str(e)
        print(f"redact: REFUSING to pass output through -- cannot load the term list "
              f"({type(e).__name__}: {why})", file=sys.stderr)
        sys.exit(2)
    out = sys.stdout.buffer
    for raw in sys.stdin.buffer:
        text = raw.decode("utf-8", "surrogateescape")
        out.write(pat.sub(repl, text).encode("utf-8", "surrogateescape"))
        out.flush()


if __name__ == "__main__":
    main()
