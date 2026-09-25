#!/usr/bin/env python3
"""Controls for scripts/redact.py and scripts/redact.sh. Run before first use.

  python3 scripts/test_redact.py            synthetic controls, then the real
                                            term list if one exists
  python3 scripts/test_redact.py --synthetic-only

POSITIVE  each term is planted into many contexts (a path, a JSON value, a
          tab-separated field, upper and lower case, its separators swapped,
          twice on one line, glued to punctuation) and must not survive the
          filter in any case or separator form.
NEGATIVE  ordinary output -- flagstat, idxstats, a VCF header, delly-style
          progress lines, JSON, the placeholders themselves -- must pass
          byte-identical. With the real list this is also the over-redaction
          check: a real term generic enough to match ordinary text fails here.
FAIL-CLOSED  a missing, empty, or group-readable term list must yield no output
          at all and a non-zero exit; and redact.sh -- CMD must keep CMD's exit
          code while filtering both of its streams.

A real term is never printed: failures name the term's label and the context
number only. Exit 0 all pass, 1 any failure.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REDACT_PY = os.path.join(HERE, "redact.py")
REDACT_SH = os.path.join(HERE, "redact.sh")
REAL = os.path.expanduser("~/patient_data/.redact_terms")

FAILS = []


def check(name, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f" -- {detail}" if detail and not ok else ""))
    if not ok:
        FAILS.append(name)


def run(stdin, terms, argv=None):
    cmd = argv or [sys.executable, REDACT_PY, "--terms", terms]
    p = subprocess.run(cmd, input=stdin, capture_output=True)
    return p.returncode, p.stdout, p.stderr


def contexts(term):
    swapped = re.sub(r"[-_. ]", lambda m: {"-": "_", "_": ".", ".": "-", " ": "_"}[m.group(0)], term)
    return [
        term,
        f"/data/incoming/{term}.bam",
        f"/data/incoming/{term}.bam.bai",
        f'{{"sample": "{term}", "n": 3}}',
        f"@RG\tID:x\tSM:{term}\tPL:ILLUMINA",
        term.upper(),
        term.lower(),
        swapped,
        f"{term} and again {term}",
        f"({term}),{term};",
        f"prefix{term}suffix",
        f"[{term}]",
    ]


def residual(term, text):
    """Does the term survive in any case or separator form?"""
    body = "".join("[-_.\\s]" if ch in "-_. " else re.escape(ch) for ch in term)
    return re.search(body, text, re.I) is not None


ORDINARY = (
    "38959428903 + 0 in total (QC-passed reads + QC-failed reads)\n"
    "631618015 + 0 primary mapped (97.83% : N/A)\n"
    "chr20\t64444167\t17681348\t0\n*\t0\t0\t12345\n"
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
    "[2026-Sep-25 12:00:00] Paired-end and split-read scanning\n"
    "[2026-Sep-25 12:00:00] Split-read clustering\n"
    "[2026-Sep-25 12:00:00] Genotyping\n"
    '{"label": "SAMPLE_A", "records": 30980, "bnd": 9187, "pass_bnd": 896}\n'
    "SAMPLE_A SAMPLE_B deid/SAMPLE_A.bam NA12878.chr20_chr21.bam IMP01 IMP10\n"
    "hs38DH GRCh38_full_analysis_set_plus_decoy_hla.fa human.hg38.excl.tsv\n"
    "Maximum resident set size (kbytes): 1212416\nElapsed (wall clock) time: 45:33.10\n"
    "Lietuviškas tekstas: ąčęėįšųūž — em dash, ±500 bp\n"
)


def synthetic(tmp):
    print("synthetic term list")
    terms = os.path.join(tmp, "terms")
    planted = [("ZQX-4481-ready", "SAMPLE_A:filename"), ("ZQX4481", "SAMPLE_A:SM"),
               ("Mx_77_Lab", "SAMPLE_B:RG-LB")]
    with open(terms, "w") as f:
        f.writelines(f"{t}\t{l}\n" for t, l in planted)
    os.chmod(terms, 0o600)
    bad = []
    for t, label in planted:
        for i, c in enumerate(contexts(t)):
            rc, out, _ = run((c + "\n").encode(), terms)
            if rc != 0 or residual(t, out.decode()) or f"[{label}]".encode() not in out:
                bad.append(f"{label} context {i}")
    check(f"positive: {len(planted)} planted terms x {len(contexts('x'))} contexts all masked, "
          f"each with its label", not bad, ", ".join(bad[:5]))
    rc, out, _ = run(ORDINARY.encode(), terms)
    check("negative: ordinary output passes byte-identical", rc == 0 and out == ORDINARY.encode(),
          f"rc={rc}")
    # longest match wins: the filename stem contains the SM value
    # Built at run time: an identifier-shaped NAME.bam literal in this source would,
    # rightly, trip scripts/identifier_gate.py.
    rc, out, _ = run(f"{planted[0][0]}.bam\n".encode(), terms)
    check("longest term wins (filename stem, not its SM part)",
          out == b"[SAMPLE_A:filename].bam\n", out.decode(errors="replace").strip())

    print("fail-closed")
    rc, out, err = run(b"ZQX4481 leaks if passed through\n", os.path.join(tmp, "absent"))
    check("missing list -> no output, non-zero exit", rc != 0 and out == b"", f"rc={rc} out={len(out)}B")
    empty = os.path.join(tmp, "empty")
    open(empty, "w").close()
    os.chmod(empty, 0o600)
    rc, out, _ = run(b"ZQX4481 leaks if passed through\n", empty)
    check("empty list -> no output, non-zero exit", rc != 0 and out == b"", f"rc={rc} out={len(out)}B")
    os.chmod(terms, 0o644)
    rc, out, err = run(b"ZQX4481 leaks if passed through\n", terms)
    check("group/world-readable list -> refused, no output", rc != 0 and out == b"",
          f"rc={rc} out={len(out)}B")
    os.chmod(terms, 0o600)
    rc, out, _ = run(b"ZQX4481\n", terms)
    check("control for the control: the same list at 0600 works again", rc == 0 and out == b"[SAMPLE_A:SM]\n")

    print("redact.sh -- CMD")
    env = dict(os.environ, HOME=tmp)          # redact.sh uses the default list under $HOME
    os.makedirs(os.path.join(tmp, "patient_data"), exist_ok=True)
    default = os.path.join(tmp, "patient_data", ".redact_terms")
    with open(default, "w") as f:
        f.writelines(f"{t}\t{l}\n" for t, l in planted)
    os.chmod(default, 0o600)
    p = subprocess.run(["bash", REDACT_SH, "--", "sh", "-c", "echo out ZQX4481; echo err Mx-77-lab >&2; exit 7"],
                       capture_output=True, env=env)
    out = p.stdout.decode()
    check("wrapper keeps the command's exit code (7)", p.returncode == 7, f"rc={p.returncode}")
    check("wrapper masks stdout and stderr", "[SAMPLE_A:SM]" in out and "[SAMPLE_B:RG-LB]" in out
          and not residual("ZQX4481", out) and not residual("Mx_77_Lab", out), out.strip())
    os.remove(default)
    p = subprocess.run(["bash", REDACT_SH, "--", "sh", "-c", "echo ZQX4481; exit 0"],
                       capture_output=True, env=env)
    check("wrapper with no list: output withheld, non-zero exit even though the command exited 0",
          p.returncode != 0 and p.stdout == b"", f"rc={p.returncode} out={len(p.stdout)}B")
    p = subprocess.run(["bash", REDACT_SH], input=b"ZQX4481\n", capture_output=True, env=env)
    check("pipe form with no list: output withheld", p.returncode != 0 and p.stdout == b"",
          f"rc={p.returncode}")


def real():
    print("real term list (terms are never printed)")
    if not os.path.exists(REAL):
        print("  NOT RUN -- no real term list yet")
        return False
    terms = []
    with open(REAL, encoding="utf-8") as f:
        for line in f:
            t, _, l = line.rstrip("\n").partition("\t")
            if t:
                terms.append((t, l))
    bad, n = [], 0
    for t, label in terms:
        for i, c in enumerate(contexts(t)):
            n += 1
            rc, out, _ = run((c + "\n").encode(), REAL)
            if rc != 0 or residual(t, out.decode("utf-8", "replace")):
                bad.append(f"{label} context {i}")
    check(f"positive: {len(terms)} real terms x {len(contexts('x'))} contexts = {n} cases, "
          f"none survives", not bad, ", ".join(bad[:5]))
    rc, out, _ = run(ORDINARY.encode(), REAL)
    changed = [i for i, (a, b) in enumerate(zip(ORDINARY.encode().splitlines(), out.splitlines())) if a != b]
    check("negative: ordinary output passes byte-identical under the REAL list "
          "(no real term is generic enough to mangle it)",
          rc == 0 and out == ORDINARY.encode(), f"rc={rc}, changed line numbers {changed}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--synthetic-only", action="store_true")
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        synthetic(tmp)
    ran_real = False if a.synthetic_only else real()
    print()
    if FAILS:
        print(f"FAILED: {len(FAILS)} check(s)")
        sys.exit(1)
    print("ALL PASSED" + ("" if ran_real else " (synthetic only)"))


if __name__ == "__main__":
    main()
