#!/usr/bin/env python3
"""Build the two private term lists from the patient BAMs, printing no term.

  ~/patient_data/.redact_terms     "term<TAB>label" per line, mode 600 -- read by
                                   scripts/redact.py, which masks every term in
                                   anything displayed
  ~/patient_data/.identifier_list  one term per line, mode 600 -- passed to
                                   scripts/identifier_gate.py --identifiers

Sources, per sample:
  SM         every @RG SM value (required by the prompt)
  filename   the BAM's file-name stem, and the stem of every other non-hidden
             top-level file whose name contains it (required by the prompt)
  RG-ID/LB/PU/CN/DS  other @RG values that can name a person, a library, a
             flowcell or a centre; redaction only, and only if they pass the
             safety filter below
  component  identifier-shaped parts of the SM values and stems (split on
             . _ - and space; 4+ characters and containing a digit)

Safety filter: a term shorter than 4 characters (3 for CN), all digits, or on
the generic-genomics stoplist is dropped, because masking it would mangle
ordinary output ("GRCh38", "L001", "ready"). Dropped terms are counted, never
shown. Written atomically with mode 600 from the start.
"""
import os
import re
import stat
import sys
import tempfile

import pysam

from common import PATIENT_DIR, TERMS, die, locate, top_level_entries

IDLIST = os.path.join(PATIENT_DIR, ".identifier_list")
STOP = re.compile(r"^(grch3[78]|hg(19|38)|b37|hs3[78]d?[h5]?|hs38dh|hs37d5|t2t|chm13|"
                  r"l\d{3}|s\d{1,3}|r[12]|v\d+|lane\d*|ready|sorted|dedup|markdup|md|recal|"
                  r"bqsr|final|merged|aligned|realigned|rmdup|bwa|mem|wgs|wes|bam|cram|"
                  r"illumina|novaseq\d*|hiseq\w*|mgiseq\w*|dnbseq\w*|unknown|none|null)$", re.I)
SPLIT = re.compile(r"[._\-\s]+")


def keep(term, min_len=4):
    return (len(term) >= min_len and not term.isdigit() and not STOP.match(term))


def write_private(path, lines):
    fd, tmp = tempfile.mkstemp(dir=PATIENT_DIR, prefix=".tmp_terms_")
    try:
        os.fchmod(fd, 0o600)
        with os.fdopen(fd, "w") as f:
            f.writelines(lines)
        os.replace(tmp, path)
    except BaseException:
        os.unlink(tmp)
        raise
    mode = stat.S_IMODE(os.stat(path).st_mode)
    if mode != 0o600:
        die(f"a private list ended up with mode {oct(mode)}, not 0o600")


def main():
    samples = locate()
    names = [n for n, st in top_level_entries() if not n.startswith(".") and stat.S_ISREG(st.st_mode)]
    terms = {}          # term -> set(labels)
    dropped = 0

    def add(term, label, min_len=4):
        nonlocal dropped
        term = term.strip()
        if not term:
            return
        if keep(term, min_len):
            terms.setdefault(term, set()).add(label)
        else:
            dropped += 1

    for label, rec in samples.items():
        with pysam.AlignmentFile(rec["bam"], "rb") as f:
            rgs = f.header.to_dict().get("RG", [])
        if not rgs:
            die(f"{label}: header has no @RG lines -- no SM value to build from")
        sms = {rg.get("SM", "") for rg in rgs} - {""}
        if not sms:
            die(f"{label}: no @RG line carries an SM value")
        for sm in sms:
            add(sm, f"{label}:SM", min_len=1)          # an SM is always masked, however short
            for part in SPLIT.split(sm):
                if re.search(r"\d", part):
                    add(part, f"{label}:component")
        stem = os.path.basename(rec["bam"])[:-4]
        add(stem, f"{label}:filename", min_len=1)
        for part in SPLIT.split(stem):
            if re.search(r"\d", part):
                add(part, f"{label}:component")
        for n in names:                                  # e.g. <stem>.bam.md5, <stem>_report.pdf
            if stem in n and n != os.path.basename(rec["bam"]):
                other = n.split(".")[0]
                if other != stem:
                    add(other, f"{label}:filename")
        for rg in rgs:
            for key in ("ID", "LB", "PU", "DS"):
                if key in rg:
                    add(str(rg[key]), f"{label}:RG-{key}")
            if "CN" in rg:
                add(str(rg["CN"]), f"{label}:RG-CN", min_len=3)

    # A term found for both samples keeps both labels.
    rows = []
    for term in sorted(terms, key=lambda t: (-len(t), t)):
        labels = sorted(terms[term])
        if "\t" in term or "\n" in term:
            die("a term contains a tab or newline; refusing to write a list that would parse wrongly")
        rows.append((term, ",".join(labels)))
    write_private(TERMS, [f"{t}\t{l}\n" for t, l in rows])
    # The gate list: identifiers only -- CN and DS name a centre or describe a run,
    # which the gate's purpose (patient identifiers) does not cover.
    gate = [t for t, l in rows if any(x.split(":")[1] in ("SM", "filename", "component", "RG-ID",
                                                          "RG-LB", "RG-PU") for x in l.split(","))]
    write_private(IDLIST, [f"{t}\n" for t in gate])

    by_source = {}
    for t, l in rows:
        for x in l.split(","):
            by_source[x] = by_source.get(x, 0) + 1
    print(f"redaction terms: {len(rows)} written (mode 600); dropped by the safety filter: {dropped}")
    for k in sorted(by_source):
        print(f"  {k:24s} {by_source[k]}")
    lens = sorted(len(t) for t, _ in rows)
    print(f"term lengths: min {lens[0]}, max {lens[-1]}")
    shared = sum(1 for _, l in rows if "SAMPLE_A" in l and "SAMPLE_B" in l)
    print(f"terms found in both samples: {shared}")
    print(f"identifier list for the gate: {len(gate)} entries (mode 600)")


if __name__ == "__main__":
    main()
