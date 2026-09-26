#!/usr/bin/env python3
"""Phase 12 Task 6(b): the supplementary-flag contradiction -- aggregate counts only.

    supplementary_flags.py count LABEL [--threads N]   one streaming pass over the records flagged 0x800
    supplementary_flags.py oldscan LABEL CHROM          the 2026-08-29 scan, reconstructed, vs direct counts
    supplementary_flags.py pg LABEL                     @PG PN values in header order, gated before display
    supplementary_flags.py record                       both samples' results -> the committed record

LABEL is SAMPLE_A or SAMPLE_B, read through ~/patient_data/deid/LABEL.bam (a link,
never named in output). Results go to ~/patient_data/rerun_2026-09/suppl_*_LABEL.json
and are printed as counts only. Display through scripts/redact.sh anyway.

THE CONTRADICTION. samtools flagstat (scripts/patient/verify_bams.py, 2026-09-25)
counts 3,694,802 and 3,726,497 supplementary records in the two BAMs, while
stage1_igv_assistant/results/REAL_PATIENT_DATA_VALIDATION.md (2026-08-29) says
"bwa mem -M produces zero 0x800 reads -- confirmed at chromosome scale over ~9-10 M
reads per chromosome per BAM". The code of that scan was never committed (the
commit that added the document, e8e3c83, holds the document only), so `oldscan`
runs reconstructions of it, each written out below, beside direct counts.

BLINDING. Only counts by contig class (from RNAME and nothing else), flag bits and
tag presence are kept. No position, read name, read group or sequence is kept.

FIGURES AND THEIR DEFINITIONS (written into the record beside each figure)
  flag_0x800          records whose FLAG has bit 0x800 set, whatever else is set, from
                      `samtools view -f 0x800` over the whole BAM, counted by the contig
                      class of their RNAME
  with_SA             of those, records carrying an SA tag
  also_0x100          of those, records that also have bit 0x100 (secondary) set
  hard_clipped        of those, records whose CIGAR contains an H operation
  sa_first_class      of those carrying SA, counts by the contig class of the first SA
                      entry's contig name
  flagstat_suppl      samtools flagstat's "supplementary" (QC-passed + QC-failed):
                      records with 0x800 and WITHOUT 0x100 -- flagstat counts a record
                      carrying both bits as secondary -- so it should equal flag_0x800
                      minus also_0x100
  contig_classes      primary = chr1-22, chrX, chrY; chrM; decoy = name ends _decoy;
                      alt = name ends _alt; HLA = name starts HLA-; EBV = chrEBV;
                      unlocalised/unplaced = name ends _random, or starts chrUn_ and is
                      not a decoy; other = anything else
  oldscan             on one chromosome, one pysam pass over fetch(CHROM), i.e. every
                      record whose RNAME is CHROM:
                      records_on_chrom -- every record returned (to set beside the
                      document's "~9-10 M reads per chromosome");
                      direct_0x800 -- records with 0x800 set, no other filter;
                      direct_0x800_samtools -- `samtools view -c -f 0x800 BAM CHROM`, an
                      independent count of the same thing;
                      recon_tool_filter -- skip unmapped, secondary and supplementary
                      records (the read filter of the split-read tool,
                      bam_tools.py line 1467), then count 0x800: zero by construction;
                      recon_secondary_skip -- skip unmapped and secondary records (the
                      filter at bam_tools.py line 1004), then count 0x800;
                      secondary -- records with 0x100;
                      sa_any -- records with an SA tag; sa_on_primary -- records with an
                      SA tag that pass the tool filter; sa_on_secondary, sa_on_supplementary
                      -- records with an SA tag and 0x100, and with 0x800 and not 0x100;
                      the document's identity "SA_any - SA_on_primary equals the secondary
                      count" is tested as sa_any - sa_on_primary == secondary
  pg_programs         the PN value of each @PG header record, in header order (a record
                      without PN is listed as "(no PN)"), shown only after the identifier
                      gate's scanner and the redaction terms find nothing in them
"""
import json
import os
import subprocess
import sys
from collections import Counter

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from common import DEID_DIR, RUN_DIR, TERMS, die  # noqa: E402

SAMTOOLS = os.path.expanduser("~/miniconda3/envs/synth-hts/bin/samtools")
PRIMARY = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
CLASSES = ["primary", "chrM", "unlocalised/unplaced", "decoy", "alt", "HLA", "EBV", "other"]
LABELS = ("SAMPLE_A", "SAMPLE_B")
IDENTIFIERS = os.path.expanduser("~/patient_data/.identifier_list")
RECORD = os.path.join(REPO, "stage1_igv_assistant", "results", "supplementary_flags_2026-09-26.json")


def definitions():
    doc = __doc__.split("FIGURES AND THEIR DEFINITIONS (written into the record beside each figure)")[1]
    out, key = {}, None
    for line in doc.splitlines():
        if line.startswith("  ") and not line.startswith("    ") and line.strip():
            key, _, rest = line.strip().partition(" ")
            out[key] = rest.strip()
        elif key and line.strip():
            out[key] += " " + line.strip()
    return {k: " ".join(v.split()) for k, v in out.items()}


def contig_class(name):
    if name in PRIMARY:
        return "primary"
    if name == "chrM":
        return "chrM"
    if name.endswith("_decoy"):
        return "decoy"
    if name.endswith("_alt"):
        return "alt"
    if name.startswith("HLA-"):
        return "HLA"
    if name == "chrEBV":
        return "EBV"
    if name.endswith("_random") or name.startswith("chrUn_"):
        return "unlocalised/unplaced"
    return "other"


def bam_path(label):
    if label not in LABELS:
        die("LABEL must be SAMPLE_A or SAMPLE_B")
    p = os.path.join(DEID_DIR, f"{label}.bam")
    if not os.path.exists(p):
        die(f"{label}: the de-identified link is missing")
    return p


def save(name, label, rec):
    os.makedirs(RUN_DIR, mode=0o700, exist_ok=True)
    with open(os.path.join(RUN_DIR, f"suppl_{name}_{label}.json"), "w") as f:
        json.dump(rec, f, indent=1)
    print(json.dumps(rec, indent=1))


def step_count(label, threads):
    path = bam_path(label)
    zero = lambda: {c: 0 for c in CLASSES}  # noqa: E731
    n, sa, sec, hard = zero(), zero(), zero(), zero()
    sa_first = {c: zero() for c in CLASSES}
    p = subprocess.Popen([SAMTOOLS, "view", "--no-PG", "-@", str(threads), "-f", "0x800", path],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.split("\t", 11)
        flag, cls, cigar = int(f[1]), contig_class(f[2]), f[5]
        n[cls] += 1
        if flag & 0x100:
            sec[cls] += 1
        if "H" in cigar:
            hard[cls] += 1
        tags = f[11] if len(f) > 11 else ""
        i = tags.find("SA:Z:")
        if i == 0 or (i > 0 and tags[i - 1] == "\t"):
            sa[cls] += 1
            sa_first[cls][contig_class(tags[i + 5:].split(",", 1)[0])] += 1
    rc = p.wait()
    if rc != 0:
        die(f"{label}: samtools view failed (exit {rc})")
    fs = json.load(open(os.path.join(RUN_DIR, f"verify_{label}_flagstat.json")))["flagstat"]
    flagstat_suppl = fs["QC-passed reads"]["supplementary"] + fs["QC-failed reads"]["supplementary"]
    total = sum(n.values())
    rec = {"label": label,
           "flag_0x800": {"total": total, "by_class": n},
           "with_SA": {"total": sum(sa.values()), "by_class": sa},
           "also_0x100": {"total": sum(sec.values()), "by_class": sec},
           "hard_clipped": {"total": sum(hard.values()), "by_class": hard},
           "sa_first_class": {c: {k: v for k, v in d.items() if v} for c, d in sa_first.items() if any(d.values())},
           "flagstat_suppl": flagstat_suppl,
           "flag_0x800_minus_also_0x100_equals_flagstat_suppl": total - sum(sec.values()) == flagstat_suppl}
    save("count", label, rec)
    return 0


def step_oldscan(label, chrom):
    path = bam_path(label)
    c = Counter()
    with pysam.AlignmentFile(path) as bam:
        if chrom not in bam.references:
            die("unknown contig")
        for r in bam.fetch(chrom):
            c["records_on_chrom"] += 1
            supp, sec, has_sa = r.is_supplementary, r.is_secondary, r.has_tag("SA")
            c["direct_0x800"] += supp
            c["secondary"] += sec
            c["sa_any"] += has_sa
            c["sa_on_secondary"] += has_sa and sec
            c["sa_on_supplementary"] += has_sa and supp and not sec
            if not (r.is_unmapped or sec):
                c["recon_secondary_skip"] += supp
                if not supp:
                    c["sa_on_primary"] += has_sa
                    c["recon_tool_filter"] += supp  # always 0: supplementary records were skipped
    p = subprocess.run([SAMTOOLS, "view", "-c", "-f", "0x800", path, chrom], capture_output=True, text=True)
    if p.returncode != 0:
        die(f"{label}: samtools view -c failed (exit {p.returncode})")
    rec = {"label": label, "chrom": chrom, **{k: c[k] for k in (
        "records_on_chrom", "direct_0x800", "recon_tool_filter", "recon_secondary_skip", "secondary",
        "sa_any", "sa_on_primary", "sa_on_secondary", "sa_on_supplementary")},
           "direct_0x800_samtools": int(p.stdout.strip())}
    rec["identity_sa_any_minus_sa_on_primary_equals_secondary"] = rec["sa_any"] - rec["sa_on_primary"] == rec["secondary"]
    save(f"oldscan_{chrom}", label, rec)
    return 0


def step_pg(label):
    from identifier_gate import Gate  # scripts/identifier_gate.py
    from redact import compile_terms, load_terms  # scripts/redact.py
    with pysam.AlignmentFile(bam_path(label)) as bam:
        pns = [str(e.get("PN", "(no PN)")) for e in bam.header.to_dict().get("PG", [])]
    if not os.path.exists(IDENTIFIERS):
        die("no identifier list: the PN values are not shown")
    idents = [l.strip() for l in open(IDENTIFIERS) if l.strip() and not l.startswith("#")]
    gate = Gate(idents, None, os.path.expanduser("~"), REPO)
    text = "\n".join(pns)
    gate.scan("pg", "PN values", text, paths=True)
    if not os.path.exists(TERMS):
        die("no redaction term list: the PN values are not shown")
    pat, _ = compile_terms(load_terms(TERMS))
    residual = sum(bool(pat.search(pn)) for pn in pns)
    if gate.fired or residual:
        die(f"{label}: the gate fired on the PN values ({len(gate.fired)} hits, {residual} term matches); "
            "not shown and not recorded")
    rec = {"label": label, "pg_records": len(pns), "pg_programs_in_header_order": pns,
           "pn_counts": dict(Counter(pns)), "gate": "identifier gate scanner: nothing fired",
           "redaction_terms_matched": residual}
    save("pg", label, rec)
    return 0


def step_record():
    out = {"what": "The supplementary-flag contradiction, Phase 12 Task 6(b), 2026-09-26: aggregate counts only",
           "contradiction": {"flagstat_2026_09_25": {"SAMPLE_A": 3694802, "SAMPLE_B": 3726497},
                             "validation_doc_2026_08_29": "REAL_PATIENT_DATA_VALIDATION.md: bwa mem -M produces "
                                                          "zero 0x800 reads -- confirmed at chromosome scale over "
                                                          "~9-10 M reads per chromosome per BAM",
                             "old_scan_code": "never committed: e8e3c83 added the document only; oldscan runs "
                                              "reconstructions"},
           "blinding": "counts by contig class, flag bit and tag presence only",
           "definitions": definitions(),
           "code": "scripts/patient/supplementary_flags.py"}
    for label in LABELS:
        got = {}
        for f in sorted(os.listdir(RUN_DIR)):
            if f.startswith("suppl_") and f.endswith(f"_{label}.json"):
                got[f[len("suppl_"):-len(f"_{label}.json")]] = json.load(open(os.path.join(RUN_DIR, f)))
        if "count" not in got or "pg" not in got or not any(k.startswith("oldscan_") for k in got):
            die(f"{label}: a step has not been run yet")
        out[label] = got
    if os.path.exists(RECORD):
        die("the record exists; a committed record is never overwritten")
    with open(RECORD, "w") as f:
        json.dump(out, f, indent=1)
    print("written: stage1_igv_assistant/results/supplementary_flags_2026-09-26.json")
    return 0


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[:1] == ["count"] and len(a) >= 2:
        sys.exit(step_count(a[1], int(a[3]) if len(a) > 3 and a[2] == "--threads" else 3))
    if a[:1] == ["oldscan"] and len(a) == 3:
        sys.exit(step_oldscan(a[1], a[2]))
    if a[:1] == ["pg"] and len(a) == 2:
        sys.exit(step_pg(a[1]))
    if a[:1] == ["record"]:
        sys.exit(step_record())
    die("usage: supplementary_flags.py count LABEL [--threads N] | oldscan LABEL CHROM | pg LABEL | record")
