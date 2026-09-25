#!/usr/bin/env python3
"""Verify that the re-downloaded patient BAMs are the Phase 0 files (Task 1).

    verify_bams.py header   LABEL          header, index and redaction checks (seconds)
    verify_bams.py flagstat LABEL THREADS  full read of the BAM (minutes)
    verify_bams.py pa       LABEL          primary reads on chr20+chr21, pa-tag count
    verify_bams.py compare                 every figure against the Phase 0 record

LABEL is SAMPLE_A or SAMPLE_B. Each step writes its JSON under
~/patient_data/rerun_2026-09/ and prints counts only -- no name, no header text,
no read. Display it through scripts/redact.sh anyway.

The @PG count is taken from the parsed header and, independently, from
`samtools head`, because `samtools view -H` adds a @PG record of its own to what
it prints (Phase 0 counted with `samtools head` for that reason).

The index is checked against the BAM itself, not just opened: the per-contig
mapped counts stored in the .bai must sum to flagstat's "mapped" total, which
comes from reading every record. A stale or foreign index fails that.
"""
import json
import os
import sys
import time

import pysam

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import (BACKGROUND, PHASE0, PHASE0_AH_CONTIGS, PHASE0_CONTIG_MD5, RUN_DIR, TERMS,
                    contig_md5, die, locate)
from redact import compile_terms, load_terms  # noqa: E402  (scripts/redact.py)

# Fixed test regions for the index check. Chosen before looking at any patient
# data and unrelated to any call; their coordinates are never printed.
TEST_REGIONS = [("chr1", 10_000_000, 10_010_000), ("chr7", 50_000_000, 50_010_000),
                ("chr20", 40_000_000, 40_010_000), ("chrX", 80_000_000, 80_010_000)]


def out_path(label, step):
    os.makedirs(RUN_DIR, mode=0o700, exist_ok=True)
    return os.path.join(RUN_DIR, f"verify_{label}_{step}.json")


def save(label, step, rec):
    with open(out_path(label, step), "w") as f:
        json.dump(rec, f, indent=1)


def ah_set(header_dict):
    return {sq["SN"] for sq in header_dict.get("SQ", []) if "AH" in sq}


def step_header(label):
    s = locate()[label]
    rec = {"label": label, "bytes": os.path.getsize(s["bam"]), "index_form": s["index_form"]}
    if not s["index"]:
        die(f"{label}: no index beside the BAM -- report and stop (never re-index)")
    rec["index_bytes"] = os.path.getsize(s["index"])
    with pysam.AlignmentFile(s["bam"], "rb") as f:
        h = f.header.to_dict()
        rec["pg_records_parsed"] = len(h.get("PG", []))
        rec["rg_records"] = len(h.get("RG", []))
        rec["sq_records"] = len(h.get("SQ", []))
        rec["hd_so"] = h.get("HD", {}).get("SO")
        rec["contig_md5"] = contig_md5(zip(f.references, f.lengths))
        ah = ah_set(h)
        rec["ah_contigs"] = len(ah)
        rec["has_index"] = f.has_index()
        # index usability: every read fetched for a region must overlap it
        ok_regions, fetched = 0, 0
        for chrom, a, b in TEST_REGIONS:
            n = bad = 0
            for r in f.fetch(chrom, a, b):
                n += 1
                if r.reference_name != chrom or r.reference_start >= b or \
                        (r.reference_end or r.reference_start + 1) <= a:
                    bad += 1
            fetched += n
            ok_regions += (n > 0 and bad == 0)
        rec["index_fetch_regions_ok"] = f"{ok_regions}/{len(TEST_REGIONS)}"
        rec["index_fetch_reads"] = fetched
    with pysam.AlignmentFile(BACKGROUND, "rb") as g:
        rec["ah_set_equals_na12878_background"] = ah == ah_set(g.header.to_dict())
    alt = os.path.expanduser("~/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.alt")
    if os.path.exists(alt):
        names = {l.split("\t")[0] for l in open(alt) if l.strip() and not l.startswith("@")}
        rec["ah_set_equals_hs38DH_alt_file"] = ah == names
    head = pysam.samtools.head(s["bam"])
    rec["pg_records_samtools_head"] = sum(1 for l in head.splitlines() if l.startswith("@PG\t"))
    # idxstats from the index alone
    mapped = unmapped = 0
    for line in pysam.idxstats(s["bam"]).splitlines():
        name, length, m, u = line.split("\t")
        mapped += int(m)
        unmapped += int(u)
    rec["idxstats_mapped_sum"] = mapped
    rec["idxstats_unmapped_sum"] = unmapped
    # the header is the most identifier-dense text there is: filter it and check
    pat, repl = compile_terms(load_terms(TERMS))
    lines = head.splitlines()
    filtered = [pat.sub(repl, l) for l in lines]
    rec["header_lines"] = len(lines)
    rec["header_lines_changed_by_redaction"] = sum(a != b for a, b in zip(lines, filtered))
    rec["header_residual_matches_after_redaction"] = sum(bool(pat.search(l)) for l in filtered)
    save(label, "header", rec)
    print(json.dumps(rec, indent=1))


def step_flagstat(label, threads):
    s = locate()[label]
    t0 = time.monotonic()
    out = pysam.flagstat("-@", str(threads), "-O", "json", s["bam"])
    rec = {"label": label, "threads": threads, "wall_seconds": round(time.monotonic() - t0, 1),
           "samtools": pysam.version.__samtools_version__, "flagstat": json.loads(out)}
    save(label, "flagstat", rec)
    qp, qf = rec["flagstat"]["QC-passed reads"], rec["flagstat"]["QC-failed reads"]
    print(json.dumps({"label": label, "wall_seconds": rec["wall_seconds"],
                      "total": [qp["total"], qf["total"]], "mapped": [qp["mapped"], qf["mapped"]],
                      "primary_mapped": [qp["primary mapped"], qf["primary mapped"]],
                      "secondary": [qp["secondary"], qf["secondary"]],
                      "supplementary": [qp["supplementary"], qf["supplementary"]],
                      "duplicates": [qp["duplicates"], qf["duplicates"]]}, indent=1))


def step_pa(label):
    s = locate()[label]
    t0 = time.monotonic()
    rec = {"label": label, "contigs": ["chr20", "chr21"]}
    n = pa = 0
    with pysam.AlignmentFile(s["bam"], "rb") as f:
        for chrom in rec["contigs"]:
            for r in f.fetch(chrom):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                n += 1
                if r.has_tag("pa"):
                    pa += 1
    rec.update(primary_reads_scanned=n, pa_tagged=pa,
               pa_fraction=round(pa / n, 6) if n else None,
               wall_seconds=round(time.monotonic() - t0, 1))
    save(label, "pa", rec)
    print(json.dumps(rec, indent=1))


def step_compare():
    rows, ok = [], True

    def row(label, what, got, want):
        nonlocal ok
        good = got == want
        ok &= good
        rows.append((label, what, got, want, "MATCH" if good else "DIFFERS"))

    for label, p0 in PHASE0.items():
        try:
            h = json.load(open(out_path(label, "header")))
            fs = json.load(open(out_path(label, "flagstat")))
            pa = json.load(open(out_path(label, "pa")))
        except FileNotFoundError as e:
            die(f"{label}: a step has not been run yet ({os.path.basename(e.filename)})")
        qp, qf = fs["flagstat"]["QC-passed reads"], fs["flagstat"]["QC-failed reads"]
        row(label, "bytes", h["bytes"], p0["bytes"])
        row(label, "primary mapped (flagstat, QC-passed)", qp["primary mapped"], p0["primary_mapped"])
        row(label, "@PG records (samtools head)", h["pg_records_samtools_head"], p0["pg_records"])
        row(label, "@PG records (parsed header)", h["pg_records_parsed"], p0["pg_records"])
        row(label, "contig (name,length) md5", h["contig_md5"], PHASE0_CONTIG_MD5)
        row(label, "AH-tagged contigs", h["ah_contigs"], PHASE0_AH_CONTIGS)
        row(label, "AH set = NA12878 background's", h["ah_set_equals_na12878_background"], True)
        if "ah_set_equals_hs38DH_alt_file" in h:
            row(label, "AH set = hs38DH .alt file", h["ah_set_equals_hs38DH_alt_file"], True)
        row(label, "index: .bai mapped sum = flagstat mapped", h["idxstats_mapped_sum"],
            qp["mapped"] + qf["mapped"])
        row(label, "index: fetch returns only overlapping reads", h["index_fetch_regions_ok"],
            f"{len(TEST_REGIONS)}/{len(TEST_REGIONS)}")
        row(label, "reads carrying pa (chr20+chr21) > 0", pa["pa_tagged"] > 0, True)
        row(label, "header: no term survives redaction", h["header_residual_matches_after_redaction"], 0)
    w = max(len(r[1]) for r in rows)
    for label, what, got, want, v in rows:
        g = f"{got:,}" if isinstance(got, int) and not isinstance(got, bool) else str(got)
        e = f"{want:,}" if isinstance(want, int) and not isinstance(want, bool) else str(want)
        print(f"{label}  {what:<{w}}  {g:>34}  {e:>34}  {v}")
    print("\nALL FIGURES MATCH" if ok else "\nAT LEAST ONE FIGURE DIFFERS -- STOP")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    step = sys.argv[1]
    if step == "compare":
        step_compare()
    elif step == "header":
        step_header(sys.argv[2])
    elif step == "flagstat":
        step_flagstat(sys.argv[2], int(sys.argv[3]))
    elif step == "pa":
        step_pa(sys.argv[2])
    else:
        die(f"unknown step {step!r}")
