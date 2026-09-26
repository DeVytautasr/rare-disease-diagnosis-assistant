#!/usr/bin/env python3
"""Phase 12 Task 6(a): properties of SAMPLE_A and SAMPLE_B -- aggregate counts only.

    properties.py LABEL [--threads N]      LABEL = SAMPLE_A | SAMPLE_B
    properties.py record                   both samples' results -> the committed record

Reads ~/patient_data/deid/LABEL.bam (a link, never named in output). Each sample's
result is written to ~/patient_data/rerun_2026-09/properties_LABEL.json and printed:
counts only. Display through scripts/redact.sh anyway.

BLINDING. No coordinate, gene, read name, read group or per-position base is kept
or printed. The concordance step looks at single positions and keeps only totals,
because a single position's base is a patient genotype. The random positions are
drawn from the reference (seeded) and are not stored.

FIGURES AND THEIR DEFINITIONS (written into the record beside each figure)
  bwa_version         the VN value of each @PG record whose PN is bwa -- nothing else
                      from the line
  read_length         the modal query length (sequence length, soft clips included)
                      of primary records (not secondary, not supplementary) met in the
                      random-position sample below; the share at the mode is given
  nominal_depth       index-mapped records on chr1-22, X, Y (the BAM index's mapped
                      count: every mapped record, secondary and supplementary
                      included -- the Phase 0 definition of "primary mapped") x the
                      modal read length / the summed length of chr1-22, X, Y
  duplicate_fraction  from one samtools flagstat pass (QC-passed columns):
                      (i) records flagged 0x400 / all records, and
                      (ii) primary records flagged 0x400 / primary records
  median_insert_size  median |TLEN| of records that are primary, properly paired
                      (0x2), first in pair (0x40), not duplicate, MAPQ >= 20, with
                      the mate on the same contig and TLEN != 0: the first 50 such
                      records starting at or after each random position
  concordance         at 1,000 positions drawn uniformly over chr1-22, X, Y (seed
                      below; a position whose reference base is not A/C/G/T is redrawn):
                      the majority base among bases with base quality >= 20 in reads
                      with MAPQ >= 20 (unmapped, secondary, supplementary, QC-fail and
                      duplicate records excluded), where at least 10 such bases cover
                      the position, compared with the reference base. Kept: positions
                      drawn, positions meeting the depth floor, agreeing, disagreeing
"""
import json
import os
import random
import statistics
import subprocess
import sys
from collections import Counter

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
from common import DEID_DIR, RUN_DIR, die  # noqa: E402

REF = os.path.expanduser("~/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa")
SAMTOOLS = os.path.expanduser("~/miniconda3/envs/synth-hts/bin/samtools")
PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
SEED = 20260926
N_POSITIONS = 1000
PER_POSITION_INSERTS = 50
MIN_MAPQ, MIN_BQ, MIN_DEPTH = 20, 20, 10
FLAG_FILTER = 0x4 | 0x100 | 0x200 | 0x400 | 0x800


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


def positions(fa, lengths):
    rng = random.Random(SEED)
    total = sum(lengths[c] for c in PRIMARY)
    out = []
    while len(out) < N_POSITIONS:
        x = rng.randrange(total)
        for c in PRIMARY:
            if x < lengths[c]:
                break
            x -= lengths[c]
        if fa.fetch(c, x, x + 1).upper() in "ACGT":
            out.append((c, x))
    return out


def sample(label, threads):
    bam_path = os.path.join(DEID_DIR, f"{label}.bam")
    if not os.path.exists(bam_path):
        die(f"{label}: the de-identified link is missing")
    bam = pysam.AlignmentFile(bam_path)
    fa = pysam.FastaFile(REF)
    lengths = dict(zip(bam.references, bam.lengths))
    # bwa VN values only
    pg = bam.header.to_dict().get("PG", [])
    bwa_vn = [e.get("VN") for e in pg if str(e.get("PN", "")).lower() == "bwa"]
    # index statistics
    idx = {s.contig: s.mapped for s in bam.get_index_statistics()}
    mapped = sum(idx.get(c, 0) for c in PRIMARY)
    primary_len = sum(lengths[c] for c in PRIMARY)
    # random positions: read length, insert size, concordance
    lens, inserts = Counter(), []
    drawn = covered = agree = disagree = 0
    for chrom, pos in positions(fa, lengths):
        drawn += 1
        got = 0
        for r in bam.fetch(chrom, pos, pos + 2000):
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if r.reference_start >= pos:
                lens[r.query_length] += 1
                if (got < PER_POSITION_INSERTS and r.is_proper_pair and r.is_read1 and not r.is_duplicate
                        and r.mapping_quality >= MIN_MAPQ and r.template_length != 0
                        and r.next_reference_name == r.reference_name):
                    inserts.append(abs(r.template_length))
                    got += 1
        bases = Counter()
        for col in bam.pileup(chrom, pos, pos + 1, truncate=True, min_mapping_quality=MIN_MAPQ,
                              min_base_quality=MIN_BQ, flag_filter=FLAG_FILTER, ignore_orphans=False,
                              max_depth=100000):
            for b in col.get_query_sequences(add_indels=False):
                b = (b or "").upper()
                if b in ("A", "C", "G", "T"):
                    bases[b] += 1
        if sum(bases.values()) >= MIN_DEPTH:
            covered += 1
            top = max(bases.items(), key=lambda kv: (kv[1], kv[0]))[0]
            if top == fa.fetch(chrom, pos, pos + 1).upper():
                agree += 1
            else:
                disagree += 1
    bam.close()
    mode_len, mode_n = lens.most_common(1)[0]
    # one flagstat pass
    p = subprocess.run([SAMTOOLS, "flagstat", "-@", str(threads), "-O", "json", bam_path], capture_output=True,
                       text=True)
    if p.returncode != 0:
        die(f"{label}: samtools flagstat failed (exit {p.returncode})")
    fs = json.loads(p.stdout)["QC-passed reads"]
    rec = {"label": label,
           "bwa_version": {"values": sorted(set(bwa_vn)), "bwa_pg_records": len(bwa_vn)},
           "read_length": {"mode": mode_len, "share_at_mode": round(mode_n / sum(lens.values()), 4),
                           "records_examined": sum(lens.values()), "min": min(lens), "max": max(lens)},
           "nominal_depth": {"value": round(mapped * mode_len / primary_len, 2), "index_mapped_chr1_22_X_Y": mapped,
                             "read_length_used": mode_len, "chr1_22_X_Y_length": primary_len},
           "duplicate_fraction": {
               "all_records": {"numerator": fs["duplicates"], "denominator": fs["total"],
                               "value": round(fs["duplicates"] / fs["total"], 6)},
               "primary_records": {"numerator": fs.get("primary duplicates"), "denominator": fs.get("primary"),
                                   "value": (round(fs["primary duplicates"] / fs["primary"], 6)
                                             if fs.get("primary") else None)}},
           "flagstat": {k: fs.get(k) for k in ("total", "primary", "secondary", "supplementary", "duplicates",
                                               "primary duplicates", "mapped", "primary mapped")},
           "median_insert_size": {"value": statistics.median(inserts), "records": len(inserts),
                                  "per_position": PER_POSITION_INSERTS},
           "concordance": {"positions_drawn": drawn, "positions_meeting_floors": covered, "agree": agree,
                           "disagree": disagree, "value": round(agree / covered, 6) if covered else None},
           "floors": {"min_mapq": MIN_MAPQ, "min_base_quality": MIN_BQ, "min_depth": MIN_DEPTH,
                      "excluded_flags": hex(FLAG_FILTER)},
           "seed": SEED, "positions": N_POSITIONS}
    os.makedirs(RUN_DIR, exist_ok=True)
    with open(os.path.join(RUN_DIR, f"properties_{label}.json"), "w") as f:
        json.dump(rec, f, indent=1)
    print(json.dumps(rec, indent=1))
    return 0


def record():
    out = {"what": "Properties of the two patient BAMs, Phase 12 Task 6(a), 2026-09-26: aggregate counts only",
           "blinding": "no coordinate, gene, read or per-position base; the concordance keeps totals only",
           "definitions": definitions(),
           "thesis_values_compared": {"bwa_version": "0.7.17-r1188", "nominal_depth": "30.7 / 32.9-fold",
                                      "duplicate_fraction": "12.96 / 12.97%", "median_insert_size": "317 / 311 bp",
                                      "concordance": "99.90%"},
           "code": "scripts/patient/properties.py"}
    for label in ("SAMPLE_A", "SAMPLE_B"):
        p = os.path.join(RUN_DIR, f"properties_{label}.json")
        if not os.path.exists(p):
            die(f"{label}: no result yet")
        out[label] = json.load(open(p))
    dest = os.path.join(REPO, "stage1_igv_assistant", "results", "patient_properties_2026-09-26.json")
    if os.path.exists(dest):
        die("the record exists; a committed record is never overwritten")
    with open(dest, "w") as f:
        json.dump(out, f, indent=1)
    print("written: stage1_igv_assistant/results/patient_properties_2026-09-26.json")
    return 0


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "record":
        sys.exit(record())
    if len(sys.argv) < 2 or sys.argv[1] not in ("SAMPLE_A", "SAMPLE_B"):
        die("usage: properties.py SAMPLE_A|SAMPLE_B [--threads N] | record")
    threads = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[2] == "--threads" else 3
    sys.exit(sample(sys.argv[1], threads))
