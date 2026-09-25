"""Shared paths, tools and helpers for the synthetic positive control
(12 heterozygous balanced t(20;21) implants in the NA12878 chr20+chr21 slice).

Every path is under the user's home and built from "~", so nothing here names a
machine. Large outputs (BAM, BCF, FASTQ) go under ~/public_data/sim and stay out
of git; the code and the JSON records that describe them are committed.
"""
import hashlib
import json
import os
import subprocess
import sys

HOME = os.path.expanduser("~")
BG = os.path.join(HOME, "public_data", "NA12878.chr20_chr21.bam")
BG_BYTES = 1_636_070_786                        # recorded when the slice was made, 2026-09-25
REF = os.path.join(HOME, "reference", "GRCh38_full_analysis_set_plus_decoy_hla.fa")
EXCL = os.path.join(HOME, "reference", "human.hg38.excl.tsv")
SIM = os.path.join(HOME, "public_data", "sim")
WORK = os.path.join(SIM, "work")
BAMS = os.path.join(SIM, "bams")                # make_demo_bundle.py reads IMP01/IMP10 here
DELLY_DIR = os.path.join(SIM, "delly")          # ... and their BCFs here

CONDA = os.path.join(HOME, "miniconda3", "envs")
BWA = os.path.join(CONDA, "synth-bwa", "bin", "bwa")
ART = os.path.join(CONDA, "synth-art", "bin", "art_illumina")
SAMTOOLS = os.path.join(CONDA, "synth-hts", "bin", "samtools")
BCFTOOLS = os.path.join(CONDA, "synth-hts", "bin", "bcftools")
DELLY = os.path.join(HOME, "tools", "delly", "delly")

# The background's own alignment: bwa 0.7.15-r1140 mem -Y -K 100000000 (its @PG).
BWA_ARGS = ["mem", "-Y", "-K", "100000000"]
READ_GROUP = r"@RG\tID:SIM\tSM:NA12878\tLB:NA12878\tPL:illumina"
READ_LEN = 150
CONTIGS = ("chr20", "chr21")


def tilde(path):
    """A path as it may appear in a committed record: ~ instead of the home."""
    return path.replace(HOME, "~", 1) if path.startswith(HOME) else path


def masks():
    """{contig: [(start, end, kind)]} from the exclude template, 0-based half-open."""
    out = {c: [] for c in CONTIGS}
    with open(EXCL) as f:
        for line in f:
            p = line.split()
            if len(p) >= 3 and p[0] in out:
                out[p[0]].append((int(p[1]), int(p[2]), p[3] if len(p) > 3 else ""))
    return out


def run(cmd, log=None, **kw):
    """Run a command, echoing it (with ~ paths) first; raise on a non-zero exit."""
    shown = " ".join(tilde(str(c)) for c in cmd)
    print(f"$ {shown}", file=log or sys.stderr, flush=True)
    return subprocess.run([str(c) for c in cmd], check=True, **kw)


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def tool_versions():
    """What each tool reports when invoked -- recorded, never assumed."""
    def first(cmd, pat):
        p = subprocess.run(cmd, capture_output=True, text=True)
        for line in (p.stdout + p.stderr).splitlines():
            if pat in line:
                return line.strip()
        return None
    return {"bwa": first([BWA], "Version"), "art_illumina": first([ART], "Version"),
            "samtools": first([SAMTOOLS, "--version"], "samtools"),
            "bcftools": first([BCFTOOLS, "--version"], "bcftools"),
            "htslib (samtools)": first([SAMTOOLS, "--version"], "htslib"),
            "delly": first([DELLY], "Version")}


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as f:
        json.dump(obj, f, indent=1, sort_keys=False)
    os.replace(path + ".tmp", path)
