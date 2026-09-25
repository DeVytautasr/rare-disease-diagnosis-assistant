"""Shared helpers for work on the two patient BAMs in ~/patient_data.

The samples are located by exact byte size against the Phase 0 record and are
referred to only as SAMPLE_A and SAMPLE_B. Nothing in this module prints or
returns a file name for display, and every error it raises is worded without
one: an exception message that quotes a path would carry the identifier out.
"""
import hashlib
import os
import stat
import sys

PATIENT_DIR = os.path.expanduser("~/patient_data")
DEID_DIR = os.path.join(PATIENT_DIR, "deid")
RUN_DIR = os.path.join(PATIENT_DIR, "rerun_2026-09")
TERMS = os.path.join(PATIENT_DIR, ".redact_terms")

# The Phase 0 record, as given in the 2026-09-25 prompt.
PHASE0 = {
    "SAMPLE_A": {"bytes": 38_959_428_903, "primary_mapped": 631_618_015, "pg_records": 34},
    "SAMPLE_B": {"bytes": 41_617_797_998, "primary_mapped": 677_604_873, "pg_records": 20},
}
# md5 of "name<TAB>length\n" for every @SQ line, in header order. The method was
# identified on 2026-09-25: it is the only one of 17 serializations tried that
# reproduces this value from the hs38DH .fai.
PHASE0_CONTIG_MD5 = "c325b75cc2827f4e6a1b39091a12f718"
PHASE0_AH_CONTIGS = 3171

BACKGROUND = os.path.expanduser("~/public_data/NA12878.chr20_chr21.bam")


def contig_md5(pairs):
    return hashlib.md5("".join(f"{n}\t{l}\n" for n, l in pairs).encode()).hexdigest()


def die(msg, code=2):
    """Exit with a message that must not contain a path or a file name."""
    print(f"STOPPED: {msg}", file=sys.stderr)
    sys.exit(code)


def top_level_entries():
    """[(name, lstat)] for the top level of PATIENT_DIR. Names stay in memory."""
    try:
        names = os.listdir(PATIENT_DIR)
    except OSError as e:
        die(f"cannot list the patient directory ({type(e).__name__}, errno {e.errno})")
    out = []
    for n in names:
        try:
            out.append((n, os.lstat(os.path.join(PATIENT_DIR, n))))
        except OSError as e:
            die(f"cannot stat an entry ({type(e).__name__}, errno {e.errno})")
    return out


def index_for(bam):
    """(index path, form) or (None, None). Form is a generic description."""
    for path, form in ((bam + ".bai", "<bam>.bai"), (bam[:-4] + ".bai", "<stem>.bai"),
                       (bam + ".csi", "<bam>.csi")):
        if os.path.isfile(path):
            return path, form
    return None, None


def locate():
    """{label: {"bam": path, "index": path, "index_form": str}} by exact size.

    Only regular, non-hidden top-level *.bam files are candidates -- never a
    symlink (deid/ links point back at these) and never an rsync temporary."""
    bams = [(n, st) for n, st in top_level_entries()
            if not n.startswith(".") and n.endswith(".bam") and stat.S_ISREG(st.st_mode)]
    found = {}
    for label, rec in PHASE0.items():
        hits = [n for n, st in bams if st.st_size == rec["bytes"]]
        if len(hits) != 1:
            die(f"{label}: {len(hits)} regular .bam files have the Phase 0 size "
                f"{rec['bytes']:,} (need exactly 1); sizes present: "
                f"{sorted(st.st_size for _, st in bams)}")
        bam = os.path.join(PATIENT_DIR, hits[0])
        idx, form = index_for(bam)
        found[label] = {"bam": bam, "index": idx, "index_form": form}
    return found


def deid_bam(label):
    return os.path.join(DEID_DIR, f"{label}.bam")
