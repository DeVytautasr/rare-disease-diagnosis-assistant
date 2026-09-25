#!/usr/bin/env python3
"""Create ~/patient_data/deid/SAMPLE_A.bam and SAMPLE_B.bam, with index links,
and verify them by idxstats THROUGH the links (Task 2).

    make_deid_links.py            create (or confirm) the links, then verify
    make_deid_links.py --control  also run the swapped-link control

Every later command names only these link paths, so no identifier reaches a
command line, a process listing, a log or a BCF header's file field. (The @RG
SM value still travels inside the BAM and into delly's BCF sample column,
which is why the BCFs stay under ~/patient_data.)

Verification of a link, for its label:
  - it resolves to a regular file of that label's Phase 0 byte size
  - htslib finds the index through it (the index link is <label>.bam.bai)
  - idxstats through the link equals idxstats on the original, line for line

The --control run builds, in a temporary directory, a SAMPLE_A link that points
at SAMPLE_B's BAM, and requires the same verification to reject it.
Prints counts only.
"""
import os
import sys
import tempfile

import pysam

from common import DEID_DIR, PHASE0, deid_bam, die, locate


def link(target, path):
    if os.path.islink(path):
        if os.readlink(path) != target:
            die("an existing deid link points somewhere else -- not overwriting it")
        return "kept"
    if os.path.exists(path):
        die("a non-link file already sits at a deid path -- not overwriting it")
    os.symlink(target, path)
    return "created"


def verify(label, bam_link, original):
    """(ok, reasons[], idxstats text) -- reasons carry no names."""
    reasons = []
    real = os.path.realpath(bam_link)
    if not os.path.isfile(real):
        return False, ["does not resolve to a regular file"], ""
    size = os.path.getsize(real)
    if size != PHASE0[label]["bytes"]:
        reasons.append(f"resolves to a file of {size:,} bytes, not {label}'s {PHASE0[label]['bytes']:,}")
    with pysam.AlignmentFile(bam_link, "rb") as f:
        if not f.has_index():
            reasons.append("htslib finds no index through the link")
    through = pysam.idxstats(bam_link)
    direct = pysam.idxstats(original)
    if through != direct:
        reasons.append("idxstats through the link differs from idxstats on the original")
    return not reasons, reasons, through


def summary(idx):
    rows = [l.split("\t") for l in idx.splitlines()]
    return (len(rows), sum(int(r[2]) for r in rows), sum(int(r[3]) for r in rows))


def main():
    samples = locate()
    os.makedirs(DEID_DIR, mode=0o700, exist_ok=True)
    all_ok = True
    for label, s in samples.items():
        if not s["index"]:
            die(f"{label}: no index beside the BAM -- report and stop (never re-index)")
        b = link(s["bam"], deid_bam(label))
        i = link(s["index"], deid_bam(label) + ".bai")
        ok, reasons, through = verify(label, deid_bam(label), s["bam"])
        n, mapped, unmapped = summary(through)
        print(f"{label}: bam link {b}, index link {i}; resolves to {label}'s Phase 0 size; "
              f"idxstats through the link: {n} rows, {mapped:,} mapped, {unmapped:,} unmapped; "
              f"identical to the original: {'yes' if ok else 'NO'}"
              + ("" if ok else f" -- {'; '.join(reasons)}"))
        all_ok &= ok
    if "--control" in sys.argv:
        with tempfile.TemporaryDirectory() as tmp:
            wrong = os.path.join(tmp, "SAMPLE_A.bam")
            os.symlink(samples["SAMPLE_B"]["bam"], wrong)
            os.symlink(samples["SAMPLE_B"]["index"], wrong + ".bai")
            ok, reasons, _ = verify("SAMPLE_A", wrong, samples["SAMPLE_A"]["bam"])
            print(f"control: a SAMPLE_A link pointing at SAMPLE_B's BAM is "
                  f"{'ACCEPTED -- the check is broken' if ok else 'rejected'}"
                  + ("" if ok else f" ({len(reasons)} reasons: {'; '.join(reasons)})"))
            all_ok &= not ok
    print("deid links VERIFIED" if all_ok else "deid links NOT verified")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
