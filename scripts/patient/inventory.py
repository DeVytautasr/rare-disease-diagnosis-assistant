#!/usr/bin/env python3
"""Name-free inventory of ~/patient_data and the size-based SAMPLE_A/SAMPLE_B
assignment (Task 1 of the 2026-09-25 recovery prompt).

Prints counts by kind, byte sizes, which Phase 0 size each BAM matches, each
BAM's index and its form, and, for hidden entries such as rsync temporaries left
by an interrupted transfer, their sizes and which sample they belong to. Never
a name. Read-only: nothing is created, moved or deleted.
"""
import os
import stat
import time

from common import PATIENT_DIR, PHASE0, index_for, top_level_entries

KNOWN_DIRS = {"deid", "rerun_2026-09"}
KNOWN_FILES = {"SAMPLE_MAP.md", ".redact_terms", ".identifier_list", "redact.sh", "redact.py"}
KNOWN_EXT = {"md", "txt", "json", "md5", "tsv", "csv", "pdf", "vcf", "gz", "bcf", "csi",
             "tbi", "crai", "cram", "bam", "bai", "log", "sh", "py"}


def main():
    entries = top_level_entries()
    print(f"top-level entries: {len(entries)}")
    bams = {n: st for n, st in entries
            if not n.startswith(".") and n.endswith(".bam") and stat.S_ISREG(st.st_mode)}
    by_size = {rec["bytes"]: label for label, rec in PHASE0.items()}
    print(f"regular .bam files: {len(bams)}")
    label_of = {}
    for n, st in sorted(bams.items(), key=lambda kv: kv[1].st_size):
        label = by_size.get(st.st_size)
        label_of[n] = label
        line = f"  {st.st_size:>16,} bytes -> " + (f"{label} (equals the Phase 0 size)" if label
                                                  else "matches NO Phase 0 size")
        idx, form = index_for(os.path.join(PATIENT_DIR, n))
        if idx:
            ist = os.stat(idx)
            age = ist.st_mtime - st.st_mtime
            line += (f"; index {form}, {ist.st_size:,} bytes, modified "
                     f"{'after' if age >= 0 else 'BEFORE'} the BAM ({age:+.0f} s)")
        else:
            line += "; NO index"
        print(line)
    for label in PHASE0:
        if label not in label_of.values():
            print(f"  {label}: no .bam of its Phase 0 size ({PHASE0[label]['bytes']:,})")
    bai = [n for n, st in entries if n.endswith(".bai") and not n.startswith(".")]
    print(f"regular index files (.bai): {len(bai)}")

    hidden = [(n, st) for n, st in entries if n.startswith(".") and n not in KNOWN_FILES]
    print(f"hidden entries (excluding this workflow's own files): {len(hidden)}")
    for n, st in sorted(hidden, key=lambda kv: kv[1].st_size):
        kind = ("dir" if stat.S_ISDIR(st.st_mode) else "symlink" if stat.S_ISLNK(st.st_mode)
                else "file")
        owner = next((label_of[b] or "an unmatched BAM" for b in bams if n.startswith("." + b)), None)
        if owner is None:
            owner = next((f"the index of {label_of[b]}" for b in bams
                          if n.startswith("." + b[:-4] + ".bai") or n.startswith("." + b + ".bai")), None)
        print(f"  {kind}, {st.st_size:,} bytes, modified {time.strftime('%Y-%m-%d %H:%M', time.localtime(st.st_mtime))}"
              + (f", name has the rsync-temporary shape for {owner}" if owner else ", name matches no BAM"))

    dirs = [n for n, st in entries if stat.S_ISDIR(st.st_mode) and not n.startswith(".")]
    print(f"directories: {len(dirs)}" + (f" ({', '.join(sorted(d for d in dirs if d in KNOWN_DIRS))} "
                                          f"+ {sum(d not in KNOWN_DIRS for d in dirs)} other, names withheld)"
                                          if dirs else ""))
    others = [n for n, st in entries if not n.startswith(".") and not stat.S_ISDIR(st.st_mode)
              and not n.endswith(".bam") and not n.endswith(".bai")]
    ext = {}
    for n in others:
        e = n.rsplit(".", 1)[-1].lower() if "." in n else ""
        key = n if n in KNOWN_FILES else ("." + e if e in KNOWN_EXT else "<other extension>")
        ext[key] = ext.get(key, 0) + 1
    print(f"other entries: {len(others)} {dict(sorted(ext.items()))}")
    print(f"SAMPLE_MAP.md present: {os.path.isfile(os.path.join(PATIENT_DIR, 'SAMPLE_MAP.md'))}")


if __name__ == "__main__":
    main()
