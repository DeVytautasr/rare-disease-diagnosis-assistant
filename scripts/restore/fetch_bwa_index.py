#!/usr/bin/env python3
"""Fetch the 1000 Genomes prebuilt BWA index for hs38DH, WITH its .alt file,
into ~/reference/, and verify every file against the md5 that 1000 Genomes
publishes in its current.tree listing. Downloaded, never built.

    fetch_bwa_index.py [--dest DIR] [--ebi-conns N] [--aws-conns N]

The three large files are fetched as parallel byte ranges from both mirrors: on
2026-09-25 one connection to EBI measured 1.6 MB/s, eight together 5.0 MB/s,
and four more to AWS added 1.1 MB/s. Ranges are appended to per-segment files,
so an interrupted run resumes instead of restarting.

A file gets its final name only after the md5 of the assembled file equals the
published md5. Nothing is trusted from the transfer itself: sizes are checked
against the servers' Content-Length, and the md5 check covers the whole file.

Exit 0 only if all six files are present under their final names with
verified md5s. The published md5 table is kept beside the index
(hs38DH_published_md5.tsv) so later runs and later sessions can re-check it.
"""
import argparse
import hashlib
import os
import queue
import re
import sys
import threading
import time

import requests

EBI = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
AWS = "https://1000genomes.s3.amazonaws.com"
DIR = "technical/reference/GRCh38_reference_genome"
FA = "GRCh38_full_analysis_set_plus_decoy_hla.fa"
EXTS = ["alt", "amb", "ann", "pac", "sa", "bwt"]
# Observed with HEAD requests to both mirrors on 2026-09-25; a mismatch stops the run.
EXPECTED_BYTES = {"amb": 20199, "ann": 448319, "bwt": 3217347004,
                  "pac": 804336731, "sa": 1608673512, "alt": 487553}
# Published md5s already verified in earlier sessions, used as a positive
# control that the current.tree parse finds the right rows.
KNOWN_MD5 = {FA: "64b32de2fc934679c16e83a2bc072064",
             FA + ".fai": "5ccc91e56dc4a05448dd5b9507ec6bc6",
             FA + ".alt": "b07e65aa4425bc365141756f5c98328c"}
SEG = 64 * 1024 * 1024
CHUNK = 1 << 20

_lock = threading.Lock()
_done_bytes = 0


def log(msg):
    print(time.strftime("%Y-%m-%dT%H:%M:%S%z"), msg, flush=True)


def published_md5s(dest):
    """{filename: md5} for the reference directory, parsed from current.tree."""
    table = os.path.join(dest, "hs38DH_published_md5.tsv")
    url = f"{EBI}/current.tree"
    log(f"published md5s: streaming {url}")
    rows = {}
    with requests.get(url, stream=True, timeout=(30, 300)) as r:
        r.raise_for_status()
        last_modified = r.headers.get("Last-Modified", "?")
        # current.tree is served as text/plain with no charset, so requests
        # would hand back bytes; decode each line explicitly.
        for raw in r.iter_lines():
            line = raw.decode("utf-8", "replace")
            if line and f"{DIR}/{FA}" in line:
                path = line.split("\t")[0]
                md5 = next((f for f in line.split("\t") if re.fullmatch(r"[0-9a-f]{32}", f)), None)
                if md5:
                    rows[os.path.basename(path)] = md5
    for name, want in KNOWN_MD5.items():
        if rows.get(name) != want:
            sys.exit(f"current.tree parse failed its control: {name} -> {rows.get(name)}, "
                     f"expected {want}")
    with open(table, "w") as f:
        f.write(f"# source: {url} (Last-Modified: {last_modified}), fetched "
                f"{time.strftime('%Y-%m-%d')}\n")
        for name in sorted(rows):
            f.write(f"{name}\t{rows[name]}\n")
    log(f"published md5s: {len(rows)} rows for {FA}*; the 3 previously verified md5s "
        f"match; table -> {table}")
    return rows


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(8 * CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def content_length(url):
    r = requests.head(url, timeout=30, allow_redirects=True)
    r.raise_for_status()
    return int(r.headers["Content-Length"])


def fetch_segment(mirror, ext, idx, start, end, segpath):
    """Append bytes [start+have, end] to segpath. Returns when complete; raises on error."""
    global _done_bytes
    need = end - start + 1
    have = os.path.getsize(segpath) if os.path.exists(segpath) else 0
    if have > need:
        os.remove(segpath)
        have = 0
    if have == need:
        return
    url = f"{mirror}/{DIR}/{FA}.{ext}"
    hdr = {"Range": f"bytes={start + have}-{end}"}
    with requests.get(url, headers=hdr, stream=True, timeout=(30, 120)) as r:
        if r.status_code != 206:
            raise RuntimeError(f"{ext} seg {idx}: HTTP {r.status_code} for a range request")
        with open(segpath, "ab") as f:
            for block in r.iter_content(CHUNK):
                f.write(block)
                with _lock:
                    _done_bytes += len(block)
    got = os.path.getsize(segpath)
    if got != need:
        raise RuntimeError(f"{ext} seg {idx}: {got} of {need} bytes after the response ended")


def worker(mirror, q, errors):
    while True:
        try:
            item = q.get_nowait()
        except queue.Empty:
            return
        ext, idx, start, end, segpath, tries = item
        try:
            fetch_segment(mirror, ext, idx, start, end, segpath)
        except Exception as e:  # requeue and let any worker (either mirror) resume it
            errors.append(f"{ext} seg {idx} via {mirror.split('/')[2]}: {type(e).__name__}: {e}")
            if tries >= 30:
                log(f"GIVING UP on {ext} seg {idx} after {tries} attempts")
                continue
            time.sleep(min(60, 5 * (tries + 1)))
            q.put((ext, idx, start, end, segpath, tries + 1))
        finally:
            q.task_done()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", default=os.path.expanduser("~/reference"))
    ap.add_argument("--ebi-conns", type=int, default=10)
    ap.add_argument("--aws-conns", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.dest, exist_ok=True)
    segdir = os.path.join(a.dest, ".bwa_index_segments")
    os.makedirs(segdir, exist_ok=True)

    md5s = published_md5s(a.dest)
    want = {ext: md5s.get(f"{FA}.{ext}") for ext in EXTS}
    missing = [e for e, m in want.items() if not m]
    if missing:
        sys.exit(f"no published md5 for: {missing}")

    todo, plan = [], {}
    for ext in EXTS:
        final = os.path.join(a.dest, f"{FA}.{ext}")
        if os.path.exists(final) and md5_of(final) == want[ext]:
            log(f".{ext}: already present with the published md5 -- skipped")
            continue
        size = content_length(f"{EBI}/{DIR}/{FA}.{ext}")
        size_aws = content_length(f"{AWS}/{DIR}/{FA}.{ext}")
        if size != EXPECTED_BYTES[ext] or size_aws != size:
            sys.exit(f".{ext}: Content-Length EBI {size} / AWS {size_aws}, "
                     f"expected {EXPECTED_BYTES[ext]} -- stopping")
        segs = []
        for i, start in enumerate(range(0, size, SEG)):
            end = min(size, start + SEG) - 1
            segs.append((ext, i, start, end, os.path.join(segdir, f"{ext}.{i:04d}")))
        plan[ext] = (final, size, segs)
        todo.extend(segs)
    total = sum(p[1] for p in plan.values())
    already = sum(min(os.path.getsize(s[4]), s[3] - s[2] + 1) for s in todo if os.path.exists(s[4]))
    log(f"to fetch: {len(plan)} files, {total:,} bytes in {len(todo)} segments "
        f"({already:,} bytes already on disk); {a.ebi_conns} EBI + {a.aws_conns} AWS connections")

    q = queue.Queue()
    # largest files first so the tail is short
    for s in sorted(todo, key=lambda s: (-plan[s[0]][1], s[1])):
        q.put(s + (0,))
    errors = []
    threads = [threading.Thread(target=worker, args=(EBI, q, errors), daemon=True)
               for _ in range(a.ebi_conns)]
    threads += [threading.Thread(target=worker, args=(AWS, q, errors), daemon=True)
                for _ in range(a.aws_conns)]
    t0 = time.monotonic()
    for t in threads:
        t.start()
    last_err = 0
    while any(t.is_alive() for t in threads):
        time.sleep(30)
        with _lock:
            done = _done_bytes
        el = time.monotonic() - t0
        rate = done / el if el else 0
        left = total - already - done
        log(f"progress: {already + done:,} / {total:,} bytes | {rate / 1e6:.2f} MB/s | "
            f"ETA {left / rate / 60 if rate else float('inf'):.1f} min | "
            f"{len(errors)} transient errors")
        for e in errors[last_err:]:
            log(f"  transient: {e}")
        last_err = len(errors)

    ok = True
    for ext, (final, size, segs) in plan.items():
        short = [s for s in segs
                 if not os.path.exists(s[4]) or os.path.getsize(s[4]) != s[3] - s[2] + 1]
        if short:
            log(f".{ext}: FAILED -- {len(short)} segment(s) incomplete")
            ok = False
            continue
        part = final + ".part"
        h = hashlib.md5()
        with open(part, "wb") as out:
            for s in segs:
                with open(s[4], "rb") as f:
                    for block in iter(lambda: f.read(8 * CHUNK), b""):
                        h.update(block)
                        out.write(block)
        got, nbytes = h.hexdigest(), os.path.getsize(part)
        if got == want[ext] and nbytes == size:
            os.replace(part, final)
            for s in segs:
                os.remove(s[4])
            log(f".{ext}: OK md5 {got} matches the published md5 | {nbytes:,} bytes")
        else:
            log(f".{ext}: FAILED md5 {got} (published {want[ext]}), {nbytes:,} bytes "
                f"(expected {size:,}); kept as {part}")
            ok = False
    if ok:
        try:
            os.rmdir(segdir)
        except OSError:
            pass
    # Final state, re-derived from disk rather than from the loop above.
    final_ok = all(os.path.exists(os.path.join(a.dest, f"{FA}.{e}")) for e in EXTS)
    log(f"=== index fetch {'COMPLETE' if ok and final_ok else 'INCOMPLETE'}: "
        f"{sum(os.path.exists(os.path.join(a.dest, f'{FA}.{e}')) for e in EXTS)}/6 files "
        f"under final names; wall {time.monotonic() - t0:.0f}s ===")
    sys.exit(0 if ok and final_ok else 1)


if __name__ == "__main__":
    main()
