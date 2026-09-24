#!/usr/bin/env python3
"""
identifier_gate.py -- run before pushing to the public remote.

Scans three scopes for what must never be published:
  commits   the messages of every commit in BASE..HEAD
  diff      the lines those commits ADD (removed lines are reported, not gated:
            they were already public)
  tree      every blob in HEAD's tree, read from git objects (the working tree
            is never touched); .docx/.xlsx/.pptx are opened and their XML scanned

Rules -- mechanical, so nothing is judged case by case at gate time:
  ID-list     any entry of --identifiers FILE (one per line; e.g. built from
              ~/patient_data/SAMPLE_MAP.md). Without it the run is INCOMPLETE.
  ID-sample   a sample-file token (NAME.bam/.cram/.bai/.crai/.csi/.vcf/.bcf/.fastq)
              whose NAME has identifier shape -- contains a digit, or ends in
              "-ready" as the real transferred BAMs do -- unless NAME is a public
              reference sample, a synthetic implant, or a documented placeholder
  ID-lt-code  an 11-digit Lithuanian personal code with a valid date and checksum
  ID-email    an email address outside the allowlist below
  KEY         the actual API key (--keyfile), or a known secret format
  PATH        (diff and commit messages only) an absolute path under a personal or
              local root: /home/<name>, /Users/<name>, /root, a Windows user folder
              under /mnt/<drive>/Users, the Claude scratch area under /tmp, anything
              naming patient_data, or this machine's home or repository path.
              Other added absolute paths are listed, not gated.

Hits are printed with location and a MASKED token -- an identifier is never
printed in full, and key material never at all.

Exit codes: 0 clean, with an identifier list
            1 something fired -- do not push
            2 nothing fired, but no identifier list was supplied (INCOMPLETE)

  python3 scripts/identifier_gate.py --base COMMIT [--identifiers FILE] [--keyfile .api/claude_api_key]
  python3 scripts/identifier_gate.py --self-test     # positive and negative controls
"""
import argparse
import gzip
import io
import os
import random
import re
import string
import struct
import subprocess
import sys
import tempfile
import zipfile
import zlib

EXT = r"(?:bam|cram|bai|crai|csi|bcf|vcf(?:\.gz)?|fastq(?:\.gz)?|fq(?:\.gz)?)"
SAMPLE_TOKEN = re.compile(rf"(?<![\w.+-])([A-Za-z0-9][\w.+-]*?)\.{EXT}(?!\w)")
PUBLIC = re.compile(r"^(NA|HG|GM)\d{5}(?!\d)|^HG00[1-7](?!\d)|^HCC\d+|^CHM13|^GRCh3[78]")
SYNTHETIC = re.compile(r"^IMP\d{2}$")
PLACEHOLDER = re.compile(r"^SAMPLE_[A-Z](-ready)?$|^(patient|sample)_\d{3}$")
LT_CODE = re.compile(r"(?<!\d)([1-6])(\d{2})(\d{2})(\d{2})(\d{3})(\d)(?!\d)")
EMAIL = re.compile(r"(?<![\w.+-])[\w.+-]+@[\w-]+(?:\.[\w-]+)+")
# The thesis author's own addresses (git author; TUTORIAL.md byline), the commit
# co-author trailer, and reserved example domains. Anything else fires.
EMAIL_OK = re.compile(r"^(noreply@anthropic\.com|vytis\.official@gmail\.com|vytautas\.rimas@mf\.stud\.vu\.lt|[\w.+-]+@example\.(com|org|invalid))$", re.I)
SECRETS = re.compile(r"sk-ant-[A-Za-z0-9_\-]{20,}|AKIA[0-9A-Z]{16}|gh[pousr]_[A-Za-z0-9]{36,}"
                     r"|github_pat_[A-Za-z0-9_]{20,}|-----BEGIN [A-Z ]*PRIVATE KEY-----|xox[baprs]-[A-Za-z0-9-]{10,}")
ABS_PATH = re.compile(r"(?<![\w.~:/@-])/[\w.+@%~-]+(?:/[\w.+@%~.*-]*)+")
PERSONAL = re.compile(r"^/home/(?!\.\.\.(/|$))[^/]+|^/Users/[^/]+|^/root(/|$)|^/mnt/[a-z]/Users/[^/]+"
                      r"|^/tmp/claude[-]|patient_data")


def mask(tok):
    return tok[:2] + "*" * max(0, len(tok) - 2) + f" (len {len(tok)})"


def lt_code_valid(m):
    digits = [int(c) for c in m.group(0)]
    century = {1: 1800, 2: 1800, 3: 1900, 4: 1900, 5: 2000, 6: 2000}[digits[0]]
    yy, mm, dd = int(m.group(2)), int(m.group(3)), int(m.group(4))
    if not (1 <= mm <= 12 and 1 <= dd <= 31) or century + yy < 1870:
        return False
    s = sum(d * w for d, w in zip(digits[:10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 1])) % 11
    if s == 10:
        s = sum(d * w for d, w in zip(digits[:10], [3, 4, 5, 6, 7, 8, 9, 1, 2, 3])) % 11
        s = 0 if s == 10 else s
    return s == digits[10]


class Gate:
    def __init__(self, identifiers, key, home, repo):
        self.identifiers = [i for i in identifiers if i]
        self.key = key
        self.local = [p for p in (home, repo) if p and len(p) > 1]
        self.fired, self.listed = [], []

    def scan(self, scope, where, text, paths=False):
        low = text.lower()
        for ident in self.identifiers:
            if ident.lower() in low:
                self.fired.append(("ID-list", scope, where, mask(ident)))
        for m in SAMPLE_TOKEN.finditer(text):
            name = m.group(1)
            shaped = bool(re.search(r"\d", name)) or name.lower().endswith("-ready")
            if shaped and not (PUBLIC.search(name) or SYNTHETIC.search(name) or PLACEHOLDER.search(name)):
                self.fired.append(("ID-sample", scope, where, mask(name)))
        for m in LT_CODE.finditer(text):
            if lt_code_valid(m):
                self.fired.append(("ID-lt-code", scope, where, mask(m.group(0))))
        for m in EMAIL.finditer(text):
            if not EMAIL_OK.match(m.group(0)):
                self.fired.append(("ID-email", scope, where, mask(m.group(0))))
        if (self.key and self.key in text) or SECRETS.search(text):
            self.fired.append(("KEY", scope, where, "(key material -- not shown)"))
        if paths:
            for m in ABS_PATH.finditer(text):
                p = m.group(0)
                if PERSONAL.search(p) or any(p.startswith(l) for l in self.local):
                    self.fired.append(("PATH", scope, where, re.sub(r"^(/[^/]+/)([^/]+)", r"\1<masked>", p)))
                else:
                    self.listed.append((scope, where, p))


def git(repo, *args):
    out = subprocess.run(["git", "-C", repo, *args], capture_output=True, check=True).stdout
    return out.decode("utf-8", "replace")


UNSCANNED = []          # binary blobs with no text container this gate can read


def blob_texts(path, data):
    """Text the blob can carry. Compressed or binary bytes are never scanned as if
    they were text -- a PNG's image data once 'contained' an email address by
    chance. OOXML (.docx/.xlsx/.pptx) -> its XML; gzip -> decompressed; PNG -> its
    tEXt/zTXt/iTXt chunks; any other binary is counted as not scanned."""
    if data[:2] == b"PK":
        try:
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                for n in z.namelist():
                    if n.endswith((".xml", ".rels", ".txt")):
                        yield f"{path}!{n}", z.read(n).decode("utf-8", "replace")
            return
        except zipfile.BadZipFile:
            pass
    if data[:2] == b"\x1f\x8b":
        yield f"{path}!gunzip", gzip.decompress(data).decode("utf-8", "replace")
        return
    if data[:8] == b"\x89PNG\r\n\x1a\n":
        i = 8
        while i + 8 <= len(data):
            n, kind = struct.unpack(">I4s", data[i:i + 8])
            body = data[i + 8:i + 8 + n]
            if kind == b"tEXt":
                yield f"{path}!tEXt", body.decode("latin-1")
            elif kind == b"zTXt":
                k, _, rest = body.partition(b"\0")
                yield f"{path}!zTXt", zlib.decompress(rest[1:]).decode("latin-1")
            elif kind == b"iTXt":
                k, _, rest = body.partition(b"\0")
                flag, rest = rest[0], rest[2:]
                _lang, _, rest = rest.partition(b"\0")
                _tk, _, text = rest.partition(b"\0")
                yield f"{path}!iTXt", (zlib.decompress(text) if flag else text).decode("utf-8", "replace")
            i += 12 + n
        return
    if b"\0" in data[:8192]:
        UNSCANNED.append(path)
        return
    yield path, data.decode("latin-1")


def run(repo, base, identifiers, key, quiet=False):
    home = os.path.expanduser("~")
    g = Gate(identifiers, key, home, os.path.abspath(repo))
    commits = git(repo, "rev-list", "--reverse", f"{base}..HEAD").split()
    for c in commits:
        g.scan("commits", f"commit {c[:10]} message", git(repo, "log", "-1", "--format=%B", c), paths=True)
    removed = []
    # Binary files (a .docx is a zip) cannot be diffed as lines: each one changed in
    # the range is scanned whole, as its new blob, the same way the tree scan does.
    for rec in git(repo, "diff", "--numstat", "-z", f"{base}..HEAD").split("\0"):
        parts = rec.split("\t")
        if len(parts) == 3 and parts[0] == "-" and parts[1] == "-" and parts[2]:
            try:
                data = subprocess.run(["git", "-C", repo, "show", f"HEAD:{parts[2]}"],
                                      capture_output=True, check=True).stdout
            except subprocess.CalledProcessError:
                continue                           # deleted in the range: nothing added
            for where, text in blob_texts(parts[2], data):
                for i, ln in enumerate(text.splitlines(), 1):
                    g.scan("diff", f"{where}:{i} (binary, added)", ln, paths=True)
    diff = git(repo, "diff", "--unified=0", "--no-color", f"{base}..HEAD")
    f, line = None, 0
    for raw in diff.splitlines():
        if raw.startswith("+++ "):
            f = raw[6:] if raw.startswith("+++ b/") else raw[4:]
        elif raw.startswith("@@"):
            line = int(re.search(r"\+(\d+)", raw).group(1))
        elif raw.startswith("+") and not raw.startswith("+++"):
            g.scan("diff", f"{f}:{line} (added)", raw[1:], paths=True)
            line += 1
        elif raw.startswith("-") and not raw.startswith("---"):
            removed.append((f, raw[1:]))
    for rec in git(repo, "ls-tree", "-r", "-z", "HEAD").split("\0"):
        if not rec:
            continue
        meta, path = rec.split("\t", 1)
        if meta.split()[1] != "blob":
            continue
        data = subprocess.run(["git", "-C", repo, "cat-file", "blob", meta.split()[2]],
                              capture_output=True, check=True).stdout
        for where, text in blob_texts(path, data):
            for i, ln in enumerate(text.splitlines(), 1):
                g.scan("tree", f"{where}:{i}", ln)
    pre = Gate(identifiers, key, home, os.path.abspath(repo))
    for fpath, text in removed:
        pre.scan("diff", f"{fpath} (removed line)", text)
    if not quiet:
        print(f"identifier gate: {len(commits)} commit(s) in {base[:10]}..HEAD; "
              f"identifier list: {'%d entr%s' % (len(identifiers), 'y' if len(identifiers) == 1 else 'ies') if identifiers else 'NOT SUPPLIED'}")
        for rule, scope, where, tok in g.fired:
            print(f"  FIRED  {rule:10s} [{scope}] {where}  {tok}")
        for rule, scope, where, tok in pre.fired:
            print(f"  (pre-existing, removed by these commits) {rule} {where} {tok}")
        print(f"  binary blobs with no readable text container, NOT scanned: {len(set(UNSCANNED))} "
              f"{sorted(set(UNSCANNED))[:5]}")
        print(f"  added absolute paths not under a personal root ({len(g.listed)}), listed not gated:")
        for p in sorted({p for _, _, p in g.listed}):
            print(f"    {p}")
        verdict = ("FIRED -- do not push" if g.fired else
                   "clean" if identifiers else
                   "INCOMPLETE -- nothing fired, but no identifier list was supplied; "
                   "identifiers that match no pattern rule cannot be excluded")
        print(f"identifier gate: {verdict}")
    return g, (1 if g.fired else (0 if identifiers else 2))


def self_test():
    """Positive controls CREATE each condition; the negative control must stay silent."""
    rnd = random.Random()
    # the ZQ prefix keeps a random draw from ever matching a public-sample rule (HG, NA, HCC, ...)
    ident = "ZQ" + rnd.choice(string.ascii_uppercase) + str(rnd.randint(1000, 9999)) + "-ready"
    ident2 = "ZQ" + rnd.choice(string.ascii_uppercase) + str(rnd.randint(10000, 99999))
    while True:                                   # a synthetic but checksum-valid personal code
        body = f"{rnd.choice('3456')}{rnd.randint(0, 99):02d}{rnd.randint(1, 12):02d}{rnd.randint(1, 28):02d}{rnd.randint(0, 999):03d}"
        for last in range(10):
            m = LT_CODE.search(body + str(last))
            if m and lt_code_valid(m):
                break
        else:
            continue
        code = body + str(last)
        break
    key = "sk-ant-" + "".join(rnd.choices(string.ascii_letters + string.digits, k=40))
    home_path = "/home/" + "".join(rnd.choices(string.ascii_lowercase, k=8)) + "/data/run.log"

    ident3 = "ZQ" + rnd.choice(string.ascii_uppercase) + str(rnd.randint(100, 999))
    ident4 = "ZQ" + rnd.choice(string.ascii_uppercase) + str(rnd.randint(100000, 999999))

    def png(comment):
        def chunk(kind, body):
            return struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body))
        ihdr = struct.pack(">IIBBBBB", 1, 1, 8, 0, 0, 0, 0)
        return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"tEXt", b"Comment\0" + comment.encode())
                + chunk(b"IDAT", zlib.compress(b"\0\0")) + chunk(b"IEND", b""))

    def docx(xml_text):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as z:
            z.writestr("[Content_Types].xml", '<?xml version="1.0"?><Types/>')
            z.writestr("word/document.xml", f'<?xml version="1.0"?><w:document><w:body><w:p><w:r>'
                                            f'<w:t>{xml_text}</w:t></w:r></w:p></w:body></w:document>')
        return buf.getvalue()

    def repo_with(planted):
        d = tempfile.mkdtemp(prefix="idgate_")
        env = dict(os.environ, GIT_AUTHOR_NAME="c", GIT_AUTHOR_EMAIL="c@example.invalid",
                   GIT_COMMITTER_NAME="c", GIT_COMMITTER_EMAIL="c@example.invalid")
        run_git = lambda *a: subprocess.run(["git", "-C", d, *a], capture_output=True, check=True, env=env)
        run_git("init", "-q")
        open(os.path.join(d, "base.txt"), "w").write("base\n")
        run_git("add", "-A"); run_git("commit", "-qm", "base")
        base = subprocess.run(["git", "-C", d, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        text = (f"the sample {ident}.bam was realigned\n" if planted else "the sample SAMPLE_A-ready.bam was realigned\n")
        text += (f"log at {home_path}\n" if planted else "log at /home/.../data/run.log\n")
        text += (f"token {key}\n" if planted else "token (none)\n")
        open(os.path.join(d, "notes.md"), "w").write(text)
        open(os.path.join(d, "chapter.docx"), "wb").write(
            docx(f"cohort {ident2}.cram, code {code}" if planted else "cohort NA12878.cram, code none"))
        open(os.path.join(d, "figure.png"), "wb").write(
            png(f"source {ident3}.bam" if planted else "source HG002.GRCh38.300x.bam"))
        open(os.path.join(d, "calls.vcf.gz"), "wb").write(
            gzip.compress((f"##sample={ident4}.bam\n" if planted else "##sample=NA12878.bam\n").encode()))
        run_git("add", "-A"); run_git("commit", "-qm", "add notes" if not planted else f"notes for {ident}.bam")
        return d, base

    ok = True
    d, base = repo_with(planted=True)
    g, rc = run(d, base, [], None, quiet=True)
    want = {("ID-sample", "diff"), ("ID-sample", "tree"), ("ID-sample", "commits"), ("PATH", "diff"),
            ("KEY", "diff"), ("KEY", "tree"), ("ID-lt-code", "tree")}
    got = {(r, s) for r, s, _, _ in g.fired}
    docx_hits = {r for r, s, w, _ in g.fired if s == "tree" and "chapter.docx!word/document.xml" in w}
    print(f"POSITIVE  planted in text, commit message and .docx XML -> exit {rc}")
    for r, s in sorted(want):
        print(f"  {'fired ' if (r, s) in got else 'MISSED'}  {r:10s} in {s}")
        ok &= (r, s) in got
    print(f"  {'fired ' if {'ID-sample', 'ID-lt-code'} <= docx_hits else 'MISSED'}  both identifiers found INSIDE the .docx XML")
    ok &= {"ID-sample", "ID-lt-code"} <= docx_hits and rc == 1
    png_hit = any(s_ == "tree" and "figure.png!tEXt" in w for _, s_, w, _ in g.fired)
    gz_hit = any(s_ == "tree" and "calls.vcf.gz!gunzip" in w for _, s_, w, _ in g.fired)
    print(f"  {'fired ' if png_hit else 'MISSED'}  identifier found in a PNG tEXt chunk")
    print(f"  {'fired ' if gz_hit else 'MISSED'}  identifier found inside a gzip member")
    ok &= png_hit and gz_hit
    d2, base2 = repo_with(planted=False)
    g2, rc2 = run(d2, base2, [], None, quiet=True)
    print(f"NEGATIVE  same files, planted items replaced by placeholders -> exit {rc2}, fired: {len(g2.fired)}")
    ok &= not g2.fired and rc2 == 2
    g3, rc3 = run(d2, base2, ["SOMETHING-NOT-PRESENT"], None, quiet=True)
    print(f"NEGATIVE  with an identifier list that matches nothing -> exit {rc3} (0 = clean)")
    ok &= rc3 == 0
    g4, rc4 = run(d, base, [ident2], None, quiet=True)
    print(f"POSITIVE  identifier list entry planted in .docx -> "
          f"{'fired' if any(r == 'ID-list' for r, *_ in g4.fired) else 'MISSED'} (exit {rc4})")
    ok &= any(r == "ID-list" for r, *_ in g4.fired)
    print(f"self-test: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".")
    ap.add_argument("--base")
    ap.add_argument("--identifiers", help="file, one identifier per line ('#' comments allowed)")
    ap.add_argument("--keyfile", default=".api/claude_api_key")
    ap.add_argument("--self-test", action="store_true")
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if not a.base:
        ap.error("--base is required")
    idents = []
    if a.identifiers:
        idents = [l.strip() for l in open(a.identifiers) if l.strip() and not l.startswith("#")]
    key = None
    if a.keyfile and os.path.exists(os.path.join(a.repo, a.keyfile)):
        key = open(os.path.join(a.repo, a.keyfile)).read().strip() or None
    return run(a.repo, a.base, idents, key)[1]


if __name__ == "__main__":
    sys.exit(main())
