#!/usr/bin/env python3
"""
Regression tests: the install checks must be able to fail.

Named for the condition that exposed them. Phase 11's `--check` hard-coded its
MINIMAL verdict (ok=True, ANDed with a count compared against >= 0), and the
FULL-tier indicator only tested that a file existed at an igv.sh path: an
empty, non-executable file with no Java anywhere showed [x]. Neither could
report a failure, so neither verified anything. The banner also resolved a
config-file IGV path that the panel tool itself never read, so it could claim an
IGV the tool could not find.

Every check below is driven both ways -- break the condition and require FAIL,
restore it and require PASS. A check that has only ever been seen to pass has
not been tested.

Run: python3 stage1_igv_assistant/tests/test_checks_can_fail.py
"""
import contextlib
import os
import subprocess
import sys
import tempfile
import zipfile

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, REPO)

from stage1_igv_assistant import ui  # noqa: E402
from stage1_igv_assistant import server as ev  # noqa: E402
from stage1_igv_assistant import candidate_server as cs  # noqa: E402
from stage1_igv_assistant.score_tiers import derive_tiers, derive_bands  # noqa: E402

FAILURES = []
LAYERS = "four evidence layers at a hand-entered coordinate"
CHAIN = "candidate set through load and the filter chain"
TIERS = "scoring tiers derivable from bam_tools' source"
FILES = "every registered file exists (each BAM indexed)"


def check(label, condition, detail=""):
    if condition:
        print(f"  PASSED ✓  {label}")
    else:
        print(f"  FAILED ✗  {label}\n             {detail}")
        FAILURES.append(label)


def minimal():
    return {name: ok for name, ok, _ in ui.verify_minimal()}


@contextlib.contextmanager
def patched(obj, attr, value):
    real = getattr(obj, attr)
    setattr(obj, attr, value)
    try:
        yield
    finally:
        setattr(obj, attr, real)


@contextlib.contextmanager
def environ(**kv):
    saved = {k: os.environ.get(k) for k in kv}
    for k, v in kv.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v
    try:
        yield
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def broken(*a, **k):
    return {"error": "control: deliberately broken", "error_type": "control"}


def write(path, text, mode):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)
    os.chmod(path, mode)


def main():
    ui.TIERS, ui.BANDS = derive_tiers(), derive_bands()

    print("MINIMAL self-test (what --check and the banner report)")
    base = minimal()
    check(f"every condition passes on a working install ({len(base)})", all(base.values()), str(base))
    with patched(ev, "count_soft_clipped_reads", broken):
        got = minimal()
    check("a broken evidence layer fails the layers condition, and only that one",
          got[LAYERS] is False and sum(not v for v in got.values()) == 1, str(got))
    with patched(cs, "_list_candidates", broken):
        got = minimal()
    check("a broken filter chain fails the candidate-set condition, and only that one",
          got[CHAIN] is False and sum(not v for v in got.values()) == 1, str(got))
    with patched(ui, "TIERS", None):
        got = minimal()
    check("underivable scoring tiers fail the tier condition", got[TIERS] is False, str(got))
    with tempfile.TemporaryDirectory() as d:
        ui.DATASETS["GHOST"] = os.path.join(d, "ghost.bam")
        try:
            got = minimal()
        finally:
            ui.DATASETS.pop("GHOST")
        check("a registered BAM that does not exist fails the files condition", got[FILES] is False)
        unindexed = os.path.join(d, "unindexed.bam")
        ui._selftest_fixture(d)
        os.rename(os.path.join(d, "selftest.bam"), unindexed)
        ui.DATASETS["NOINDEX"] = unindexed
        try:
            got = minimal()
        finally:
            ui.DATASETS.pop("NOINDEX")
        check("a registered BAM without an index fails the files condition", got[FILES] is False)
    got = minimal()
    check("restored: every condition passes again", all(got.values()), str(got))
    check("the self-test's calls stay out of the session log", not any(
        "__selftest__" in str(c.get("params")) for c in ui.RECORDER.calls))

    print("\nFULL-tier indicator (IGV panels)")
    with tempfile.TemporaryDirectory() as d:
        igv = os.path.join(d, "IGV", "igv.sh")
        nobin = os.path.join(d, "nobin")
        os.makedirs(nobin)
        write(igv, "", 0o644)                        # the original defect: empty, not executable
        with environ(IGV_PATH=igv, PATH=nobin, DISPLAY=None):
            check("empty non-executable igv.sh, no Java, no display -> NOT available",
                  not ui.verify_full()["ok"])
        write(igv, "#!/bin/sh\nexit 3\n", 0o755)
        jar = os.path.join(d, "IGV", "lib", "igv.jar")
        os.makedirs(os.path.dirname(jar))
        with zipfile.ZipFile(jar, "w") as z:         # class-file major 65 = compiled for Java 21
            z.writestr("org/broad/igv/ui/Main.class", b"\xca\xfe\xba\xbe\x00\x00\x00\x41" + b"\x00" * 16)
        bindir = os.path.join(d, "bin")
        java = os.path.join(bindir, "java")

        def set_java(version):
            write(java, f"#!/bin/sh\necho 'openjdk version \"{version}\"' >&2\n", 0o755)

        set_java("25.0.1")
        full = dict(IGV_PATH=igv, PATH=bindir + os.pathsep + "/usr/bin:/bin", DISPLAY="localhost:10.0")
        with environ(**full):
            check("every prerequisite present -> available", ui.verify_full()["ok"],
                  str(ui.verify_full()["conditions"]))
        os.chmod(igv, 0o644)
        with environ(**full):
            check("igv.sh not executable -> NOT available", not ui.verify_full()["ok"])
        os.chmod(igv, 0o755)
        os.rename(jar, jar + ".hidden")
        with environ(**full):
            check("no lib/igv.jar -> NOT available", not ui.verify_full()["ok"])
        os.rename(jar + ".hidden", jar)
        set_java("17.0.2")
        with environ(**full):
            check("Java older than IGV was compiled for -> NOT available", not ui.verify_full()["ok"])
        set_java("25.0.1")
        with environ(**dict(full, PATH=nobin)):
            check("no Java anywhere -> NOT available", not ui.verify_full()["ok"])
        with environ(**dict(full, DISPLAY=None)):
            check("no DISPLAY -> NOT available", not ui.verify_full()["ok"])
        with environ(**dict(full, DISPLAY=":987654")):
            check("a local DISPLAY without its X socket -> NOT available", not ui.verify_full()["ok"])
        with environ(**full):
            check("restored -> available again", ui.verify_full()["ok"])

        print("\nA config-file IGV path reaches the panel tool, not only the banner")
        log = os.path.join(d, "launched.log")
        write(igv, f"#!/bin/sh\necho \"$0\" >> '{log}'\nexit 3\n", 0o755)
        conf = os.path.join(d, "sv.conf")
        with open(conf, "w") as f:
            f.write(f"[paths]\nigv = {igv}\n")
        child = (
            "import sys\n"
            f"sys.path.insert(0, {os.path.abspath(REPO)!r})\n"
            "from stage1_igv_assistant import ui\n"
            "from stage1_igv_assistant.tools import bam_tools\n"
            "sys.argv = ['ui', '--check', '--no-autodiscover']\n"
            "try:\n    ui.main()\nexcept SystemExit:\n    pass\n"
            "bam_tools.run_igv_screenshot(bam_paths=['x.bam'], chromosome='chr1', start=1, end=2,\n"
            f"                             output_path={os.path.join(d, 'out.png')!r}, timeout_sec=10)\n")
        env = {k: v for k, v in os.environ.items() if k != "IGV_PATH"}
        env.update(SV_CONFIG=conf, DISPLAY="localhost:10.0")
        subprocess.run([sys.executable, "-c", child], env=env, capture_output=True, text=True, timeout=120)
        launched = open(log).read().split() if os.path.exists(log) else []
        check("the panel tool launched the igv.sh named only in the config file",
              launched and all(p == igv for p in launched), str(launched))

    print("\n" + "=" * 60)
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S): {FAILURES}")
        sys.exit(1)
    print("ALL CHECK-FALSIFIABILITY TESTS PASSED")


if __name__ == "__main__":
    main()
