#!/usr/bin/env python3
"""
PreToolUse guard for restricted subagents.

Wired into individual subagents through their frontmatter `hooks:` block, so
it fires only while that subagent is running. Two modes:

  guard.py no-git       patient-data   — no route to git / a commit
  guard.py no-mutation  verifier,      — read and run, but do not change
                        thesis-editor    anything in the working tree

Contract (docs/en/hooks): read a JSON object on stdin, print a
hookSpecificOutput block with permissionDecision "deny" to block the call,
or exit 0 silently to fall through to the normal permission flow.

FAIL CLOSED. Any parse failure, unexpected shape, or internal error denies.
A guard that allows on error is worse than no guard, because it looks like
one. Read LIMITS.md in this directory before trusting any of this: the holes
are real, enumerated, and some of them cannot be closed at this layer.
"""
import json
import os
import re
import sys

MODE = sys.argv[1] if len(sys.argv) > 1 else "no-mutation"

# Paths a restricted agent may still write to: scratch space only, never the
# repository and never the patient data directory.
ALLOWED_WRITE_ROOTS = ("/tmp/", "/var/tmp/", "/dev/null", "/dev/stdout", "/dev/stderr")

# Substring identifying the governed patient directory.
PATIENT_MARKER = "patient_data"

# The ONLY location a restricted agent may write BAM-derived output to.
# It sits *inside* patient_data on purpose: every rule that matches a patient
# path still matches a slice written here, so the slice stays governed. A
# slice written to /tmp would leave the namespace and become ungoverned bytes
# that these rules no longer see -- the laundering path this route exists to
# prevent.
DERIVED_MARKER = "patient_data/derived"

# Used to spot patient-derived files being written into the repository.
# CLAUDE_PROJECT_DIR is set by Claude Code when it runs a hook.
REPO_ROOT = os.path.abspath(os.environ.get("CLAUDE_PROJECT_DIR", "")) if os.environ.get("CLAUDE_PROJECT_DIR") else ""

# Wrappers that take the real command as their argument. Stripped before the
# leading token is identified, so `sudo git commit` resolves to `git`.
WRAPPERS = {
    "env", "sudo", "doas", "nohup", "nice", "ionice", "time", "command",
    "builtin", "exec", "stdbuf", "setsid", "xargs", "watch", "script",
}
# Wrappers whose own first argument is a value, not the command.
WRAPPERS_WITH_ARG = {"timeout", "flock", "chrt", "taskset"}

# git subcommands that only read. Everything else is denied in both modes:
# an allowlist cannot silently miss a new mutating subcommand the way a
# blocklist can.
GIT_READONLY = {
    "log", "show", "diff", "status", "blame", "annotate", "ls-files",
    "ls-tree", "ls-remote", "rev-parse", "rev-list", "cat-file", "describe",
    "shortlog", "count-objects", "check-ignore", "check-attr", "grep",
    "whatchanged", "reflog", "for-each-ref", "show-ref", "symbolic-ref",
    "var", "version", "help", "merge-base", "name-rev", "verify-commit",
    "diff-tree", "diff-index", "difftool", "fsck", "bisect",
}

# Anything here is a general-purpose write primitive regardless of arguments.
BLOCKED_ALWAYS = {
    "gh", "hub", "glab", "dd", "shred", "mkfs", "fdisk", "parted",
    "crontab", "at", "systemctl", "service",
    # Interactive/scriptable editors: `vim -c wq`, `ex -sc wq` write files.
    "vim", "vi", "ex", "nano", "emacs", "ed", "pico", "joe", "micro",
    # Build and package tooling: all of it writes, none of it is needed to
    # check a claim.
    "make", "npm", "yarn", "pnpm", "cargo", "go", "pip", "pip3", "conda",
    "gradle", "mvn", "cmake", "meson", "ninja", "setup.py",
}

# Network and transfer tools. Denied for patient-data: each is both a route
# to a remote commit (the GitHub API over curl) and a route to moving
# sequencing data off this machine.
NETWORK_TOOLS = {
    "curl", "wget", "nc", "ncat", "socat", "ssh", "scp", "sftp", "rsync",
    "ftp", "telnet", "aws", "gcloud", "az", "s3cmd", "rclone",
}

# Mutating file commands: denied unless every path argument is under an
# allowed write root.
PATH_MUTATORS = {
    "rm", "mv", "cp", "ln", "install", "truncate", "mkdir", "rmdir",
    "touch", "chmod", "chown", "chgrp", "patch", "unlink", "tee", "rsync",
}

# Interpreters. Running a script file is allowed (the test suite is a script
# and must run); inline code and heredocs are denied, because `python -c` is
# an unrestricted write primitive.
INTERPRETERS = {
    "python", "python3", "python2", "perl", "ruby", "node", "nodejs",
    "php", "lua", "Rscript", "osascript", "deno", "bun",
}
SHELLS = {"bash", "sh", "zsh", "ksh", "dash", "fish", "csh", "tcsh"}
INLINE_CODE_FLAGS = {"-c", "-e", "-E", "--command", "--eval", "--exec"}

# Archiving/compression. Denied on a patient path in BOTH modes: an archive
# is the classic staging step before data leaves a host, and nothing in this
# project needs one. Denied even into derived/, unlike samtools output.
ARCHIVE_TOOLS = {
    "tar", "zip", "gzip", "bgzip", "bzip2", "xz", "zstd", "7z", "compress",
}

# Mutators whose every patient path argument is itself a write target.
SOURCE_DESTRUCTIVE = {
    "rm", "unlink", "mv", "truncate", "chmod", "chown", "chgrp", "touch",
    "mkdir", "rmdir", "patch",
}
# Mutators whose LAST path argument is the destination; earlier ones are
# sources and may legitimately be the read-only BAMs.
DEST_LAST = {"cp", "install", "ln", "tee"}

# Commands whose whole purpose is to read a file to stdout in a form that
# survives copy-paste. Denied on a patient path in both modes.
#
# DEFENCE IN DEPTH ONLY -- read LIMITS.md hole 2 before trusting this. The
# transcript is the exfiltration channel, and `samtools view` reaches it
# through exactly the same route while having to stay open, because reading
# the BAM is the patient-data agent's entire job. Denying these removes an
# obvious accidental route at zero cost (no workflow here needs them); it
# does not close the class.
CONTENT_DUMPERS = {
    "base64", "base32", "xxd", "od", "hexdump", "strings", "uuencode",
}
# `head`/`tail` are denied only in byte mode (-c). Line mode stays open so
# the transfer log can be read, which is ordinary work.
BYTE_MODE_READERS = {"head", "tail"}


def deny(reason: str) -> None:
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "deny",
            "permissionDecisionReason": reason,
        }
    }))
    sys.exit(0)


def allow() -> None:
    """Emit no decision; the normal permission flow continues."""
    sys.exit(0)


def under_allowed_root(path: str) -> bool:
    p = path.strip().strip("'\"")
    if not p:
        return False
    if not p.startswith("/"):
        return False  # relative path -> resolves against cwd, i.e. the repo
    return any(p == r.rstrip("/") or p.startswith(r) for r in ALLOWED_WRITE_ROOTS)


def _quote_mask(text: str):
    """
    Per-character flags: True where the character sits inside a quoted span.

    Returns (mask, balanced). If quoting is unbalanced the mask is unusable,
    and every caller falls back to the quote-blind behaviour -- inspecting
    too much rather than too little. Fail closed.
    """
    mask = [False] * len(text)
    quote = None
    i = 0
    while i < len(text):
        ch = text[i]
        if quote is None and ch == "\\":
            if i + 1 < len(text):
                mask[i] = mask[i + 1] = False
            i += 2
            continue
        if quote is None and ch in "'\"":
            quote = ch
            mask[i] = True
            i += 1
            continue
        if quote is not None:
            mask[i] = True
            if ch == "\\" and quote == '"' and i + 1 < len(text):
                mask[i + 1] = True
                i += 2
                continue
            if ch == quote:
                quote = None
            i += 1
            continue
        i += 1
    return mask, quote is None


def _split_unquoted(text: str) -> list:
    """Split on ; && || | & and newlines, ignoring separators inside quotes."""
    mask, balanced = _quote_mask(text)
    if not balanced:
        # Unbalanced quoting: fall back to the blunt split so nothing hides.
        return re.split(r"(?:\|\||&&|[;|\n&])", text)
    parts, start, i = [], 0, 0
    while i < len(text):
        if mask[i]:
            i += 1
            continue
        if text[i:i + 2] in ("&&", "||"):
            parts.append(text[start:i])
            i += 2
            start = i
            continue
        if text[i] in ";|&\n":
            parts.append(text[start:i])
            i += 1
            start = i
            continue
        i += 1
    parts.append(text[start:])
    return parts


def split_segments(command: str) -> list:
    """
    Break a command line into individually-executed segments.

    Handles ; && || | newlines, and pulls the bodies out of $(...) and
    `...` substitutions so a command hidden inside one is still inspected.
    This is a best-effort tokenizer, not a shell parser -- see LIMITS.md.
    """
    segments = []
    # Command substitutions, innermost first.
    for m in re.finditer(r"\$\(([^()]*)\)", command):
        segments.extend(split_segments(m.group(1)))
    for m in re.finditer(r"`([^`]*)`", command):
        segments.extend(split_segments(m.group(1)))
    # Replace substitutions with a marker that still carries a '$', so a
    # command name built from one ( `$(which git) commit` ) is still caught
    # by the expansion check in check_segment.
    stripped = re.sub(r"\$\([^()]*\)", " $__SUBST__ ", command)
    stripped = re.sub(r"`[^`]*`", " $__SUBST__ ", stripped)
    segments.extend(_split_unquoted(stripped))
    return [s.strip() for s in segments if s.strip()]


def resolve_tokens(segment: str) -> list:
    """Tokens of a segment with env assignments and wrappers stripped."""
    try:
        import shlex
        tokens = shlex.split(segment, comments=True)
    except ValueError:
        tokens = segment.split()
    while tokens:
        head = tokens[0]
        if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*=.*", head):
            tokens = tokens[1:]
            continue
        base = os.path.basename(head)
        if base in WRAPPERS:
            tokens = tokens[1:]
            continue
        if base in WRAPPERS_WITH_ARG:
            tokens = tokens[2:] if len(tokens) > 1 else []
            continue
        break
    return tokens


def check_redirections(command: str) -> None:
    """Deny > and >> to anywhere outside an allowed write root."""
    # Skip fd duplications (2>&1, >&2) and here-strings.
    mask, balanced = _quote_mask(command)
    for m in re.finditer(r"(?<![0-9<>])>>?\s*([^\s;|&<>]+)", command):
        target = m.group(1)
        if target.startswith("&"):
            continue
        # A '>' inside a quoted argument (grep -n 'a>b') is data, not a
        # redirection. On unbalanced quoting the mask is untrusted, so the
        # check still runs.
        if balanced and mask[m.start()]:
            continue
        if not under_allowed_root(target):
            deny(
                f"Blocked: shell redirection writes to {target!r}, outside the "
                f"scratch roots {ALLOWED_WRITE_ROOTS}. This agent reports; it "
                f"does not modify the working tree."
            )


def check_git(tokens: list) -> None:
    sub = None
    for t in tokens[1:]:
        if t.startswith("-"):
            continue
        # `git -C <dir> commit` — the -C value is consumed as a non-flag, so
        # skip one token after the flags that take a value.
        sub = t
        break
    # Re-walk properly, honouring value-taking global flags.
    i, sub = 1, None
    while i < len(tokens):
        t = tokens[i]
        if t in ("-C", "-c", "--git-dir", "--work-tree", "--namespace", "--exec-path"):
            i += 2
            continue
        if t.startswith("-"):
            i += 1
            continue
        sub = t
        break
    if sub is None:
        allow()
    if sub not in GIT_READONLY:
        deny(
            f"Blocked: `git {sub}` is not in this agent's read-only git "
            f"allowlist. Permitted: {', '.join(sorted(GIT_READONLY))}. "
            f"Committing, staging, or altering refs is not available here."
        )


def check_interpreter(tokens: list, base: str, why: str) -> None:
    """
    Interpreters are write/exec primitives. Running a script *file* is
    allowed (the test suite must run); inline code, stdin scripts, and
    module execution are not.
    """
    if base not in INTERPRETERS | SHELLS:
        return
    args = tokens[1:]
    for t in args:
        if t in INLINE_CODE_FLAGS:
            deny(f"Blocked: `{base} {t}` runs inline code. {why}")
        if t == "-m":
            deny(f"Blocked: `{base} -m` executes a module. {why}")
    non_flag = [t for t in args if not t.startswith("-")]
    if not non_flag or "-" in args:
        deny(
            f"Blocked: `{base}` with no script file reads its program from "
            f"stdin (e.g. a pipeline into a shell). {why}"
        )


def _clean(tok: str) -> str:
    return tok.strip().strip("'\"")


def _is_pathlike(tok: str) -> bool:
    t = _clean(tok)
    if not t or t.startswith("-"):
        return False
    return "/" in t or t.startswith("~") or t.startswith(".")


def _is_patient(tok: str) -> bool:
    return PATIENT_MARKER in _clean(tok)


def _is_derived(tok: str) -> bool:
    t = _clean(tok).rstrip("/")
    return t.endswith(DERIVED_MARKER) or (DERIVED_MARKER + "/") in t + "/"


def _output_value(tokens: list):
    """Value of -o / -O / --output / --output=X, if present."""
    for i, t in enumerate(tokens):
        if t in ("-o", "-O", "--output"):
            return tokens[i + 1] if i + 1 < len(tokens) else None
        if t.startswith("--output="):
            return t.split("=", 1)[1]
    return None


def check_patient_data(tokens: list, base: str):
    """
    Rules protecting ~/patient_data, applied in BOTH modes.

    These used to live inside the `no-git` branch, which meant a `verifier`
    or `thesis-editor` run could `tar` and `curl` the very BAMs the
    `patient-data` agent exists to protect (LIMITS.md, hole 3 -- found by
    probing, not by review). The protection belongs to the DATA, not to the
    agent, so it runs before any mode-specific logic.

    Returns None to fall through to the mode rules.
    """
    involved = [t for t in tokens[1:] if _is_patient(t)]
    if not involved:
        return None

    if base in NETWORK_TOOLS:
        deny(
            f"Blocked: `{base}` with a patient_data path in its arguments. "
            f"Patient sequencing data does not leave this host, from any "
            f"agent. This rule applies in every guard mode."
        )

    if base in ARCHIVE_TOOLS:
        deny(
            f"Blocked: `{base}` archiving a patient_data path. Archiving is "
            f"the staging step before data leaves a host, and no workflow "
            f"here needs it. This rule applies in every guard mode."
        )

    if base in CONTENT_DUMPERS:
        deny(
            f"Blocked: `{base}` on a patient_data path dumps file content to "
            f"stdout, and stdout here is the conversation transcript. Note "
            f"this is defence in depth, not a boundary: `samtools view` "
            f"reaches the same channel and must stay open. See LIMITS.md "
            f"hole 2."
        )

    if base in BYTE_MODE_READERS and any(
        t == "-c" or (t.startswith("-c") and t[2:].lstrip("0123456789") == "")
        for t in tokens[1:]
    ):
        deny(
            f"Blocked: `{base} -c` on a patient_data path dumps raw bytes to "
            f"the transcript. Line mode (`{base} -n`, `{base} -3`) is still "
            f"allowed for reading the transfer log. See LIMITS.md hole 2."
        )

    paths = [t for t in tokens[1:] if _is_pathlike(t)]
    escaping = [_clean(p) for p in paths if not _is_patient(p)]

    if base in PATH_MUTATORS | SOURCE_DESTRUCTIVE | DEST_LAST:
        if escaping:
            deny(
                f"Blocked: `{base}` would move patient data to {escaping}, "
                f"outside ~/patient_data. Once patient bytes leave that "
                f"directory every path rule here stops matching them. Write "
                f"derived output to ~/{DERIVED_MARKER}/ instead, where it "
                f"stays governed. This rule applies in every guard mode."
            )
        # The derived/ route. Destination-style commands need only their
        # destination there; the sources may be the read-only BAMs.
        # Everything else (rm, mv, chmod, truncate ...) writes to every path
        # it names, so all of them must already be derived output -- that is
        # what keeps the source BAMs unmovable and undeletable.
        if base in DEST_LAST:
            permitted = bool(paths) and _is_derived(paths[-1])
        else:
            permitted = all(_is_derived(pp) for pp in paths if _is_patient(pp))
        if permitted:
            return "allow"
        deny(
            f"Blocked: `{base}` targeting a path under patient_data. The "
            f"patient BAMs are read-only: never copied, moved, re-indexed, "
            f"or archived. Derived output belongs in ~/{DERIVED_MARKER}/, "
            f"which is the only writable location here. This rule applies in "
            f"every guard mode."
        )

    if base in ("samtools", "bcftools", "bgzip", "tabix"):
        non_flag = [_clean(t) for t in tokens[1:] if not t.startswith("-")]
        sub = non_flag[0] if non_flag else None
        writers = {
            "index", "sort", "merge", "reheader", "addreplacerg", "calmd",
            "markdup", "fixmate", "collate", "depad", "split", "faidx",
            "dict", "cat", "ampliconclip", "consensus",
        }
        has_output = any(
            t == "-o" or t == "-O" or t.startswith("--output")
            for t in tokens[1:]
        )
        if sub in writers or has_output:
            if escaping:
                deny(
                    f"Blocked: `{base} {sub}` writes to {escaping}, outside "
                    f"~/patient_data, with a patient_data path among its "
                    f"inputs. A slice written to /tmp stops being governed by "
                    f"these rules -- write it to ~/{DERIVED_MARKER}/ instead."
                )
            out = _output_value(tokens)
            if out is not None and _is_derived(out):
                return "allow"
            # Operating entirely on already-derived output, e.g.
            # `samtools index ~/patient_data/derived/slice.bam`. A slice needs
            # an index to be useful, and this never touches a source BAM.
            if all(_is_derived(pp) for pp in paths if _is_patient(pp)):
                return "allow"
            deny(
                f"Blocked: `{base} {sub}` writes, and a patient_data path is "
                f"among its arguments. Reading (view, head, flagstat, "
                f"idxstats, stats, depth, coverage, quickcheck) is allowed; "
                f"writing beside the BAMs is not. Send derived output to "
                f"~/{DERIVED_MARKER}/ with -o."
            )
    return None


def check_segment(segment: str) -> None:
    tokens = resolve_tokens(segment)
    if not tokens:
        return
    base = os.path.basename(tokens[0])

    # A command name built from a variable ($G, $(which git)) defeats every
    # name-based check below. Legitimate commands here never need one.
    if "$" in tokens[0]:
        deny(
            f"Blocked: command name {tokens[0]!r} is built from a shell "
            f"expansion, which defeats this guard's command-name checks. "
            f"Write the command literally."
        )

    # Patient-data rules run in BOTH modes -- see check_patient_data. A
    # permitted write into patient_data/derived/ short-circuits: the mode
    # rules below would otherwise deny it for being outside /tmp.
    if check_patient_data(tokens, base) == "allow":
        return

    if MODE == "no-git":
        if base in ("git", "gh", "hub", "glab", "tig"):
            deny(
                f"Blocked: `{base}` is unavailable to the patient-data agent. "
                f"This agent handles real patient sequencing data and has no "
                f"route to a commit, by design."
            )
        if base in NETWORK_TOOLS:
            deny(
                f"Blocked: `{base}` can move data off this machine and can "
                f"reach a remote git API. Patient sequencing data does not "
                f"leave this host from this agent."
            )
        args = [t.strip("'\"") for t in tokens[1:] if not t.startswith("-")]

        # Patient-path handling now lives in check_patient_data, which runs
        # in both modes. What remains here is repo-specific: the patient-data
        # agent must not write into the repository even from a non-patient
        # source.
        if base in PATH_MUTATORS | ARCHIVE_TOOLS:
            if REPO_ROOT and any(
                os.path.abspath(os.path.join(os.getcwd(), a)).startswith(REPO_ROOT)
                for a in args
            ):
                deny(
                    f"Blocked: `{base}` writing inside the repository from the "
                    f"patient-data agent. Patient-derived files never enter "
                    f"the repo tree."
                )

        # samtools writes on a patient path are handled by
        # check_patient_data, in both modes.
        check_interpreter(tokens, base, "It can invoke git or move data.")
        return

    # MODE == "no-mutation"
    if base == "git":
        check_git(tokens)
        return
    if base in BLOCKED_ALWAYS:
        deny(
            f"Blocked: `{base}` can modify state outside this agent's remit. "
            f"This agent checks and reports; it does not change things."
        )
    if base in PATH_MUTATORS:
        paths = [t for t in tokens[1:] if not t.startswith("-")]
        if base == "tee":
            paths = [p for p in paths]
        offending = [p for p in paths if not under_allowed_root(p)]
        if offending or not paths:
            deny(
                f"Blocked: `{base}` targeting {offending or '(unspecified)'} "
                f"outside the scratch roots {ALLOWED_WRITE_ROOTS}. This agent "
                f"reports discrepancies; it does not reconcile them."
            )
    check_interpreter(
        tokens, base,
        "Running a script file by path (e.g. the test suite) is allowed.",
    )
    if base in ("sed", "perl") and any(
        t.startswith("-i") or t == "--in-place" for t in tokens[1:]
    ):
        deny(f"Blocked: `{base} -i` edits files in place.")


def main() -> None:
    try:
        raw = sys.stdin.read()
        payload = json.loads(raw)
    except Exception as exc:  # noqa: BLE001 - fail closed on anything
        deny(f"Blocked: guard could not parse the hook payload ({exc!r}). "
             f"Failing closed.")

    try:
        if payload.get("tool_name") != "Bash":
            allow()
        command = payload.get("tool_input", {}).get("command")
        if not isinstance(command, str) or not command.strip():
            deny("Blocked: no inspectable command string in tool_input. "
                 "Failing closed.")

        # Heredocs feed a body to whatever is on the left; treat as inline code.
        if "<<" in command:
            head = resolve_tokens(split_segments(command)[0] if split_segments(command) else command)
            if head and os.path.basename(head[0]) in INTERPRETERS | SHELLS:
                deny("Blocked: heredoc into an interpreter is inline code. "
                     "Run a script file by path instead.")

        if MODE == "no-mutation":
            check_redirections(command)

        for segment in split_segments(command):
            check_segment(segment)

        allow()
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001 - fail closed
        deny(f"Blocked: guard raised {exc!r} while inspecting the command. "
             f"Failing closed.")


if __name__ == "__main__":
    main()
