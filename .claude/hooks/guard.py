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
    segments.extend(re.split(r"(?:\|\||&&|[;|\n&])", stripped))
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
    for m in re.finditer(r"(?<![0-9<>])>>?\s*([^\s;|&<>]+)", command):
        target = m.group(1)
        if target.startswith("&"):
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
        touches_patient = any("patient_data" in a for a in args)

        # Copying, moving, linking or archiving the patient files -- anywhere,
        # including into the repository.
        if base in PATH_MUTATORS | {"tar", "zip", "gzip", "bgzip", "bzip2"}:
            if touches_patient:
                deny(
                    f"Blocked: `{base}` targeting a path under patient_data. "
                    f"The patient BAMs are read-only: never copied, moved, "
                    f"re-indexed, or archived, and never placed in the repo."
                )
            if REPO_ROOT and any(
                os.path.abspath(os.path.join(os.getcwd(), a)).startswith(REPO_ROOT)
                for a in args
            ):
                deny(
                    f"Blocked: `{base}` writing inside the repository from the "
                    f"patient-data agent. Patient-derived files never enter "
                    f"the repo tree."
                )

        # samtools reads are the agent's core job; samtools writes are not.
        if base in ("samtools", "bcftools", "bgzip", "tabix"):
            sub = args[0] if args else None
            writers = {
                "index", "sort", "merge", "reheader", "addreplacerg", "calmd",
                "markdup", "fixmate", "collate", "depad", "split", "faidx",
                "dict", "cat", "ampliconclip", "consensus",
            }
            has_output = any(
                t == "-o" or t.startswith("--output") or t == "-O"
                for t in tokens[1:]
            )
            if touches_patient and (sub in writers or has_output):
                deny(
                    f"Blocked: `{base} {sub}` writes, and a patient_data path "
                    f"is among its arguments. Reading (view, head, flagstat, "
                    f"idxstats, stats, depth, coverage, quickcheck) is "
                    f"allowed; writing beside the BAMs is not."
                )
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
