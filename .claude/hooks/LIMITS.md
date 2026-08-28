# What the subagent guard actually blocks, and what gets through

`guard.py` is a `PreToolUse` hook wired into three subagents through their
frontmatter, so it fires only while that subagent is running. It inspects the
`Bash` command string and returns `permissionDecision: "deny"` to block.

This document exists because **an incomplete deny list that looks complete is
worse than an honest advisory note.** Everything below was tested, not
assumed — see `test_guard.py`, which asserts both directions and asserts that
the documented holes are still holes.

---

## Enforcement model

| Layer | Mechanism | Scope | Strength |
|---|---|---|---|
| Tool removal | `tools:` / `disallowedTools:` frontmatter | per-agent | **Structural.** The tool is not in the agent's schema. It cannot be called. |
| Bash filtering | `hooks: PreToolUse` → `guard.py` | per-agent | **Structural but incomplete.** Blocks by command name; a determined interpreter call gets past it. |
| Prompt rules | system prompt text | per-agent | **Advisory.** The agent can decline to follow it. |

The middle row is the one to be careful about. It is real enforcement — the
call genuinely does not run — but it filters *command text*, and command text
is not a sandbox.

---

## Structural: cannot happen

**verifier, thesis-editor** — `disallowedTools: Write, Edit`. Neither tool is
in the agent's schema. There is no `Write` or `Edit` to call, regardless of
what the agent decides to do. This is the same class of fix as FIX C in the
benchmark: the parameter was removed, so it cannot be supplied.

**patient-data** — `Write`, `Edit` and `NotebookEdit` are absent from its
`tools:` list for the same reason.

---

## Structural but incomplete: the Bash guard

### Verified blocked (see `test_guard.py`)

**`no-git` mode — patient-data**

- Every `git` invocation, including `git -C`, `sudo git`, `env git`,
  `/usr/bin/git`, `$(which git) commit`, and `` `which git` commit ``
- `gh`, `hub`, `glab`, `tig`
- Inline code that could shell out: `python3 -c`, `bash -c`, `sh -c`,
  `python3 -m`, a bare interpreter reading stdin (`… | bash`,
  `printf x | python3 -`), heredocs into an interpreter
- Network and transfer tools — `curl`, `wget`, `ssh`, `scp`, `sftp`, `rsync`,
  `nc`, `socat`, `aws`, `gcloud`, `rclone`. These close the GitHub-API route
  to a commit *and* the route to moving sequencing data off this host.
- Any `cp`/`mv`/`ln`/`rm`/`tar`/`zip` touching a `patient_data` path, or
  writing anywhere inside the repository
- `samtools`/`bcftools` **writing** subcommands (`index`, `sort`, `merge`,
  `reheader`, …) or `-o` output when a `patient_data` path is among the
  arguments

Reading is deliberately untouched: `samtools view -H`, `view -c`, `idxstats`,
`flagstat`, `quickcheck`, `ls`, `du` on the patient files all pass. That is
the agent's job, and a guard that blocked it would just get switched off.

**`no-mutation` mode — verifier, thesis-editor**

- `git` outside a read-only allowlist. Permitted: `log`, `show`, `diff`,
  `status`, `blame`, `ls-files`, `rev-parse`, `cat-file`, `check-ignore`,
  `reflog`, and similar. Everything else denied — including `commit`, `add`,
  `push`, `reset`, `checkout`, `stash`, `rebase`, `rm`. An allowlist is used
  rather than a blocklist so a mutating subcommand cannot be missed by
  omission.
- `rm`, `mv`, `cp`, `ln`, `tee`, `truncate`, `chmod`, `patch`, `install`
  targeting anything outside `/tmp`
- Shell redirection `>` / `>>` to any path outside `/tmp` (`2>&1` and
  `> /dev/null` still work)
- `sed -i`, `perl -i`
- Editors: `vim -c wq`, `ex -sc wq`, `nano`, `emacs`, `ed`
- Build/package tools: `make`, `npm`, `pip`, `cargo`, `go`, `cmake`, `conda`
- `gh`, `dd`, `shred`, `crontab`, `systemctl`
- Inline code and stdin scripts, as above

Reading and running are untouched: the three test suites, `git log`, `grep`,
`cat`, `find`, `wc`, `diff`, and redirection into `/tmp` all pass.

### Patient data is protected in BOTH modes (fixed 2026-08-28)

Originally every `patient_data` rule sat inside `if MODE == "no-git"`, so a
`verifier` or `thesis-editor` run could `tar` and `curl` the patient BAMs
that the `patient-data` agent exists to protect. **The boundary was drawn
around the agent, not around the data.** That was hole 3; it was found by
probing the guard, not by reading it.

`check_patient_data()` now runs before any mode-specific logic, in every
mode. With a `patient_data` path among the arguments, these deny regardless
of which agent is running:

- network egress — `curl`, `wget`, `scp`, `rsync`, `ssh`, `nc`, `sftp`, …
- archiving — `tar`, `zip`, `gzip`, `bgzip`, `bzip2`, `xz`, `zstd`, `7z`
- copy/move/link/delete/chmod — every `PATH_MUTATORS` command
- `samtools`/`bcftools` writing subcommands and `-o` output

Asserted for both modes in `test_guard.py`. `no-mutation` keeps its own
freedoms where no patient path is involved: `tar` into `/tmp`, `curl` to a
URL, `rm /tmp/scratch` all still pass.

### The derived/ write route (added 2026-08-28)

`~/patient_data/derived/` is the only location a guarded agent may write
BAM-derived output to. It sits **inside** `patient_data` on purpose: every
rule above still matches what lands there, so a slice stays governed.

Permitted in both modes:

```
samtools view -b -o ~/patient_data/derived/slice.bam <source.bam> <region>
samtools index ~/patient_data/derived/slice.bam
samtools sort  -o ~/patient_data/derived/sorted.bam ~/patient_data/derived/slice.bam
mkdir / rm     ~/patient_data/derived/...
cp <source.bam> ~/patient_data/derived/
```

Still denied, and asserted in `test_guard.py`:

- taking derived output anywhere else — `cp`/`mv` to `/tmp` or the repo,
  `curl`, `scp`, `tar`, `base64` on a derived path
- writing to, moving, or deleting the source BAMs, including
  `mv ~/patient_data/<source>.bam ~/patient_data/derived/`
- copying *into* `derived/` from outside `patient_data` — produce output
  there directly rather than staging it in `/tmp` first

The reasoning is the laundering argument. `samtools view -b -o /tmp/slice.bam`
looks harmless and is not: `/tmp` is a permitted scratch root in
`no-mutation`, so the moment patient bytes land there, every `patient_data`
path rule stops matching them and the slice becomes freely copyable. Keeping
derived output inside the governed namespace is what makes slicing safe to
allow at all.

Commands that write to `derived/` short-circuit the mode rules, which would
otherwise deny them in `no-mutation` for being outside `/tmp`. That is
deliberate: `derived/` is not the working tree, so nothing `no-mutation`
protects is weakened.

### Verification status: live-fired vs. asserted by the suite

These are not the same claim, and this file should not blur them.

`test_guard.py` (267 assertions) runs `guard.py` as a subprocess with the real
hook payload shape, in both modes. That proves the **decision function**. It
does not prove the hook is registered with the agent and firing — a guard that
is correct but unwired blocks nothing.

**Live-fired** means a real subagent issued the command and the `PreToolUse`
hook actually intercepted it, or actually let it through.

| Rule | Live-fired | Where |
|---|---|---|
| git allowlist, incl. `-C`, `/usr/bin/git`, env prefix, `$(which git)`, `&&`/`;`, extra whitespace, `--allow-empty` | yes | verifier |
| `bash -c` / `python3 -c` inline code | yes | verifier, patient-data |
| `sed -i`, `tee`, `cp`, `mv`, redirection outside `/tmp` | yes | verifier |
| `git` and `gh` wholly unavailable | yes | patient-data |
| network egress with a patient path (`curl`, `scp`, `rsync`, `wget`) | yes | patient-data, verifier |
| archiving a patient path (`tar`) | yes | patient-data, verifier |
| copy/move/delete/index on a patient path | yes | patient-data, verifier |
| content dumps (`base64`, `xxd`, `head -c`) on a patient path | yes | patient-data, verifier, thesis-editor |
| `derived/` write route allowed, and works end to end | yes | patient-data |
| `-o /tmp/...` denied while `-o derived/...` allowed | yes | patient-data |
| fail-closed on malformed payloads | **suite only** | — |
| quote-awareness and unbalanced-quote fallback | **suite only** | — |
| the `no-mutation` half of the full deny matrix | **suite only** | see below |

Two details worth keeping straight:

- The exhaustive deny matrix was fired through the **patient-data** agent, so
  it exercised `no-git`. In `no-mutation` the individually live-fired rules are
  the ones listed above; the rest of that matrix rests on the suite, which
  asserts every case in both modes.
- Discrimination was confirmed, not just denial: the same `curl` to the same
  host **without** a patient path returned 200, and the same `tar` into `/tmp`
  succeeded. The guard keys on the path, not on the binary name.

### Why the deny side needed a decoy

Two subagents in succession refused to probe the deny rules against the real
BAMs, and both were right. Probing a guard by performing the prohibited act is
only sound when the guard is known-good; the premise each time was that it had
just changed. With no cleanup permitted, a miss would have rewritten a source
BAM's index, separated a 39 GB BAM from its index, or deleted the transfer log.

`~/patient_data_probe/` resolves this: a synthetic BAM with no provenance, its
own log and `derived/`, where destruction is free. Because `patient_data` is a
substring of `patient_data_probe`, every rule matches it identically; the suite
asserts the decoy verdict-for-verdict against the real tree in both modes, so a
decoy that silently diverges fails the suite rather than quietly invalidating
the verification. Egress probes point at a listener on `127.0.0.1`, so a guard
failure lands in a local log instead of leaving the host.

All twelve deny shapes fired live against the decoy on 2026-08-29 and were
denied; the `derived/` route was allowed on the same decoy in the same run,
confirming the guard was not simply blanket-denying; zero bytes reached the
listener.

**The residue this leaves.** A subagent also declined to delete real
patient-derived files from `derived/` — an operation the guard *permits*. The
hook was not the backstop there; the agent's own judgment was. That is the
advisory layer at the bottom of this file doing the work, and it is the honest
shape of the whole system: the guard makes the common accident impossible, and
the things it cannot reach are still declinable instructions.

### Fail-closed

A malformed payload, a missing command string, or any internal exception
returns `deny`. Tested. A guard that allowed on error would be worse than
none, because it would look like one.

---

## What gets through — the honest list

### 1. Any script file. This is the fundamental hole and it cannot be closed here.

`python3 anything.py` is allowed, because `verifier` **must** run
`test_bam_tools.py`, and that suite legitimately writes synthetic BAMs to
temporary paths. Once an interpreter runs a file, the file can do anything —
write into the repo, call git, move data.

The same applies to `./helper.sh` and any executable on `PATH`.

Claude Code's own documentation states this limit directly: deny rules
*"don't apply to arbitrary subprocesses that read or write files indirectly,
like a Python or Node script that opens files itself."*

**This is not fixable at the hook layer.** Blocking script execution would
stop verifier verifying, which is its entire purpose. Closing it requires
OS-level enforcement — Claude Code's sandbox, a container, or a read-only
bind mount. Nothing in this repository does that today.

What remains true even so: the agent has no `Write`/`Edit` tool, and would
have to *author* a script to get a write primitive — which itself requires a
write. The realistic residual risk is an agent invoking a script that already
exists and happens to mutate state, not an agent bootstrapping one.

### 2. Reading patient data into the transcript — narrowed, NOT closed

`base64`, `base32`, `xxd`, `od`, `hexdump`, `strings` and `head -c`/`tail -c`
on a `patient_data` path now deny in both modes (fixed 2026-08-28). That also
settles an incoherence found by the same probe run: `gzip -c patient.bam`
denied while `base64 patient.bam` allowed, same channel, opposite verdicts.
Both deny now.

**Do not read that as the class being closed. It is not, and it cannot be.**

The exfiltration channel here is the conversation transcript, and no hook can
inspect what an agent writes into its own report. Two things still reach it:

- **`samtools view` emits read-level content to stdout**, and it must stay
  open — reading the BAM is the patient-data agent's entire job. This is the
  same channel `base64` used.
- **`cat` is deliberately left open**, because `cat ~/patient_data/download.log`
  is ordinary work. `cat` on a BAM therefore also passes. Line-mode `head`/`tail`
  are open for the same reason; only `-c` byte mode is denied.

So the denials above raise the cost of an *accident* — no bioinformatics
workflow needs `base64` on a BAM, so nothing legitimate is lost — while
leaving the deliberate path wide open. The rule "never output read sequences,
sample identifiers, or patient-identifying data" remains **advisory**, in the
section at the bottom of this file, and it is still the most important
patient-data rule and the least enforceable one.

### 3. Obfuscation the tokenizer does not model

Segment splitting and redirection scanning are quote-aware as of 2026-08-28:
`grep -E 'tar|curl|cp ' LIMITS.md` and `grep -n 'a>b' f` used to deny, because
the `|` and `>` inside a quoted argument were read as a pipe and a
redirection. Unbalanced quoting falls back to the old blunt split, so the
failure direction is still "inspect too much".

Beyond that, `guard.py` is a best-effort tokenizer, not a shell parser. Untested and
likely to get through: nested `$(…)` beyond one level, arithmetic expansion,
`eval` of an assembled string, `alias`, shell functions defined earlier in
the same command, `xargs` reading command names from stdin, unusual quoting.

Variable-built command names (`export G=git; $G commit`) *are* caught, by
denying any command whose first token contains `$`.

### 4. Scope: the guard only fires for these three agents

`bam-tools-dev`, `test-writer`, `benchmark-runner`, `docs-writer` and the
main session are unaffected and can do everything above. That is deliberate —
they need to. But see hole 3: the two agents that *are* guarded
still differ in what they protect, so "guarded" does not imply "patient-safe".

### 5. The agents must actually be loaded

Hooks live in agent frontmatter. If the agent is invoked some other way, or
the frontmatter is edited, the guard does not run. It is configuration, not a
kernel.

### 6. Not a defence against a hostile agent

This raises the cost of an accident and makes the common mistake impossible.
It is not an adversarial boundary. A model actively trying to commit patient
data has routes — most simply, asking the main session to do it.

---

## What is still advisory, and stays advisory

These live in the agents' system prompts and cannot be enforced by any
mechanism available here:

- **Never output read sequences, sample identifiers, or patient-identifying
  data.** No hook can inspect what an agent writes in its report. This is the
  most important patient-data rule and it is the least enforceable one.
- **Never present a finding as a clinical result.**
- **Report, do not reconcile** (verifier) — the agent could still describe a
  fix instead of reporting a discrepancy.
- **Do not write thesis prose** (thesis-editor) — it has no `Write`/`Edit`,
  so it cannot edit the chapter, but it could still emit prose in its report.

Per this project's own finding: *an instruction a model can decline is not a
constraint.* These four are declinable. They are documented here so that is
visible rather than assumed.

---

## Re-running the tests

```bash
python3 .claude/hooks/test_guard.py
```

Both directions, plus fail-closed behaviour, plus assertions that the known
holes are still open. If a `KNOWN HOLES` line reports `CLOSED`, the guard got
stronger and this document is now overstating the gap — update it.
