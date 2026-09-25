#!/usr/bin/env bash
# Install the synthetic-control toolchain in user space, pinned, from bioconda
# via Miniconda:
#   bwa 0.7.15        must report 0.7.15-r1140, the version the NA12878
#                     background was aligned with (its @PG line)
#   art 2016.06.05    ART_Illumina 2.5.8
#   samtools 1.21, bcftools 1.21
#
# Each tool lives in its own environment, so no solver compromise between an
# old bwa build and a current htslib can change a pinned version. Each tool is
# verified by INVOKING it and reading the version it prints -- the installer's
# exit code is not evidence. Exit 0 only if all four report the pinned version.
#
#   install_synth_tools.sh [PREFIX]      (default ~/miniconda3)
set -uo pipefail
PREFIX=${1:-$HOME/miniconda3}
INSTALLER=Miniconda3-latest-Linux-x86_64.sh
BASEURL=https://repo.anaconda.com/miniconda
CH=(--override-channels -c conda-forge -c bioconda)
log() { echo "$(date -Is) $*"; }

if [ -x "$PREFIX/bin/conda" ]; then
  log "miniconda: already at $PREFIX -- installer skipped ($("$PREFIX/bin/conda" --version))"
else
  tmp=$(mktemp -d)
  log "miniconda: downloading $BASEURL/$INSTALLER"
  curl -fsSL --retry 5 -o "$tmp/$INSTALLER" "$BASEURL/$INSTALLER" || { log "FAILED to download the installer"; exit 1; }
  # The sha256 published in the repository's own index page, row for this file.
  want=$(curl -fsSL --retry 5 "$BASEURL/" | tr -d '\n' | grep -oE "href=\"$INSTALLER\".{0,400}" | grep -oE '[0-9a-f]{64}' | head -1)
  got=$(sha256sum "$tmp/$INSTALLER" | cut -d' ' -f1)
  if [ -z "$want" ] || [ "$want" != "$got" ]; then
    log "FAILED: installer sha256 $got does not match the published ${want:-<not found>}"; exit 1
  fi
  log "miniconda: installer sha256 matches the published value ($got)"
  bash "$tmp/$INSTALLER" -b -p "$PREFIX" || { log "FAILED: installer exited $?"; exit 1; }
  rm -rf "$tmp"
  log "miniconda: installed ($("$PREFIX/bin/conda" --version))"
fi
CONDA="$PREFIX/bin/conda"

make_env() {  # name spec...
  local name=$1; shift
  if [ -d "$PREFIX/envs/$name" ]; then
    log "env $name: exists -- create skipped"
  else
    log "env $name: creating with $*"
    "$CONDA" create -y -q -n "$name" "${CH[@]}" --strict-channel-priority "$@" \
      || log "env $name: conda create exited $? (the invocation checks below decide)"
  fi
}
make_env synth-bwa bwa=0.7.15
make_env synth-art art=2016.06.05
make_env synth-hts samtools=1.21 bcftools=1.21

fail=0
check() {  # label expected actual
  if [ "$3" = "$2" ]; then log "VERIFIED $1: invocation reports '$3'"
  else log "NOT VERIFIED $1: invocation reports '${3:-<nothing>}', expected '$2'"; fail=1; fi
}
E=$PREFIX/envs
check bwa "Version: 0.7.15-r1140" \
  "$("$E/synth-bwa/bin/bwa" 2>&1 | grep -m1 '^Version:' | tr -d '\r')"
check art_illumina "Version 2.5.8" \
  "$("$E/synth-art/bin/art_illumina" 2>&1 | grep -m1 -oE 'Version [0-9.]+')"
check samtools "samtools 1.21" "$("$E/synth-hts/bin/samtools" --version 2>&1 | head -1)"
check bcftools "bcftools 1.21" "$("$E/synth-hts/bin/bcftools" --version 2>&1 | head -1)"

# The exact package set of each environment, for the run record.
for env in synth-bwa synth-art synth-hts; do
  "$CONDA" list -n "$env" --explicit --md5 > "$E/$env.explicit.txt" 2>/dev/null \
    && log "env $env: explicit package list -> $E/$env.explicit.txt"
done
[ $fail -eq 0 ] && log "=== all four tools verified by invocation ===" || log "=== NOT all tools verified ==="
exit $fail
