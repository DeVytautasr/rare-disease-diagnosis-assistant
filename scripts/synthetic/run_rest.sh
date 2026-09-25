#!/usr/bin/env bash
# After the (revised) IMP01 gate passed: build the remaining implants and run
# delly on every implant and on the untouched background, concurrently.
#
#   run_rest.sh [IMPxx ...]      (default: IMP02..IMP12 still to build)
#
# Every step is its OWN job -- own tmux window, own log, own exit-code file --
# launched through scripts/restore/job.sh into ~/public_data/sim/logs:
#   build_IMPxx   make_implants.py build IMPxx   (at most MAXBUILD at once, default 3)
#   delly_IMPxx   run_delly_synth.sh IMPxx       (starts as soon as its build exits 0)
#   delly_background, delly_IMP01                (start at once: their BAMs exist)
# A delly job waits for at least MINFREE_MIB of available memory before starting.
# The hs38DH FASTA's md5 is verified once here, for the whole batch.
# Exit 0 only if every build and every delly job exited 0.
set -u
R=$(cd "$(dirname "$(readlink -f "$0")")/../.." && pwd)
L=$HOME/public_data/sim/logs
SESSION=${SESSION:-claude-rerun}
MAXBUILD=${MAXBUILD:-3}
MINFREE_MIB=${MINFREE_MIB:-3000}
PY=$R/.venv/bin/python
BG=$HOME/public_data/NA12878.chr20_chr21.bam
BAMS=$HOME/public_data/sim/bams
REF=$HOME/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
mkdir -p "$L"
log() { echo "$(date -Is) $*"; }

got=$(md5sum "$REF" | cut -d' ' -f1)
[ "$got" = 64b32de2fc934679c16e83a2bc072064 ] || { log "FASTA md5 $got is not the published md5 -- refusing"; exit 3; }
log "hs38DH FASTA md5 verified once for the batch"

job() {  # name command...
  local name=$1; shift
  tmux new-window -d -t "$SESSION" -n "$name" "bash $R/scripts/restore/job.sh $L $name -- $*; read"
  log "launched $name"
}
delly_job() {  # label bam
  while [ "$(awk '/^MemAvailable/{print int($2/1024)}' /proc/meminfo)" -lt "$MINFREE_MIB" ]; do sleep 10; done
  job "delly_$1" bash "$R/scripts/synthetic/run_delly_synth.sh" "$1" "$2"
}

delly_job background "$BG"
delly_job IMP01 "$BAMS/IMP01.bam"

pending=("$@")
[ ${#pending[@]} -eq 0 ] && pending=(IMP02 IMP03 IMP04 IMP05 IMP06 IMP07 IMP08 IMP09 IMP10 IMP11 IMP12)
running=()
failed=0
while [ ${#pending[@]} -gt 0 ] || [ ${#running[@]} -gt 0 ]; do
  still=()
  for i in "${running[@]}"; do
    if [ -f "$L/build_$i.rc" ]; then
      rc=$(cat "$L/build_$i.rc")
      if [ "$rc" = 0 ]; then delly_job "$i" "$BAMS/$i.bam"; else log "build_$i FAILED (exit $rc) -- no delly for it"; failed=1; fi
    else
      still+=("$i")
    fi
  done
  running=("${still[@]}")
  while [ ${#running[@]} -lt "$MAXBUILD" ] && [ ${#pending[@]} -gt 0 ]; do
    i=${pending[0]}; pending=("${pending[@]:1}")
    job "build_$i" "$PY" "$R/scripts/synthetic/make_implants.py" build "$i"
    running+=("$i")
  done
  sleep 5
done
log "all builds finished; waiting for the delly jobs"
labels=(background IMP01 IMP02 IMP03 IMP04 IMP05 IMP06 IMP07 IMP08 IMP09 IMP10 IMP11 IMP12)
for l in "${labels[@]}"; do
  [ -f "$L/delly_$l.pid" ] || continue
  until [ -f "$L/delly_$l.rc" ]; do sleep 10; done
  rc=$(cat "$L/delly_$l.rc"); log "delly_$l exit $rc"
  [ "$rc" = 0 ] || failed=1
done
log "=== pipeline finished: $([ $failed = 0 ] && echo 'every build and delly job exited 0' || echo 'AT LEAST ONE JOB FAILED') ==="
exit $failed
