#!/usr/bin/env bash
# delly call on ONE patient sample, whole genome, one thread (Task 3a / Task 6).
#
#   run_delly.sh SAMPLE_A|SAMPLE_B
#
# The two samples run as two SEPARATE processes: a joint two-sample call
# produces a joint call set, which Phase 3 rejected. The only input path on the
# command line is the deid symlink, so no identifier reaches a command line, a
# process listing or this log. The BCF still carries the @RG SM value in its
# sample column and stays under ~/patient_data.
#
# Pre-flight -- every input is re-verified against a published or recorded
# value before delly starts, and the run is refused if any check fails:
#   delly v2.6.0 sha256 (GitHub release digest), the exclude template's git
#   blob at tag v2.6.0 (unmodified), the hs38DH FASTA and .fai md5s (published
#   by 1000 Genomes), and the deid link resolving to the Phase 0 byte size.
#
# Output: ~/patient_data/rerun_2026-09/delly/<label>.bcf (+ .csi). Stdout and
# stderr -- delly's timestamped phase markers and /usr/bin/time -v's report,
# including peak RSS -- go to the job log (run it under scripts/restore/job.sh).
set -uo pipefail
label=${1:?usage: run_delly.sh SAMPLE_A|SAMPLE_B}
case $label in SAMPLE_A) want_bytes=38959428903 ;; SAMPLE_B) want_bytes=41617797998 ;;
  *) echo "unknown label"; exit 2 ;; esac
REF=$HOME/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
EXCL=$HOME/reference/human.hg38.excl.tsv
DELLY=$HOME/tools/delly/delly
BAM=$HOME/patient_data/deid/$label.bam
OUT=$HOME/patient_data/rerun_2026-09/delly
fail=0
pre() {  # name got want
  if [ "$2" = "$3" ]; then echo "pre-flight OK   $1"; else echo "pre-flight FAIL $1: got '$2', want '$3'"; fail=1; fi
}
pre "delly sha256" "$(sha256sum "$DELLY" | cut -d' ' -f1)" 85ecf4d64e23672a51c71f6e3a2dfda997753266157ade79ef598c8f96ecb469
pre "exclude template git blob (v2.6.0, unmodified)" "$(git hash-object --no-filters "$EXCL")" 3125a61491cb1c62ce46b365897eb4cd4e9fdb66
pre "hs38DH .fai md5" "$(md5sum "$REF.fai" | cut -d' ' -f1)" 5ccc91e56dc4a05448dd5b9507ec6bc6
pre "hs38DH FASTA md5" "$(md5sum "$REF" | cut -d' ' -f1)" 64b32de2fc934679c16e83a2bc072064
pre "deid link resolves to the $label Phase 0 size" "$(stat -L -c %s "$BAM")" "$want_bytes"
pre "index reachable through the link" "$([ -e "$BAM.bai" ] && echo yes)" yes
[ $fail -eq 0 ] || { echo "REFUSING to run delly: a pre-flight check failed"; exit 3; }
mkdir -p -m 700 "$OUT"
echo "delly: $("$DELLY" 2>&1 | grep -m1 'Version' | tr -s ' ')"
export OMP_NUM_THREADS=1
# delly v2.6.0 has no "call" subcommand: short-read discovery is "sr", and it
# defaults to FOUR threads (-h 4), so one thread has to be asked for explicitly.
cmd=("$DELLY" sr -h 1 -g "$REF" -x "$EXCL" -o "$OUT/$label.bcf" "$BAM")
echo "command: OMP_NUM_THREADS=1 /usr/bin/time -v ${cmd[*]/#$HOME/~}"
echo "started: $(date -Is)"
/usr/bin/time -v "${cmd[@]}"
rc=$?
echo "ended: $(date -Is) delly exit=$rc"
exit $rc
