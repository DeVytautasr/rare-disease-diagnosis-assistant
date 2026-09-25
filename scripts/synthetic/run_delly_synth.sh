#!/usr/bin/env bash
# delly on ONE synthetic-control BAM -- an implant or the untouched background --
# with the same command form and exclude template as the patient runs:
# delly v2.6.0 sr, one thread (-h 1; its default is four), the unmodified v2.6.0
# exclude template, hs38DH.
#
#   run_delly_synth.sh LABEL BAM        -> ~/public_data/sim/delly/LABEL.bcf
#
# Pre-flight re-verifies delly's sha256, the template's git blob and the .fai md5
# (the FASTA's md5 is checked once by whoever launches a batch, not 13 times at
# once). Afterwards it reports how many samples the BCF names -- each implant BAM
# must present as ONE sample -- and how many records it holds.
set -uo pipefail
label=${1:?usage: run_delly_synth.sh LABEL BAM}
bam=${2:?usage: run_delly_synth.sh LABEL BAM}
REF=$HOME/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
EXCL=$HOME/reference/human.hg38.excl.tsv
DELLY=$HOME/tools/delly/delly
BCFTOOLS=$HOME/miniconda3/envs/synth-hts/bin/bcftools
OUT=$HOME/public_data/sim/delly
fail=0
pre() { if [ "$2" = "$3" ]; then echo "pre-flight OK   $1"; else echo "pre-flight FAIL $1: got '$2', want '$3'"; fail=1; fi; }
pre "delly sha256" "$(sha256sum "$DELLY" | cut -d' ' -f1)" 85ecf4d64e23672a51c71f6e3a2dfda997753266157ade79ef598c8f96ecb469
pre "exclude template git blob (v2.6.0, unmodified)" "$(git hash-object --no-filters "$EXCL")" 3125a61491cb1c62ce46b365897eb4cd4e9fdb66
pre "hs38DH .fai md5" "$(md5sum "$REF.fai" | cut -d' ' -f1)" 5ccc91e56dc4a05448dd5b9507ec6bc6
pre "BAM and its index present" "$([ -s "$bam" ] && [ -s "$bam.bai" ] && echo yes)" yes
[ $fail -eq 0 ] || { echo "REFUSING to run delly: a pre-flight check failed"; exit 3; }
mkdir -p "$OUT"
export OMP_NUM_THREADS=1
cmd=("$DELLY" sr -h 1 -g "$REF" -x "$EXCL" -o "$OUT/$label.bcf" "$bam")
echo "command: OMP_NUM_THREADS=1 /usr/bin/time -v ${cmd[*]}" | sed "s|$HOME|~|g"
echo "started: $(date -Is)"
/usr/bin/time -v "${cmd[@]}"
rc=$?
echo "ended: $(date -Is) delly exit=$rc"
if [ $rc -eq 0 ]; then
  echo "samples in BCF: $("$BCFTOOLS" query -l "$OUT/$label.bcf" | wc -l) ($("$BCFTOOLS" query -l "$OUT/$label.bcf" | tr '\n' ' '))"
  echo "records in BCF: $("$BCFTOOLS" view -H "$OUT/$label.bcf" | wc -l)"
fi
exit $rc
