#!/usr/bin/env bash
# Route displayed output through the patient-identifier redaction (redact.py).
#
#   some_command 2>&1 | scripts/redact.sh
#   scripts/redact.sh -- some_command args ...
#
# The second form filters the command's stdout and stderr together and exits
# with the COMMAND's exit code -- unless the filter itself refused to run (it
# fails closed when it cannot load its term list), in which case it exits with
# the filter's code, because the command's output was withheld.
here=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
if [ "${1:-}" = "--" ]; then
  shift
  "$@" 2>&1 | python3 "$here/redact.py"
  st=("${PIPESTATUS[@]}")
  [ "${st[1]}" -ne 0 ] && exit "${st[1]}"
  exit "${st[0]}"
fi
exec python3 "$here/redact.py" "$@"
