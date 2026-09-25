#!/usr/bin/env bash
# Run one background job with its own log and its own exit-code file.
#
#   job.sh DIR NAME -- command [args ...]
#
# DIR/NAME.log      stdout and stderr of the command
# DIR/NAME.pid      the command's own PID (kill by this PID, never by pattern)
# DIR/NAME.started  ISO time the command started
# DIR/NAME.ended    ISO time it ended
# DIR/NAME.rc       its exit code -- written LAST, so a waiter polls for this file
#
# An existing NAME.rc is removed first: a stale exit code must never satisfy a
# waiter for a new run.
set -u
dir=$1 name=$2
shift 2
[ "${1:-}" = "--" ] && shift
mkdir -p "$dir"
rm -f "$dir/$name.rc" "$dir/$name.ended"
date -Is > "$dir/$name.started"
"$@" > "$dir/$name.log" 2>&1 &
child=$!
echo "$child" > "$dir/$name.pid"
wait "$child"
rc=$?
date -Is > "$dir/$name.ended"
echo "$rc" > "$dir/$name.rc"
exit "$rc"
