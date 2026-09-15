#!/usr/bin/env bash
# Run under nohup; one immutable state directory per submitted job.
set -u
if (( $# < 2 )); then
  echo 'usage: hpc_run.sh NEW_STATE_DIRECTORY COMMAND [ARG ...]' >&2
  exit 2
fi
state=$1
shift
mkdir -- "$state" || exit 2
printf '%s\n' "$$" > "$state/runner.pid"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$state/started.utc"
printf '%q ' "$@" > "$state/command.sh"
printf '\n' >> "$state/command.sh"
printf 'RUNNING\n' > "$state/status"
"$@" > "$state/stdout.log" 2> "$state/stderr.log" &
child=$!
printf '%s\n' "$child" > "$state/child.pid"
wait "$child"
code=$?
printf '%s\n' "$code" > "$state/exit_code"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$state/finished.utc"
if (( code == 0 )); then
  printf 'SUCCEEDED\n' > "$state/status"
else
  printf 'FAILED\n' > "$state/status"
fi
exit "$code"
