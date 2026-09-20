#!/usr/bin/env bash
# Retry transient CADD transport failures while preserving verified chunks.
set -u
root=${1:?usage: hpc_cadd_recovery.sh ROOT BUILD NEW_LOG_DIRECTORY}
build=${2:?}
log_dir=${3:?}
case "$build" in hg19|hg38) ;; *) exit 2 ;; esac
mkdir "$log_dir" || exit 2
downloader="$root/tmp/hpc_cadd_parallel_download-bounded-queue.py"
test -f "$downloader" || exit 2
for attempt in $(seq 1 48); do
  printf '%s attempt %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$attempt"
  python3 "$downloader" --root "$root" --build "$build" --workers 4 \
    > "$log_dir/attempt-$(printf '%02d' "$attempt").stdout" \
    2> "$log_dir/attempt-$(printf '%02d' "$attempt").stderr"
  code=$?
  if (( code == 0 )); then
    printf 'official CADD download and MD5 verification complete\n'
    exit 0
  fi
  last="$log_dir/attempt-$(printf '%02d' "$attempt").stderr"
  # Retry only connection failures. Identity, size, and checksum errors must
  # stop rather than reusing potentially inconsistent data.
  if grep -Eq 'ValueError:|RuntimeError:' "$last"; then
    tail -n 12 "$last" >&2
    exit "$code"
  fi
  if ! grep -Eq 'URLError|TimeoutError|timed out|Temporary failure in name resolution|truncated CADD range|ConnectionResetError|RemoteDisconnected|Connection aborted' "$last"; then
    tail -n 12 "$last" >&2
    exit "$code"
  fi
  tail -n 2 "$last" >&2
  if (( attempt < 48 )); then sleep 120; fi
done
echo 'transient CADD download retries exhausted; chunks preserved' >&2
exit 1
