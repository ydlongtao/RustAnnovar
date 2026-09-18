#!/usr/bin/env bash
# Resumable official score-only downloads; checksum before exposing final files.
set -euo pipefail
root=${1:?usage: hpc_cadd_download.sh ROOT [hg38|hg19]}
build=${2:-hg38}
case "$build" in hg38) assembly=GRCh38 ;; hg19) assembly=GRCh37 ;; *) exit 2 ;; esac
dir="$root/humandb/cadd-v1.7/$assembly"
mkdir -p "$dir"
exec 9>"$dir/download.lock"
flock -n 9 || { echo 'another CADD download is live' >&2; exit 3; }
base="https://krishna.gs.washington.edu/download/CADD/v1.7/$assembly"
curl --fail --location --retry 8 --connect-timeout 30 "$base/MD5SUMs" -o "$dir/MD5SUMs.new"
mv "$dir/MD5SUMs.new" "$dir/MD5SUMs"
for name in whole_genome_SNVs.tsv.gz.tbi whole_genome_SNVs.tsv.gz; do
  checksum=$(awk -v name="$name" '$2 == name {print $1}' "$dir/MD5SUMs")
  [[ "$checksum" =~ ^[0-9a-f]{32}$ ]] || { echo "missing official checksum for $name" >&2; exit 4; }
  if [[ -f "$dir/$name" ]]; then
    printf '%s  %s\n' "$checksum" "$dir/$name" | md5sum -c -
    continue
  fi
  curl --fail --location --continue-at - --retry 20 --retry-delay 10 \
    --connect-timeout 30 --speed-limit 32 --speed-time 90 \
    "$base/$name" -o "$dir/$name.part"
  printf '%s  %s\n' "$checksum" "$dir/$name.part" | md5sum -c -
  mv "$dir/$name.part" "$dir/$name"
done
printf '%s\n' "CADD $build v1.7: both files verified against official MD5SUMs"
