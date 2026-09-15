#!/usr/bin/env bash
set -euo pipefail
root=${RUSTANNOVAR_WORKSPACE:-/DATABANK/users/hflt/RustAnnovar}
mkdir -p "$root/datasets/giab" "$root/runs"
exec 9>"$root/runs/datasets.lock"
flock -n 9 || exit 0
: > "$root/runs/datasets.exit"
trap 'echo $? > "$root/runs/datasets.exit"' EXIT
base=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1
for build in GRCh37 GRCh38; do
  name="HG002_${build}_1_22_v4.2.1_benchmark.vcf.gz"
  if [ ! -f "$root/datasets/giab/$name" ]; then
    curl -fL --connect-timeout 20 --max-time 7200 --retry 3 --retry-all-errors --continue-at - "$base/$build/$name" -o "$root/datasets/giab/$name.part"
    gzip -t "$root/datasets/giab/$name.part"
    mv "$root/datasets/giab/$name.part" "$root/datasets/giab/$name"
  fi
  sha256sum "$root/datasets/giab/$name" >> "$root/datasets/giab/SHA256SUMS"
  printf '%s\t%s\n' "$name" "$base/$build/$name" >> "$root/datasets/giab/sources.tsv"
done
