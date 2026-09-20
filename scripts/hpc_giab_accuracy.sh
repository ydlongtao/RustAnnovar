#!/usr/bin/env bash
# Run a shared-AVinput comparison on a selected public GIAB VCF.
set -euo pipefail
root=${1:?usage: hpc_giab_accuracy.sh ROOT SAMPLE_DIRECTORY NEW_RUN_DIRECTORY BUILD PROTOCOL}
sample=${2:?}
run=${3:?}
build=${4:?}
protocol=${5:?}
case "$build" in
  hg19) assembly=GRCh37 ;;
  hg38) assembly=GRCh38 ;;
  *) echo 'build must be hg19 or hg38' >&2; exit 2 ;;
esac
annovar="$root/baseline/annovar"
database="$annovar/humandb"
binary="$root/builds/cadd-5499531/rust-annovar"
selection=${6:-"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/giab_accuracy.py"}
input="$sample/HG002_${assembly}_accuracy.vcf"
mkdir "$run"
test -s "$input" && test -s "$sample/selection.json"
test -x "$binary" && test -f "$database/${build}_${protocol}.txt"
test -f "$database/${build}_${protocol}Mrna.fa"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$run/started.utc"
uname -a > "$run/host.txt"
perl -v > "$run/perl-version.txt"
"$binary" --version > "$run/rust-version.txt"
sha256sum "$input" "$sample/selection.json" "$binary" \
  "$annovar/convert2annovar.pl" "$annovar/table_annovar.pl" \
  "$annovar/annotate_variation.pl" "$database/${build}_${protocol}.txt" \
  "$database/${build}_${protocol}Mrna.fa" "$selection" > "$run/SHA256SUMS"

/usr/bin/time -f 'elapsed_seconds=%e\npeak_rss_kib=%M\nexit_code=%x' \
  -o "$run/convert.metrics" \
  perl "$annovar/convert2annovar.pl" --format vcf4 \
  --outfile "$run/shared.avinput" "$input" \
  > "$run/convert.stdout" 2> "$run/convert.stderr"
test -s "$run/shared.avinput"
sha256sum "$run/shared.avinput" >> "$run/SHA256SUMS"

/usr/bin/time -f 'elapsed_seconds=%e\npeak_rss_kib=%M\nexit_code=%x' \
  -o "$run/perl.metrics" \
  perl "$annovar/table_annovar.pl" "$run/shared.avinput" "$database" \
  -buildver "$build" -out "$run/perl" -protocol "$protocol" -operation g \
  -nastring . -remove > "$run/perl.stdout" 2> "$run/perl.stderr"
test -s "$run/perl.${build}_multianno.txt"
sha256sum "$run/perl.${build}_multianno.txt" >> "$run/SHA256SUMS"

/usr/bin/time -f 'elapsed_seconds=%e\npeak_rss_kib=%M\nexit_code=%x' \
  -o "$run/rust.metrics" \
  "$binary" table "$run/shared.avinput" "$database" --build "$build" \
  --protocol "$protocol" --operation g --threads 16 --batch-size 10000 \
  --memory-budget 16 --output "$run/rust.tsv" \
  --report-json "$run/rust-report.json" \
  > "$run/rust.stdout" 2> "$run/rust.stderr"
test -s "$run/rust.tsv"
sha256sum "$run/rust.tsv" >> "$run/SHA256SUMS"

# Measure VCF conversion independently from shared-input annotation.
"$binary" convert "$input" --output "$run/rust-convert.avinput" \
  > "$run/rust-convert.stdout" 2> "$run/rust-convert.stderr"
sha256sum "$run/rust-convert.avinput" >> "$run/SHA256SUMS"
if python3 "$selection" compare-conversion "$run/rust-convert.avinput" \
  "$run/shared.avinput" --report "$run/conversion.json"; then
  printf 'CONVERSION_KEYS_AGREE\n' > "$run/conversion-status"
else
  printf 'CONVERSION_DIFFERENCES_RECORDED\n' > "$run/conversion-status"
fi
test -s "$run/conversion.json"

# Nonzero means a measured disagreement, not failure of the comparison job.
if python3 "$selection" compare "$run/rust.tsv" \
  "$run/perl.${build}_multianno.txt" --protocol "$protocol" \
  --report "$run/comparison.json"; then
  printf 'ALL_COMPARED_FIELDS_AGREE\n' > "$run/comparison-status"
else
  printf 'DIFFERENCES_RECORDED\n' > "$run/comparison-status"
fi
test -s "$run/comparison.json"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$run/finished.utc"
