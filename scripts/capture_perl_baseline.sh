#!/usr/bin/env bash
set -euo pipefail

: "${ANNOVAR_HOME:?Set ANNOVAR_HOME to the registered ANNOVAR installation directory}"

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
baseline_dir="$project_root/tests/fixtures/perl_expected"

test -f "$ANNOVAR_HOME/table_annovar.pl"
mkdir -p "$baseline_dir"

(
  cd "$ANNOVAR_HOME"
  shasum -a 256 ./*.pl
) > "$baseline_dir/SHA256SUMS"

perl "$ANNOVAR_HOME/table_annovar.pl" \
  "$ANNOVAR_HOME/example/ex1.avinput" \
  "$ANNOVAR_HOME/humandb" \
  -buildver hg19 \
  -out "$baseline_dir/ex1" \
  -protocol refGene \
  -operation g \
  -nastring . \
  -remove

echo "Captured baseline under $baseline_dir"
