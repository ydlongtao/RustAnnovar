#!/usr/bin/env bash
# Run after the matching CADD/GIAB downloads. Never clears global caches.
set -euo pipefail
root=${1:?usage: hpc_cadd_acceptance.sh ROOT BINARY BUILD NEW_REPORT_DIR}
binary=${2:?}
build=${3:?}
report=${4:?}
case "$build" in hg38) assembly=GRCh38 ;; hg19) assembly=GRCh37 ;; *) exit 2 ;; esac
cadd="$root/humandb/cadd-v1.7/$assembly"
# Wait for existing dependency tasks, then verify final files, not lock presence.
exec 8>"$cadd/download.lock"
flock 8
exec 7>"$root/runs/datasets.lock"
flock 7
input="$root/datasets/giab/HG002_${assembly}_1_22_v4.2.1_benchmark.vcf.gz"
test -f "$cadd/whole_genome_SNVs.tsv.gz" && test -f "$input" && test -x "$binary" || { echo 'verified final dependency files or binary missing' >&2; exit 42; }
exec 9>"$root/runs/benchmark.lock"
flock 9
mkdir "$report"
sha256sum "$binary" "$input" "$cadd/whole_genome_SNVs.tsv.gz.tbi" > "$report/SHA256SUMS"
cp "$cadd/MD5SUMs" "$report/CADD-official-MD5SUMs"
uname -a > "$report/host.txt"
lscpu > "$report/cpu.txt"
df -T "$cadd" > "$report/filesystem.txt"
# Select up to 16 physical cores on one NUMA node, within our permitted cpuset.
affinity=$(python3 - <<'PY'
import os, subprocess
allowed = os.sched_getaffinity(0)
seen = set()
selected = []
node = None
for row in subprocess.check_output(['lscpu','-p=CPU,CORE,SOCKET,NODE'],text=True).splitlines():
    if row.startswith('#'): continue
    cpu,core,socket,n = map(int,row.split(','))
    if cpu not in allowed: continue
    if node is None: node=n
    if n != node or (socket,core) in seen: continue
    seen.add((socket,core)); selected.append(cpu)
    if len(selected)==16: break
if len(selected)!=16: raise SystemExit('fewer than 16 available physical cores on one NUMA node')
print(','.join(map(str,selected)))
PY
)
printf '%s\n' "$affinity" > "$report/affinity.txt"
vmstat 1 > "$report/vmstat.log" &
monitor=$!
trap 'kill "$monitor" 2>/dev/null || true' EXIT
for run in first_run hot1 hot2 hot3 hot4 hot5; do
  cat /proc/loadavg > "$report/$run.load-before"
  cat /proc/vmstat > "$report/$run.vm-before"
  /usr/bin/time -f 'elapsed_seconds=%e\npeak_rss_kib=%M\nexit_code=%x' -o "$report/$run.metrics" \
    taskset -c "$affinity" "$binary" annotate "$input" "$cadd/whole_genome_SNVs.tsv.gz" \
      --operation cadd --protocol cadd17 --vcf-input --cadd-build "$build" --cadd-version 1.7 \
      --cadd-require-all --threads 16 --batch-size 10000 --memory-budget 16 \
      --output "$report/$run.tsv" --report-json "$report/$run.json"
  sha256sum "$report/$run.tsv" >> "$report/output-SHA256SUMS"
  if [[ "$run" != first_run ]]; then cmp "$report/first_run.tsv" "$report/$run.tsv"; fi
  cat /proc/loadavg > "$report/$run.load-after"
  cat /proc/vmstat > "$report/$run.vm-after"
done
python3 - "$report" <<'PY'
import json, statistics, sys
from pathlib import Path
p=Path(sys.argv[1]); runs=[]
for name in ['first_run','hot1','hot2','hot3','hot4','hot5']:
    record=json.loads((p/f'{name}.json').read_text())
    counts=record['cadd_counts']['cadd17']
    if counts.get('scored',0)==0: raise SystemExit('no SNVs scored')
    if any(n for state,n in counts.items() if state not in ('scored','not_snv')): raise SystemExit('unscored SNVs present')
    metric=dict(line.split('=',1) for line in (p/f'{name}.metrics').read_text().splitlines())
    runs.append({'run':name,'elapsed_seconds':float(metric['elapsed_seconds']),'peak_rss_kib':int(metric['peak_rss_kib']),'cadd_counts':counts})
hot=[r['elapsed_seconds'] for r in runs[1:]]
(p/'summary.json').write_text(json.dumps({'scope':'GIAB HG002 v4.2.1 chromosomes 1-22; strict SNV coverage, other alleles excluded from scoring','cache':'first run uncontrolled; five subsequent runs; global cache untouched','threads':16,'runs':runs,'subsequent_median_seconds':statistics.median(hot),'subsequent_range_seconds':[min(hot),max(hot)]},indent=2)+'\n')
PY
