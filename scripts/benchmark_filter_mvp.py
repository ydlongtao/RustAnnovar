#!/usr/bin/env python3
"""Serial portable release comparison. Synthetic AVinput, not a WGS claim."""
import argparse
import csv
import hashlib
import itertools
import json
import os
import pathlib
import platform
import re
import statistics
import subprocess
import time


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--directory', type=pathlib.Path, required=True)
    parser.add_argument('--native', type=pathlib.Path, default=pathlib.Path('target/release/rustannovar'))
    parser.add_argument('--baseline', type=pathlib.Path, default=pathlib.Path('target/release/rust-annovar'))
    parser.add_argument('--sizes', type=int, nargs='+', default=[10000, 100000, 1000000])
    parser.add_argument('--threads', type=int, nargs='+', default=[1, 16])
    parser.add_argument('--repeats', type=int, default=5)
    args = parser.parse_args()
    if min(args.sizes + args.threads + [args.repeats]) <= 0:
        parser.error('sizes, threads and repeats must be positive')
    root = args.directory.resolve()
    root.mkdir(parents=True, exist_ok=False)
    bins = {'baseline': args.baseline.resolve(), 'native': args.native.resolve()}
    env = {'platform': platform.platform(), 'machine': platform.machine(),
           'cpu': subprocess.check_output(['sysctl', '-n', 'machdep.cpu.brand_string'], text=True).strip() if platform.system() == 'Darwin' else platform.processor(),
           'logical_cpus': os.cpu_count(), 'compiler': subprocess.check_output(['rustc', '-Vv'], text=True),
           'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
           'source_dirty': bool(subprocess.check_output(['git', 'status', '--porcelain'], text=True)),
           'worktree_diff_sha256': hashlib.sha256(subprocess.check_output(['git', 'diff', 'HEAD'])).hexdigest(),
           'binaries': {k: {'path': str(v), 'sha256': sha(v)} for k, v in bins.items()},
           'settings': {'profile': 'portable release, no target-cpu=native', 'batch_size': 100000,
                        'cache': 'first run recorded separately, five warm runs by default; caches not cleared',
                        'database': '100K in-memory text rows; no disk index for either path',
                        'thread_affinity': 'not pinned; workstation scheduler',
                        'repeats': args.repeats}}
    db = root / 'synthetic.txt'
    with db.open('w') as f:
        f.write('#Chr\tStart\tEnd\tRef\tAlt\tvalue\n')
        for pos in range(100000, 0, -1):
            f.write(f'{pos % 22 + 1}\t{pos}\t{pos}\tA\tC\tvalue{pos}\n')
    env['database_sha256'] = sha(db)
    (root / 'environment.json').write_text(json.dumps(env, indent=2) + '\n')
    rows = []
    for size in args.sizes:
        data = root / f'{size}.avinput'
        with data.open('w') as f:
            for i in range(size):
                pos = i * 7919 % 100000 + 1
                alt = 'G' if i % 10 == 0 else 'C'
                f.write(f'{pos % 22 + 1}\t{pos}\t{pos}\tA\t{alt}\n')
        for threads in args.threads:
            # Alternate order between repetitions to reduce ordering bias.
            for repeat in range(args.repeats + 1):
                for name in (['baseline', 'native'] if repeat % 2 == 0 else ['native', 'baseline']):
                    out = root / f'{size}-{threads}-{name}.tsv'
                    report = root / f'{size}-{threads}-{name}-{repeat}.json'
                    cmd = [str(bins[name]), 'annotate']
                    cmd += ['--input', str(data), '--database', str(db)] if name == 'native' else [str(data), str(db), '--protocol', 'value']
                    cmd += ['--output', str(out), '--threads', str(threads), '--batch-size', '100000', '--report-json', str(report)]
                    timer = ['/usr/bin/time', '-lp'] if platform.system() == 'Darwin' else ['/usr/bin/time', '-f', 'BENCH %e %M %U %S']
                    before = os.getloadavg(); start = time.perf_counter()
                    run = subprocess.run(timer + cmd, text=True, capture_output=True)
                    elapsed = time.perf_counter() - start
                    prefix = root / f'{size}-{threads}-{name}-{repeat}'
                    prefix.with_suffix('.stderr').write_text(run.stderr)
                    if run.returncode:
                        raise RuntimeError(f'{cmd}: {run.stderr}')
                    if platform.system() == 'Darwin':
                        m = re.search(r'(\d+)\s+maximum resident set size', run.stderr)
                        rss = int(m[1]) / 1024 if m else None
                        user = re.search(r'^user\s+([\d.]+)', run.stderr, re.M)
                        system = re.search(r'^sys\s+([\d.]+)', run.stderr, re.M)
                        cpu = float(user[1]) + float(system[1]) if user and system else None
                    else:
                        m = re.search(r'BENCH ([\d.]+) (\d+) ([\d.]+) ([\d.]+)', run.stderr)
                        rss = int(m[2]) if m else None
                        cpu = float(m[3]) + float(m[4]) if m else None
                    row = {'engine': name, 'size': size, 'threads': threads, 'repeat': repeat,
                           'wall_seconds': elapsed, 'cpu_seconds': cpu, 'rss_kib': rss,
                           'load_before': before, 'load_after': os.getloadavg(), 'command': cmd,
                           'input_sha256': sha(data), 'stages': json.loads(report.read_text())}
                    rows.append(row)
                    with (root / 'measurements.jsonl').open('a') as f:
                        f.write(json.dumps(row) + '\n')
            with (root / f'{size}-{threads}-baseline.tsv').open() as a, (root / f'{size}-{threads}-native.tsv').open() as b:
                for number, (x, y) in enumerate(itertools.zip_longest(csv.reader(a, delimiter='\t'), csv.reader(b, delimiter='\t')), 1):
                    if x is None or y is None or x[:-1] != y:
                        raise RuntimeError(f'Output mismatch at line {number}')
    summary = []
    for size, threads, name in itertools.product(args.sizes, args.threads, bins):
        group = [r for r in rows if (r['size'], r['threads'], r['engine']) == (size, threads, name) and r['repeat'] > 0]
        median = statistics.median(r['wall_seconds'] for r in group)
        summary.append({'engine': name, 'size': size, 'threads': threads, 'median_seconds': median,
                        'min_seconds': min(r['wall_seconds'] for r in group), 'max_seconds': max(r['wall_seconds'] for r in group),
                        'peak_rss_kib': max((r['rss_kib'] for r in group if r['rss_kib'] is not None), default=None),
                        'variants_per_second': size / median})
    (root / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (root / 'SUCCESS').write_text('Output equality checked excluding baseline AnnotationStatus. Synthetic filter scope only.\n')
    print(root)


if __name__ == '__main__':
    main()
