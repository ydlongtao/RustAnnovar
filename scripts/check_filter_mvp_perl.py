#!/usr/bin/env python3
"""Black-box generic filter comparison using original synthetic fixtures only."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def canonical(values):
    chrom = values[0].removeprefix('chr')
    return (chrom, *values[1:])


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--annovar', type=Path, required=True)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--native', type=Path, default=Path('target/release/rustannovar'))
    args = p.parse_args()
    root = args.directory.resolve()
    root.mkdir(parents=True, exist_ok=False)
    dbdir = root / 'humandb'
    dbdir.mkdir()
    fixtures = Path(__file__).resolve().parents[1] / 'crates/rustannovar-cli/tests/fixtures'
    data = root / 'input.avinput'
    db = dbdir / 'hg19_mvp.txt'
    shutil.copyfile(fixtures / 'exact.avinput', data)
    shutil.copyfile(fixtures / 'exact-db.txt', db)
    prefix = root / 'perl'
    commands = [
        ['perl', str(args.annovar.resolve() / 'annotate_variation.pl'), '-filter', '-dbtype', 'generic', '-buildver', 'hg19', '-genericdbfile', db.name, '-outfile', str(prefix), str(data), str(dbdir)],
        [str(args.native.resolve()), 'annotate', '--input', str(data), '--database', str(db), '--output', str(root / 'native.tsv'), '--threads', '4'],
    ]
    for i, cmd in enumerate(commands):
        r = subprocess.run(cmd, text=True, capture_output=True)
        (root / f'command-{i}.stderr').write_text(r.stderr)
        (root / f'command-{i}.stdout').write_text(r.stdout)
        if r.returncode:
            raise RuntimeError(f'Command exited {r.returncode}: {r.stderr}')
    with (root / 'native.tsv').open() as f:
        rows = list(csv.reader(f, delimiter='\t'))[1:]
    native_hits = Counter(canonical(row[:5]) + (row[5],) for row in rows if row[5] != '.')
    native_missing = Counter(canonical(row[:5]) for row in rows if row[5] == '.')
    with Path(str(prefix) + '.hg19_generic_dropped').open() as f:
        expected_hits = Counter(canonical(row[2:7]) + (row[1],) for row in csv.reader(f, delimiter='\t'))
    with Path(str(prefix) + '.hg19_generic_filtered').open() as f:
        expected_missing = Counter(canonical(row[:5]) for row in csv.reader(f, delimiter='\t'))
    input_rows = [line.split('\t') for line in data.read_text().splitlines()]
    order_ok = [row[:5] for row in rows] == input_rows
    differences = {name: [{'record': key, 'count': count} for key, count in values.items()] for name, values in {
        'native_only_hits': native_hits - expected_hits,
        'perl_only_hits': expected_hits - native_hits,
        'native_only_missing': native_missing - expected_missing,
        'perl_only_missing': expected_missing - native_missing,
    }.items()}
    report = {'passed': order_ok and not any(differences.values()), 'records': len(rows),
              'matched_records': sum(native_hits.values()), 'unmatched_records': sum(native_missing.values()),
              'input_order_preserved': order_ok, 'differences': differences, 'commands': commands,
              'sha256': {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in [data, db, args.native, args.annovar / 'annotate_variation.pl']},
              'scope': 'Synthetic generic filter: SNV, MNV, insertion, deletion, X/Y, chr prefix, repeated input, alternate mismatch. No real database release certification or conflicting duplicate database rows.'}
    (root / 'comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ['passed', 'records', 'matched_records', 'unmatched_records']}))
    raise SystemExit(0 if report['passed'] else 1)


if __name__ == '__main__':
    main()
