#!/usr/bin/env python3
"""Deterministic GIAB VCF sampling and unfiltered ANNOVAR table comparison."""

import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path


FIELD_NAMES = ("Func", "Gene", "GeneDetail", "ExonicFunc", "AAChange")
KEY_FIELDS = ("Chr", "Start", "End", "Ref", "Alt")
RATES = {"snv": 100, "indel": 300, "mnv": 1000, "multiallelic": 1000, "other": 1000}


def input_file(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def variant_type(reference, alternate):
    alleles = alternate.split(",")
    if len(alleles) > 1:
        return "multiallelic"
    alt = alleles[0]
    if not set(reference + alt) <= set("ACGT"):
        return "other"
    if len(reference) == len(alt) == 1:
        return "snv"
    if reference and alt and len(reference) == len(alt):
        return "mnv"
    return "indel"


def sample(args):
    counts = collections.Counter()
    selected = collections.Counter()
    chroms = collections.Counter()
    with input_file(args.input) as source, open(args.output, "w") as target:
        for line in source:
            if line.startswith("#"):
                target.write(line)
                continue
            parts = line.split("\t", 5)
            if len(parts) < 5:
                raise ValueError("malformed VCF row")
            chrom, pos, _, ref, alt = parts[:5]
            category = variant_type(ref.upper(), alt.upper())
            counts[category] += 1
            identity = "\t".join((args.seed, chrom, pos, ref, alt)).encode()
            value = int.from_bytes(hashlib.blake2b(identity, digest_size=8).digest(), "big")
            if value % 10000 < RATES[category]:
                target.write(line)
                selected[category] += 1
                chroms[chrom] += 1
    result = {
        "input": str(args.input),
        "source_sha256": file_sha256(args.input),
        "output_sha256": file_sha256(args.output),
        "seed": args.seed,
        "denominator": 10000,
        "rates": RATES,
        "source_records_by_type": dict(counts),
        "selected_records_by_type": dict(selected),
        "selected_records_by_chromosome": dict(sorted(chroms.items())),
    }
    write_json(args.report, result)


def file_sha256(path):
    value = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def read_table(path, fields):
    rows = {}
    occurrences = collections.Counter()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"missing header: {path}")
        missing = set(KEY_FIELDS + fields) - set(reader.fieldnames)
        if missing:
            raise ValueError(f"missing columns in {path}: {sorted(missing)}")
        for number, row in enumerate(reader, 1):
            key = tuple(row[field] for field in KEY_FIELDS)
            ordinal = occurrences[key]
            occurrences[key] += 1
            identity = key + (ordinal,)
            if identity in rows:
                raise ValueError(f"duplicate identity in {path}: {identity}")
            rows[identity] = (number, row)
    return rows


def compare(args):
    protocol = getattr(args, "protocol", "refGene")
    fields = tuple(f"{name}.{protocol}" for name in FIELD_NAMES)
    rust = read_table(args.rust, fields)
    perl = read_table(args.perl, fields)
    common = rust.keys() & perl.keys()
    field_counts = {field: {"exact": 0, "order_insensitive": 0, "different": 0} for field in fields}
    nonmissing_perl = {field: {"total": 0, "exact": 0, "order_insensitive": 0, "different": 0} for field in fields}
    category_counts = collections.Counter()
    category_mismatches = collections.Counter()
    category_rows_exact = collections.Counter()
    category_rows_equivalent = collections.Counter()
    examples = {field: [] for field in fields}
    rows_all_exact = 0
    rows_all_order_insensitive = 0
    for identity in sorted(common):
        _, a = rust[identity]
        _, b = perl[identity]
        category = variant_type(identity[3].replace("-", ""), identity[4].replace("-", ""))
        category_counts[category] += 1
        exact_row = True
        equivalent_row = True
        for field in fields:
            av, bv = a[field], b[field]
            nonmissing = bv not in ("", ".")
            if nonmissing:
                nonmissing_perl[field]["total"] += 1
            if av == bv:
                field_counts[field]["exact"] += 1
                if nonmissing:
                    nonmissing_perl[field]["exact"] += 1
                continue
            exact_row = False
            if field in (f"Gene.{protocol}", f"AAChange.{protocol}") and sorted(av.split(",")) == sorted(bv.split(",")):
                field_counts[field]["order_insensitive"] += 1
                if nonmissing:
                    nonmissing_perl[field]["order_insensitive"] += 1
                continue
            equivalent_row = False
            field_counts[field]["different"] += 1
            if nonmissing:
                nonmissing_perl[field]["different"] += 1
            category_mismatches[(category, field)] += 1
            if len(examples[field]) < args.max_examples:
                examples[field].append({"variant": list(identity), "rust": av, "perl": bv})
        rows_all_exact += exact_row
        rows_all_order_insensitive += equivalent_row
        category_rows_exact[category] += exact_row
        category_rows_equivalent[category] += equivalent_row
    missing_in_rust = perl.keys() - rust.keys()
    missing_in_perl = rust.keys() - perl.keys()
    result = {
        "rust_table_sha256": file_sha256(args.rust),
        "perl_table_sha256": file_sha256(args.perl),
        "protocol": protocol,
        "rust_rows": len(rust),
        "perl_rows": len(perl),
        "matched_variant_occurrences": len(common),
        "rows_all_fields_exact": rows_all_exact,
        "rows_all_fields_order_insensitive": rows_all_order_insensitive,
        "missing_in_rust": len(missing_in_rust),
        "missing_in_perl": len(missing_in_perl),
        "missing_in_rust_examples": [list(v) for v in list(sorted(missing_in_rust))[:args.max_examples]],
        "missing_in_perl_examples": [list(v) for v in list(sorted(missing_in_perl))[:args.max_examples]],
        "field_counts": field_counts,
        "field_counts_when_perl_nonmissing": nonmissing_perl,
        "alleles_by_type": dict(sorted(category_counts.items())),
        "rows_all_fields_exact_by_type": dict(sorted(category_rows_exact.items())),
        "rows_all_fields_order_insensitive_by_type": dict(sorted(category_rows_equivalent.items())),
        "different_by_type_and_field": {
            f"{category}/{field}": count for (category, field), count in sorted(category_mismatches.items())
        },
        "examples": examples,
        "comparison_policy": "Match Chr/Start/End/Ref/Alt plus duplicate occurrence; report every shared field. Gene and AAChange comma-list order is counted separately; no semantic differences are suppressed.",
    }
    write_json(args.report, result)
    return 0 if not missing_in_rust and not missing_in_perl and not any(x["different"] for x in field_counts.values()) else 1


def compare_conversion(args):
    def read_keys(path):
        keys = []
        with open(path) as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise ValueError(f"{path}:{number}: fewer than five AVinput columns")
                keys.append(tuple(fields[:5]))
        return keys

    rust = read_keys(args.rust)
    perl = read_keys(args.perl)
    a, b = collections.Counter(rust), collections.Counter(perl)
    only_rust, only_perl = a - b, b - a

    def normalize_chrom(key):
        chrom = key[0]
        if chrom.lower().startswith("chr"):
            chrom = chrom[3:]
        return (chrom,) + key[1:]

    alias_rust = collections.Counter(normalize_chrom(key) for key in rust)
    alias_perl = collections.Counter(normalize_chrom(key) for key in perl)
    alias_only_rust, alias_only_perl = alias_rust - alias_perl, alias_perl - alias_rust

    def av_type(key):
        ref, alt = key[3:5]
        if ref == "-" or alt == "-":
            return "indel"
        if set(ref + alt) <= set("ACGT"):
            if len(ref) == len(alt) == 1:
                return "snv"
            if len(ref) == len(alt):
                return "mnv"
            return "complex"
        return "symbolic_or_ambiguous"

    def types(counter):
        total = collections.Counter()
        for key, n in counter.items():
            total[av_type(key)] += n
        return dict(sorted(total.items()))

    result = {
        "rust_sha256": file_sha256(args.rust),
        "perl_sha256": file_sha256(args.perl),
        "rust_rows": len(rust),
        "perl_rows": len(perl),
        "ordered_equal": rust == perl,
        "multiset_equal": a == b,
        "multiset_equal_after_chr_prefix_normalization": alias_rust == alias_perl,
        "only_rust": sum(only_rust.values()),
        "only_perl": sum(only_perl.values()),
        "only_rust_after_chr_prefix_normalization": sum(alias_only_rust.values()),
        "only_perl_after_chr_prefix_normalization": sum(alias_only_perl.values()),
        "shared_by_type": types(a & b),
        "only_rust_by_type": types(only_rust),
        "only_perl_by_type": types(only_perl),
        "only_rust_after_chr_prefix_normalization_by_type": types(alias_only_rust),
        "only_perl_after_chr_prefix_normalization_by_type": types(alias_only_perl),
        "only_rust_examples": [list(v) + [n] for v, n in only_rust.most_common(args.max_examples)],
        "only_perl_examples": [list(v) + [n] for v, n in only_perl.most_common(args.max_examples)],
        "only_rust_after_chr_prefix_normalization_examples": [list(v) + [n] for v, n in alias_only_rust.most_common(args.max_examples)],
        "only_perl_after_chr_prefix_normalization_examples": [list(v) + [n] for v, n in alias_only_perl.most_common(args.max_examples)],
        "policy": "Compare all five AVinput columns and duplicate multiplicity. Report exact and chr-prefix-normalized multisets separately; no position or allele differences are normalized.",
    }
    write_json(args.report, result)
    return 0 if a == b else 1


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    sampler = sub.add_parser("sample")
    sampler.add_argument("input", type=Path)
    sampler.add_argument("output", type=Path)
    sampler.add_argument("--report", type=Path, required=True)
    sampler.add_argument("--seed", default="RustAnnovar-GIAB-v1")
    comparator = sub.add_parser("compare")
    comparator.add_argument("rust", type=Path)
    comparator.add_argument("perl", type=Path)
    comparator.add_argument("--report", type=Path, required=True)
    comparator.add_argument("--protocol", default="refGene")
    comparator.add_argument("--max-examples", type=int, default=10)
    conversion = sub.add_parser("compare-conversion")
    conversion.add_argument("rust", type=Path)
    conversion.add_argument("perl", type=Path)
    conversion.add_argument("--report", type=Path, required=True)
    conversion.add_argument("--max-examples", type=int, default=10)
    args = parser.parse_args()
    if args.command == "sample":
        sample(args)
        return 0
    if args.command == "compare":
        return compare(args)
    return compare_conversion(args)


if __name__ == "__main__":
    raise SystemExit(main())
