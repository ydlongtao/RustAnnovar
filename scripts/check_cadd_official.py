#!/usr/bin/env python3
"""Check public probe SNVs against independent official CADD API scores.

This sends only the fixed public probes below, never a user's variant file.
No genomic data or large score databases are distributed with this script.
"""
import argparse
import csv
import json
import subprocess
import urllib.request
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    results = []
    for build, assembly in [("hg19", "GRCh37"), ("hg38", "GRCh38")]:
        expected = []
        urls = []
        for locus in ["1:10001", "2:1000000", "7:117603632", "22:20000000", "X:10000000", "Y:10000000"]:
            url = f"https://cadd.kircherlab.bihealth.org/api/v1.0/{assembly}-v1.7/{locus}"
            with urllib.request.urlopen(url, timeout=60) as response:
                scores = json.load(response)
            if not scores:
                raise RuntimeError(f"official probe returned no scores: {url}")
            for row in scores:
                expected.append(row)
            urls.append(url)
        # Deliberate unsorted duplicate inputs, including all three ALT alleles.
        expected = list(reversed(expected)) + expected[:3]
        input_path = args.output_dir / f"{build}.avinput"
        input_path.write_text("".join(f'{r["Chrom"]}\t{r["Pos"]}\t{r["Pos"]}\t{r["Ref"]}\t{r["Alt"]}\n' for r in expected))
        source = f"https://krishna.gs.washington.edu/download/CADD/v1.7/{assembly}/whole_genome_SNVs.tsv.gz"
        output = args.output_dir / f"{build}.tsv"
        report = args.output_dir / f"{build}.json"
        command = [str(args.binary.resolve()), "annotate", str(input_path), source, "--operation", "cadd", "--protocol", "cadd17", "--cadd-build", build, "--cadd-require-all", "--output", str(output), "--report-json", str(report)]
        subprocess.run(command, check=True)
        with output.open() as handle:
            observed = list(csv.DictReader(handle, delimiter="\t"))
        assert len(observed) == len(expected)
        for actual, oracle in zip(observed, expected):
            assert actual["CADD_raw.cadd17"] == oracle["RawScore"], (actual, oracle)
            assert actual["CADD_phred.cadd17"] == oracle["PHRED"], (actual, oracle)
            assert actual["CADD_status.cadd17"] == "scored"
        evidence = {"build": build, "alleles": len(expected), "api_urls": urls, "source": source, "expected": expected, "command": command}
        (args.output_dir / f"{build}-oracle.json").write_text(json.dumps(evidence, indent=2) + "\n")
        results.append({"build": build, "alleles": len(expected), "matched": True})
        print(f"{build}: verified {len(expected)} public probe scores", flush=True)
    (args.output_dir / "summary.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
