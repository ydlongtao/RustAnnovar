# Full-file CADD validation on GIAB HG002

The native local-Tabix CADD path passed a strict, autosomal GIAB HG002
GRCh38 run using the complete official CADD v1.7 score-only SNV database.
**GRCh37/hg19 validation remains pending.** This is an open-beta research
validation, not clinical validation or a comparison with Perl ANNOVAR.

## GRCh38 result

The official `whole_genome_SNVs.tsv.gz` download was published locally only
after its MD5 matched the official `MD5SUMs`. The adjacent official Tabix index
was also checksum-verified. The workload was GIAB HG002 v4.2.1 GRCh38
chromosomes 1–22, not all possible genomic positions or a representative
unfiltered WGS call set. The immutable Rust executable was built from
`5499531`. The command used `--operation cadd --protocol cadd17`,
`--cadd-build hg38 --cadd-version 1.7 --cadd-require-all`, 16 threads,
10,000-record batches, and a 16 GiB managed memory budget.

| Check | Result |
| --- | ---: |
| Original VCF records processed | 4,048,342 |
| Output alleles / TSV rows | 4,096,123 |
| Alleles with `CADD_status.cadd17=scored` | 3,463,511 |
| Alleles with `CADD_status.cadd17=not_snv` | 632,612 |
| Unscored SNVs in the JSON report or TSV | 0 |
| First run wall time | 226.88 s |
| Five subsequent runs, median (range) | 229.72 s (228.61–230.22 s) |
| Peak RSS across six runs | 395,944–404,260 KiB |

All six runs exited successfully and produced byte-identical TSV output,
SHA-256 `7487d7717ae462d3488307b9b40a034cc5089a82659a5800643932b2c0389e96`.
An independent streaming read of all 4,096,123 output rows confirmed the
status counts above, populated RawScore and PHRED fields on every `scored`
row, and empty score fields on every `not_snv` row. Strict mode returned
success because no SNV was unscored. Non-SNV alleles are outside CADD's
score-only SNV scope. Input alleles with ambiguous or unsupported bases
cannot be assigned a CADD SNV score and are reported separately by the tool.

The five subsequent runs are **not controlled cold-cache measurements**:
global caches were left intact on the shared HPC host. These timings include
input parsing, local score queries, and TSV writing; they exclude database
download and its MD5 validation. No speedup over the Perl implementation is
implied. Reproducibility files, including the exact command, binary/input/index
SHA-256 checksums, six JSON reports, `/usr/bin/time` metrics, output hashes,
host data, and load snapshots, are retained at
`/DATABANK/users/hflt/RustAnnovar/reports/cadd-hg38-giab-5499531-bounded/`.
The original database and GIAB VCF are not redistributed.

## Remaining gate

The matching GRCh37/hg19 official score file must finish downloading and
pass its official MD5 before the queued full GIAB HG002 v4.2.1 GRCh37
strict-SNV acceptance can run. A connection timeout stopped one transfer
after preserving its completed chunks; resumable downloading and a fresh
acceptance job were started from those chunks. Do not use the GRCh38 result
as evidence of GRCh37 coverage.
