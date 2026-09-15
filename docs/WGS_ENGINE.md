# WGS engine development: 0.2.0-beta.1

This candidate adds streaming CLI annotation, compact batch ownership, interval indexes, and an immutable disk index for filter databases. It remains an open beta. Public WES/WGS compatibility and performance claims require the remote acceptance report; passing synthetic tests alone is insufficient.

See the [interim acceptance report](WGS_ACCEPTANCE.md) for measured regressions and outstanding release gates.

## Running the candidate

```bash
cargo build --release --locked
rust-annovar db index humandb/hg38_clinvar.txt --kind filter --tmp-dir /path/to/scratch
rust-annovar db validate humandb/hg38_clinvar.txt --full
rust-annovar table sample.vcf.gz humandb \
  --build hg38 --protocol refGene,cytoBand,clinvar --operation g,r,f \
  --vcf-input --output sample.tsv --vcf-output sample.annotated.vcf \
  --threads 16 --batch-size 10000 --memory-budget 16 \
  --tmp-dir /path/to/scratch --report-json sample.report.json
```

Use the actual installed database protocol names (including release suffixes). `.txt.gz` databases are also resolved. Filters larger than 64 MiB on disk require a new index. The new index builder consumes plain or gzip text, externally sorts records in bounded runs, and preserves duplicate source rows. It creates a `.rai.json` manifest and a content-addressed `*.rai` data file next to the database. Keep both files. Data generations are published before the manifest; interrupted rebuilds leave the previous generation usable. Obsolete generations are not automatically deleted.

The index uses zero-based internal coordinates, an explicit normalization policy, source SHA-256, per-block checksums, and a data SHA-256. Normal opens check source path, size, nanosecond modification time and the first MiB hash. Full validation hashes both complete files. Normal checks cannot detect an intentionally preserved timestamp plus a same-size modification outside the prefix: use `db validate --full` after moving or modifying data. Rebuild on another machine because source paths and modification times are validated.

Gene and region indexes are in-memory static interval trees. Transcript exon offsets and CDS bounds are computed at load time. Filter blocks have a per-database 64 MiB decoded-data cache; database blocks are read in order for each batch. `--memory-budget` is in GiB and limits input batch bytes; it is not a hard RSS ceiling. Gene models/transcript FASTA, result strings and aggregate database caches still need separate measurement. Default batch size is 10,000 complete input records and the byte cap is at most 64 MiB. A record is never split across batches. AVinput gets a schema prepass; stdin AVinput is temporarily spooled to disk.

## Coordinates, statuses and output changes

- `--normalize annovar` is the default. It retains existing ANNOVAR-compatible coordinate conversion and VCF anchor trimming.
- `--normalize left-align --reference genome.fa` validates REF and left-aligns sequence alleles using a samtools-compatible `genome.fa.fai`. Build filter indexes with the same two options. Queries reject a different reference hash or normalization policy. Reference hashing happens when opening this mode and is included in total time.
- CLI tables add `AnnotationStatus` before `Otherinfo` columns. Values currently include `ok`, `unsupported_alt`, and `coding_unknown`. GeneDetail records reasons such as missing sequence, reference mismatch or unsupported coding boundary.
- VCF output retains original alleles, GT, FORMAT and samples. `FA_<protocol>` and `FA_STATUS` use one entry per original ALT. Annotation fields are separated by `|`; reserved characters inside field values are percent-escaped. Unsupported ALT slots remain present. Existing `FA_` INFO headers are rejected to avoid collisions.
- `--unsupported error` fails on unsupported VCF ALT; the default `preserve` emits the record and status. This option does not promise structural-variant annotation.
- File outputs use temporary sibling files and rename only after successful processing. Each file is atomic individually; multiple output files are not a cross-file transaction. stdout cannot be rolled back.
- Existing Rust full-document helpers remain available. The bounded-memory and new VCF serialization guarantees apply to the streaming CLI path; library callers should use that path for large jobs.

## Gene consequences and compatibility boundaries

The candidate fixes ncRNA intronic/splice labels, uses intronic splice windows, and produces cDNA details for coding-transcript UTR/splice SNVs. Reference mismatches and incomplete CDS information are reported rather than assigned a confident protein consequence.

MNV and small insertion/deletion/replacement events completely inside one coding exon now calculate mutated coding sequence and protein changes on both strands. These outputs are preliminary consequence descriptions, not full HGVS certification or complete ANNOVAR indel-text compatibility. Cross-exon/CDS boundary events retain positional classification with `coding_unknown`. Haplotypes across VCF records are not combined.

Nearest-gene lookup preserves the baseline's deterministic first tied transcript/gene, rather than expanding the main Gene column to all tied names. This is an intentional compatibility adjustment to the proposed plan.

ClinVar and gnomAD ANNOVAR-format files use the generic exact-match engine and retain input fields. Dedicated version-specific protocol certification, numeric threshold options, and real-file validation are still pending. Raw ClinVar/gnomAD VCF is not interchangeable with an ANNOVAR text database. Parquet, complete HGVS and additional specialized gene-model adapters are deferred.

## Remote acceptance

Workspace: `/DATABANK/users/hflt/RustAnnovar/` on `gse234129-hpc`. Large data never enter Git. `scripts/hpc_datasets.sh` downloads fixed GIAB HG002 v4.2.1 GRCh37/GRCh38 inputs, checks gzip integrity and records checksums/URLs. The download script uses a lock and records an exit code.

`HPC benchmark binaries` builds baseline `deb7a25f1ae178ad8c6060629b1f75f19bec834d` and the candidate with the same Linux runner/toolchain. Artifacts include executable SHA-256, commit and compiler information. They can be staged under `builds/baseline/` and `builds/candidate/` without installing system-wide Rust.

```bash
python3 scripts/hpc_benchmark.py --mode smoke
python3 scripts/hpc_benchmark.py --mode scale
python3 scripts/hpc_benchmark.py --mode giab
```

The harness serializes jobs with a lock, selects distinct physical cores on one NUMA node, records machine/storage details and load/swap counters, and runs one first-run plus five warm repeats. Smoke tests use 100,000 synthetic variants at 1/2/4/8/16/32 threads. Scale tests use 1M/10M records at 16 threads. Synthetic inputs repeat 100,000 positions and are explicitly not representative WGS samples. Output row equality is checked after excluding the new status column.

GIAB mode records missing inputs as not validated. It times refGene annotation and Perl separately and writes a field-level comparison report using `compare_annotations.py`. It does not silently count a missing database as passed. No global caches are cleared. A large machine does not waive the 20 GiB peak-RSS target or the <=20% growth target from 1M to 10M records.

No 0.2 WGS speedup is claimed until these acceptance checks complete. No new public registry release should be inferred from the development version number.
