# CADD SNV scores (integration in progress)

## Native indexed scoring

The `cadd` operation queries an official six-column score-only BGZF file
and its adjacent `.tbi` directly in Rust, without conversion or full loading.
Use `--cadd-build hg19|hg38` with `annotate` to verify the official assembly
header. `table --build` performs this verification for the `cadd` operation.
A missing or mismatching assembly declaration fails. This checks metadata,
not the integrity or completeness of the bulk download. `--cadd-version 1.7`
checks the official version declaration. Reports include `cadd_databases`
with source, assembly, version and local/remote mode.

```bash
./target/release/rust-annovar annotate input.vcf /data/whole_genome_SNVs.tsv.gz \
  --operation cadd --protocol cadd17 --vcf-input --output scores.tsv \
  --cadd-build hg38 --cadd-version 1.7 --cadd-require-all --report-json coverage.json
```

This path adds `CADD_status.cadd17`: `scored`, `reference_mismatch`,
`score_not_found`, `contig_not_covered`, `not_snv`, or `unsupported_snv` (ambiguous N). Every processed allele
receives a state. Missing scores are not invented. `--cadd-require-all` fails
if any SNV is unscored, including ambiguous N alleles and preserves existing file outputs; stdout
cannot be rolled back. The JSON report contains per-protocol `cadd_counts`
by state. Non-SNV alleles are reported but do not fail this SNV-only policy. Each touched 16,384-base Tabix bin is queried once per batch. Local bins
use the configured Rayon pool; remote queries remain sequential with a bounded
reader cache. Full WGS throughput remains unmeasured.


## Explicit remote mode

A HTTPS URL may replace the local score-file path. This downloads the small
Tabix index and requested compressed byte ranges, using a bounded 4 MiB
range cache per reader. HTTP 206 responses, file length, and ETag or
Last-Modified are checked; network failures return errors. This accesses
the official full-genome file, not a limited preselected variant database.

```bash
./target/release/rust-annovar annotate input.vcf \
  https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz \
  --operation cadd --protocol cadd17 --vcf-input --output scores.tsv \
  --cadd-build hg38 --cadd-version 1.7 --cadd-require-all --report-json coverage.json
```

Remote mode sends byte-range requests to the selected data host. It does
not upload the VCF or sample columns, but requested ranges can disclose
approximate genomic locations to that host. Use local mode for private
workflows. Remote mode is suited to small jobs and evaluation; full WGS
should use a local copy to avoid network latency. It is explicitly selected
by passing a URL and is never an automatic fallback for missing scores.

For combined native `table` annotation, symlink the official local files
as `humandb/hg38_cadd17.txt.gz` and `humandb/hg38_cadd17.txt.gz.tbi`, and use
`--protocol refGene,cadd17 --operation g,cadd --build hg38`. This includes
CADD status and coverage statistics; the converted `f` path does not.

On HPC, `bash scripts/hpc_cadd_download.sh /DATABANK/users/hflt/RustAnnovar hg38`
(or `hg19`) downloads resumably under `humandb/cadd-v1.7/GRCh38` (or GRCh37),
verifies both files against the official MD5SUMs, and publishes final filenames
only after verification. Existing files are checked rather than overwritten.

For faster HPC transfers, the bounded four-connection downloader preserves
and reuses existing single-stream prefixes, resumes individual chunks, checks
HTTP ranges and source identity, and verifies the merged file against official
MD5 before exposing it:

```bash
python3 scripts/hpc_cadd_parallel_download.py \
  --root /DATABANK/users/hflt/RustAnnovar --build hg38 --workers 4
python3 scripts/test_cadd_download.py
```

Chunks and the original partial file are retained for recovery. Their presence
does not mean a completed or verified database exists. This method temporarily
needs about twice the compressed database size for chunks and the merge.

## Converted filter scoring

RustAnnovar can import official precomputed CADD tables and annotate exact
chromosome/position/reference/alternate matches with `CADD_raw` and
`CADD_phred`. This reads published scores; it does not reproduce or retrain
the CADD prediction model.

Download the appropriate **score-only all-possible-SNV** file from the
[official CADD downloads](https://cadd.kircherlab.bihealth.org/download).
Use GRCh37 for hg19 or GRCh38 for hg38. The importer accepts plain TSV or
gzip TSV, locates columns by their official names, and also accepts
annotation-inclusive tables. The build and version flags are provenance
declarations by the operator, not automatic proof of the downloaded assembly.

```bash
cargo build --release --bin rust-annovar --locked
mkdir -p humandb
./target/release/rust-annovar db import-cadd whole_genome_SNVs.tsv.gz \
  humandb/hg38_cadd17.txt --build hg38 --version 1.7
./target/release/rust-annovar db index humandb/hg38_cadd17.txt --tmp-dir /path/to/tmp
./target/release/rust-annovar db validate humandb/hg38_cadd17.txt --full
./target/release/rust-annovar table input.vcf humandb --build hg38 \
  --protocol cadd17 --operation f --vcf-input --threads 16 \
  --output output.tsv --vcf-output output.vcf
```

TSV score columns are `CADD_raw.cadd17` and `CADD_phred.cadd17`.
VCF `FA_cadd17` contains RawScore followed by PHRED, separated by `|`,
with values grouped in original ALT order. Original sample fields are retained.
Import refuses existing destinations, malformed coordinates, noncanonical SNVs,
and nonfinite scores. A failed import leaves no new destination file.

**Current limitations:** conversion writes an expanded text file and the
subsequent external sort needs substantial temporary disk space. Full-genome
acceptance is not certified yet.
On this converted generic-filter path, missing matches use the normal filter missing value (`.` by default),
including reference mismatches, uncovered contigs, and non-SNVs. A successful
run therefore does not prove that every input SNV received a score. Dedicated
coverage/error reporting remains required for that goal. General
`AnnotationStatus=ok` means the row was processed, not that CADD was found.

Obtain and use CADD data under its own terms; the RustAnnovar software license
does not license the data. See the [CADD site](https://cadd.gs.washington.edu/)
and [CADD v1.7 publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC10767851/).
No CADD bulk data are distributed in this repository.

## Real-data verification

A range download of the first BGZF block of the official GRCh38 v1.7
score-only file yielded 2,409 complete allele records. Native lookup recovered
all original RawScore/PHRED strings exactly. A separate official GRCh37 v1.7
fragment with 2,397 complete alleles also passed. This is a real-data fragment check,
not a full-genome coverage or performance certificate. Reproduce with an
independently downloaded official score-only TSV:

```bash
CADD_TEST_TSV=/path/to/official.tsv cargo test --test cadd official_records \
  --locked -- --ignored --nocapture
```

Cross-chromosome native HTTP-range lookup also matched the independent
official API for 21 inputs per build (hg19 and hg38): chromosome 1, 2, 7,
22, X and Y, all three ALT bases, reversed input order, and duplicates.
This validates query semantics across the full-file address space; it does
not establish whole-WGS speed or exhaustive coverage of every genomic site.

```bash
python3 scripts/check_cadd_official.py --binary target/release/rust-annovar \
  --output-dir /path/to/new-verification-directory
```

## Query optimization measurements

See the [regional query profile](CADD-query-profile.md) for five-run timings and their limitations. Full WGS acceptance remains pending; the regional speedup must not be presented as a WGS speedup.

The parallel downloader stores 64 MiB resumable pieces but limits each HTTP Range request to 4 MiB. This reduces the duration of individual connections; truncated responses retain downloaded bytes and retry from the saved offset. Source identity checks and the final official MD5 check still apply. If a job exhausts its retries, inspect its terminal exit status before resubmitting; existing pieces are reused.

Only up to the configured worker count of piece tasks is submitted at once. A fatal piece error stops submission of further pieces; already active transfers finish before the error is reported. This avoids waiting for an entire genome of queued transfers after a failure.
