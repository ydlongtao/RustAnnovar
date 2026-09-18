# CADD query optimization profile

This is a small regional query profile, not a whole-genome benchmark or a claim of complete SNV coverage.

## Workload and method

On 18 September 2026, we scored 310 canonical SNV alleles from the first megabase of chromosome 1 in GIAB HG002 v4.2.1 GRCh38. The database was a 64 MiB prefix of the official CADD v1.7 GRCh38 score-only BGZF file, queried with the official full Tabix index. All 310 requested alleles were present and scored. This deliberately partial database must not be used to certify full-genome coverage.

Both Rust release builds ran on the same HPC host and storage, with 16 threads, strict SNV coverage, and the same input/database paths. Five subsequent runs alternated baseline and optimized builds under the exclusive benchmark lock. Global caches were not cleared; these are subsequent-run measurements, not controlled cold-cache measurements. Elapsed time and peak RSS were recorded with `/usr/bin/time`. All ten output files were compared byte for byte against the baseline output and matched.

| Build | Elapsed seconds (five runs) | Median | Peak RSS range (KiB) |
| --- | --- | --- | --- |
| `72f40c3`, repeated position queries | 3.74, 3.74, 3.74, 3.71, 3.72 | 3.74 s | 47,012–47,300 |
| `5499531`, batched Tabix bins | 0.25, 0.24, 0.25, 0.24, 0.25 | 0.25 s | 63,688–63,928 |

The median ratio is 14.96× for this workload. The optimized implementation queries each touched 16,384-base Tabix bin once per batch and uses the configured Rayon thread pool for local files. Remote HTTPS queries remain serial to avoid request fan-out. The extra parallel readers increase RSS in this small workload.

## Correctness and remaining acceptance

The optimized build also matched the independent official CADD API for 21 fixed public probe records in each of hg19 and hg38, including unsorted inputs and repeated records. Unit/integration tests cover bin boundaries, duplicate input preservation, exact REF/ALT matching, missing-score statuses, and one-thread/four-thread consistency.

Full CADD files for both assemblies are downloading with official MD5 verification. Full GIAB autosomal strict-SNV acceptance jobs are queued and serialized. Until those jobs finish, full-genome coverage, WGS throughput, and large-input memory behavior remain unverified. These timings compare two Rust CADD implementations; they are not a Perl ANNOVAR comparison.

Raw profiling files are retained on the HPC under `/DATABANK/users/hflt/RustAnnovar/runs/cadd-query-profile-20260918/`. Large databases and variant files are not distributed in this repository.
