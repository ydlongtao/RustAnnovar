# Filter MVP performance

These measurements cover the experimental `rustannovar` AVinput/filter frontend only. They do not establish a speedup for VCF, large disk-indexed databases, gene annotation, real WES/WGS, or Perl ANNOVAR. The previous WGS candidate's failed performance and compatibility gates remain documented in [WGS_ACCEPTANCE.md](docs/WGS_ACCEPTANCE.md).

## Method and provenance

Hardware: Apple M1, 8 logical CPUs; macOS-26.5-arm64-arm-64bit-Mach-O. Compiler: `rustc 1.96.0 (ac68faa20 2026-05-25)`. Both executables use portable `cargo build --workspace --release`, without native CPU tuning. 16 requested threads oversubscribe this workstation; scheduling is not pinned.

The baseline is the existing `rust-annovar` frontend at source revision `354456b`, rebuilt with the same toolchain and release settings. The new workspace implementation was measured from the working tree based on that commit; all recorded native Rust source hashes match implementation commit `1091dff8cbee57b1dfbdbe87b776e6a32a14edae`. Exact native source-file hashes and executable hashes are retained; the metadata explicitly marks that working tree as dirty. Do not confuse these results with binaries from the earlier remote WGS campaign.

The text database contains 100,000 rows across 22 chromosomes, deliberately unsorted, with one annotation field. Inputs contain 10K, 100K or 1M records with permuted positions, repeated records at 1M and 10% absent ALT matches. Both frontends use a 100,000-record batch size. No persistent disk index is used. Startup, database load/index construction, parsing, querying and file output are included. The baseline includes its AVinput schema prepass and compatibility status field; the native MVP validates extra-column width while streaming and does not generate that status field. This is a narrower frontend comparison, not an isolated comparison of hash maps versus binary search.

Each combination runs once as a first-run observation and five further times for reported warm statistics. Baseline/native execution order alternates. Caches are not cleared. Peak RSS comes from `/usr/bin/time -l`; wall time is measured by the parent process, including the timer wrapper. CPU time, load snapshots and per-stage timings are retained. Every output record and header is compared after excluding only the baseline AnnotationStatus column; the baseline protocol is explicitly named `value` to align headers. No missing-value or annotation differences are ignored.

## Results

| Records | Threads | Baseline seconds | MVP seconds (min–max) | Speedup | Baseline / MVP peak MiB | MVP variants/s |
|---|---|---|---|---|---|---|
| 10,000 | 1 | 0.1080 | 0.0368 (0.0358–0.0371) | 2.93× | 54.9 / 30.5 | 271,623 |
| 10,000 | 16 | 0.1085 | 0.0359 (0.0354–0.0373) | 3.02× | 55.4 / 31.3 | 278,286 |
| 100,000 | 1 | 0.3009 | 0.0829 (0.0799–0.0866) | 3.63× | 117.1 / 46.4 | 1,206,110 |
| 100,000 | 16 | 0.2490 | 0.0734 (0.0723–0.0767) | 3.39× | 117.9 / 48.4 | 1,361,884 |
| 1,000,000 | 1 | 2.1261 | 0.5633 (0.5360–0.5957) | 3.77× | 201.7 / 46.4 | 1,775,192 |
| 1,000,000 | 16 | 1.6206 | 0.4636 (0.4319–0.6111) | 3.50× | 202.0 / 50.3 | 2,157,221 |

All comparisons passed. Throughput is records divided by median wall time. The MVP's 16-thread versus one-thread improvement is modest; no linear-scaling claim is made. Its database must still fit in memory, so these small-database RSS measurements are not evidence of bounded memory for arbitrary database size.

## Measured writer change

The initial sorted-index MVP took 0.5602 seconds at 1M records / 16 threads. Stage timing showed output formatting was the dominant remaining cost. The retained optimization writes validated unquoted TSV coordinate/extra-column slices directly, while retaining an escaping fallback for quotes and whitespace-delimited AVinput. The final median is 0.4636 seconds, a 17.2% reduction versus that initial MVP. New tests cover empty trailing fields, quoted values and whitespace input. Initial and final source hashes and measurements distinguish these implementations.

Final stage medians at 1M records / 16 threads (these medians need not sum exactly to a separately measured total):

| Stage | Seconds |
|---|---|
| database_load_index_seconds | 0.0209 |
| parse_normalize_seconds | 0.1578 |
| index_lookup_seconds | 0.0405 |
| merge_write_seconds | 0.2408 |
| total_seconds | 0.4546 |

Parsing and representation normalization share one measured stage. Lookup includes worker scheduling; duplicate merging/escaping is included in writing. External peak RSS and CPU time are in the raw measurements. This instrumentation does not substitute for a sampling profiler or isolate each allocation cost.

## Correctness and reproduction

The workspace regression suite covers coordinate conversion, chromosome aliases, alternate-specific hits, simple indels/MNVs, randomized sorted-index versus full-scan comparisons, duplicate source-order preservation, strict malformed-input/schema errors, quoted/extra fields, atomic failure behavior, and thread/batch output invariance. A separate black-box Perl generic-filter comparison passes 10 synthetic inputs (eight matches, two misses), including repeated input and indels. It does not certify real ClinVar/dbSNP releases or conflicting duplicate database rows.

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --locked
cargo build --workspace --release --locked
python3 scripts/benchmark_filter_mvp.py --directory /tmp/NEW-benchmark-directory
python3 scripts/check_filter_mvp_perl.py --annovar /path/to/registered/annovar \
  --directory /tmp/NEW-oracle-directory
```

The scripts reject an existing run directory. `--sizes`, `--threads` and `--repeats` configure the benchmark. Large generated data and output files stay outside Git.

[Final summary](docs/benchmarks/filter-mvp-2026-09-15/final-summary.json) · [Environment](docs/benchmarks/filter-mvp-2026-09-15/final-environment.json) · [All timing records](docs/benchmarks/filter-mvp-2026-09-15/measurements.jsonl) · [Source hashes](docs/benchmarks/filter-mvp-2026-09-15/final-source.json) · [Perl comparison](docs/benchmarks/filter-mvp-2026-09-15/perl-comparison.json)

The MVP milestone ends here. Next work should evaluate integration with existing persistent indexes and real versioned filter databases before considering mmap/RADB or expanding biological consequences. No new registry release or full ANNOVAR compatibility is asserted.
