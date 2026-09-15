# WGS engine acceptance: interim report

**Status: development candidate, not an accepted WGS replacement.** The 0.2.0-beta.1 source version does not imply a public release. Synthetic memory scaling passes in the initial candidate; the requested speed target fails. GIAB, real ClinVar/gnomAD, and hg38 gene-model acceptance remain pending. Existing 0.1 benchmark claims must not be reused for this engine.

## Fixed environment and method

Host: `gse234129-hpc` (hostname `dell`), Intel Xeon Platinum 8368Q, 76 physical cores / 152 logical CPUs, approximately 1 TiB RAM. Workspace: `/DATABANK/users/hflt/RustAnnovar/`, ext4 shared filesystem. Jobs are serialized with an advisory lock. Thread tests bind distinct physical cores on one NUMA node. The requested application budget is 16 GiB; this option currently caps batch input bytes, not all application allocations.

Rust baseline: `deb7a25f1ae178ad8c6060629b1f75f19bec834d`. Baseline and candidate binaries are built with the same GitHub Actions Linux/toolchain setup. Remote run directories retain executable hashes, compiler/commit artifact metadata, topology, commands, exit codes, wall time, maximum RSS, load averages and swap counters. Warm results below are medians of five runs, preceded by one separately recorded first run. Global caches are never cleared; first runs are not cold-cache measurements.

Synthetic input repeats 100,000 unique A>C sites. The filter database contains those sites with one text value per site. All outputs are written to the same storage. The candidate uses its disk index, whose build is timed separately; the baseline uses its original database path. Table equality is checked after excluding the new status column. This intentionally simple workload measures streaming behavior and must not be labeled a WGS biological benchmark.

## Initial candidate results

Candidate `ef68553`, run `runs/20260915T233408-scale`, 16 threads:

| Records | Implementation | Median seconds (range) | Maximum RSS |
|---|---|---|---|
| 1,000,000 | Rust baseline | 4.94 (4.91–5.00) | 1,522.9 MiB |
| 1,000,000 | Candidate | 7.64 (7.53–7.66) | 39.2 MiB |
| 10,000,000 | Rust baseline | 19.34 (19.10–19.67) | 14,843.7 MiB |
| 10,000,000 | Candidate | 75.56 (74.83–76.55) | 40.9 MiB |

Candidate RSS growth is 4.3%, below the 20% target, and peak RSS is below 20 GiB for this workload. Candidate runtime is 3.91 times baseline at 10M records, failing the target of at most 0.5 times baseline. Output equality passed. See [raw summary](benchmarks/2026-09-15/ef68553-scale.json).

Candidate `5a7ea58` caches decoded filter blocks and reduces small-batch task overhead. Its 100,000-record smoke test at 16 threads measures 0.88 seconds (0.88–0.89), versus baseline 0.66 seconds (0.64–0.67), also failing the speed target. This short workload has substantial fixed costs. Results for 1/2/4/8/16/32 threads are in the [raw summary](benchmarks/2026-09-15/5a7ea58-smoke.json). Its 1M/10M run is recorded separately when complete. These timings do not include later gene-classification fixes.

The initial scale campaign recorded one-minute system load between 2.52 and 8.70 and zero swap-in/out page increments during measured commands. The revised smoke campaign recorded load 2.67–3.03 and zero swap increments. These system-wide observations do not prove absence of per-core or storage contention. The [environment audit](benchmarks/2026-09-15/environment.json) retains CPU affinity and executable hashes.

Perl speed on matching real GIAB inputs has not yet been measured in this campaign. No Perl acceleration factor is asserted.

## Correctness evidence and remaining gates

Local validation passes 22 regular tests and two explicitly enabled Perl oracle tests. Coverage includes gzip external sorting, duplicate hits, corrupt index blocks, stale sources, nested intervals, mixed supported/symbolic ALT, preserved samples, batch/thread invariance, atomic failure handling, left alignment policy agreement, and positive/negative-strand small coding changes. Clippy with warnings denied passes.

A separate 50,000-site chr1 positional check found and fixed exclusion of noncoding isoforms for genes with coding isoforms, ncRNA category precedence, combined function separators and flank boundaries. After correction, Func.refGene and Gene.refGene match Perl at every tested site. These arbitrary A>G sites are not reference-consistent and therefore provide **no evidence of protein-consequence agreement**. The input recipe and checksum are recorded in [the positional summary](benchmarks/2026-09-15/positional-check.json). Gene model and Perl checksums are in [the baseline manifest](WGS_BASELINE.json).

| Gate | State |
|---|---|
| Streaming, compact batch representation, interval/disk indexes | Implemented; compatibility adapters materialize variants per batch |
| Fixed example / synthetic regression | Passed within the stated test scope |
| 1M→10M RSS growth and 20 GiB peak | Passed for initial synthetic filter workload only |
| Runtime ≤50% of previous Rust baseline | Failed on measured synthetic workloads |
| GIAB HG002 v4.2.1 GRCh37/GRCh38 | Downloading with resume; not validated |
| Matching hg38 model/FASTA and versioned ClinVar/gnomAD | Not available in the fixed local oracle; not validated |
| Complete field equivalence including MNV/Indel text | Pending; preliminary protein notation is not fully compatible |
| Application-wide memory accounting and bounded result expansion | Pending; batch budget is not an RSS guarantee |
| Public GitHub/RubyGems/GitHub Packages release | Deferred until acceptance evidence supports the release |

## Reproduction

Build and invocation examples are in [WGS_ENGINE.md](WGS_ENGINE.md). Use `scripts/hpc_benchmark.py --candidate-dir <immutable-build-directory> --mode smoke|scale|giab --wait`. Keep `scripts/compare_annotations.py` beside the harness. GIAB runs write per-build field differences and mark missing inputs explicitly. A completed timing run is not a compatibility pass. Inspect load/swap records before accepting measurements affected by shared-host contention.

Large inputs, licensed ANNOVAR scripts and databases remain in the dedicated workspace, outside Git. Run logs and resumable download state stay on the server; only source, data manifests and compact reports are published.
