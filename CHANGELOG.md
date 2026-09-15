# Changelog

## Unreleased — isolated filter MVP

- Added a five-package Rust workspace and experimental `rustannovar annotate` AVinput frontend.
- Added compact chromosome/allele keys, chromosome-sorted binary-search indexes and borrowed annotation-pool matches.
- Added bounded-record streaming batches, ordered parallel chunks, strict schema/coordinate validation and buffered atomic TSV output.
- Added stage timing, portable 10K/100K/1M benchmark tooling and a synthetic black-box Perl filter comparison.
- Kept the existing `rust-annovar` commands and package version intact. Workspace MVP package version 0.1.0 is an internal development version, not a public release or a rollback of the root 0.2.0-beta.1 candidate.
- Updated workspace CI and source-gem file inclusion. RubyGems continues to install the existing `rust-annovar` binary.
