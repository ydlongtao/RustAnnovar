# Changelog

## Unreleased — CADD SNV integration

- Added reproducible public GIAB HG002 v4.2.1 hg19/hg38 ANNOVAR gene-annotation concordance analysis, including exact and list-order-aware field comparisons, independent VCF conversion checks, compact evidence JSON, and an English report. The results expose unresolved differences; they do not establish full compatibility.
- Added native local Tabix and explicit HTTPS byte-range queries of official score-only CADD tables.
- Added RawScore, PHRED, per-allele status, strict SNV coverage mode, assembly/version checks and JSON provenance.
- Added streaming import, checksum-verified resumable HPC downloads and independent official API probe verification.
- Real-data fragment and cross-chromosome checks pass. Complete official GRCh38 v1.7 score-only data passed strict autosomal GIAB HG002 validation; GRCh37 full-file acceptance remains pending. See `docs/CADD-GIAB-validation.md` for scope and measured runtime.

## 0.2.0-beta.2 — main source update

- Merged streaming annotation and the isolated sorted-filter MVP into `main`.
- Unified all six Cargo package versions at `0.2.0-beta.2`; source-gem metadata uses `0.2.0.beta.2`.
- Retained open-beta labeling, measured performance scope and unresolved WGS/Indel compatibility gates. This source update does not publish registry packages.

## Development history — isolated filter MVP

- Added a five-package Rust workspace and experimental `rustannovar annotate` AVinput frontend.
- Added compact chromosome/allele keys, chromosome-sorted binary-search indexes and borrowed annotation-pool matches.
- Added bounded-record streaming batches, ordered parallel chunks, strict schema/coordinate validation and buffered atomic TSV output.
- Added stage timing, portable 10K/100K/1M benchmark tooling and a synthetic black-box Perl filter comparison.
- Kept the existing `rust-annovar` commands and package version intact. Workspace MVP package version 0.1.0 is an internal development version, not a public release or a rollback of the root 0.2.0-beta.1 candidate.
- Updated workspace CI and source-gem file inclusion. RubyGems continues to install the existing `rust-annovar` binary.
