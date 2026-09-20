# Public-data annotation concordance: RustAnnovar versus ANNOVAR

**Result (19 September 2026): the two programs are not fully concordant.** On stratified samples of the public GIAB HG002 v4.2.1 small-variant benchmark, all five compared gene-annotation fields agreed for **59,369/60,389 alleles (98.311%)** with hg19/refGene and **59,683/60,687 (98.346%)** with hg38/refGeneWithVer after disregarding only comma-list order in `Gene` and `AAChange`. There were **1,020** and **1,004** substantive row-level disagreements, respectively. These are *agreement rates against a particular ANNOVAR installation*, not estimates of the biological accuracy of either annotator. GIAB establishes variant calls, not ground-truth transcript consequences.

RustAnnovar remains in **open beta**. Users should review discordant annotations and independently validate results before consequential use. This report does not certify complete ANNOVAR compatibility.

## Public input and controlled comparison

The source is the [NIST Genome in a Bottle HG002 v4.2.1 benchmark](https://www.nist.gov/programs-projects/genome-bottle), using its [GRCh37 VCF](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz) and [GRCh38 VCF](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz). Only chromosomes 1–22 in those files were evaluated. Source and sample SHA-256 values, selection rates, and counts are in the [hg19](benchmarks/2026-09-19/hg19-selection.json) and [hg38](benchmarks/2026-09-19/hg38-selection.json) manifests.

Sampling was deterministic: BLAKE2b of seed `RustAnnovar-GIAB-v1` plus `CHROM`, `POS`, `REF`, `ALT` selected 1% of simple SNVs, 3% of indels, and 10% of multiallelic/other records. This deliberately enriches difficult alleles; the overall percentage **must not be interpreted as a population-wide concordance estimate**. It selected 55,586 VCF records for GRCh37 and 55,882 for GRCh38. Neither source VCF contained a separately classified MNV under this sampler, so MNV accuracy is not assessed here.

For each assembly, the registered ANNOVAR `convert2annovar.pl --format vcf4` produced **one shared AVinput**. The Perl `table_annovar.pl` and the Rust `table` command then annotated that identical input against the same fixed `humandb` gene model and transcript FASTA. The hg19 protocol was `refGene`; hg38 used `refGeneWithVer`. Perl's default table polishing was retained. The Rust binary was version `0.2.0-beta.2`, built from source commit `5499531`; this documentation update does not change that executable. The Perl installation used Perl 5.34.0. The registered ANNOVAR package and database files are not redistributed.

For exact provenance, the Rust binary SHA-256 was `69d47f65f6c8c96c2162b3f6859b9de9e7418aa12029fbbde6e46ca2fb83c0fa`; the final comparator SHA-256 was `30ddf7327b74f8799b22447ac3375f138c2e47484f452deb8b54f31284499cc4`. The full per-run checksum lists remain on the HPC.

Rows were paired by `Chr/Start/End/Ref/Alt` plus duplicate occurrence. The comparison covered `Func`, `Gene`, `GeneDetail`, `ExonicFunc`, and `AAChange` for the requested protocol. The comparator reports literal equality and, separately, equality after sorting only comma-separated `Gene` and `AAChange` entries. It does not suppress any other difference. Both programs emitted the same number of alleles, with no missing row or duplicate-key mismatch on either assembly. [Comparator source](../scripts/giab_accuracy.py) and its [unit tests](../scripts/test_giab_accuracy.py) are public.

## Annotation results

| Measure | hg19 / refGene | hg38 / refGeneWithVer |
|---|---:|---:|
| Shared AVinput alleles; rows from each program | 60,389 | 60,687 |
| All five fields literally identical | 59,288 (98.177%) | 59,605 (98.217%) |
| All five fields equal allowing only list order | 59,369 (98.311%) | 59,683 (98.346%) |
| Substantive disagreements after list-order allowance | 1,020 | 1,004 |
| SNV rows equal after list-order allowance | 35,010 / 35,427 (98.823%) | 34,775 / 35,195 (98.807%) |
| Indel rows equal after list-order allowance | 24,306 / 24,907 (97.587%) | 24,854 / 25,437 (97.708%) |

| Field | hg19 substantive differences / all alleles | hg38 substantive differences / all alleles |
|---|---:|---:|
| `Func` | 4 / 60,389 | 1 / 60,687 |
| `Gene` | 507 / 60,389 | 528 / 60,687 |
| `GeneDetail` | 1,015 / 60,389 | 998 / 60,687 |
| `ExonicFunc` | 9 / 60,389 | 7 / 60,687 |
| `AAChange` | 12 / 60,389; another 81 list-order-only | 10 / 60,687; another 78 list-order-only |

Most `ExonicFunc` and `AAChange` cells are missing because most sampled alleles do not receive coding annotation. Conditional on **Perl producing a nonmissing value**, `ExonicFunc` agrees for **243/250 (97.2%)** on hg19 and **220/227 (96.9%)** on hg38. `AAChange`, allowing list order, agrees for **238/250 (95.2%)** and **217/227 (95.6%)**, respectively. These conditional values are more informative for coding consequences than an all-allele percentage inflated by shared missing cells. `GeneDetail` agrees for 32,667/33,673 and 32,794/33,787 nonmissing Perl values.

Examples of unresolved differences include an hg19 insertion at `20:31108088` classified `intronic` by Rust and `ncRNA_exonic` by Perl; an insertion at `6:30558477` classified `frameshift insertion` versus `nonframeshift insertion`; and an hg38 deletion at `chr1:248626808-248626809` classified `frameshift deletion` versus `stopgain`. These are **observed disagreements**, not adjudicated errors. Complete field counts and up to ten concrete examples per field are in [hg19 results](benchmarks/2026-09-19/hg19-annotation-comparison.json) and [hg38 results](benchmarks/2026-09-19/hg38-annotation-comparison.json).

## Separate VCF conversion check

The shared-AVinput test isolates the annotation engine. A second test converted each selected VCF independently with Rust and Perl, then compared AVinput `Chr/Start/End/Ref/Alt` multisets, retaining duplicate multiplicity. Its [hg19](benchmarks/2026-09-19/hg19-conversion-comparison.json) and [hg38](benchmarks/2026-09-19/hg38-conversion-comparison.json) results are **not equivalent**:

| Measure | hg19 | hg38 |
|---|---:|---:|
| Rust / Perl output rows | 60,390 / 60,389 | 60,687 / 60,687 |
| Unmatched Rust / Perl rows, exact five columns | 2,966 / 2,965 | 60,687 / 60,687 |
| Unmatched Rust / Perl rows, after stripping only `chr` prefix | 2,966 / 2,965 | 2,981 / 2,981 |
| Of normalized Rust unmatched rows: indel / symbolic or ambiguous | 2,911 / 55 | 2,926 / 55 |

The GRCh38 exact-key result is dominated by a naming difference: Perl retains `chr1`, whereas Rust emits `1`. After that prefix alone is harmonized, **all simple SNV keys match** (35,427 hg19; 35,195 hg38), but thousands of indel positions/representations and symbolic alleles still differ. For example, one hg19 deletion is emitted at `1:1701409-1701414` by Rust and `1:1701415-1701420` by Perl. The comparison did not use a reference genome to decide whether such shifted indels are biologically equivalent. One extra hg19 Rust AVinput row also needs investigation. Thus direct VCF-to-output equivalence is **not established**, even where annotation of a shared AVinput agrees.

## Reproduction and provenance

All large files and registered ANNOVAR assets remain under `/DATABANK/users/hflt/RustAnnovar/` on the remote HPC. The fixed VCFs are in `datasets/giab/`; selected inputs and manifests are in `datasets/giab/accuracy-sample-20260919/` and `datasets/giab/accuracy-sample-hg38-20260919/`. Full logs, time/memory measurements, SHA-256 lists, converted AVinput, and complete outputs remain in `reports/giab-accuracy-refgene-20260919/` and `reports/giab-accuracy-refgenewithver-hg38-20260919/`. Both annotation runners completed successfully; their `comparison-status` is `DIFFERENCES_RECORDED`, **not a compatibility pass**. The public JSONs are compact copies of final reports. Their SHA-256 values and the [runner](../scripts/hpc_giab_accuracy.sh) make the evidence auditable without publishing licensed files.

From the repository checkout on that HPC, with the same local ANNOVAR package, databases, and fixed Rust binary:

```bash
ROOT=/DATABANK/users/hflt/RustAnnovar
mkdir -p "$ROOT/datasets/giab/accuracy-sample-new"
python3 scripts/giab_accuracy.py sample \
  "$ROOT/datasets/giab/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz" \
  "$ROOT/datasets/giab/accuracy-sample-new/HG002_GRCh37_accuracy.vcf" \
  --report "$ROOT/datasets/giab/accuracy-sample-new/selection.json"
bash scripts/hpc_giab_accuracy.sh "$ROOT" \
  "$ROOT/datasets/giab/accuracy-sample-new" \
  "$ROOT/reports/giab-accuracy-new" hg19 refGene
```

Substitute `GRCh38`, `hg38`, and `refGeneWithVer` for the other assembly. The runner refuses to overwrite an existing run directory and produces both annotation and conversion comparisons. A comparator's nonzero exit code intentionally signals a measured difference; the runner records that status without treating the measurement itself as a failed job. The checked-in manifests record the exact source and selected-VCF checksums required to verify a rerun.

**Scope:** this is a stratified, sampled **gene-annotation concordance** study on hg19/hg38, not full WGS annotation, a performance benchmark, or a clinical validation. It does not assess ClinVar, gnomAD, cytoBand, other annotation protocols, structural variants, MNV consequences, or whether a reported difference reflects a defect in ANNOVAR. Fixes should be followed by rerunning both the shared-AVinput and independent conversion comparisons.
