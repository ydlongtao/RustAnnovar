use flate2::Compression;
use flate2::write::GzEncoder;
use rust_annovar::database::{FilterDatabase, IndexManifest, RegionDatabase};
use rust_annovar::gene::GeneDatabase;
use rust_annovar::io::{normalize_vcf_alleles, read_avinput, read_vcf, write_avinput};
use rust_annovar::model::Variant;
use rust_annovar::pipeline::{Operation, Protocol, annotate_table};
use std::fs;
use std::io::Write;
use tempfile::tempdir;

#[test]
fn converts_multiallelic_vcf_and_normalizes_indels() {
    let dir = tempdir().unwrap();
    let vcf = dir.path().join("input.vcf");
    fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC,ATG\t.\tPASS\t.\n",
    ).unwrap();
    let document = read_vcf(&vcf).unwrap();
    assert_eq!(document.variants.len(), 2);
    assert_eq!(
        document.variants[0].avinput_fields(),
        ["1", "10", "10", "A", "C"]
    );
    assert_eq!(
        document.variants[1].avinput_fields(),
        ["1", "10", "10", "-", "TG"]
    );
}

#[test]
fn avinput_round_trip_preserves_event_coordinates() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("input.avinput");
    fs::write(
        &input,
        "1\t10\t10\tA\tC\tnote\n1\t20\t20\t-\tTG\n1\t30\t31\tAC\t-\n",
    )
    .unwrap();
    let variants = read_avinput(&input).unwrap();
    let output = dir.path().join("output.avinput");
    write_avinput(&variants, &output, true).unwrap();
    assert_eq!(
        fs::read_to_string(output).unwrap(),
        "1\t10\t10\tA\tC\tnote\n1\t20\t20\t-\tTG\n1\t30\t31\tAC\t-\n"
    );
}

#[test]
fn filter_requires_exact_allele_match_and_joins_duplicates() {
    let dir = tempdir().unwrap();
    let db = dir.path().join("filter.txt");
    fs::write(
        &db,
        "#Chr\tStart\tEnd\tRef\tAlt\tAF\n1\t10\t10\tA\tC\t0.1\n1\t10\t10\tA\tC\t0.2\n",
    )
    .unwrap();
    let database = FilterDatabase::load(&db).unwrap();
    let hit = database
        .annotate(&Variant::new("chr1", 9, 10, "A", "C").unwrap(), "freq")
        .unwrap();
    assert_eq!(hit.values, ["0.1;0.2"]);
    assert!(
        database
            .annotate(&Variant::new("1", 9, 10, "A", "G").unwrap(), "freq")
            .is_none()
    );
}

#[test]
fn region_annotation_uses_overlap_not_alleles() {
    let dir = tempdir().unwrap();
    let db = dir.path().join("region.txt");
    fs::write(&db, "1\t8\t12\tband1\n1\t20\t30\tband2\n").unwrap();
    let database = RegionDatabase::load(&db).unwrap();
    let hit = database
        .annotate(&Variant::new("1", 9, 10, "A", "C").unwrap(), "band", 0.0)
        .unwrap();
    assert_eq!(hit.values, ["band1"]);
}

#[test]
fn region_annotation_reads_gff3_coordinates() {
    let dir = tempdir().unwrap();
    let db = dir.path().join("region.gff3");
    fs::write(
        &db,
        "##gff-version 3\nchr1\tsource\tTF_binding_site\t10\t20\t849\t-\t.\tID=x;Name=TFBS\n",
    )
    .unwrap();
    let database = RegionDatabase::load(&db).unwrap();
    let hit = database
        .annotate(&Variant::new("1", 9, 10, "A", "C").unwrap(), "gff3", 0.0)
        .unwrap();
    assert_eq!(hit.values, ["Score=849;Name=x"]);
}

#[test]
fn gene_annotation_calculates_snv_coding_change() {
    let dir = tempdir().unwrap();
    let model = dir.path().join("refGene.txt");
    let fasta = dir.path().join("refGeneMrna.fa");
    fs::write(
        &model,
        "NM_1\tchr1\t+\t0\t12\t0\t12\t1\t0,\t12,\t0\tGENE1\n",
    )
    .unwrap();
    fs::write(&fasta, ">NM_1\nATGGAATTTTAA\n").unwrap();
    let database = GeneDatabase::load(&model, Some(&fasta)).unwrap();
    let annotation = database.annotate(
        &Variant::new("1", 4, 5, "A", "G").unwrap(),
        "refGene",
        2,
        1000,
    );
    assert_eq!(annotation.values[0], "exonic");
    assert_eq!(annotation.values[1], "GENE1");
    assert_eq!(annotation.values[3], "nonsynonymous SNV");
    assert_eq!(annotation.values[4], "GENE1:NM_1:exon1:c.A5G:p.E2G");
}

#[test]
fn table_keeps_protocol_order_and_missing_values() {
    let dir = tempdir().unwrap();
    let db = dir.path().join("db.txt");
    fs::write(
        &db,
        "#Chr\tStart\tEnd\tRef\tAlt\tAF\n1\t10\t10\tA\tC\t0.1\n",
    )
    .unwrap();
    let protocols = vec![Protocol {
        name: "freq".into(),
        operation: Operation::Filter,
        database: db,
        fasta: None,
    }];
    let variants = vec![
        Variant::new("1", 9, 10, "A", "C").unwrap(),
        Variant::new("1", 10, 11, "A", "G").unwrap(),
    ];
    let result = annotate_table(&variants, &protocols, ".").unwrap();
    assert_eq!(
        result.headers,
        ["Chr", "Start", "End", "Ref", "Alt", "AF.freq"]
    );
    assert_eq!(result.rows[0][5], "0.1");
    assert_eq!(result.rows[1][5], ".");
}

#[test]
fn manifest_detects_source_changes() {
    let dir = tempdir().unwrap();
    let source = dir.path().join("db.txt");
    fs::write(&source, "one").unwrap();
    let manifest = IndexManifest::build(&source, "filter").unwrap();
    assert!(manifest.is_current().unwrap());
    fs::write(&source, "changed content").unwrap();
    assert!(!manifest.is_current().unwrap());
}

#[test]
fn filter_sidecar_loads_only_requested_bins() {
    let dir = tempdir().unwrap();
    let source = dir.path().join("filter.txt");
    fs::write(
        &source,
        "#Chr\tStart\tEnd\tRef\tAlt\tAF\n1\t10\t10\tA\tC\t0.1\n1\t2000010\t2000010\tG\tT\t0.2\n",
    )
    .unwrap();
    let index = IndexManifest::build(&source, "filter").unwrap();
    index.write(&source.with_extension("fai.json")).unwrap();
    let query = Variant::new("1", 9, 10, "A", "C").unwrap();
    let database =
        FilterDatabase::load_for_variants(&source, std::slice::from_ref(&query)).unwrap();
    assert_eq!(database.annotate(&query, "freq").unwrap().values, ["0.1"]);
    assert!(
        database
            .annotate(
                &Variant::new("1", 2_000_009, 2_000_010, "G", "T").unwrap(),
                "freq"
            )
            .is_none()
    );
}

#[test]
fn normalization_handles_deletion_anchor() {
    assert_eq!(
        normalize_vcf_alleles(9, "ATC", "A"),
        (10, 12, "TC".into(), "".into())
    );
}

#[test]
fn reads_gzip_vcf() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("input.vcf.gz");
    let file = fs::File::create(&path).unwrap();
    let mut gzip = GzEncoder::new(file, Compression::fast());
    gzip
        .write_all(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t7\t.\tG\tA\t.\t.\t.\n")
        .unwrap();
    gzip.finish().unwrap();
    let document = read_vcf(&path).unwrap();
    assert_eq!(
        document.variants[0].avinput_fields(),
        ["1", "7", "7", "G", "A"]
    );
}

#[test]
fn avinput_zero_reference_retains_large_deletion_span() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("large.avinput");
    fs::write(&input, "13\t20797176\t21105944\t0\t-\tlarge deletion\n").unwrap();
    let variants = read_avinput(&input).unwrap();
    assert_eq!(variants[0].start, 20_797_175);
    assert_eq!(variants[0].end, 21_105_944);
    assert_eq!(variants[0].reference, "0");
    assert_eq!(variants[0].extra, ["large deletion"]);
    assert_eq!(
        variants[0].avinput_fields(),
        ["13", "20797176", "21105944", "0", "-"]
    );
}
