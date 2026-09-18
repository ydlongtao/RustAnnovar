use rust_annovar::{
    cadd,
    disk_index::DiskFilter,
    model::Variant,
    pipeline::{AnnotationEngine, Operation, Protocol},
};
use std::fs;

/// Run against independently downloaded official records, without distributing data.
#[test]
#[ignore = "requires CADD_TEST_TSV pointing to an official score-only TSV"]
fn official_records_are_all_scored_by_native_lookup() {
    use noodles_core::Position;
    use noodles_csi::binning_index::index::{header::Builder, reference_sequence::bin::Chunk};
    use std::io::{BufRead, Write};
    let input = std::env::var("CADD_TEST_TSV").unwrap();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("official.tsv.gz");
    let mut writer = noodles_bgzf::Writer::new(fs::File::create(&path).unwrap());
    let mut indexer = noodles_tabix::index::Indexer::default();
    indexer.set_header(
        Builder::gff()
            .set_start_position_index(1)
            .set_end_position_index(None)
            .build(),
    );
    let mut variants = Vec::new();
    let mut expected = Vec::new();
    for line in std::io::BufReader::new(fs::File::open(input).unwrap()).lines() {
        let line = line.unwrap();
        let begin = writer.virtual_position();
        writeln!(writer, "{line}").unwrap();
        if line.starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let pos: u64 = fields[1].parse().unwrap();
        let position = Position::try_from(pos as usize).unwrap();
        indexer
            .add_record(
                fields[0],
                position,
                position,
                Chunk::new(begin, writer.virtual_position()),
            )
            .unwrap();
        variants.push(Variant::new(fields[0], pos - 1, pos, fields[2], fields[3]).unwrap());
        expected.push(vec![
            fields[4].to_string(),
            fields[5].to_string(),
            "scored".to_string(),
        ]);
    }
    writer.finish().unwrap();
    noodles_tabix::io::Writer::new(fs::File::create(format!("{}.tbi", path.display())).unwrap())
        .write_index(&indexer.build())
        .unwrap();
    let db = rust_annovar::cadd::Database::open(&path).unwrap();
    assert!(!variants.is_empty());
    for (input, expected) in variants.chunks(10000).zip(expected.chunks(10000)) {
        assert_eq!(db.batch(input, ".").unwrap(), expected);
    }
    println!("verified {} official CADD alleles", variants.len());
}

#[test]
fn native_tabix_scores_and_reports_all_missing_categories() {
    use noodles_core::Position;
    use noodles_csi::binning_index::index::{header::Builder, reference_sequence::bin::Chunk};
    use std::io::Write;
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("scores.tsv.gz");
    let mut writer = noodles_bgzf::Writer::new(fs::File::create(&path).unwrap());
    writeln!(writer, "##CADD GRCh37-v1.7").unwrap();
    writeln!(writer, "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED").unwrap();
    let mut indexer = noodles_tabix::index::Indexer::default();
    indexer.set_header(
        Builder::gff()
            .set_start_position_index(1)
            .set_end_position_index(None)
            .build(),
    );
    let begin = writer.virtual_position();
    writeln!(writer, "11\t118373578\tA\tG\t3.477832\t22.3").unwrap();
    let end = writer.virtual_position();
    let pos = Position::try_from(118373578).unwrap();
    indexer
        .add_record("11", pos, pos, Chunk::new(begin, end))
        .unwrap();
    writer.finish().unwrap();
    noodles_tabix::io::Writer::new(fs::File::create(format!("{}.tbi", path.display())).unwrap())
        .write_index(&indexer.build())
        .unwrap();
    let db = rust_annovar::cadd::Database::open(&path).unwrap();
    db.check_build("hg19").unwrap();
    db.check_version("1.7").unwrap();
    assert!(db.check_version("1.6").is_err());
    assert!(db.check_build("hg38").is_err());
    let variants = [
        Variant::new("chr11", 118373577, 118373578, "A", "G").unwrap(),
        Variant::new("11", 118373577, 118373578, "C", "G").unwrap(),
        Variant::new("11", 118373578, 118373579, "A", "C").unwrap(),
        Variant::new("MT", 9, 10, "A", "C").unwrap(),
        Variant::new("11", 9, 11, "AA", "CC").unwrap(),
        Variant::new("11", 9, 10, "N", "C").unwrap(),
    ];
    let rows = db.batch(&variants, ".").unwrap();
    assert_eq!(rows[0], ["3.477832", "22.3", "scored"]);
    let vcf = dir.path().join("sample.vcf");
    fs::write(&vcf, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n11\t118373578\t.\tA\tG,C\t.\tPASS\t.\tGT\t1/2\n").unwrap();
    let output = dir.path().join("result.tsv");
    let report = dir.path().join("report.json");
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_rust-annovar"))
        .args([
            "annotate",
            vcf.to_str().unwrap(),
            path.to_str().unwrap(),
            "--operation",
            "cadd",
            "--protocol",
            "cadd17",
            "--vcf-input",
            "--output",
            output.to_str().unwrap(),
            "--report-json",
            report.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let text = fs::read_to_string(&output).unwrap();
    assert!(text.contains("3.477832\t22.3\tscored"));
    let counts: serde_json::Value =
        serde_json::from_reader(fs::File::open(report).unwrap()).unwrap();
    assert_eq!(counts["cadd_counts"]["cadd17"]["scored"], 1);
    assert_eq!(counts["cadd_counts"]["cadd17"]["score_not_found"], 1);
    let table_db = dir.path().join("hg19_cadd17.txt.gz");
    fs::copy(&path, &table_db).unwrap();
    fs::copy(
        format!("{}.tbi", path.display()),
        format!("{}.tbi", table_db.display()),
    )
    .unwrap();
    let vcf_output = dir.path().join("annotated.vcf");
    let status = std::process::Command::new(env!("CARGO_BIN_EXE_rust-annovar"))
        .args([
            "table",
            vcf.to_str().unwrap(),
            dir.path().to_str().unwrap(),
            "--build",
            "hg19",
            "--protocol",
            "cadd17",
            "--operation",
            "cadd",
            "--vcf-input",
            "--output",
            output.to_str().unwrap(),
            "--vcf-output",
            vcf_output.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let annotated = fs::read_to_string(vcf_output).unwrap();
    assert!(annotated.contains("FA_cadd17=3.477832|22.3|scored,.|.|score_not_found"));
    assert!(annotated.contains("\tGT\t1/2"));
    fs::write(&output, "existing destination").unwrap();
    let strict = std::process::Command::new(env!("CARGO_BIN_EXE_rust-annovar"))
        .args([
            "annotate",
            vcf.to_str().unwrap(),
            path.to_str().unwrap(),
            "--operation",
            "cadd",
            "--vcf-input",
            "--output",
            output.to_str().unwrap(),
            "--cadd-require-all",
            "--batch-size",
            "1",
        ])
        .status()
        .unwrap();
    assert!(!strict.success());
    assert_eq!(fs::read_to_string(output).unwrap(), "existing destination");
    assert_eq!(
        rows.iter().map(|r| r[2].as_str()).collect::<Vec<_>>(),
        [
            "scored",
            "reference_mismatch",
            "score_not_found",
            "contig_not_covered",
            "not_snv",
            "unsupported_snv"
        ]
    );
}

#[test]
fn cadd_known_score_survives_index_and_exact_allele_matching() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("official.tsv");
    let output = dir.path().join("hg19_cadd17.txt");
    // Public CADD GRCh37-v1.7 lookup: chr11:118373578 A>G.
    fs::write(&input, "## CADD v1.7\n#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n11\t118373578\tA\tG\t3.477832\t22.3\n").unwrap();
    assert_eq!(cadd::import(&input, &output, "hg19", "1.7").unwrap(), 1);
    let protocols = [Protocol {
        name: "cadd17".into(),
        operation: Operation::Filter,
        database: output.clone(),
        fasta: None,
    }];
    let variants = [
        Variant::new("chr11", 118373577, 118373578, "A", "G").unwrap(),
        Variant::new("11", 118373577, 118373578, "C", "G").unwrap(),
    ];
    let scan = AnnotationEngine::new(&protocols)
        .unwrap()
        .annotate(&variants, ".")
        .unwrap();
    DiskFilter::build_normalized(&output, None, None).unwrap();
    let indexed = AnnotationEngine::new(&protocols)
        .unwrap()
        .annotate(&variants, ".")
        .unwrap();
    assert_eq!(scan.rows, indexed.rows);
    assert_eq!(&indexed.rows[0][5..7], ["3.477832", "22.3"]);
    assert_eq!(&indexed.rows[1][5..7], [".", "."]);
}

#[test]
fn cadd_bad_scores_fail_atomically() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("bad.tsv");
    let output = dir.path().join("out.txt");
    fs::write(
        &input,
        "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n1\t10\tA\tG\tNaN\t22\n",
    )
    .unwrap();
    assert!(cadd::import(&input, &output, "hg38", "1.7").is_err());
    assert!(!output.exists());
}
