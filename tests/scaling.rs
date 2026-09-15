use rust_annovar::{
    Variant, database::FilterDatabase, disk_index::DiskFilter, interval::IntervalIndex,
};
use std::{fs, process::Command};
use tempfile::tempdir;
#[test]
fn interval_index_matches_full_scan_with_nested_and_touching_intervals() {
    let spans: Vec<_> = (0..200).map(|i| (i * 7, i * 7 + (i % 13) * 19)).collect();
    let index = IntervalIndex::new(spans.iter().copied());
    for pos in (0..1600).step_by(3) {
        let expected: Vec<_> = spans
            .iter()
            .enumerate()
            .filter(|(_, r)| r.0 <= pos + 2 && r.1 >= pos)
            .map(|(i, _)| i)
            .collect();
        assert_eq!(index.query(pos, pos + 2), expected);
    }
}
#[test]
fn disk_index_matches_unsorted_duplicate_records_and_rejects_mutation() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("db.txt");
    let mut text = "#Chr\tStart\tEnd\tRef\tAlt\tAF\n".to_string();
    for i in (1..6000).rev() {
        text.push_str(&format!("1\t{i}\t{i}\tA\tC\t{i}\n"));
    }
    text.push_str("1\t10\t10\tA\tC\tduplicate\n");
    fs::write(&path, &text).unwrap();
    let scan = FilterDatabase::load(&path).unwrap();
    DiskFilter::build(&path, Some(dir.path())).unwrap();
    let disk = DiskFilter::open(&path, true).unwrap();
    let vs = vec![
        Variant::new("1", 9, 10, "A", "C").unwrap(),
        Variant::new("1", 3000, 3001, "A", "C").unwrap(),
        Variant::new("2", 9, 10, "A", "C").unwrap(),
    ];
    for _ in 0..2 {
        let batch = disk.load_batch(&vs).unwrap();
        for v in &vs {
            assert_eq!(
                scan.annotate(v, "x").map(|a| a.values),
                batch.annotate(v, "x").map(|a| a.values)
            );
        }
    }
    fs::write(&path, text + "1\t1\t1\tA\tG\tother\n").unwrap();
    assert!(DiskFilter::open(&path, false).is_err());
}
#[test]
fn streaming_preserves_unsupported_alt_slots_and_is_batch_thread_invariant() {
    let d = tempdir().unwrap();
    let db = d.path().join("db.txt");
    fs::write(
        &db,
        "#Chr\tStart\tEnd\tRef\tAlt\tvalue\n1\t10\t10\tA\tC\tfirst\n1\t10\t10\tA\tG\tsecond\n",
    )
    .unwrap();
    fs::create_dir(d.path().join("humandb")).unwrap();
    fs::copy(&db, d.path().join("humandb/hg38_test.txt")).unwrap();
    let input = d.path().join("input.vcf");
    fs::write(&input,"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n1\t10\t.\tA\tC,<DEL>,G\t.\tPASS\t.\tGT:AD\t1/3:1,2,0,3\n1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\n").unwrap();
    let mut outputs = Vec::new();
    for (batch, threads) in [(1, 1), (3, 4)] {
        let out = d.path().join(format!("{batch}.tsv"));
        let vout = d.path().join(format!("{batch}.vcf"));
        let status = Command::new(env!("CARGO_BIN_EXE_rust-annovar"))
            .arg("table")
            .arg(&input)
            .arg(d.path().join("humandb"))
            .args([
                "--build",
                "hg38",
                "--protocol",
                "test",
                "--operation",
                "f",
                "--vcf-input",
                "--output",
            ])
            .arg(&out)
            .arg("--vcf-output")
            .arg(&vout)
            .args([
                "--batch-size",
                &batch.to_string(),
                "--threads",
                &threads.to_string(),
            ])
            .status()
            .unwrap();
        assert!(status.success());
        outputs.push((
            fs::read_to_string(out).unwrap(),
            fs::read_to_string(vout).unwrap(),
        ));
    }
    assert_eq!(outputs[0], outputs[1]);
    assert!(outputs[0].0.contains("unsupported_alt"));
    assert!(
        outputs[0]
            .1
            .contains("FA_test=first,.,second;FA_STATUS=ok,unsupported_alt,ok\tGT:AD\t1/3:1,2,0,3")
    );
}
#[test]
fn malformed_late_input_does_not_replace_existing_output() {
    let d = tempdir().unwrap();
    let input = d.path().join("input.avinput");
    let db = d.path().join("db.txt");
    let out = d.path().join("out.tsv");
    fs::write(&input, "1\t1\t1\tA\tC\nbad\n").unwrap();
    fs::write(&db, "1\t1\t1\tA\tC\tvalue\n").unwrap();
    fs::write(&out, "existing").unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_rust-annovar"))
        .arg("annotate")
        .arg(input)
        .arg(db)
        .arg("--output")
        .arg(&out)
        .args(["--batch-size", "1"])
        .status()
        .unwrap();
    assert!(!status.success());
    assert_eq!(fs::read_to_string(out).unwrap(), "existing");
}

#[test]
fn reference_normalization_and_index_policy_agree() {
    use rust_annovar::reference::ReferenceGenome;
    let d = tempdir().unwrap();
    let fasta = d.path().join("ref.fa");
    fs::write(&fasta, ">1\nCAAAAG\n").unwrap();
    fs::write(d.path().join("ref.fa.fai"), "1\t6\t3\t6\t7\n").unwrap();
    let reference = ReferenceGenome::open(&fasta).unwrap();
    let mut a = Variant::new("1", 4, 4, "", "A").unwrap();
    reference.normalize(&mut a).unwrap();
    assert_eq!((a.start, a.end), (1, 1));
    let source = d.path().join("db.txt");
    fs::write(&source, "1\t4\t4\t-\tA\tmatch\n").unwrap();
    DiskFilter::build_normalized(&source, Some(d.path()), Some(&reference)).unwrap();
    let index = DiskFilter::open(&source, true).unwrap();
    assert!(index.check_normalization(None).is_err());
    index.check_normalization(Some(&reference)).unwrap();
    assert_eq!(
        index
            .load_batch(&[a.clone()])
            .unwrap()
            .annotate(&a, "test")
            .unwrap()
            .values,
        ["match"]
    );
}

#[test]
fn coding_small_changes_validate_sequence_and_classify_both_strands() {
    use rust_annovar::gene::GeneDatabase;
    let d = tempdir().unwrap();
    let model = d.path().join("gene.txt");
    let fa = d.path().join("tx.fa");
    fs::write(&model,"NM_P\t1\t+\t0\t12\t0\t12\t1\t0,\t12,\t0\tPLUS\nNM_M\t2\t-\t0\t12\t0\t12\t1\t0,\t12,\t0\tMINUS\n").unwrap();
    fs::write(&fa, ">NM_P\nATGGAATTTTAA\n>NM_M\nATGGAATTTTAA\n").unwrap();
    let db = GeneDatabase::load(&model, Some(&fa)).unwrap();
    for (v, expected) in [
        (
            Variant::new("1", 3, 6, "GAA", "GGC").unwrap(),
            "nonsynonymous block substitution",
        ),
        (
            Variant::new("1", 3, 6, "GAA", "").unwrap(),
            "nonframeshift deletion",
        ),
        (
            Variant::new("1", 6, 6, "", "A").unwrap(),
            "frameshift insertion",
        ),
        (
            Variant::new("2", 6, 9, "TTC", "").unwrap(),
            "nonframeshift deletion",
        ),
    ] {
        let hit = db.annotate(&v, "refGene", 2, 1000);
        assert_eq!(hit.values[3], expected);
        assert!(hit.values[4].contains(":c."));
        assert!(hit.values[4].contains(":p."));
        assert!(!hit.values[4].contains('?'));
    }
    let hit = db.annotate(
        &Variant::new("1", 3, 6, "CCC", "").unwrap(),
        "refGene",
        2,
        1000,
    );
    assert_eq!(hit.values[3], "unknown");
    assert_eq!(hit.values[2], "reference_mismatch");
}

#[test]
fn disk_index_external_merge_handles_gzip_and_corrupt_data() {
    use flate2::{Compression, write::GzEncoder};
    use std::io::Write;
    let d = tempdir().unwrap();
    let source = d.path().join("db.txt.gz");
    let mut gzip = GzEncoder::new(fs::File::create(&source).unwrap(), Compression::fast());
    writeln!(gzip, "#Chr\tStart\tEnd\tRef\tAlt\tvalue").unwrap();
    for i in (1..22000).rev() {
        writeln!(gzip, "1\t{i}\t{i}\tA\tC\t{}", "x".repeat(500)).unwrap();
    }
    gzip.finish().unwrap();
    let manifest = DiskFilter::build(&source, Some(d.path())).unwrap();
    let disk = DiskFilter::open(&source, true).unwrap();
    let v = Variant::new("1", 999, 1000, "A", "C").unwrap();
    assert_eq!(
        disk.load_batch(std::slice::from_ref(&v))
            .unwrap()
            .annotate(&v, "x")
            .unwrap()
            .values,
        ["x".repeat(500)]
    );
    let json: serde_json::Value =
        serde_json::from_reader(fs::File::open(manifest).unwrap()).unwrap();
    let data = d.path().join(json["data_file"].as_str().unwrap());
    let mut bytes = fs::read(&data).unwrap();
    bytes[0] ^= 1;
    fs::write(data, bytes).unwrap();
    assert!(DiskFilter::open(&source, true).is_err());
    let fresh = DiskFilter::open(&source, false).unwrap();
    let first = Variant::new("1", 0, 1, "A", "C").unwrap();
    assert!(fresh.load_batch(&[first]).is_err());
}

#[test]
fn coding_genes_exclude_noncoding_isoforms_and_ncrna_exons_outrank_introns() {
    use rust_annovar::gene::GeneDatabase;
    let dir = tempdir().unwrap();
    let model = dir.path().join("genes.txt");
    fs::write(
        &model,
        concat!(
            "0\tNM_TEST\tchr1\t+\t0\t100\t0\t100\t2\t0,80,\t20,100,\t0\tCODING\n",
            "0\tNR_TEST\tchr1\t+\t0\t100\t100\t100\t1\t0,\t100,\t0\tCODING\n",
            "0\tNR_INTRON\tchr2\t+\t0\t100\t100\t100\t2\t0,80,\t20,100,\t0\tRNA_INTRON\n",
            "0\tNR_EXON\tchr2\t+\t0\t100\t100\t100\t1\t0,\t100,\t0\tRNA_EXON\n"
        ),
    )
    .unwrap();
    let db = GeneDatabase::load(&model, None).unwrap();
    for (chrom, function, gene) in [
        ("1", "intronic", "CODING"),
        ("2", "ncRNA_exonic", "RNA_EXON"),
    ] {
        let hit = db.annotate(
            &Variant::new(chrom, 40, 41, "A", "G").unwrap(),
            "refGene",
            2,
            1000,
        );
        assert_eq!(&hit.values[..2], &[function.to_string(), gene.to_string()]);
    }
}
