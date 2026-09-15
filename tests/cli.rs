use std::fs;
use std::process::Command;
use tempfile::tempdir;

#[test]
fn table_cli_annotates_vcf_and_preserves_sample_columns() {
    let dir = tempdir().unwrap();
    let db_dir = dir.path().join("humandb");
    fs::create_dir(&db_dir).unwrap();
    fs::write(
        db_dir.join("hg38_clinvar.txt"),
        "#Chr\tStart\tEnd\tRef\tAlt\tCLNSIG\n1\t10\t10\tA\tC\tPathogenic\n",
    )
    .unwrap();
    let input = dir.path().join("sample.vcf");
    fs::write(
        &input,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t10\trs1\tA\tC\t60\tPASS\tDP=20\tGT\t0/1\n",
    )
    .unwrap();
    let table = dir.path().join("output.tsv");
    let vcf = dir.path().join("output.vcf");
    let status = Command::new(env!("CARGO_BIN_EXE_rust-annovar"))
        .args([
            "table",
            input.to_str().unwrap(),
            db_dir.to_str().unwrap(),
            "--build",
            "hg38",
            "--protocol",
            "clinvar",
            "--operation",
            "f",
            "--vcf-input",
            "--output",
            table.to_str().unwrap(),
            "--vcf-output",
            vcf.to_str().unwrap(),
        ])
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(
        fs::read_to_string(table).unwrap(),
        "Chr\tStart\tEnd\tRef\tAlt\tCLNSIG.clinvar\tAnnotationStatus\n1\t10\t10\tA\tC\tPathogenic\tok\n"
    );
    let annotated_vcf = fs::read_to_string(vcf).unwrap();
    assert!(annotated_vcf.contains("##INFO=<ID=FA_clinvar"));
    assert!(annotated_vcf.contains("DP=20;FA_clinvar=Pathogenic;FA_STATUS=ok\tGT\t0/1"));
}

#[test]
fn invalid_protocol_lists_fail_without_partial_output() {
    let dir = tempdir().unwrap();
    let input = dir.path().join("input.avinput");
    fs::write(&input, "1\t10\t10\tA\tC\n").unwrap();
    let output = dir.path().join("output.tsv");
    let result = Command::new(env!("CARGO_BIN_EXE_rust-annovar"))
        .args([
            "table",
            input.to_str().unwrap(),
            dir.path().to_str().unwrap(),
            "--build",
            "hg38",
            "--protocol",
            "one,two",
            "--operation",
            "f",
            "--output",
            output.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(!output.exists());
    assert!(String::from_utf8_lossy(&result.stderr).contains("counts differ"));
}
