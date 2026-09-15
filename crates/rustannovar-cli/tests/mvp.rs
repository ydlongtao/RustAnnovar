use std::{fs, process::Command};
use tempfile::tempdir;
fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rustannovar"))
}
#[test]
fn output_retains_order_extra_fields_missing_and_duplicate_annotations() {
    let dir = tempdir().unwrap();
    let p = dir.path();
    fs::write(p.join("db"),"#Chr\tStart\tEnd\tRef\tAlt\tCLNSIG\n1\t10\t10\tA\tC\tfirst\n1\t10\t10\tA\tG\tother-alt\n1\t10\t10\tA\tC\tsecond\nMT\t5\t5\t-\tAC\tinsertion\nY\t7\t8\tAC\t-\tdeletion\nX\t3\t4\tAC\tGT\tmnv\n").unwrap();
    let input = "chrM\t5\t5\t-\tAC\tsample1\nchr1\t10\t10\tA\tC\tsample2\nY\t7\t8\tAC\t-\tsample3\n1\t10\t10\tA\tT\tsample4\nX\t3\t4\tAC\tGT\tsample5\n";
    fs::write(p.join("in"), input.repeat(2000)).unwrap();
    for (threads, batch, out) in [(1, 1, "a"), (4, 10000, "b")] {
        let r = binary()
            .args(["annotate", "--input"])
            .arg(p.join("in"))
            .arg("--database")
            .arg(p.join("db"))
            .arg("--output")
            .arg(p.join(out))
            .args([
                "--threads",
                &threads.to_string(),
                "--batch-size",
                &batch.to_string(),
                "--nastring",
                "NA",
                "--quiet",
            ])
            .output()
            .unwrap();
        assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
        assert!(r.stderr.is_empty());
    }
    let text = fs::read_to_string(p.join("a")).unwrap();
    assert_eq!(text, fs::read_to_string(p.join("b")).unwrap());
    assert_eq!(
        text.lines().take(6).collect::<Vec<_>>(),
        vec![
            "Chr\tStart\tEnd\tRef\tAlt\tCLNSIG\tOtherinfo1",
            "chrM\t5\t5\t-\tAC\tinsertion\tsample1",
            "chr1\t10\t10\tA\tC\tfirst;second\tsample2",
            "Y\t7\t8\tAC\t-\tdeletion\tsample3",
            "1\t10\t10\tA\tT\tNA\tsample4",
            "X\t3\t4\tAC\tGT\tmnv\tsample5"
        ]
    );
}
#[test]
fn late_parse_error_preserves_destination_and_reports_location() {
    let d = tempdir().unwrap();
    let p = d.path();
    fs::write(p.join("db"), "1\t1\t1\tA\tC\tx\n").unwrap();
    fs::write(p.join("in"), "1\t1\t1\tA\tC\n1\t12A\t12\tA\tC\n").unwrap();
    fs::write(p.join("out"), "keep me").unwrap();
    let r = binary()
        .args(["annotate", "--input"])
        .arg(p.join("in"))
        .arg("--database")
        .arg(p.join("db"))
        .arg("--output")
        .arg(p.join("out"))
        .args(["--batch-size", "1"])
        .output()
        .unwrap();
    assert!(!r.status.success());
    assert!(String::from_utf8_lossy(&r.stderr).contains("input line 2"));
    assert_eq!(fs::read_to_string(p.join("out")).unwrap(), "keep me");
}
#[test]
fn pipeline_stdout_and_match_only_database_work() {
    let d = tempdir().unwrap();
    let p = d.path();
    fs::write(p.join("db"), "1\t1\t1\tA\tC\n").unwrap();
    fs::write(p.join("in"), "1\t1\t1\tA\tC\n").unwrap();
    let r = binary()
        .args(["annotate", "--input"])
        .arg(p.join("in"))
        .arg("--database")
        .arg(p.join("db"))
        .args(["--output", "-", "--quiet"])
        .output()
        .unwrap();
    assert!(r.status.success());
    assert_eq!(
        String::from_utf8(r.stdout).unwrap(),
        "Chr\tStart\tEnd\tRef\tAlt\tmatch\n1\t1\t1\tA\tC\t1\n"
    );
}

#[test]
fn raw_tsv_fast_path_and_escaping_fallback_have_correct_columns() {
    let d = tempdir().unwrap();
    let p = d.path();
    fs::write(p.join("db"), "1\t1\t1\tA\tC\tvalue\n").unwrap();
    for (i, input, expected) in [
        (0, "1\t1\t1\tA\tC\t\n", "1\t1\t1\tA\tC\tvalue\t\n"),
        (1, "1 1 1 A C sample\n", "1\t1\t1\tA\tC\tvalue\tsample\n"),
        (
            2,
            "1\t1\t1\tA\tC\t\"sample\"\n",
            "1\t1\t1\tA\tC\tvalue\t\"\"\"sample\"\"\"\n",
        ),
    ] {
        fs::write(p.join("in"), input).unwrap();
        let r = binary()
            .args(["annotate", "--input"])
            .arg(p.join("in"))
            .arg("--database")
            .arg(p.join("db"))
            .args(["--quiet"])
            .output()
            .unwrap();
        assert!(r.status.success(), "case {i}");
        assert_eq!(
            String::from_utf8(r.stdout).unwrap(),
            format!("Chr\tStart\tEnd\tRef\tAlt\tvalue\tOtherinfo1\n{expected}")
        );
    }
}
