use rustannovar_core::{AnnotationResult, Variant};
use rustannovar_db::Database;
use std::io::Cursor;
#[test]
fn sorted_index_matches_brute_force_and_preserves_duplicate_source_order() {
    let mut lines = Vec::new();
    for i in (1..=2000).rev() {
        lines.push(format!(
            "{}\t{}\t{}\tA\t{}\trow{i}",
            i % 3 + 1,
            i % 71 + 1,
            i % 71 + 1,
            if i % 2 == 0 { "C" } else { "G" }
        ));
    }
    lines.extend([
        "MT\t7\t7\tA\tC\tmitochondrial".into(),
        "X\t10\t10\t-\tAC\tins".into(),
        "Y\t12\t13\tAC\t-\tdel".into(),
    ]);
    let db = Database::load(Cursor::new(lines.join("\n"))).unwrap();
    let variants: Vec<_> = lines
        .iter()
        .enumerate()
        .map(|(i, l)| Variant::parse_avinput(l, i as u64).unwrap())
        .collect();
    for q in variants
        .iter()
        .step_by(7)
        .chain(variants.iter().rev().take(3))
    {
        let expected: Vec<_> = variants
            .iter()
            .enumerate()
            .filter(|(_, v)| v.chrom == q.chrom && v.locus == q.locus)
            .map(|(i, _)| i)
            .collect();
        match db.lookup(q) {
            AnnotationResult::Matches(actual) => assert_eq!(actual, expected),
            _ => panic!("missing exact match"),
        }
    }
    for q in ["1 7 7 A T", "1 7 8 AC C", "2 10 10 - AC"] {
        assert!(matches!(
            db.lookup(&Variant::parse_avinput(q, 0).unwrap()),
            AnnotationResult::Missing
        ));
    }
    assert_eq!(db.records(), lines.len());
}
#[test]
fn database_schema_is_strict_and_error_identifies_line() {
    for text in [
        "#Chr\tStart\tEnd\tRef\tAlt\tAF\n1\t1\t1\tA\tC\t.1\textra",
        "1\t12A\t12\tA\tC\tx",
        "#Chr\tStart\tEnd\tRef\tAlt\tAF\tAF\n1\t1\t1\tA\tC\t.1\t.2",
    ] {
        assert!(Database::load(Cursor::new(text)).is_err());
    }
    let e = Database::load(Cursor::new("##metadata\n1\t12A\t12\tA\tC\tx")).unwrap_err();
    assert!(format!("{e:#}").contains("database line 2"));
}
