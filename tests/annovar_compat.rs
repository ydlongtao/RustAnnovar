use rust_annovar::database::{FilterDatabase, RegionDatabase};
use rust_annovar::io::read_avinput;
use rust_annovar::pipeline::{Operation, Protocol, annotate_table};
use std::collections::BTreeSet;
use std::fs;
use std::path::PathBuf;

/// This test uses the locally registered ANNOVAR package and its 270 MB
/// transcript FASTA. It is ignored in routine CI but is the release oracle for
/// gene annotation compatibility.
#[test]
#[ignore = "requires the registered local annovar/ package"]
fn refgene_example_matches_perl_for_core_snv_consequences() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let annovar = root.join("annovar");
    assert!(
        annovar.join("annotate_variation.pl").exists(),
        "place the registered ANNOVAR package at annovar/"
    );
    let variants = read_avinput(&annovar.join("example/ex1.avinput")).unwrap();
    let result = annotate_table(
        &variants,
        &[Protocol {
            name: "refGene".into(),
            operation: Operation::Gene,
            database: annovar.join("humandb/hg19_refGene.txt"),
            fasta: Some(annovar.join("humandb/hg19_refGeneMrna.fa")),
        }],
        ".",
    )
    .unwrap();
    let expected_text =
        fs::read_to_string(root.join("tests/fixtures/perl_expected/ex1.hg19_multianno.txt"))
            .unwrap();
    let expected = expected_text
        .lines()
        .skip(1)
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .collect::<Vec<_>>();
    assert_eq!(result.rows.len(), expected.len());

    for (variant, (actual, expected)) in variants.iter().zip(result.rows.iter().zip(expected)) {
        assert_eq!(&actual[..5], &expected[..5], "coordinate mismatch");
        assert_eq!(
            actual[5], expected[5],
            "function mismatch at source line {}",
            variant.source_line
        );
        assert_eq!(
            as_set(&actual[6]),
            as_set(expected[6]),
            "gene mismatch at source line {}",
            variant.source_line
        );
        if variant.reference.len() == 1
            && variant.reference != "0"
            && variant.alternate.len() == 1
            && actual[5] == "exonic"
        {
            assert_eq!(
                actual[8], expected[8],
                "SNV class mismatch at source line {}",
                variant.source_line
            );
            assert_eq!(
                as_set(&actual[9]),
                as_set(expected[9]),
                "AAChange mismatch at source line {}",
                variant.source_line
            );
        }
    }
}

#[test]
#[ignore = "requires the registered local annovar/ package"]
fn bundled_filter_and_gff3_examples_match_perl() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let annovar = root.join("annovar");
    let variants = read_avinput(&annovar.join("example/ex1.avinput")).unwrap();

    let filter =
        FilterDatabase::load(&annovar.join("humandb/hg19_example_db_generic.txt")).unwrap();
    let filter_hit = filter.annotate(&variants[11], "generic").unwrap();
    assert_eq!(filter_hit.values, ["0.05"]);

    let region = RegionDatabase::load(&annovar.join("humandb/hg19_example_db_gff3.txt")).unwrap();
    let region_hit = region.annotate(&variants[14], "gff3", 0.0).unwrap();
    assert_eq!(region_hit.values, ["Score=843;Name=20975115"]);
}

fn as_set(value: &str) -> BTreeSet<&str> {
    value.split(',').collect()
}
