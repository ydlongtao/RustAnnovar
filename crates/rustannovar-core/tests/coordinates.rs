use rustannovar_core::{Allele, Chromosome, Variant};
#[test]
fn chromosomes_and_alleles_are_normalized_without_conflating_contigs() {
    for (a, b) in [
        ("1", "chr1"),
        ("X", "chrX"),
        ("Y", "chrY"),
        ("M", "chrM"),
        ("chrM", "MT"),
        ("chrMT", "MT"),
    ] {
        assert_eq!(Chromosome::parse(a).unwrap(), Chromosome::parse(b).unwrap());
    }
    assert_ne!(
        Chromosome::parse("01").unwrap(),
        Chromosome::parse("1").unwrap()
    );
    assert_ne!(
        Chromosome::parse("GL0001").unwrap(),
        Chromosome::parse("GL0002").unwrap()
    );
    assert_eq!(Allele::parse("a", true).unwrap(), Allele::A);
    assert_eq!(
        Allele::parse("acn", true).unwrap(),
        Allele::parse("ACN", true).unwrap()
    );
    assert!(Allele::parse("<DEL>", false).is_err());
}
#[test]
fn avinput_coordinates_and_invalid_records() {
    let snv = Variant::parse_avinput("chr1\t10\t10\tA\tC", 7).unwrap();
    assert_eq!((snv.locus.start, snv.locus.end, snv.row_id), (9, 10, 7));
    let ins = Variant::parse_avinput("1 10 10 - AT", 0).unwrap();
    assert_eq!((ins.locus.start, ins.locus.end), (10, 10));
    assert!(Variant::parse_avinput("1 10 11 AC GT", 0).is_ok());
    assert!(Variant::parse_avinput("1 10 100 0 -", 0).is_ok());
    for line in [
        "1 1",
        "1 12A 12 A C",
        "1 0 0 A C",
        "1 5 4 A C",
        "1 5 6 A C",
        "1 5 6 - A",
        "1 5 5 - -",
        "1 5 5 A C,G",
    ] {
        assert!(
            Variant::parse_avinput(line, 0).is_err(),
            "accepted malformed record {line}"
        );
    }
}
