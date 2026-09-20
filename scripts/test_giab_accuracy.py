#!/usr/bin/env python3
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from giab_accuracy import compare, compare_conversion, variant_type


HEADER = (
    "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\t"
    "GeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\n"
)


class AccuracyComparatorTest(unittest.TestCase):
    def test_duplicate_occurrences_and_list_order(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rust = root / "rust.tsv"
            perl = root / "perl.tsv"
            report = root / "comparison.json"
            prefix = "1\t100\t100\tA\tC\texonic\t"
            rust.write_text(
                HEADER
                + prefix + "B,A\t.\tnonsynonymous SNV\tB:x,A:y\n"
                + prefix + "C\t.\tsynonymous SNV\tC:z\n"
            )
            perl.write_text(
                HEADER
                + prefix + "A,B\t.\tnonsynonymous SNV\tA:y,B:x\n"
                + prefix + "C\t.\tstopgain\tC:z\n"
            )
            self.assertEqual(compare(SimpleNamespace(rust=rust, perl=perl, report=report, max_examples=5)), 1)
            result = json.loads(report.read_text())
            self.assertEqual(result["matched_variant_occurrences"], 2)
            self.assertEqual(result["missing_in_perl"], 0)
            self.assertEqual(result["field_counts"]["Gene.refGene"]["order_insensitive"], 1)
            self.assertEqual(result["field_counts"]["AAChange.refGene"]["order_insensitive"], 1)
            self.assertEqual(result["field_counts"]["ExonicFunc.refGene"]["different"], 1)
            self.assertEqual(result["rows_all_fields_order_insensitive_by_type"]["snv"], 1)

    def test_variant_type(self):
        self.assertEqual(variant_type("A", "C"), "snv")
        self.assertEqual(variant_type("AC", "GT"), "mnv")
        self.assertEqual(variant_type("A", "AG"), "indel")
        self.assertEqual(variant_type("A", "C,G"), "multiallelic")
        self.assertEqual(variant_type("A", "<DEL>"), "other")

    def test_conversion_compares_duplicate_multiplicity(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rust = root / "rust.avinput"
            perl = root / "perl.avinput"
            report = root / "conversion.json"
            rust.write_text("1\t10\t10\tA\tC\n1\t10\t10\tA\tC\n")
            perl.write_text("1\t10\t10\tA\tC\n")
            self.assertEqual(
                compare_conversion(SimpleNamespace(rust=rust, perl=perl, report=report, max_examples=5)), 1
            )
            result = json.loads(report.read_text())
            self.assertEqual(result["only_rust"], 1)
            self.assertFalse(result["multiset_equal"])

    def test_conversion_reports_chr_prefix_separately(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            rust = root / "rust.avinput"
            perl = root / "perl.avinput"
            report = root / "conversion.json"
            rust.write_text("1\t10\t10\tA\tC\n")
            perl.write_text("chr1\t10\t10\tA\tC\n")
            self.assertEqual(
                compare_conversion(SimpleNamespace(rust=rust, perl=perl, report=report, max_examples=5)), 1
            )
            result = json.loads(report.read_text())
            self.assertEqual(result["only_rust_after_chr_prefix_normalization"], 0)
            self.assertTrue(result["multiset_equal_after_chr_prefix_normalization"])


if __name__ == "__main__":
    unittest.main()
