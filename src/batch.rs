//! Compact ownership while input batches are being accumulated.
use crate::model::Variant;
use std::collections::HashMap;
#[derive(Debug, Clone, Copy)]
pub struct ChromId(u32);
#[derive(Debug, Clone, Copy)]
pub struct VariantId {
    pub source_line: usize,
    pub allele_index: usize,
}
#[derive(Debug)]
struct PackedVariant {
    id: VariantId,
    chrom: ChromId,
    output_chrom: ChromId,
    start: u64,
    end: u64,
    reference: std::ops::Range<usize>,
    alternate: std::ops::Range<usize>,
    extra: Vec<String>,
    source_record: Option<usize>,
}
#[derive(Debug, Default)]
pub struct VariantBatch {
    names: Vec<String>,
    dictionary: HashMap<String, ChromId>,
    alleles: Vec<u8>,
    variants: Vec<PackedVariant>,
}
impl VariantBatch {
    fn intern(&mut self, s: String) -> ChromId {
        if let Some(id) = self.dictionary.get(&s) {
            return *id;
        }
        let id = ChromId(
            self.names
                .len()
                .try_into()
                .expect("batch chromosome count exceeds u32"),
        );
        self.names.push(s.clone());
        self.dictionary.insert(s, id);
        id
    }
    fn allele(&mut self, s: String) -> std::ops::Range<usize> {
        let start = self.alleles.len();
        self.alleles.extend_from_slice(s.as_bytes());
        start..self.alleles.len()
    }
    pub fn push(&mut self, v: Variant) {
        let chrom = self.intern(v.chrom);
        let output_chrom = self.intern(v.output_chrom);
        let reference = self.allele(v.reference);
        let alternate = self.allele(v.alternate);
        self.variants.push(PackedVariant {
            id: VariantId {
                source_line: v.source_line,
                allele_index: v.allele_index,
            },
            chrom,
            output_chrom,
            start: v.start,
            end: v.end,
            reference,
            alternate,
            extra: v.extra,
            source_record: v.source_record,
        });
    }
    /// Compatibility adapter for the public annotation API. The packed storage
    /// is released before annotation results accumulate.
    pub fn into_variants(self) -> Vec<Variant> {
        self.variants
            .into_iter()
            .map(|v| Variant {
                chrom: self.names[v.chrom.0 as usize].clone(),
                output_chrom: self.names[v.output_chrom.0 as usize].clone(),
                start: v.start,
                end: v.end,
                reference: String::from_utf8(self.alleles[v.reference].to_vec())
                    .expect("UTF-8 allele"),
                alternate: String::from_utf8(self.alleles[v.alternate].to_vec())
                    .expect("UTF-8 allele"),
                extra: v.extra,
                source_line: v.id.source_line,
                source_record: v.source_record,
                allele_index: v.id.allele_index,
            })
            .collect()
    }
}
