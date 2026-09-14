use anyhow::{Result, bail};
use serde::{Deserialize, Serialize};

/// Canonical internal representation. `start == end` represents an insertion
/// between bases; all other records use a zero-based half-open interval.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Variant {
    pub chrom: String,
    pub output_chrom: String,
    pub start: u64,
    pub end: u64,
    pub reference: String,
    pub alternate: String,
    pub extra: Vec<String>,
    pub source_line: usize,
    pub source_record: Option<usize>,
    pub allele_index: usize,
}

impl Variant {
    pub fn new(
        chrom: impl Into<String>,
        start: u64,
        end: u64,
        reference: impl Into<String>,
        alternate: impl Into<String>,
    ) -> Result<Self> {
        if end < start {
            bail!("variant end ({end}) precedes start ({start})");
        }
        let chrom = normalize_chrom(&chrom.into());
        Ok(Self {
            output_chrom: chrom.clone(),
            chrom,
            start,
            end,
            reference: reference.into().to_ascii_uppercase(),
            alternate: alternate.into().to_ascii_uppercase(),
            extra: Vec::new(),
            source_line: 0,
            source_record: None,
            allele_index: 0,
        })
    }

    pub fn key(&self) -> VariantKey<'_> {
        VariantKey {
            chrom: &self.chrom,
            start: self.start,
            end: self.end,
            reference: &self.reference,
            alternate: &self.alternate,
        }
    }

    pub fn avinput_fields(&self) -> [String; 5] {
        let (start, end, reference, alternate) = if self.start == self.end {
            (
                self.start,
                self.start,
                "-".to_string(),
                display_allele(&self.alternate),
            )
        } else {
            (
                self.start + 1,
                self.end,
                display_allele(&self.reference),
                display_allele(&self.alternate),
            )
        };
        [
            self.output_chrom.clone(),
            start.to_string(),
            end.to_string(),
            reference,
            alternate,
        ]
    }

    pub fn overlaps(&self, start: u64, end: u64) -> bool {
        if self.start == self.end {
            self.start >= start && self.start <= end
        } else {
            self.start < end && start < self.end
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct VariantKey<'a> {
    pub chrom: &'a str,
    pub start: u64,
    pub end: u64,
    pub reference: &'a str,
    pub alternate: &'a str,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnnotationKind {
    Filter,
    Region,
    Gene,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Annotation {
    pub protocol: String,
    pub kind: AnnotationKind,
    pub values: Vec<String>,
}

pub fn normalize_chrom(chrom: &str) -> String {
    chrom.strip_prefix("chr").unwrap_or(chrom).to_string()
}

fn display_allele(allele: &str) -> String {
    if allele.is_empty() {
        "-".to_string()
    } else {
        allele.to_string()
    }
}
