//! Compact exact-match keys. Coordinates are zero-based and half-open;
//! empty reference alleles denote interbase insertion sites.
use anyhow::{Result, bail};
use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Chromosome {
    Numeric(u8),
    X,
    Y,
    MT,
    Other(Box<str>),
}
impl Chromosome {
    pub fn parse(s: &str) -> Result<Self> {
        let s = s.strip_prefix("chr").unwrap_or(s);
        if s.is_empty() {
            bail!("chromosome must not be empty");
        }
        Ok(match s {
            "X" => Self::X,
            "Y" => Self::Y,
            "M" | "MT" => Self::MT,
            _ => match s.parse::<u8>() {
                Ok(n) if n > 0 && !s.starts_with('0') => Self::Numeric(n),
                _ => Self::Other(s.into()),
            },
        })
    }
}
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Allele {
    Empty,
    A,
    C,
    G,
    T,
    UnknownReference,
    Sequence(Box<[u8]>),
}
impl Allele {
    pub fn parse(s: &str, reference: bool) -> Result<Self> {
        Ok(match s {
            "-" => Self::Empty,
            "A" | "a" => Self::A,
            "C" | "c" => Self::C,
            "G" | "g" => Self::G,
            "T" | "t" => Self::T,
            "0" if reference => Self::UnknownReference,
            _ if !s.is_empty() && s.bytes().all(|b| b"ACGTNacgtn".contains(&b)) => {
                Self::Sequence(s.to_ascii_uppercase().into_bytes().into_boxed_slice())
            }
            _ => bail!("unsupported allele {s:?}; expected A/C/G/T/N sequence or '-'"),
        })
    }
    pub fn bases(&self) -> &[u8] {
        match self {
            Self::Empty => b"",
            Self::A => b"A",
            Self::C => b"C",
            Self::G => b"G",
            Self::T => b"T",
            Self::UnknownReference => b"0",
            Self::Sequence(s) => s,
        }
    }
}
impl fmt::Display for Allele {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if matches!(self, Self::Empty) {
            return f.write_str("-");
        }
        f.write_str(std::str::from_utf8(self.bases()).map_err(|_| fmt::Error)?)
    }
}
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Locus {
    pub start: u64,
    pub end: u64,
    pub reference: Allele,
    pub alternate: Allele,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Variant {
    pub chrom: Chromosome,
    pub locus: Locus,
    pub row_id: u64,
}

pub fn fields(line: &str) -> impl Iterator<Item = &str> {
    let tabs = line.contains('\t');
    line.split(move |c: char| {
        if tabs {
            c == '\t'
        } else {
            c.is_ascii_whitespace()
        }
    })
    .filter(move |s| tabs || !s.is_empty())
}
impl Variant {
    pub fn parse_avinput(line: &str, row_id: u64) -> Result<Self> {
        let mut f = fields(line);
        let mut next = || {
            f.next()
                .ok_or_else(|| anyhow::anyhow!("expected at least five AVinput columns"))
        };
        let chrom = Chromosome::parse(next()?)?;
        let start_text = next()?;
        let start: u64 = start_text
            .parse()
            .map_err(|_| anyhow::anyhow!("invalid start position {start_text:?}"))?;
        let end_text = next()?;
        let end: u64 = end_text
            .parse()
            .map_err(|_| anyhow::anyhow!("invalid end position {end_text:?}"))?;
        let reference = Allele::parse(next()?, true)?;
        let alternate = Allele::parse(next()?, false)?;
        let locus = if reference == Allele::Empty {
            if start != end || alternate == Allele::Empty {
                bail!("insertions require Start=End and a nonempty ALT");
            }
            Locus {
                start,
                end,
                reference,
                alternate,
            }
        } else {
            if start == 0 || end < start {
                bail!("invalid one-based interval {start}..{end}");
            }
            if reference != Allele::UnknownReference
                && end - start + 1 != reference.bases().len() as u64
            {
                bail!("REF length does not match interval {start}..{end}");
            }
            Locus {
                start: start - 1,
                end,
                reference,
                alternate,
            }
        };
        Ok(Self {
            chrom,
            locus,
            row_id,
        })
    }
}

/// Borrowed result rows; formatting is deferred to the output layer.
#[derive(Debug)]
pub enum AnnotationResult<'a> {
    Missing,
    Matches(&'a [usize]),
}
