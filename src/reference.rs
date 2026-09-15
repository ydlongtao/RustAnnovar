//! Reference-backed minimal representation and left alignment.
use crate::model::{Variant, normalize_chrom};
use anyhow::{Context, Result, bail};
use sha2::{Digest, Sha256};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Read, Seek, SeekFrom},
    path::Path,
    sync::Mutex,
};
#[derive(Debug)]
struct Entry {
    len: u64,
    offset: u64,
    bases: u64,
    bytes: u64,
}
#[derive(Debug)]
pub struct ReferenceGenome {
    file: Mutex<File>,
    entries: HashMap<String, Entry>,
    pub sha256: String,
}
impl ReferenceGenome {
    pub fn open(path: &Path) -> Result<Self> {
        let index = std::path::PathBuf::from(format!("{}.fai", path.display()));
        let mut entries = HashMap::new();
        for line in BufReader::new(
            File::open(&index).context("reference requires a samtools-compatible .fai index")?,
        )
        .lines()
        {
            let line = line?;
            let f: Vec<_> = line.split('\t').collect();
            if f.len() < 5 {
                bail!("invalid FASTA index");
            }
            let e = Entry {
                len: f[1].parse()?,
                offset: f[2].parse()?,
                bases: f[3].parse()?,
                bytes: f[4].parse()?,
            };
            if e.bases == 0 || e.bytes < e.bases {
                bail!("invalid FASTA line width");
            }
            if entries.insert(normalize_chrom(f[0]), e).is_some() {
                bail!("ambiguous reference chromosome aliases");
            }
        }
        let mut h = Sha256::new();
        std::io::copy(&mut File::open(path)?, &mut h)?;
        Ok(Self {
            file: Mutex::new(File::open(path)?),
            entries,
            sha256: format!("{:x}", h.finalize()),
        })
    }
    pub fn sequence(&self, chrom: &str, start: u64, end: u64) -> Result<Vec<u8>> {
        let e = self
            .entries
            .get(chrom)
            .context("chromosome absent from reference")?;
        if end < start || end > e.len {
            bail!("reference query out of bounds");
        }
        let mut file = self
            .file
            .lock()
            .map_err(|_| anyhow::anyhow!("reference lock poisoned"))?;
        let mut out = Vec::with_capacity(usize::try_from(end - start)?);
        let mut pos = start;
        while pos < end {
            let n = (e.bases - pos % e.bases).min(end - pos);
            let offset = e.offset + pos / e.bases * e.bytes + pos % e.bases;
            file.seek(SeekFrom::Start(offset))?;
            let old = out.len();
            out.resize(old + n as usize, 0);
            file.read_exact(&mut out[old..])?;
            pos += n;
        }
        out.make_ascii_uppercase();
        Ok(out)
    }
    pub fn normalize(&self, v: &mut Variant) -> Result<()> {
        if !v
            .reference
            .bytes()
            .chain(v.alternate.bytes())
            .all(|b| b"ACGTN".contains(&b))
        {
            bail!("reference normalization requires explicit sequence alleles");
        }
        if self.sequence(&v.chrom, v.start, v.end)? != v.reference.as_bytes() {
            bail!("reference allele mismatch at {}:{}", v.chrom, v.start + 1);
        }
        let mut r = v.reference.as_bytes().to_vec();
        let mut a = v.alternate.as_bytes().to_vec();
        let mut pos = v.start;
        // Extend on the left whenever suffix trimming exhausts one allele.
        // This handles insertions, deletions and unequal-length replacements.
        loop {
            while !r.is_empty() && !a.is_empty() && r.last() == a.last() {
                r.pop();
                a.pop();
            }
            if (!r.is_empty() && !a.is_empty()) || pos == 0 || (r.is_empty() && a.is_empty()) {
                break;
            }
            let base = self.sequence(&v.chrom, pos - 1, pos)?[0];
            r.insert(0, base);
            a.insert(0, base);
            pos -= 1;
        }
        let prefix = r.iter().zip(&a).take_while(|(x, y)| x == y).count();
        pos += prefix as u64;
        r.drain(..prefix);
        a.drain(..prefix);
        v.start = pos;
        v.end = pos + r.len() as u64;
        v.reference = String::from_utf8(r)?;
        v.alternate = String::from_utf8(a)?;
        Ok(())
    }
}
