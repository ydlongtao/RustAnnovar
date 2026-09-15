//! ANNOVAR text ingestion and a chromosome-partitioned sorted exact-match index.
use anyhow::{Context, Result, bail};
use rustannovar_core::{AnnotationResult, Chromosome, Locus, Variant, fields};
use std::{collections::BTreeMap, io::BufRead};

#[derive(Debug)]
struct Entry {
    locus: Locus,
    annotations: std::ops::Range<usize>,
}
#[derive(Debug)]
pub struct Database {
    pub headers: Vec<String>,
    by_chrom: BTreeMap<Chromosome, Vec<Entry>>,
    // Keys carry offsets, not cloned annotation strings. Pool rows retain source order.
    pool: Vec<Vec<Box<str>>>,
    row_ids: Vec<usize>,
}
impl Database {
    pub fn load(mut reader: impl BufRead) -> Result<Self> {
        let mut raw: BTreeMap<Chromosome, Vec<(Locus, usize)>> = BTreeMap::new();
        let mut pool = Vec::new();
        let mut headers = Vec::new();
        let mut width = None;
        let mut buffer = String::new();
        let mut line_no = 0;
        loop {
            buffer.clear();
            if reader.read_line(&mut buffer)? == 0 {
                break;
            }
            line_no += 1;
            let line = buffer.trim_end_matches(['\r', '\n']);
            if line.trim().is_empty() || line.starts_with("##") {
                continue;
            }
            if let Some(line) = line.strip_prefix('#') {
                if width.is_some() || !pool.is_empty() {
                    bail!("database line {line_no}: unexpected repeated/late schema header");
                }
                let f: Vec<_> = fields(line).collect();
                if f.len() < 5 {
                    bail!("database line {line_no}: expected five coordinate header columns");
                }
                headers = f[5..].iter().map(|s| s.to_string()).collect();
                width = Some(headers.len());
                continue;
            }
            let variant = Variant::parse_avinput(line, pool.len() as u64)
                .with_context(|| format!("database line {line_no}"))?;
            let values: Vec<Box<str>> = fields(line).skip(5).map(Box::from).collect();
            let expected = *width.get_or_insert(values.len());
            if expected != values.len() {
                bail!(
                    "database line {line_no}: expected {expected} annotation fields, found {}",
                    values.len()
                );
            }
            raw.entry(variant.chrom)
                .or_default()
                .push((variant.locus, pool.len()));
            pool.push(values);
        }
        let width = width.unwrap_or(0);
        if headers.is_empty() {
            headers = (0..width)
                .map(|i| {
                    if width == 1 {
                        "value".into()
                    } else {
                        format!("value{}", i + 1)
                    }
                })
                .collect();
        }
        if width == 0 {
            headers.push("match".into());
        }
        let mut seen = std::collections::HashSet::new();
        for h in &headers {
            if h.is_empty()
                || !seen.insert(h)
                || ["Chr", "Start", "End", "Ref", "Alt"].contains(&h.as_str())
            {
                bail!("invalid or duplicate database annotation header {h:?}");
            }
        }
        let mut row_ids = Vec::with_capacity(pool.len());
        let mut by_chrom = BTreeMap::new();
        for (chrom, mut rows) in raw {
            rows.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            let mut entries: Vec<Entry> = Vec::new();
            for (locus, id) in rows {
                row_ids.push(id);
                if let Some(last) = entries.last_mut().filter(|e| e.locus == locus) {
                    last.annotations.end += 1;
                } else {
                    entries.push(Entry {
                        locus,
                        annotations: row_ids.len() - 1..row_ids.len(),
                    });
                }
            }
            by_chrom.insert(chrom, entries);
        }
        Ok(Self {
            headers,
            by_chrom,
            pool,
            row_ids,
        })
    }
    pub fn lookup(&self, variant: &Variant) -> AnnotationResult<'_> {
        let Some(entries) = self.by_chrom.get(&variant.chrom) else {
            return AnnotationResult::Missing;
        };
        match entries.binary_search_by(|e| e.locus.cmp(&variant.locus)) {
            Ok(i) => AnnotationResult::Matches(&self.row_ids[entries[i].annotations.clone()]),
            Err(_) => AnnotationResult::Missing,
        }
    }
    pub fn row(&self, id: usize) -> &[Box<str>] {
        &self.pool[id]
    }
    pub fn records(&self) -> usize {
        self.pool.len()
    }
    pub fn chromosomes(&self) -> usize {
        self.by_chrom.len()
    }
}
