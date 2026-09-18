use crate::database::{FilterDatabase, RegionDatabase};
use crate::gene::GeneDatabase;
use crate::io::{VcfDocument, open_writer};
use crate::model::{Annotation, Variant};
use anyhow::{Context, Result, bail};
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operation {
    Gene,
    Region,
    Filter,
    Cadd,
}

impl Operation {
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "g" | "gx" | "gene" => Ok(Self::Gene),
            "r" | "region" => Ok(Self::Region),
            "f" | "filter" => Ok(Self::Filter),
            "cadd" => Ok(Self::Cadd),
            _ => bail!("unknown operation {value:?}; use g, r, f, or cadd"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Protocol {
    pub name: String,
    pub operation: Operation,
    pub database: PathBuf,
    pub fasta: Option<PathBuf>,
}

#[derive(Debug)]
pub struct TableResult {
    pub headers: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

pub fn resolve_protocols(
    db_dir: &Path,
    build: &str,
    names: &[String],
    operations: &[Operation],
) -> Result<Vec<Protocol>> {
    if names.len() != operations.len() {
        bail!("protocol and operation counts differ");
    }
    names
        .iter()
        .zip(operations)
        .map(|(name, operation)| {
            let mut database = db_dir.join(format!("{build}_{name}.txt"));
            if !database.exists() {
                database = db_dir.join(format!("{build}_{name}.txt.gz"));
            }
            if !database.exists() {
                bail!("database does not exist: {}", database.display());
            }
            if *operation == Operation::Cadd {
                crate::cadd::Database::open(&database)?.check_build(build)?;
            }
            let fasta = (*operation == Operation::Gene)
                .then(|| db_dir.join(format!("{build}_{name}Mrna.fa")))
                .filter(|p| p.exists());
            Ok(Protocol {
                name: name.clone(),
                operation: *operation,
                database,
                fasta,
            })
        })
        .collect()
}

enum Loaded {
    Cadd(crate::cadd::Database),
    Disk(crate::disk_index::DiskFilter),
    Filter(FilterDatabase),
    Region(RegionDatabase),
    Gene(GeneDatabase),
}

pub struct AnnotationEngine {
    protocols: Vec<Protocol>,
    loaded: Vec<Loaded>,
    headers: Vec<String>,
}
impl AnnotationEngine {
    pub fn new(protocols: &[Protocol]) -> Result<Self> {
        Self::new_normalized(protocols, None)
    }
    pub fn new_normalized(
        protocols: &[Protocol],
        reference: Option<&crate::reference::ReferenceGenome>,
    ) -> Result<Self> {
        let mut loaded = Vec::new();
        let mut headers = vec!["Chr", "Start", "End", "Ref", "Alt"]
            .into_iter()
            .map(str::to_string)
            .collect::<Vec<_>>();
        for protocol in protocols {
            let database = match protocol.operation {
                Operation::Cadd => {
                    let db = crate::cadd::Database::open(&protocol.database)?;
                    headers.extend(
                        crate::cadd::HEADERS
                            .iter()
                            .map(|h| qualified_header(h, &protocol.name)),
                    );
                    Loaded::Cadd(db)
                }
                Operation::Filter => {
                    if crate::disk_index::DiskFilter::exists(&protocol.database) {
                        let db = crate::disk_index::DiskFilter::open(&protocol.database, false)?;
                        db.check_normalization(reference)?;
                        headers.extend(
                            db.headers()
                                .iter()
                                .map(|h| qualified_header(h, &protocol.name)),
                        );
                        loaded.push(Loaded::Disk(db));
                        continue;
                    }
                    if reference.is_some() {
                        bail!(
                            "left-align filtering requires an index built with the same reference"
                        );
                    }
                    if std::fs::metadata(&protocol.database)?.len() > 64 * 1024 * 1024 {
                        bail!(
                            "filter database exceeds 64 MiB; build a .rai index with db index first"
                        );
                    }
                    let db = FilterDatabase::load(&protocol.database)?;
                    headers.extend(
                        db.headers
                            .iter()
                            .map(|header| qualified_header(header, &protocol.name)),
                    );
                    Loaded::Filter(db)
                }
                Operation::Region => {
                    let db = RegionDatabase::load(&protocol.database)?;
                    headers.extend(
                        db.headers
                            .iter()
                            .map(|header| qualified_header(header, &protocol.name)),
                    );
                    Loaded::Region(db)
                }
                Operation::Gene => {
                    headers.extend(
                        ["Func", "Gene", "GeneDetail", "ExonicFunc", "AAChange"]
                            .map(|header| format!("{header}.{}", protocol.name)),
                    );
                    Loaded::Gene(GeneDatabase::load(
                        &protocol.database,
                        protocol.fasta.as_deref(),
                    )?)
                }
            };
            loaded.push(database);
        }
        let mut seen = std::collections::HashSet::new();
        for name in &headers {
            if !seen.insert(name) {
                bail!("duplicate output column {name}");
            }
        }
        let mut info_ids = std::collections::HashSet::new();
        for protocol in protocols {
            if sanitize_info(&protocol.name) == "STATUS"
                || !info_ids.insert(sanitize_info(&protocol.name))
            {
                bail!("duplicate VCF protocol ID");
            }
        }
        Ok(Self {
            protocols: protocols.to_vec(),
            loaded,
            headers,
        })
    }
    pub fn ranges(&self) -> Vec<std::ops::Range<usize>> {
        let mut start = 5;
        self.loaded
            .iter()
            .map(|db| {
                let width = match db {
                    Loaded::Cadd(_) => 3,
                    Loaded::Disk(db) => db.headers().len(),
                    Loaded::Filter(db) => db.headers.len(),
                    Loaded::Region(db) => db.headers.len(),
                    Loaded::Gene(_) => 5,
                };
                let r = start..start + width;
                start += width;
                r
            })
            .collect()
    }
    pub fn headers(&self) -> &[String] {
        &self.headers
    }
    pub fn check_cadd_build(&self, build: &str) -> Result<()> {
        for db in &self.loaded {
            if let Loaded::Cadd(db) = db {
                db.check_build(build)?;
            }
        }
        Ok(())
    }
    pub fn check_cadd_version(&self, version: &str) -> Result<()> {
        for db in &self.loaded {
            if let Loaded::Cadd(db) = db {
                db.check_version(version)?;
            }
        }
        Ok(())
    }
    pub fn cadd_metadata(&self) -> std::collections::BTreeMap<String, serde_json::Value> {
        self.protocols
            .iter()
            .zip(&self.loaded)
            .filter_map(|(protocol, db)| {
                if let Loaded::Cadd(db) = db {
                    Some((protocol.name.clone(), db.metadata()))
                } else {
                    None
                }
            })
            .collect()
    }
    pub fn annotate(&self, variants: &[Variant], nastring: &str) -> Result<TableResult> {
        let mut headers = self.headers.clone();
        let protocols = &self.protocols;
        let loaded = &self.loaded;
        let prepared: Vec<Option<FilterDatabase>> = loaded
            .iter()
            .map(|db| match db {
                Loaded::Disk(db) => db.load_batch(variants).map(Some),
                _ => Ok(None),
            })
            .collect::<Result<_>>()?;
        let cadd_rows = loaded
            .iter()
            .map(|db| match db {
                Loaded::Cadd(db) => db.batch(variants, nastring).map(Some),
                _ => Ok(None),
            })
            .collect::<Result<Vec<_>>>()?;
        let extra_width = variants
            .iter()
            .map(|variant| variant.extra.len())
            .max()
            .unwrap_or(0);
        headers.extend((1..=extra_width).map(|index| format!("Otherinfo{index}")));
        let min_len = if protocols.iter().any(|p| p.operation == Operation::Gene) {
            128
        } else {
            8192
        };
        let rows = variants
            .par_iter()
            .with_min_len(min_len)
            .enumerate()
            .map(|(variant_index, variant)| {
                let mut row = variant.avinput_fields().to_vec();
                for (index, (protocol, database)) in protocols.iter().zip(loaded).enumerate() {
                    if let Loaded::Cadd(_) = database {
                        row.extend(cadd_rows[index].as_ref().unwrap()[variant_index].clone());
                        continue;
                    }
                    let unsupported = variant.source_record.is_some()
                        && !variant.alternate.is_empty()
                        && !crate::io::supported_alt(&variant.alternate);
                    let (annotation, width) = match database {
                        Loaded::Cadd(_) => unreachable!(),
                        Loaded::Disk(_) => {
                            let db = prepared[index].as_ref().unwrap();
                            (db.annotate(variant, &protocol.name), db.headers.len())
                        }
                        Loaded::Filter(db) => {
                            (db.annotate(variant, &protocol.name), db.headers.len())
                        }
                        Loaded::Region(db) => {
                            (db.annotate(variant, &protocol.name, 0.0), db.headers.len())
                        }
                        Loaded::Gene(db) => {
                            (Some(db.annotate(variant, &protocol.name, 2, 1000)), 5)
                        }
                    };
                    append_annotation(
                        &mut row,
                        if unsupported { None } else { annotation },
                        width,
                        nastring,
                    );
                }
                row.extend(variant.extra.clone());
                while row.len() < headers.len() {
                    row.push(nastring.to_string());
                }
                row
            })
            .collect();
        Ok(TableResult { headers, rows })
    }
}

pub fn annotate_table(
    variants: &[Variant],
    protocols: &[Protocol],
    nastring: &str,
) -> Result<TableResult> {
    AnnotationEngine::new(protocols)?.annotate(variants, nastring)
}

pub fn write_table(result: &TableResult, path: &Path, csv: bool) -> Result<()> {
    let mut writer = open_writer(path)?;
    let delimiter = if csv { ',' } else { '\t' };
    write_record(&mut writer, &result.headers, delimiter)?;
    for row in &result.rows {
        write_record(&mut writer, row, delimiter)?;
    }
    Ok(())
}

pub fn write_annotated_vcf(
    document: &VcfDocument,
    result: &TableResult,
    protocols: &[Protocol],
    path: &Path,
    nastring: &str,
) -> Result<()> {
    let mut writer = open_writer(path)?;
    for header in &document.headers {
        if header.starts_with("#CHROM") {
            for protocol in protocols {
                writeln!(
                    writer,
                    "##INFO=<ID=FA_{},Number=.,Type=String,Description=\"RustAnnovar {} annotation\">",
                    sanitize_info(&protocol.name),
                    protocol.name
                )?;
            }
        }
        writeln!(writer, "{header}")?;
    }
    let mut by_record: HashMap<usize, Vec<&Vec<String>>> = HashMap::new();
    for (variant, row) in document.variants.iter().zip(&result.rows) {
        if let Some(record) = variant.source_record {
            by_record.entry(record).or_default().push(row);
        }
    }
    let mut offsets = Vec::new();
    let mut cursor = 5usize;
    for protocol in protocols {
        let width = match protocol.operation {
            Operation::Gene => 5,
            _ => result.headers[cursor..]
                .iter()
                .take_while(|h| h.ends_with(&format!(".{}", protocol.name)))
                .count()
                .max(1),
        };
        offsets.push((cursor, width));
        cursor += width;
    }
    for (record_index, original) in document.records.iter().enumerate() {
        let mut record = original.clone();
        if let Some(rows) = by_record.get(&record_index) {
            let mut additions = Vec::new();
            for (protocol, (start, width)) in protocols.iter().zip(&offsets) {
                let alleles = rows
                    .iter()
                    .map(|row| row[*start..*start + *width].join("|"))
                    .collect::<Vec<_>>()
                    .join(",");
                if alleles != nastring {
                    additions.push(format!(
                        "FA_{}={}",
                        sanitize_info(&protocol.name),
                        escape_info(&alleles)
                    ));
                }
            }
            if !additions.is_empty() {
                if record[7] == "." {
                    record[7].clear();
                } else {
                    record[7].push(';');
                }
                record[7].push_str(&additions.join(";"));
            }
        }
        writeln!(writer, "{}", record.join("\t"))?;
    }
    Ok(())
}

fn append_annotation(
    row: &mut Vec<String>,
    annotation: Option<Annotation>,
    width: usize,
    nastring: &str,
) {
    if let Some(annotation) = annotation {
        let before = row.len();
        row.extend(annotation.values);
        while row.len() < before + width {
            row.push(nastring.to_string());
        }
    } else {
        row.extend(std::iter::repeat_n(nastring.to_string(), width));
    }
}

fn qualified_header(header: &str, protocol: &str) -> String {
    if header == "value" {
        protocol.to_string()
    } else {
        format!("{header}.{protocol}")
    }
}

pub(crate) fn write_record(
    writer: &mut dyn Write,
    fields: &[String],
    delimiter: char,
) -> Result<()> {
    for (i, value) in fields.iter().enumerate() {
        if i > 0 {
            writer.write_all(&[delimiter as u8])?;
        }
        if delimiter == ',' && value.contains([',', '"', '\n', '\r']) {
            writer.write_all(b"\"")?;
            for (j, part) in value.split('"').enumerate() {
                if j > 0 {
                    writer.write_all(b"\"\"")?;
                }
                writer.write_all(part.as_bytes())?;
            }
            writer.write_all(b"\"")?;
        } else {
            writer.write_all(value.as_bytes())?;
        }
    }
    writer.write_all(b"\n").context("failed writing output")
}

fn sanitize_info(value: &str) -> String {
    value
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect()
}
fn escape_info(value: &str) -> String {
    value.replace([' ', ';', '='], "_")
}
