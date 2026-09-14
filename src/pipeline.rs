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
}

impl Operation {
    pub fn parse(value: &str) -> Result<Self> {
        match value {
            "g" | "gx" | "gene" => Ok(Self::Gene),
            "r" | "region" => Ok(Self::Region),
            "f" | "filter" => Ok(Self::Filter),
            _ => bail!("unknown operation {value:?}; use g, r, or f"),
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
            let database = db_dir.join(format!("{build}_{name}.txt"));
            if !database.exists() {
                bail!("database does not exist: {}", database.display());
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

pub fn annotate_table(
    variants: &[Variant],
    protocols: &[Protocol],
    nastring: &str,
) -> Result<TableResult> {
    enum Loaded {
        Filter(FilterDatabase),
        Region(RegionDatabase),
        Gene(GeneDatabase),
    }
    let mut loaded = Vec::new();
    let mut headers = vec!["Chr", "Start", "End", "Ref", "Alt"]
        .into_iter()
        .map(str::to_string)
        .collect::<Vec<_>>();
    for protocol in protocols {
        let database = match protocol.operation {
            Operation::Filter => {
                let db = FilterDatabase::load_for_variants(&protocol.database, variants)?;
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
    let extra_width = variants
        .iter()
        .map(|variant| variant.extra.len())
        .max()
        .unwrap_or(0);
    headers.extend((1..=extra_width).map(|index| format!("Otherinfo{index}")));
    let rows = variants
        .par_iter()
        .map(|variant| {
            let mut row = variant.avinput_fields().to_vec();
            for (protocol, database) in protocols.iter().zip(&loaded) {
                let (annotation, width) = match database {
                    Loaded::Filter(db) => (db.annotate(variant, &protocol.name), db.headers.len()),
                    Loaded::Region(db) => {
                        (db.annotate(variant, &protocol.name, 0.0), db.headers.len())
                    }
                    Loaded::Gene(db) => (Some(db.annotate(variant, &protocol.name, 2, 1000)), 5),
                };
                append_annotation(&mut row, annotation, width, nastring);
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

fn write_record(writer: &mut dyn Write, fields: &[String], delimiter: char) -> Result<()> {
    let encoded = fields
        .iter()
        .map(|value| {
            if delimiter == ','
                && (value.contains(',') || value.contains('"') || value.contains('\n'))
            {
                format!("\"{}\"", value.replace('"', "\"\""))
            } else {
                value.clone()
            }
        })
        .collect::<Vec<_>>();
    writeln!(writer, "{}", encoded.join(&delimiter.to_string())).context("failed writing output")
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
