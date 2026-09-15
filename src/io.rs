use crate::model::Variant;
use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub fn open_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    if path == Path::new("-") {
        return Ok(Box::new(BufReader::new(io::stdin())));
    }
    let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    if path.extension().is_some_and(|ext| ext == "gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

pub fn open_writer(path: &Path) -> Result<Box<dyn Write>> {
    if path == Path::new("-") {
        Ok(Box::new(BufWriter::new(io::stdout())))
    } else {
        Ok(Box::new(BufWriter::new(File::create(path).with_context(
            || format!("cannot create {}", path.display()),
        )?)))
    }
}

pub fn read_avinput(path: &Path) -> Result<Vec<Variant>> {
    open_reader(path)?
        .lines()
        .enumerate()
        .filter_map(|(i, line)| match line {
            Err(e) => Some(Err(e.into())),
            Ok(line) if line.trim().is_empty() || line.starts_with('#') => None,
            Ok(line) => Some(parse_avinput_record(&line, i + 1)),
        })
        .collect()
}

pub fn parse_avinput_record(line: &str, line_no: usize) -> Result<Variant> {
    let (fields, remainder) = split_avinput_fields(line);
    if fields.len() < 5 {
        bail!("line {line_no}: expected at least five AVinput columns");
    }
    let start: u64 = fields[1].parse().context("invalid AVinput start")?;
    let end: u64 = fields[2].parse().context("invalid AVinput end")?;
    let reference = allele_from_avinput(fields[3]);
    let alternate = allele_from_avinput(fields[4]);
    if end < start {
        bail!("line {line_no}: end precedes start");
    }
    let (start, end) = if reference.is_empty() {
        (start, start)
    } else {
        (
            start
                .checked_sub(1)
                .context("AVinput coordinates are one-based")?,
            end,
        )
    };
    let mut v = Variant::new(fields[0], start, end, reference, alternate)?;
    v.output_chrom = fields[0].into();
    v.source_line = line_no;
    v.extra = remainder
        .filter(|v| !v.is_empty())
        .map(|v| v.split('\t').map(str::to_string).collect())
        .unwrap_or_default();
    Ok(v)
}

pub fn parse_vcf_record(
    line: &str,
    line_no: usize,
    record: usize,
) -> Result<(Vec<String>, Vec<Variant>)> {
    let fields: Vec<String> = line.split('\t').map(str::to_string).collect();
    if fields.len() < 8 {
        bail!("line {line_no}: VCF requires at least eight columns");
    }
    let pos: u64 = fields[1].parse().context("invalid VCF POS")?;
    if pos == 0 {
        bail!("VCF POS is one-based");
    }
    let reference = fields[3].to_ascii_uppercase();
    if reference.is_empty() || !reference.bytes().all(|b| b"ACGTN".contains(&b)) {
        bail!("line {line_no}: invalid VCF REF");
    }
    let mut variants = Vec::new();
    for (allele_index, alt) in fields[4].split(',').enumerate() {
        if alt.is_empty() {
            bail!("line {line_no}: empty ALT");
        }
        let (start, end, r, a) = if supported_alt(alt) {
            normalize_vcf_alleles(pos - 1, &reference, &alt.to_ascii_uppercase())
        } else {
            (
                pos - 1,
                (pos - 1)
                    .checked_add(reference.len() as u64)
                    .context("coordinate overflow")?,
                reference.clone(),
                alt.into(),
            )
        };
        let mut v = Variant::new(&fields[0], start, end, r, a)?;
        if !supported_alt(alt) {
            v.alternate = alt.to_string();
        }
        v.source_line = line_no;
        v.source_record = Some(record);
        v.allele_index = allele_index;
        variants.push(v);
    }
    Ok((fields, variants))
}

pub fn supported_alt(alt: &str) -> bool {
    !alt.is_empty() && alt.bytes().all(|b| b"ACGTNacgtn".contains(&b))
}

pub fn read_vcf(path: &Path) -> Result<VcfDocument> {
    let mut document = VcfDocument {
        headers: Vec::new(),
        records: Vec::new(),
        variants: Vec::new(),
    };
    for (i, line) in open_reader(path)?.lines().enumerate() {
        let line = line?;
        if line.starts_with('#') {
            document.headers.push(line);
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let (fields, variants) = parse_vcf_record(&line, i + 1, document.records.len())?;
        document.records.push(fields);
        document.variants.extend(variants);
    }
    Ok(document)
}

/// Remove shared sequence without reference-based left alignment.
pub fn normalize_vcf_alleles(
    mut start: u64,
    reference: &str,
    alternate: &str,
) -> (u64, u64, String, String) {
    let mut r = reference.as_bytes().to_vec();
    let mut a = alternate.as_bytes().to_vec();
    while r.len() > 1 && a.len() > 1 && r.last() == a.last() {
        r.pop();
        a.pop();
    }
    while r.len() > 1 && a.len() > 1 && r.first() == a.first() {
        r.remove(0);
        a.remove(0);
        start += 1;
    }
    if r.len() == 1 && a.len() > 1 && r[0] == a[0] {
        start += 1;
        a.remove(0);
        r.clear();
    } else if a.len() == 1 && r.len() > 1 && a[0] == r[0] {
        start += 1;
        r.remove(0);
        a.clear();
    }
    let end = start + r.len() as u64;
    (
        start,
        end,
        String::from_utf8(r).unwrap(),
        String::from_utf8(a).unwrap(),
    )
}

pub fn write_avinput(variants: &[Variant], path: &Path, include_info: bool) -> Result<()> {
    let mut writer = open_writer(path)?;
    for variant in variants {
        let mut fields = variant.avinput_fields().to_vec();
        if include_info {
            fields.extend(variant.extra.clone());
        }
        writeln!(writer, "{}", fields.join("\t"))?;
    }
    Ok(())
}

#[derive(Debug, Clone)]
pub struct VcfDocument {
    pub headers: Vec<String>,
    pub records: Vec<Vec<String>>,
    pub variants: Vec<Variant>,
}

fn allele_from_avinput(value: &str) -> String {
    if value == "-" {
        String::new()
    } else {
        value.to_ascii_uppercase()
    }
}

fn split_avinput_fields(line: &str) -> (Vec<&str>, Option<&str>) {
    let bytes = line.as_bytes();
    let mut fields = Vec::with_capacity(5);
    let mut cursor = 0usize;
    while fields.len() < 5 {
        while cursor < bytes.len() && bytes[cursor].is_ascii_whitespace() {
            cursor += 1;
        }
        if cursor == bytes.len() {
            break;
        }
        let start = cursor;
        while cursor < bytes.len() && !bytes[cursor].is_ascii_whitespace() {
            cursor += 1;
        }
        fields.push(&line[start..cursor]);
    }
    while cursor < bytes.len() && bytes[cursor].is_ascii_whitespace() {
        cursor += 1;
    }
    (fields, (cursor < line.len()).then(|| &line[cursor..]))
}
