use crate::model::Variant;
use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
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
    let reader = open_reader(path)?;
    let mut variants = Vec::new();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line.with_context(|| format!("failed reading line {}", line_no + 1))?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let (fields, remainder) = split_avinput_fields(&line);
        if fields.len() < 5 {
            bail!(
                "{}:{}: expected at least five AVinput columns",
                path.display(),
                line_no + 1
            );
        }
        let av_start: u64 = fields[1]
            .parse()
            .with_context(|| format!("invalid start at line {}", line_no + 1))?;
        let av_end: u64 = fields[2]
            .parse()
            .with_context(|| format!("invalid end at line {}", line_no + 1))?;
        let reference = allele_from_avinput(fields[3]);
        let alternate = allele_from_avinput(fields[4]);
        let (start, end) = if reference.is_empty() {
            (av_start, av_start)
        } else {
            if av_start == 0 {
                bail!("AVinput coordinates are one-based at line {}", line_no + 1);
            }
            (av_start - 1, av_end)
        };
        let mut variant = Variant::new(fields[0], start, end, reference, alternate)?;
        variant.output_chrom = fields[0].to_string();
        variant.extra = remainder
            .filter(|value| !value.is_empty())
            .map(|value| value.split('\t').map(str::to_string).collect())
            .unwrap_or_default();
        variant.source_line = line_no + 1;
        variants.push(variant);
    }
    Ok(variants)
}

pub fn read_vcf(path: &Path) -> Result<VcfDocument> {
    let mut reader = open_reader(path)?;
    let mut text = String::new();
    reader.read_to_string(&mut text)?;
    let mut headers = Vec::new();
    let mut records = Vec::new();
    let mut variants = Vec::new();
    for (line_no, line) in text.lines().enumerate() {
        if line.starts_with('#') {
            headers.push(line.to_string());
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<String> = line.split('\t').map(str::to_string).collect();
        if fields.len() < 8 {
            bail!(
                "{}:{}: VCF requires at least eight columns",
                path.display(),
                line_no + 1
            );
        }
        let pos: u64 = fields[1]
            .parse()
            .with_context(|| format!("invalid VCF POS at line {}", line_no + 1))?;
        if pos == 0 {
            bail!("VCF POS is one-based at line {}", line_no + 1);
        }
        let reference = fields[3].to_ascii_uppercase();
        let alts: Vec<&str> = fields[4].split(',').collect();
        let record_index = records.len();
        for (allele_index, alt) in alts.iter().enumerate() {
            if alt.starts_with('<')
                || *alt == "*"
                || *alt == "."
                || alt.contains('[')
                || alt.contains(']')
            {
                continue;
            }
            let (start, end, normalized_ref, normalized_alt) =
                normalize_vcf_alleles(pos - 1, &reference, alt);
            let mut variant = Variant::new(&fields[0], start, end, normalized_ref, normalized_alt)?;
            variant.source_line = line_no + 1;
            variant.source_record = Some(record_index);
            variant.allele_index = allele_index;
            variants.push(variant);
        }
        records.push(fields);
    }
    Ok(VcfDocument {
        headers,
        records,
        variants,
    })
}

/// Trim identical suffixes and prefixes while preserving at least one side of
/// the event. Left alignment against a reference genome is intentionally a
/// separate operation.
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
