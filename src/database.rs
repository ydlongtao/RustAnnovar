use crate::model::{Annotation, AnnotationKind, Variant, normalize_chrom};
use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;

#[derive(Debug, Clone)]
pub struct FilterRecord {
    pub variant: Variant,
    pub values: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct RegionRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub values: Vec<String>,
    pub score: Option<f64>,
    pub name: Option<String>,
}

#[derive(Debug)]
pub struct FilterDatabase {
    pub headers: Vec<String>,
    pub(crate) records: HashMap<OwnedKey, Vec<Vec<String>>>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub(crate) struct OwnedKey {
    pub(crate) chrom: String,
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) reference: String,
    pub(crate) alternate: String,
}

impl FilterDatabase {
    pub fn load(path: &Path) -> Result<Self> {
        let reader = open_text(path)?;
        let mut headers = Vec::new();
        let mut records: HashMap<OwnedKey, Vec<Vec<String>>> = HashMap::new();
        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            if line.starts_with('#') {
                if headers.is_empty() && !line.starts_with("##") {
                    headers = line
                        .trim_start_matches('#')
                        .split('\t')
                        .skip(5)
                        .map(str::to_string)
                        .collect();
                }
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                bail!(
                    "{}:{}: filter database needs Chr,Start,End,Ref,Alt",
                    path.display(),
                    line_no + 1
                );
            }
            let start: u64 = fields[1]
                .parse()
                .with_context(|| format!("invalid database start at line {}", line_no + 1))?;
            let end: u64 = fields[2]
                .parse()
                .with_context(|| format!("invalid database end at line {}", line_no + 1))?;
            let reference = db_allele(fields[3]);
            let alternate = db_allele(fields[4]);
            let internal_start = if reference.is_empty() {
                start
            } else {
                start
                    .checked_sub(1)
                    .context("database coordinates must be one-based")?
            };
            let internal_end = if reference.is_empty() { start } else { end };
            let key = OwnedKey {
                chrom: normalize_chrom(fields[0]),
                start: internal_start,
                end: internal_end,
                reference,
                alternate,
            };
            records
                .entry(key)
                .or_default()
                .push(fields[5..].iter().map(|v| (*v).to_string()).collect());
        }
        if headers.is_empty() {
            let width = records
                .values()
                .next()
                .and_then(|v| v.first())
                .map_or(1, Vec::len);
            headers = (1..=width)
                .map(|i| {
                    if width == 1 {
                        "value".into()
                    } else {
                        format!("value{i}")
                    }
                })
                .collect();
            if headers.is_empty() {
                headers.push("match".to_string());
            }
        }
        Ok(Self { headers, records })
    }

    /// Load only source-file bins needed by the query when a current
    /// `.fai.json` sidecar is available. Falls back to a full load otherwise.
    pub fn load_for_variants(path: &Path, variants: &[Variant]) -> Result<Self> {
        let index_path = path.with_extension("fai.json");
        let manifest = File::open(&index_path)
            .ok()
            .and_then(|file| serde_json::from_reader::<_, IndexManifest>(file).ok());
        let Some(manifest) = manifest else {
            return Self::load(path);
        };
        if manifest.bins.is_empty()
            || !manifest.is_current()?
            || path.extension().is_some_and(|ext| ext == "gz")
        {
            return Self::load(path);
        }

        let mut requested = HashSet::new();
        for variant in variants {
            requested.insert((variant.chrom.clone(), variant.start / manifest.bin_size));
            requested.insert((
                variant.chrom.clone(),
                variant.start.saturating_sub(1) / manifest.bin_size,
            ));
        }
        let mut ranges = requested
            .into_iter()
            .filter_map(|(chrom, bin)| manifest.bins.get(&chrom)?.get(&bin).copied())
            .collect::<Vec<_>>();
        ranges.sort_by_key(|range| range.start);
        ranges.dedup();

        let headers = read_filter_headers(path)?;
        let mut records: HashMap<OwnedKey, Vec<Vec<String>>> = HashMap::new();
        let mut file = File::open(path)?;
        for range in ranges {
            file.seek(SeekFrom::Start(range.start))?;
            let mut bytes = vec![0; usize::try_from(range.end - range.start)?];
            file.read_exact(&mut bytes)?;
            for line in String::from_utf8(bytes)?.lines() {
                if line.is_empty() || line.starts_with('#') {
                    continue;
                }
                let (key, values) = parse_filter_record(line)
                    .with_context(|| format!("invalid indexed record in {}", path.display()))?;
                records.entry(key).or_default().push(values);
            }
        }
        let headers = infer_filter_headers(headers, &records);
        Ok(Self { headers, records })
    }

    pub fn annotate(&self, variant: &Variant, protocol: &str) -> Option<Annotation> {
        let key = OwnedKey {
            chrom: variant.chrom.clone(),
            start: variant.start,
            end: variant.end,
            reference: variant.reference.clone(),
            alternate: variant.alternate.clone(),
        };
        self.records.get(&key).map(|matches| Annotation {
            protocol: protocol.to_string(),
            kind: AnnotationKind::Filter,
            values: transpose_join(matches),
        })
    }
}

#[derive(Debug)]
pub struct RegionDatabase {
    indexes: HashMap<String, crate::interval::IntervalIndex>,
    pub headers: Vec<String>,
    by_chrom: HashMap<String, Vec<RegionRecord>>,
    gff3: bool,
}

impl RegionDatabase {
    pub fn load(path: &Path) -> Result<Self> {
        let reader = open_text(path)?;
        let mut headers = Vec::new();
        let mut by_chrom: HashMap<String, Vec<RegionRecord>> = HashMap::new();
        let mut value_width = 0usize;
        let mut gff3 = false;
        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            if line.starts_with('#') {
                if headers.is_empty() {
                    let columns = line.trim_start_matches('#').split('\t').collect::<Vec<_>>();
                    let coordinate_columns = if columns.first().is_some_and(|value| *value == "bin")
                    {
                        4
                    } else {
                        3
                    };
                    headers = columns
                        .into_iter()
                        .skip(coordinate_columns)
                        .map(str::to_string)
                        .collect();
                }
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                bail!(
                    "{}:{}: region database needs at least three columns",
                    path.display(),
                    line_no + 1
                );
            }
            let (chrom, start, end, values, score, name) = if is_gff3(&fields) {
                gff3 = true;
                let start = fields[3]
                    .parse::<u64>()?
                    .checked_sub(1)
                    .context("GFF3 coordinates must be one-based")?;
                (
                    normalize_chrom(fields[0]),
                    start,
                    fields[4].parse()?,
                    Vec::new(),
                    fields[5].parse::<f64>().ok(),
                    Some(gff3_id(fields[8])),
                )
            } else {
                let (chrom_col, start_col, end_col, value_col) = detect_region_columns(&fields)?;
                (
                    normalize_chrom(fields[chrom_col]),
                    fields[start_col].parse()?,
                    fields[end_col].parse()?,
                    fields[value_col..]
                        .iter()
                        .map(|v| (*v).to_string())
                        .collect(),
                    None,
                    None,
                )
            };
            value_width = value_width.max(values.len());
            by_chrom
                .entry(chrom.clone())
                .or_default()
                .push(RegionRecord {
                    chrom,
                    start,
                    end,
                    values,
                    score,
                    name,
                });
        }
        for records in by_chrom.values_mut() {
            records.sort_by_key(|record| record.start);
        }
        if gff3 {
            headers = vec!["value".to_string()];
        } else if headers.is_empty() {
            headers = (1..=value_width.max(1))
                .map(|index| {
                    if value_width <= 1 {
                        "value".to_string()
                    } else {
                        format!("value{index}")
                    }
                })
                .collect();
        } else {
            while headers.len() < value_width {
                headers.push(format!("value{}", headers.len() + 1));
            }
        }
        let indexes = by_chrom
            .iter()
            .map(|(chrom, records)| {
                (
                    chrom.clone(),
                    crate::interval::IntervalIndex::new(records.iter().map(|r| (r.start, r.end))),
                )
            })
            .collect();
        Ok(Self {
            headers,
            by_chrom,
            gff3,
            indexes,
        })
    }

    pub fn annotate(
        &self,
        variant: &Variant,
        protocol: &str,
        min_overlap: f64,
    ) -> Option<Annotation> {
        let records = self.by_chrom.get(&variant.chrom)?;
        let query_len = (variant.end.saturating_sub(variant.start)).max(1) as f64;
        let mut hits: Vec<&RegionRecord> = Vec::new();
        for index in self
            .indexes
            .get(&variant.chrom)?
            .query(variant.start, variant.end)
        {
            let record = &records[index];
            if variant.overlaps(record.start, record.end) {
                let overlap = variant
                    .end
                    .min(record.end)
                    .saturating_sub(variant.start.max(record.start))
                    .max(1) as f64;
                if overlap / query_len >= min_overlap {
                    hits.push(record);
                }
            }
        }
        if hits.is_empty() {
            return None;
        }
        let values = if self.gff3 {
            let max_score = hits
                .iter()
                .filter_map(|record| record.score)
                .max_by(f64::total_cmp);
            let names = hits
                .iter()
                .filter(|record| max_score.is_none() || record.score == max_score)
                .filter_map(|record| record.name.as_deref())
                .collect::<Vec<_>>()
                .join(",");
            let mut value = String::new();
            if let Some(score) = max_score {
                value.push_str(&format!("Score={score};"));
            }
            if !names.is_empty() {
                value.push_str(&format!("Name={names}"));
            }
            vec![value]
        } else {
            transpose_join(
                &hits
                    .iter()
                    .map(|record| record.values.clone())
                    .collect::<Vec<_>>(),
            )
        };
        Some(Annotation {
            protocol: protocol.to_string(),
            kind: AnnotationKind::Region,
            values,
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct IndexManifest {
    pub format_version: u32,
    pub source: PathBuf,
    pub source_size: u64,
    pub modified_unix_seconds: u64,
    pub prefix_sha256: String,
    pub kind: String,
    pub bin_size: u64,
    pub bins: HashMap<String, BTreeMap<u64, ByteRange>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct ByteRange {
    pub start: u64,
    pub end: u64,
}

impl IndexManifest {
    pub fn build(source: &Path, kind: &str) -> Result<Self> {
        let (source_size, modified_unix_seconds, prefix_sha256) = source_fingerprint(source)?;
        let bin_size = 1_000_000;
        let bins = if matches!(kind, "filter" | "generic")
            && source.extension().is_none_or(|ext| ext != "gz")
        {
            build_filter_bins(source, bin_size)?
        } else {
            HashMap::new()
        };
        Ok(Self {
            format_version: 1,
            source: source.to_path_buf(),
            source_size,
            modified_unix_seconds,
            prefix_sha256,
            kind: kind.to_string(),
            bin_size,
            bins,
        })
    }

    pub fn write(&self, path: &Path) -> Result<()> {
        let temp = path.with_extension("tmp");
        let mut writer = File::create(&temp)?;
        serde_json::to_writer_pretty(&mut writer, self)?;
        writeln!(writer)?;
        fs::rename(temp, path)?;
        Ok(())
    }

    pub fn is_current(&self) -> Result<bool> {
        let (size, modified, sha256) = source_fingerprint(&self.source)?;
        Ok(size == self.source_size
            && modified == self.modified_unix_seconds
            && sha256 == self.prefix_sha256)
    }
}

fn source_fingerprint(source: &Path) -> Result<(u64, u64, String)> {
    let metadata = fs::metadata(source)?;
    let modified = metadata.modified()?.duration_since(UNIX_EPOCH)?.as_secs();
    let mut reader = File::open(source)?;
    let mut limited = std::io::Read::by_ref(&mut reader).take(1024 * 1024);
    let mut hasher = Sha256::new();
    std::io::copy(&mut limited, &mut hasher)?;
    Ok((metadata.len(), modified, format!("{:x}", hasher.finalize())))
}

fn build_filter_bins(
    source: &Path,
    bin_size: u64,
) -> Result<HashMap<String, BTreeMap<u64, ByteRange>>> {
    let mut reader = BufReader::new(File::open(source)?);
    let mut bins: HashMap<String, BTreeMap<u64, ByteRange>> = HashMap::new();
    let mut line = String::new();
    loop {
        let start = reader.stream_position()?;
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        let end = reader.stream_position()?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 2 {
            continue;
        }
        let Ok(position) = fields[1].parse::<u64>() else {
            continue;
        };
        let chrom = normalize_chrom(fields[0]);
        let bin = position.saturating_sub(1) / bin_size;
        bins.entry(chrom)
            .or_default()
            .entry(bin)
            .and_modify(|range| range.end = end)
            .or_insert(ByteRange { start, end });
    }
    Ok(bins)
}

fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    crate::io::open_reader(path)
}

pub(crate) fn parse_filter_record(line: &str) -> Result<(OwnedKey, Vec<String>)> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 5 {
        bail!("filter database needs Chr,Start,End,Ref,Alt");
    }
    let start: u64 = fields[1].parse()?;
    let end: u64 = fields[2].parse()?;
    let reference = db_allele(fields[3]);
    let alternate = db_allele(fields[4]);
    let internal_start = if reference.is_empty() {
        start
    } else {
        start
            .checked_sub(1)
            .context("database coordinates must be one-based")?
    };
    let internal_end = if reference.is_empty() { start } else { end };
    Ok((
        OwnedKey {
            chrom: normalize_chrom(fields[0]),
            start: internal_start,
            end: internal_end,
            reference,
            alternate,
        },
        fields[5..]
            .iter()
            .map(|value| (*value).to_string())
            .collect(),
    ))
}

pub(crate) fn read_filter_headers(path: &Path) -> Result<Vec<String>> {
    for line in open_text(path)?.lines() {
        let line = line?;
        if line.starts_with('#') {
            return Ok(line
                .trim_start_matches('#')
                .split('\t')
                .skip(5)
                .map(str::to_string)
                .collect());
        }
        if !line.trim().is_empty() {
            break;
        }
    }
    Ok(Vec::new())
}

fn infer_filter_headers(
    mut headers: Vec<String>,
    records: &HashMap<OwnedKey, Vec<Vec<String>>>,
) -> Vec<String> {
    if headers.is_empty() {
        let width = records
            .values()
            .next()
            .and_then(|values| values.first())
            .map_or(1, Vec::len)
            .max(1);
        headers = (1..=width)
            .map(|index| {
                if width == 1 {
                    "value".to_string()
                } else {
                    format!("value{index}")
                }
            })
            .collect();
    }
    headers
}

fn db_allele(value: &str) -> String {
    if value == "-" || value == "0" {
        String::new()
    } else {
        value.to_ascii_uppercase()
    }
}

fn transpose_join(rows: &[Vec<String>]) -> Vec<String> {
    let width = rows.iter().map(Vec::len).max().unwrap_or(0);
    (0..width)
        .map(|column| {
            rows.iter()
                .filter_map(|row| row.get(column))
                .cloned()
                .collect::<Vec<_>>()
                .join(";")
        })
        .collect()
}

fn detect_region_columns(fields: &[&str]) -> Result<(usize, usize, usize, usize)> {
    if fields.len() >= 4
        && fields[0].parse::<u64>().is_ok()
        && fields[2].parse::<u64>().is_ok()
        && fields[3].parse::<u64>().is_ok()
    {
        Ok((1, 2, 3, 4)) // UCSC bin, chrom, chromStart, chromEnd
    } else if fields[1].parse::<u64>().is_ok() && fields[2].parse::<u64>().is_ok() {
        Ok((0, 1, 2, 3))
    } else {
        bail!("unable to detect region database coordinate columns")
    }
}

fn is_gff3(fields: &[&str]) -> bool {
    fields.len() >= 9
        && fields[3].parse::<u64>().is_ok()
        && fields[4].parse::<u64>().is_ok()
        && matches!(fields[6], "+" | "-" | "." | "?")
}

fn gff3_id(attributes: &str) -> String {
    attributes
        .split(';')
        .find_map(|field| field.trim().strip_prefix("ID="))
        .unwrap_or("NA")
        .to_string()
}
