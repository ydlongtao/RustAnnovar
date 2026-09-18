//! Streaming import of official CADD precomputed SNV tables.
use anyhow::{Context, Result, bail, ensure};
use std::{
    io::{BufRead, Write},
    path::Path,
};

pub const HEADERS: [&str; 3] = ["CADD_raw", "CADD_phred", "CADD_status"];

pub struct Database {
    path: std::path::PathBuf,
    index: noodles_tabix::Index,
    chroms: std::collections::HashMap<String, String>,
    pub build: Option<String>,
    pub version: Option<String>,
    remote: Option<crate::http_range::Source>,
}
trait ReadSeek: std::io::Read + std::io::Seek + Send {}
impl<T: std::io::Read + std::io::Seek + Send> ReadSeek for T {}
impl Database {
    pub fn open(path: &Path) -> Result<Self> {
        use noodles_csi::BinningIndex;
        let url = path.to_str().filter(|s| s.starts_with("https://"));
        let (index, source, remote): (_, Box<dyn ReadSeek>, _) = if let Some(url) = url {
            let response = ureq::AgentBuilder::new()
                .https_only(true)
                .timeout(std::time::Duration::from_secs(60))
                .build()
                .get(&format!("{url}.tbi"))
                .call()?;
            let mut bytes = Vec::new();
            std::io::Read::read_to_end(
                &mut std::io::Read::take(response.into_reader(), 16 * 1024 * 1024 + 1),
                &mut bytes,
            )?;
            ensure!(
                bytes.len() <= 16 * 1024 * 1024,
                "oversized CADD tabix index"
            );
            let index = noodles_tabix::io::Reader::new(std::io::Cursor::new(bytes)).read_index()?;
            let reader = crate::http_range::Reader::open(url)?;
            let remote = Some(reader.source.clone());
            (index, Box::new(reader), remote)
        } else {
            (
                noodles_tabix::fs::read(format!("{}.tbi", path.display()))?,
                Box::new(std::fs::File::open(path)?),
                None,
            )
        };
        let header = index.header().context("CADD tabix index has no header")?;
        ensure!(
            header.reference_sequence_name_index() == 0
                && header.start_position_index() == 1
                && matches!(header.end_position_index(), None | Some(1)),
            "CADD tabix must index chromosome column 1 and one-based position column 2"
        );
        ensure!(
            header.format()
                == noodles_csi::binning_index::index::header::Format::Generic(
                    noodles_csi::binning_index::index::header::format::CoordinateSystem::Gff
                ),
            "CADD tabix requires one-based generic coordinates"
        );
        let mut file = noodles_bgzf::Reader::new(source);
        let mut line = String::new();
        let mut found = false;
        let mut build = None;
        let mut version = None;
        for _ in 0..100 {
            line.clear();
            if std::io::Read::take(&mut file, 65537).read_line(&mut line)? == 0 {
                break;
            }
            ensure!(line.len() <= 65536, "oversized CADD header");
            if line.starts_with("##CADD GRCh37-") {
                build = Some("hg19".to_string());
            }
            if line.starts_with("##CADD GRCh38-") {
                build = Some("hg38".to_string());
            }
            if let Some(value) = line
                .strip_prefix("##CADD GRCh37-v")
                .or_else(|| line.strip_prefix("##CADD GRCh38-v"))
            {
                version = value.split_whitespace().next().map(str::to_string);
            }
            if line.starts_with("#Chrom\t") {
                ensure!(
                    line.trim_end() == "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED",
                    "CADD native lookup requires official six-column score-only header"
                );
                found = true;
                break;
            }
            ensure!(line.starts_with('#'), "CADD score file has no header");
        }
        ensure!(found, "CADD score header not found");
        let mut chroms = std::collections::HashMap::new();
        for name in header.reference_sequence_names() {
            let name = std::str::from_utf8(name.as_ref())?.to_string();
            let normalized = crate::model::normalize_chrom(&name);
            ensure!(
                chroms.insert(normalized, name).is_none(),
                "ambiguous CADD chromosome aliases"
            );
        }
        Ok(Self {
            path: path.into(),
            index,
            chroms,
            build,
            version,
            remote,
        })
    }
    pub fn check_build(&self, expected: &str) -> Result<()> {
        ensure!(
            matches!(expected, "hg19" | "hg38"),
            "CADD build must be hg19 or hg38"
        );
        ensure!(
            self.build.as_deref() == Some(expected),
            "CADD assembly header missing or mismatches {expected}"
        );
        Ok(())
    }
    pub fn check_version(&self, expected: &str) -> Result<()> {
        ensure!(
            self.version.as_deref() == Some(expected),
            "CADD version header missing or mismatches {expected}"
        );
        Ok(())
    }
    pub fn metadata(&self) -> serde_json::Value {
        serde_json::json!({"source": self.path, "build": self.build, "version": self.version, "remote": self.remote.is_some()})
    }
    pub fn batch(&self, variants: &[crate::Variant], missing: &str) -> Result<Vec<Vec<String>>> {
        use std::collections::{BTreeSet, HashMap};
        let source: Box<dyn ReadSeek> = if let Some(remote) = &self.remote {
            Box::new(crate::http_range::Reader::from_source(remote.clone()))
        } else {
            Box::new(std::fs::File::open(&self.path)?)
        };
        let mut reader = noodles_csi::io::IndexedReader::new(source, self.index.clone());
        let mut positions = BTreeSet::new();
        for v in variants {
            if is_snv(v) && self.chroms.contains_key(&v.chrom) {
                positions.insert((v.chrom.clone(), v.end));
            }
        }
        type ScoreRow = (String, String, String, String);
        let mut records: HashMap<(String, u64), Vec<ScoreRow>> = HashMap::new();
        // Coalesce nearby requested positions, but bound decoding work per query.
        let mut windows: Vec<(String, u64, u64)> = Vec::new();
        for (chrom, pos) in &positions {
            if let Some((last_chrom, start, end)) = windows.last_mut() {
                if last_chrom == chrom && pos - *end <= 32 && pos - *start < 4096 {
                    *end = *pos;
                    continue;
                }
            }
            windows.push((chrom.clone(), *pos, *pos));
        }
        for (chrom, start, end) in windows {
            let start = noodles_core::Position::try_from(usize::try_from(start)?)?;
            let end = noodles_core::Position::try_from(usize::try_from(end)?)?;
            let region = noodles_core::Region::new(self.chroms[&chrom].clone(), start..=end);
            for result in reader.query(&region)? {
                let record = result?;
                let fields = record.as_ref().split('\t').collect::<Vec<_>>();
                ensure!(
                    fields.len() == 6,
                    "CADD native lookup requires score-only six-column table"
                );
                let actual: u64 = fields[1].parse()?;
                ensure!(
                    crate::model::normalize_chrom(fields[0]) == chrom,
                    "CADD index/data chromosome mismatch"
                );
                if !positions.contains(&(chrom.clone(), actual)) {
                    continue;
                }
                for (i, text) in fields[4..].iter().enumerate() {
                    let score: f64 = text.parse()?;
                    ensure!(
                        score.is_finite() && (i == 0 || score >= 0.0),
                        "invalid CADD score"
                    );
                }
                let rows = records.entry((chrom.clone(), actual)).or_default();
                ensure!(
                    !rows.iter().any(|r| r.0 == fields[2] && r.1 == fields[3]),
                    "duplicate CADD allele score"
                );
                rows.push((
                    fields[2].into(),
                    fields[3].into(),
                    fields[4].into(),
                    fields[5].into(),
                ));
            }
        }
        Ok(variants
            .iter()
            .map(|v| {
                let mut status = if is_snv_shape(v) && !is_snv(v) {
                    "unsupported_snv"
                } else if !is_snv(v) {
                    "not_snv"
                } else if !self.chroms.contains_key(&v.chrom) {
                    "contig_not_covered"
                } else {
                    "score_not_found"
                };
                if is_snv(v) {
                    if let Some(rows) = records.get(&(v.chrom.clone(), v.end)) {
                        if !rows.iter().any(|r| r.0 == v.reference) {
                            status = "reference_mismatch";
                        }
                        if let Some(r) = rows
                            .iter()
                            .find(|r| r.0 == v.reference && r.1 == v.alternate)
                        {
                            return vec![r.2.clone(), r.3.clone(), "scored".into()];
                        }
                    }
                }
                vec![missing.into(), missing.into(), status.into()]
            })
            .collect())
    }
}
pub fn is_snv(v: &crate::Variant) -> bool {
    let base = |s: &str| s.len() == 1 && matches!(s.as_bytes()[0], b'A' | b'C' | b'G' | b'T');
    v.start.checked_add(1) == Some(v.end)
        && base(&v.reference)
        && base(&v.alternate)
        && v.reference != v.alternate
}
pub fn is_snv_shape(v: &crate::Variant) -> bool {
    let base =
        |s: &str| s.len() == 1 && matches!(s.as_bytes()[0], b'A' | b'C' | b'G' | b'T' | b'N');
    v.start.checked_add(1) == Some(v.end)
        && base(&v.reference)
        && base(&v.alternate)
        && v.reference != v.alternate
}

pub fn import(input: &Path, output: &Path, build: &str, version: &str) -> Result<u64> {
    ensure!(
        matches!(build, "hg19" | "hg38"),
        "CADD build must be hg19 or hg38"
    );
    ensure!(
        !version.is_empty() && version.chars().all(|c| c.is_ascii_digit() || c == '.'),
        "invalid CADD version"
    );
    ensure!(
        !output.exists(),
        "refusing to overwrite {}",
        output.display()
    );
    ensure!(
        output != Path::new("-"),
        "CADD import requires a file destination"
    );
    let mut writer = crate::stream::AtomicOutput::new(output)?;
    writeln!(
        writer,
        "##CADD_build={build};version={version};source={}",
        input.display()
    )?;
    writeln!(writer, "#Chr\tStart\tEnd\tRef\tAlt\tCADD_raw\tCADD_phred")?;
    let mut columns = None;
    let mut count = 0;
    for (number, line) in crate::io::open_reader(input)?.lines().enumerate() {
        let line = line.with_context(|| format!("CADD line {}", number + 1))?;
        if line.is_empty() || line.starts_with("##") {
            continue;
        }
        if line.starts_with('#') {
            let fields = line.trim_start_matches('#').split('\t').collect::<Vec<_>>();
            let find = |name: &str| {
                fields
                    .iter()
                    .position(|v| *v == name)
                    .with_context(|| format!("CADD header missing {name}"))
            };
            columns = Some([
                find("Chrom")?,
                find("Pos")?,
                find("Ref")?,
                find("Alt")?,
                find("RawScore")?,
                find("PHRED")?,
            ]);
            continue;
        }
        let indices =
            columns.context("CADD table requires #Chrom/Pos/Ref/Alt/RawScore/PHRED header")?;
        let fields = line.split('\t').collect::<Vec<_>>();
        let mut values = Vec::with_capacity(6);
        for index in indices {
            values.push(
                *fields
                    .get(index)
                    .with_context(|| format!("truncated CADD line {}", number + 1))?,
            );
        }
        let pos: u64 = values[1]
            .parse()
            .with_context(|| format!("invalid CADD position at line {}", number + 1))?;
        ensure!(
            pos > 0,
            "CADD position must be one-based at line {}",
            number + 1
        );
        let base = |v: &str| v.len() == 1 && matches!(v.as_bytes()[0], b'A' | b'C' | b'G' | b'T');
        ensure!(
            base(values[2]) && base(values[3]) && values[2] != values[3],
            "not a canonical SNV at CADD line {}",
            number + 1
        );
        for (i, v) in values[4..].iter().enumerate() {
            let score: f64 = v
                .parse()
                .with_context(|| format!("invalid CADD score at line {}", number + 1))?;
            ensure!(
                score.is_finite() && (i == 0 || score >= 0.0),
                "invalid CADD score at line {}",
                number + 1
            );
        }
        writeln!(
            writer,
            "{}\t{pos}\t{pos}\t{}\t{}\t{}\t{}",
            values[0], values[2], values[3], values[4], values[5]
        )?;
        count += 1;
    }
    if count == 0 {
        bail!("CADD table contains no SNVs");
    }
    writer.commit()?;
    Ok(count)
}
