//! Immutable sorted filter records, external merge sorting and bounded block cache.
use crate::database::{FilterDatabase, OwnedKey, parse_filter_record, read_filter_headers};
use crate::io::open_reader;
use crate::model::Variant;
use anyhow::{Context, Result, bail};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::time::UNIX_EPOCH;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
struct Row {
    key: OwnedKey,
    source_row: u64,
    values: Vec<String>,
}
#[derive(Debug, Serialize, Deserialize)]
struct Block {
    checksum: String,
    first: OwnedKey,
    last: OwnedKey,
    start: u64,
    len: u64,
}
#[derive(Debug, Serialize, Deserialize)]
pub struct Manifest {
    format_version: u32,
    source: PathBuf,
    source_size: u64,
    modified_ns: u128,
    prefix_sha256: String,
    sha256: String,
    data_sha256: String,
    data_file: PathBuf,
    normalization: String,
    reference_sha256: Option<String>,
    headers: Vec<String>,
    blocks: Vec<Block>,
}
#[derive(Debug, Default)]
struct Cache {
    bytes: usize,
    tick: u64,
    blocks: HashMap<usize, (u64, Vec<Row>, usize)>,
}
#[derive(Debug)]
pub struct DiskFilter {
    manifest: Manifest,
    path: PathBuf,
    cache: Mutex<Cache>,
}
fn paths(source: &Path) -> (PathBuf, PathBuf) {
    (
        source.with_extension("rai"),
        source.with_extension("rai.json"),
    )
}
fn hash(path: &Path, limit: u64) -> Result<String> {
    let mut h = Sha256::new();
    std::io::copy(&mut File::open(path)?.take(limit), &mut h)?;
    Ok(format!("{:x}", h.finalize()))
}
fn stamp(path: &Path) -> Result<(u64, u128)> {
    let m = fs::metadata(path)?;
    Ok((
        m.len(),
        m.modified()?.duration_since(UNIX_EPOCH)?.as_nanos(),
    ))
}
fn next_row(r: &mut impl BufRead) -> Result<Option<Row>> {
    let mut line = String::new();
    if r.read_line(&mut line)? == 0 {
        Ok(None)
    } else {
        Ok(Some(serde_json::from_str(&line)?))
    }
}
fn write_row(w: &mut impl Write, row: &Row) -> Result<()> {
    serde_json::to_writer(&mut *w, row)?;
    w.write_all(b"\n")?;
    Ok(())
}
fn merge(paths: &[PathBuf], destination: &Path) -> Result<()> {
    let mut readers = paths
        .iter()
        .map(|p| File::open(p).map(BufReader::new))
        .collect::<std::io::Result<Vec<_>>>()?;
    let mut heap = BinaryHeap::new();
    for (i, r) in readers.iter_mut().enumerate() {
        if let Some(row) = next_row(r)? {
            heap.push(Reverse((row, i)));
        }
    }
    let mut out = BufWriter::new(File::create(destination)?);
    while let Some(Reverse((row, i))) = heap.pop() {
        write_row(&mut out, &row)?;
        if let Some(row) = next_row(&mut readers[i])? {
            heap.push(Reverse((row, i)));
        }
    }
    out.flush()?;
    Ok(())
}
impl DiskFilter {
    pub fn exists(source: &Path) -> bool {
        paths(source).1.exists()
    }
    pub fn check_normalization(
        &self,
        reference: Option<&crate::reference::ReferenceGenome>,
    ) -> Result<()> {
        if self.manifest.reference_sha256.as_deref() != reference.map(|r| r.sha256.as_str())
            || (self.manifest.normalization == "left-align") != reference.is_some()
        {
            bail!("index/query normalization or reference mismatch");
        }
        Ok(())
    }
    pub fn headers(&self) -> &[String] {
        &self.manifest.headers
    }
    pub fn build(source: &Path, tmp: Option<&Path>) -> Result<PathBuf> {
        Self::build_normalized(source, tmp, None)
    }
    pub fn build_normalized(
        source: &Path,
        tmp: Option<&Path>,
        reference: Option<&crate::reference::ReferenceGenome>,
    ) -> Result<PathBuf> {
        let (data, manifest_path) = paths(source);
        let (size, modified) = stamp(source)?;
        let tmp = if let Some(p) = tmp {
            tempfile::tempdir_in(p)?
        } else {
            tempfile::tempdir()?
        };
        let mut runs = Vec::new();
        let mut rows = Vec::new();
        let mut bytes = 0;
        let mut headers = read_filter_headers(source)?;
        let mut width = None;
        for (i, line) in open_reader(source)?.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let (mut key, values) =
                parse_filter_record(&line).with_context(|| format!("database line {}", i + 1))?;
            if let Some(reference) = reference {
                let mut v = Variant::new(
                    &key.chrom,
                    key.start,
                    key.end,
                    &key.reference,
                    &key.alternate,
                )?;
                reference.normalize(&mut v)?;
                key = OwnedKey {
                    chrom: v.chrom,
                    start: v.start,
                    end: v.end,
                    reference: v.reference,
                    alternate: v.alternate,
                };
            }
            if key.end < key.start {
                bail!("database end before start at line {}", i + 1);
            }
            if let Some(w) = width {
                if values.len() != w {
                    bail!("inconsistent database column count at line {}", i + 1);
                }
            } else {
                width = Some(values.len());
            }
            bytes += line.len() * 3;
            rows.push(Row {
                key,
                source_row: i as u64,
                values,
            });
            if bytes >= 16 * 1024 * 1024 {
                flush_run(&mut rows, &mut runs, tmp.path())?;
                bytes = 0;
            }
        }
        if !rows.is_empty() {
            flush_run(&mut rows, &mut runs, tmp.path())?;
        }
        if headers.is_empty() {
            headers = (0..width.unwrap_or(1).max(1))
                .map(|i| {
                    if width.unwrap_or(1) <= 1 {
                        "value".into()
                    } else {
                        format!("value{}", i + 1)
                    }
                })
                .collect();
        }
        if width.is_some_and(|w| w > 0 && headers.len() != w) {
            bail!("database header width differs from records");
        }
        let mut generation = 0;
        while runs.len() > 1 {
            let mut merged = Vec::new();
            for (i, chunk) in runs.chunks(32).enumerate() {
                let path = tmp.path().join(format!("merge-{generation}-{i}"));
                merge(chunk, &path)?;
                merged.push(path);
            }
            for p in runs {
                fs::remove_file(p)?;
            }
            runs = merged;
            generation += 1;
        }
        let parent = data.parent().unwrap_or(Path::new("."));
        let mut file = tempfile::NamedTempFile::new_in(parent)?;
        let mut blocks = Vec::new();
        let mut start = 0;
        let mut buffer = Vec::new();
        let mut first = None;
        let mut last = None;
        if let Some(run) = runs.first() {
            let mut reader = BufReader::new(File::open(run)?);
            while let Some(row) = next_row(&mut reader)? {
                if first.is_none() {
                    first = Some(row.key.clone());
                }
                last = Some(row.key.clone());
                write_row(&mut buffer, &row)?;
                if buffer.len() >= 256 * 1024 {
                    let len = buffer.len() as u64;
                    file.write_all(&buffer)?;
                    blocks.push(Block {
                        checksum: format!("{:x}", Sha256::digest(&buffer)),
                        first: first.take().unwrap(),
                        last: last.take().unwrap(),
                        start,
                        len,
                    });
                    start += len;
                    buffer.clear();
                }
            }
        }
        if !buffer.is_empty() {
            file.write_all(&buffer)?;
            blocks.push(Block {
                checksum: format!("{:x}", Sha256::digest(&buffer)),
                first: first.unwrap(),
                last: last.unwrap(),
                start,
                len: buffer.len() as u64,
            });
        }
        file.flush()?;
        file.as_file().sync_all()?;
        if stamp(source)? != (size, modified) {
            bail!("source changed during indexing");
        }
        let data_sha256 = hash(file.path(), u64::MAX)?;
        let data_file = PathBuf::from(format!(
            "{}.{}.rai",
            source
                .file_name()
                .context("missing source name")?
                .to_string_lossy(),
            data_sha256
        ));
        let manifest = Manifest {
            format_version: 2,
            source: fs::canonicalize(source)?,
            source_size: size,
            modified_ns: modified,
            prefix_sha256: hash(source, 1024 * 1024)?,
            sha256: hash(source, u64::MAX)?,
            data_sha256,
            data_file: data_file.clone(),
            normalization: if reference.is_some() {
                "left-align"
            } else {
                "annovar"
            }
            .into(),
            reference_sha256: reference.map(|r| r.sha256.clone()),
            headers,
            blocks,
        };
        let mut mf = tempfile::NamedTempFile::new_in(parent)?;
        serde_json::to_writer(&mut mf, &manifest)?;
        mf.flush()?;
        // Publish data first; readers seeing an old manifest reject changed size/hash in full validation.
        file.persist(parent.join(data_file)).map_err(|e| e.error)?;
        mf.persist(&manifest_path).map_err(|e| e.error)?;
        Ok(manifest_path)
    }
    pub fn open(source: &Path, full: bool) -> Result<Self> {
        let (_, mp) = paths(source);
        let m: Manifest = serde_json::from_reader(File::open(mp)?)?;
        if m.data_file.components().count() != 1 || m.data_file.is_absolute() {
            bail!("invalid index data filename");
        }
        let path = source.parent().unwrap_or(Path::new(".")).join(&m.data_file);
        if m.format_version != 2 || !matches!(m.normalization.as_str(), "annovar" | "left-align") {
            bail!("unsupported index version/normalization; rebuild index");
        }
        if fs::canonicalize(source)? != m.source
            || stamp(source)? != (m.source_size, m.modified_ns)
            || hash(source, 1024 * 1024)? != m.prefix_sha256
        {
            bail!(
                "stale index for {}; rebuild with db index",
                source.display()
            );
        }
        let len = fs::metadata(&path)?.len();
        let mut end = 0;
        let mut previous = None;
        for b in &m.blocks {
            if previous.is_some_and(|last| last > &b.first) {
                bail!("unsorted index block directory");
            }
            previous = Some(&b.last);
            if b.start != end || b.len == 0 || b.first > b.last {
                bail!("corrupt index block directory");
            }
            end = end.checked_add(b.len).context("index overflow")?;
        }
        if end != len {
            bail!("index data length mismatch");
        }
        if full && (hash(source, u64::MAX)? != m.sha256 || hash(&path, u64::MAX)? != m.data_sha256)
        {
            bail!("index SHA-256 mismatch");
        }
        Ok(Self {
            manifest: m,
            path,
            cache: Mutex::new(Cache::default()),
        })
    }
    pub fn load_batch(&self, variants: &[Variant]) -> Result<FilterDatabase> {
        let keys: HashSet<_> = variants
            .iter()
            .map(|v| OwnedKey {
                chrom: v.chrom.clone(),
                start: v.start,
                end: v.end,
                reference: v.reference.clone(),
                alternate: v.alternate.clone(),
            })
            .collect();
        let mut selected = HashSet::new();
        for key in &keys {
            let mut i = self.manifest.blocks.partition_point(|b| b.last < *key);
            while i < self.manifest.blocks.len() && self.manifest.blocks[i].first <= *key {
                selected.insert(i);
                i += 1;
            }
        }
        let mut selected: Vec<_> = selected.into_iter().collect();
        selected.sort_unstable();
        let mut records: HashMap<OwnedKey, Vec<Vec<String>>> = HashMap::new();
        let mut file = File::open(&self.path)?;
        let mut cache = self
            .cache
            .lock()
            .map_err(|_| anyhow::anyhow!("index cache poisoned"))?;
        for i in selected {
            cache.tick += 1;
            let tick = cache.tick;
            if !cache.blocks.contains_key(&i) {
                let b = &self.manifest.blocks[i];
                let size = usize::try_from(b.len)?;
                if size > 64 * 1024 * 1024 {
                    bail!("index block exceeds 64 MiB limit");
                }
                while cache.bytes + size > 64 * 1024 * 1024 && !cache.blocks.is_empty() {
                    let id = *cache.blocks.iter().min_by_key(|(_, v)| v.0).unwrap().0;
                    let (_, _, bytes) = cache.blocks.remove(&id).unwrap();
                    cache.bytes -= bytes;
                }
                file.seek(SeekFrom::Start(b.start))?;
                let mut bytes = vec![0; size];
                file.read_exact(&mut bytes)?;
                if format!("{:x}", Sha256::digest(&bytes)) != b.checksum {
                    bail!("index block checksum mismatch");
                }
                let mut reader = std::io::Cursor::new(bytes);
                let mut decoded = Vec::new();
                let mut retained = 0;
                while let Some(row) = next_row(&mut reader)? {
                    retained += std::mem::size_of::<Row>()
                        + row.key.chrom.capacity()
                        + row.key.reference.capacity()
                        + row.key.alternate.capacity()
                        + row.values.capacity() * std::mem::size_of::<String>()
                        + row.values.iter().map(String::capacity).sum::<usize>();
                    decoded.push(row);
                }
                retained += (decoded.capacity() - decoded.len()) * std::mem::size_of::<Row>();
                if retained > 64 * 1024 * 1024 {
                    bail!("decoded index block exceeds cache budget");
                }
                while cache.bytes + retained > 64 * 1024 * 1024 && !cache.blocks.is_empty() {
                    let id = *cache.blocks.iter().min_by_key(|(_, v)| v.0).unwrap().0;
                    let (_, _, n) = cache.blocks.remove(&id).unwrap();
                    cache.bytes -= n;
                }
                cache.bytes += retained;
                cache.blocks.insert(i, (tick, decoded, retained));
            }
            let (used, rows, _) = cache.blocks.get_mut(&i).unwrap();
            *used = tick;
            for row in rows {
                if keys.contains(&row.key) {
                    records
                        .entry(row.key.clone())
                        .or_default()
                        .push(row.values.clone());
                }
            }
        }
        Ok(FilterDatabase {
            headers: self.manifest.headers.clone(),
            records,
        })
    }
}
fn flush_run(rows: &mut Vec<Row>, runs: &mut Vec<PathBuf>, dir: &Path) -> Result<()> {
    rows.sort();
    let path = dir.join(format!("run-{}", runs.len()));
    let mut out = BufWriter::new(File::create(&path)?);
    for row in rows.drain(..) {
        write_row(&mut out, &row)?;
    }
    out.flush()?;
    runs.push(path);
    Ok(())
}
