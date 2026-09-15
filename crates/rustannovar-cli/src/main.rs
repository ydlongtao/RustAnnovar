use anyhow::{Context, Result, bail, ensure};
use clap::{Parser, Subcommand};
use rustannovar_core::{AnnotationResult, fields};
use rustannovar_db::Database;
use rustannovar_filter::Engine;
use rustannovar_vcf::{AvinputReader, Record};
use std::{
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

#[derive(Parser)]
#[command(
    version,
    about = "Experimental sorted-index ANNOVAR filter engine (AVinput MVP)"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}
#[derive(Subcommand)]
enum Command {
    Annotate {
        #[arg(long)]
        input: PathBuf,
        #[arg(long)]
        database: PathBuf,
        #[arg(long, default_value = "-")]
        output: PathBuf,
        #[arg(long, default_value_t = 1)]
        threads: usize,
        #[arg(long, default_value_t = 100_000)]
        batch_size: usize,
        /// Approximate retained input bytes per batch, excluding database and parser buffer.
        #[arg(long, default_value_t = 64 * 1024 * 1024)]
        batch_bytes: usize,
        #[arg(long, default_value = ".")]
        nastring: String,
        #[arg(long)]
        report_json: Option<PathBuf>,
        #[arg(long)]
        quiet: bool,
    },
}
struct Output {
    writer: BufWriter<Box<dyn Write>>,
    temp: Option<tempfile::TempPath>,
    destination: PathBuf,
}
impl Output {
    fn new(path: &Path) -> Result<Self> {
        let (writer, temp): (Box<dyn Write>, _) = if path == Path::new("-") {
            (Box::new(std::io::stdout()), None)
        } else {
            let parent = path
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or(Path::new("."));
            let (file, temp) = tempfile::NamedTempFile::new_in(parent)?.into_parts();
            (Box::new(file), Some(temp))
        };
        Ok(Self {
            writer: BufWriter::with_capacity(256 * 1024, writer),
            temp,
            destination: path.into(),
        })
    }
    fn commit(mut self) -> Result<()> {
        self.writer.flush()?;
        drop(self.writer);
        if let Some(temp) = self.temp {
            temp.persist(self.destination).map_err(|e| e.error)?;
        }
        Ok(())
    }
}
// All values are already TSV fields. Quotes are CSV-escaped for readers such
// as Python's csv module; normal fields take a single buffered write.
fn field(w: &mut impl Write, s: &str) -> Result<()> {
    if s.contains(['"', '\t', '\n', '\r']) {
        w.write_all(b"\"")?;
        let mut first = true;
        for p in s.split('"') {
            if !first {
                w.write_all(b"\"\"")?;
            }
            first = false;
            w.write_all(p.as_bytes())?;
        }
        w.write_all(b"\"")?;
    } else {
        w.write_all(s.as_bytes())?;
    }
    Ok(())
}
fn write_row(
    w: &mut impl Write,
    db: &Database,
    record: &Record,
    result: &AnnotationResult<'_>,
    missing: &str,
    scratch: &mut String,
) -> Result<()> {
    if let Some(end) = record.tsv_core_end {
        w.write_all(&record.raw.as_bytes()[..end])?;
    } else {
        for (i, value) in fields(&record.raw).take(5).enumerate() {
            if i != 0 {
                w.write_all(b"\t")?;
            }
            field(w, value)?;
        }
    }
    for col in 0..db.headers.len() {
        w.write_all(b"\t")?;
        match result {
            AnnotationResult::Missing => field(w, missing)?,
            AnnotationResult::Matches(ids) => {
                if ids.len() == 1 {
                    field(w, db.row(ids[0]).get(col).map_or("1", AsRef::as_ref))?;
                } else {
                    scratch.clear();
                    for (i, &id) in ids.iter().enumerate() {
                        if i != 0 {
                            scratch.push(';');
                        }
                        scratch.push_str(db.row(id).get(col).map_or("1", AsRef::as_ref));
                    }
                    field(w, scratch)?;
                }
            }
        }
    }
    if let Some(end) = record.tsv_core_end {
        w.write_all(&record.raw.as_bytes()[end..])?;
    } else {
        for value in fields(&record.raw).skip(5) {
            w.write_all(b"\t")?;
            field(w, value)?;
        }
    }
    w.write_all(b"\n")?;
    Ok(())
}
fn resolved(path: &Path) -> Result<PathBuf> {
    if path.exists() {
        return Ok(path.canonicalize()?);
    }
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    Ok(parent
        .canonicalize()?
        .join(path.file_name().context("missing filename")?))
}
fn run() -> Result<()> {
    let Command::Annotate {
        input,
        database,
        output,
        threads,
        batch_size,
        batch_bytes,
        nastring,
        report_json,
        quiet,
    } = Cli::parse().command;
    ensure!(
        batch_size > 0 && batch_bytes > 0,
        "batch-size and batch-bytes must be positive"
    );
    ensure!(
        report_json.as_deref() != Some(Path::new("-")),
        "report-json requires a file path"
    );
    // Reject aliases as well as identical strings before opening any output.
    for destination in std::iter::once(&output)
        .chain(report_json.iter())
        .filter(|p| p.as_path() != Path::new("-"))
    {
        for source in [&input, &database]
            .into_iter()
            .filter(|p| p.as_path() != Path::new("-"))
        {
            ensure!(
                resolved(destination)? != resolved(source)?,
                "output/report must not replace an input database or variant file"
            );
        }
    }
    if let Some(report) = &report_json {
        if output != Path::new("-") {
            ensure!(
                resolved(&output)? != resolved(report)?,
                "report and output paths must differ"
            );
        }
    }
    let total = Instant::now();
    let engine = Engine::new(threads)?;
    let load = Instant::now();
    let db = Database::load(BufReader::with_capacity(
        256 * 1024,
        File::open(&database).with_context(|| format!("opening {}", database.display()))?,
    ))
    .with_context(|| format!("loading {}", database.display()))?;
    let load_seconds = load.elapsed().as_secs_f64();
    let source: Box<dyn std::io::BufRead> = if input == Path::new("-") {
        Box::new(BufReader::new(std::io::stdin()))
    } else {
        Box::new(BufReader::with_capacity(
            256 * 1024,
            File::open(&input).with_context(|| format!("opening {}", input.display()))?,
        ))
    };
    let mut reader = AvinputReader::new(source);
    let mut out = Output::new(&output)?;
    let mut records = Vec::with_capacity(batch_size.min(100_000));
    let mut schema = None;
    let mut wrote_header = false;
    let mut count = 0_u64;
    let mut hits = 0_u64;
    let mut batches = 0_u64;
    let mut parse_seconds = 0.0;
    let mut lookup_seconds = 0.0;
    let mut write_seconds = 0.0;
    let mut scratch = String::new();
    loop {
        let parse = Instant::now();
        records.clear();
        let mut bytes = 0_usize;
        while records.len() < batch_size && bytes < batch_bytes {
            let Some(record) = reader
                .next_record()
                .with_context(|| format!("parsing {}", input.display()))?
            else {
                break;
            };
            let extra = fields(&record.raw).skip(5).count();
            if *schema.get_or_insert(extra) != extra {
                bail!(
                    "input record {}: inconsistent extra-column count",
                    record.variant.row_id + 1
                );
            }
            let retained = record
                .raw
                .len()
                .saturating_mul(2)
                .saturating_add(std::mem::size_of::<Record>());
            ensure!(
                retained <= batch_bytes,
                "input record {} exceeds batch-bytes; increase the budget",
                record.variant.row_id + 1
            );
            bytes = bytes.saturating_add(retained);
            records.push(record);
        }
        parse_seconds += parse.elapsed().as_secs_f64();
        if !wrote_header {
            let mut headers: Vec<_> = ["Chr", "Start", "End", "Ref", "Alt"]
                .into_iter()
                .map(String::from)
                .collect();
            headers.extend(db.headers.iter().cloned());
            headers.extend((1..=schema.unwrap_or(0)).map(|i| format!("Otherinfo{i}")));
            let unique: std::collections::HashSet<_> = headers.iter().collect();
            ensure!(
                unique.len() == headers.len(),
                "database and input extra-column names collide"
            );
            for (i, h) in headers.iter().enumerate() {
                if i != 0 {
                    out.writer.write_all(b"\t")?;
                }
                field(&mut out.writer, h)?;
            }
            out.writer.write_all(b"\n")?;
            wrote_header = true;
        }
        if records.is_empty() {
            break;
        }
        let lookup = Instant::now();
        let annotations = engine.annotate(&db, &records, |r| &r.variant);
        lookup_seconds += lookup.elapsed().as_secs_f64();
        let write = Instant::now();
        for (record, annotation) in records.iter().zip(&annotations) {
            hits += u64::from(matches!(annotation, AnnotationResult::Matches(_)));
            write_row(
                &mut out.writer,
                &db,
                record,
                annotation,
                &nastring,
                &mut scratch,
            )?;
        }
        write_seconds += write.elapsed().as_secs_f64();
        count += records.len() as u64;
        batches += 1;
    }
    let flush = Instant::now();
    out.commit()?;
    write_seconds += flush.elapsed().as_secs_f64();
    let elapsed = total.elapsed().as_secs_f64();
    if let Some(path) = report_json {
        let mut report = Output::new(&path)?;
        serde_json::to_writer_pretty(
            &mut report.writer,
            &serde_json::json!({
                "engine":"sorted-filter-mvp", "version":env!("CARGO_PKG_VERSION"), "records":count, "matched_records":hits,
                "database_records":db.records(), "database_chromosomes":db.chromosomes(), "threads":threads,
                "batch_size":batch_size,"batch_bytes":batch_bytes,"batches":batches,
                "database_load_index_seconds":load_seconds,"parse_normalize_seconds":parse_seconds,
                "index_lookup_seconds":lookup_seconds,"merge_write_seconds":write_seconds,
                "total_seconds":elapsed,"variants_per_second":count as f64 / elapsed.max(f64::EPSILON)
            }),
        )?;
        report.writer.write_all(b"\n")?;
        report.commit()?;
    }
    if !quiet {
        eprintln!(
            "{count} records, {hits} matched; {elapsed:.3}s; {:.0} variants/s",
            count as f64 / elapsed.max(f64::EPSILON)
        );
    }
    Ok(())
}
fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error:#}");
        std::process::exit(1);
    }
}
