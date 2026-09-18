//! Bounded, ordered annotation with atomic file publication.
use crate::batch::VariantBatch;
use crate::io::{open_reader, parse_avinput_record, parse_vcf_record, supported_alt};
use crate::pipeline::{AnnotationEngine, Protocol, write_record};
use anyhow::{Context, Result, bail};
use serde::Serialize;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Debug, Clone, clap::Args)]
pub struct ExecutionOptions {
    #[arg(long,value_enum,default_value_t=Normalization::Annovar)]
    pub normalize: Normalization,
    #[arg(long)]
    pub reference: Option<PathBuf>,
    #[arg(long, default_value_t = 16)]
    pub threads: usize,
    #[arg(long, default_value_t = 10_000)]
    pub batch_size: usize,
    /// Managed memory budget in GiB (not an operating-system RSS limit).
    #[arg(long, default_value_t = 16)]
    pub memory_budget: usize,
    #[arg(long)]
    pub tmp_dir: Option<PathBuf>,
    #[arg(long)]
    pub report_json: Option<PathBuf>,
    #[arg(long, value_enum, default_value_t = Unsupported::Preserve)]
    pub unsupported: Unsupported,
    /// Fail atomically if any SNV lacks a CADD score (including ambiguous N).
    #[arg(long)]
    pub cadd_require_all: bool,
    /// Verify CADD's official assembly header (hg19 or hg38).
    #[arg(long)]
    pub cadd_build: Option<String>,
    /// Verify CADD version in the official score-file header.
    #[arg(long)]
    pub cadd_version: Option<String>,
}
#[derive(Debug, Clone, Copy, clap::ValueEnum, PartialEq, Eq)]
pub enum Normalization {
    Annovar,
    LeftAlign,
}
#[derive(Debug, Clone, Copy, clap::ValueEnum, PartialEq, Eq)]
pub enum Unsupported {
    Preserve,
    Error,
}

#[derive(Default, Serialize)]
struct Report {
    version: &'static str,
    records: u64,
    alleles: u64,
    unsupported: u64,
    batches: u64,
    threads: usize,
    database_seconds: f64,
    annotation_seconds: f64,
    write_seconds: f64,
    total_seconds: f64,
    cadd_counts: std::collections::BTreeMap<String, std::collections::BTreeMap<String, u64>>,
    cadd_databases: std::collections::BTreeMap<String, serde_json::Value>,
}

pub struct AtomicOutput {
    writer: Box<dyn Write>,
    temp: Option<tempfile::TempPath>,
    destination: PathBuf,
}
impl AtomicOutput {
    pub fn new(path: &Path) -> Result<Self> {
        if path == Path::new("-") {
            return Ok(Self {
                writer: Box::new(std::io::BufWriter::new(std::io::stdout())),
                temp: None,
                destination: path.into(),
            });
        }
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let file = tempfile::NamedTempFile::new_in(parent)?;
        let (file, temp) = file.into_parts();
        Ok(Self {
            writer: Box::new(std::io::BufWriter::new(file)),
            temp: Some(temp),
            destination: path.into(),
        })
    }
    pub fn commit(mut self) -> Result<()> {
        self.writer.flush()?;
        drop(self.writer);
        if let Some(temp) = self.temp {
            temp.persist(&self.destination).map_err(|e| e.error)?;
        }
        Ok(())
    }
}
impl Write for AtomicOutput {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.writer.write(b)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.writer.flush()
    }
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    input: &Path,
    vcf: bool,
    protocols: &[Protocol],
    output: &Path,
    vcf_output: Option<&Path>,
    csv: bool,
    nastring: &str,
    opts: &ExecutionOptions,
) -> Result<()> {
    if (opts.cadd_require_all || opts.cadd_build.is_some() || opts.cadd_version.is_some())
        && !protocols
            .iter()
            .any(|p| p.operation == crate::pipeline::Operation::Cadd)
    {
        bail!("--cadd-require-all requires --operation cadd");
    }
    if opts.threads == 0 || opts.batch_size == 0 || opts.memory_budget == 0 {
        bail!("threads, batch-size and memory-budget must be positive");
    }
    if vcf_output.is_some() && !vcf {
        bail!("--vcf-output requires --vcf-input");
    }
    if vcf_output == Some(output) {
        bail!("table and VCF outputs must differ");
    }
    if opts.report_json.as_deref() == Some(output)
        || (opts.report_json.is_some() && opts.report_json.as_deref() == vcf_output)
    {
        bail!("report path conflicts with output");
    }
    let start = Instant::now();
    let reference = if opts.normalize == Normalization::LeftAlign {
        Some(crate::reference::ReferenceGenome::open(
            opts.reference
                .as_deref()
                .context("--normalize left-align requires --reference")?,
        )?)
    } else {
        None
    };
    let engine = AnnotationEngine::new_normalized(protocols, reference.as_ref())?;
    if let Some(build) = &opts.cadd_build {
        engine.check_cadd_build(build)?;
    }
    if let Some(version) = &opts.cadd_version {
        engine.check_cadd_version(version)?;
    }
    let mut report = Report {
        cadd_databases: engine.cadd_metadata(),
        version: env!("CARGO_PKG_VERSION"),
        threads: opts.threads,
        database_seconds: start.elapsed().as_secs_f64(),
        ..Default::default()
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.threads)
        .build()?;
    let mut table = AtomicOutput::new(output)?;
    let mut vout = vcf_output.map(AtomicOutput::new).transpose()?;
    let mut headers = engine.headers().to_vec();
    headers.push("AnnotationStatus".into());
    // AVinput extra columns have unknown width. A prepass keeps a fixed schema;
    // stdin is spooled into a bounded-memory temporary file.
    let mut spool = None;
    let mut extra_width = 0;
    if !vcf {
        if input == Path::new("-") {
            let mut file = if let Some(dir) = &opts.tmp_dir {
                tempfile::NamedTempFile::new_in(dir)?
            } else {
                tempfile::NamedTempFile::new()?
            };
            std::io::copy(&mut std::io::stdin(), &mut file)?;
            spool = Some(file);
        }
        let path = spool.as_ref().map_or(input, |f| f.path());
        for (i, line) in open_reader(path)?.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            extra_width = extra_width.max(parse_avinput_record(&line, i + 1)?.extra.len());
        }
    }
    headers.extend((1..=extra_width).map(|i| format!("Otherinfo{i}")));
    write_record(&mut table, &headers, if csv { ',' } else { '\t' })?;
    let path = spool.as_ref().map_or(input, |f| f.path());
    let mut variants = VariantBatch::default();
    let mut records = Vec::new();
    let mut batch_bytes = 0usize;
    let max_bytes =
        (opts.memory_budget.saturating_mul(1024 * 1024 * 1024) / 16).min(64 * 1024 * 1024);
    let mut seen_chrom = false;
    let mut vcf_columns = 0;
    for (i, line) in open_reader(path)?.lines().enumerate() {
        let line = line.with_context(|| format!("reading input line {}", i + 1))?;
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('#') {
            if vcf {
                if report.records > 0 {
                    bail!("VCF header after data at line {}", i + 1);
                }
                if line.starts_with("##INFO=<ID=FA_") {
                    bail!("input already contains RustAnnovar INFO fields");
                }
                if line.starts_with("#CHROM\t") {
                    if seen_chrom {
                        bail!("duplicate VCF column header");
                    }
                    seen_chrom = true;
                    vcf_columns = line.split('\t').count();
                    if let Some(w) = &mut vout {
                        for (p, range) in protocols.iter().zip(engine.ranges()) {
                            writeln!(
                                w,
                                "##INFO=<ID=FA_{},Number=A,Type=String,Description=\"Fields: {}. Values percent-escaped; fields separated by |\">",
                                info_id(&p.name),
                                engine.headers()[range].join("|")
                            )?;
                        }
                        writeln!(
                            w,
                            "##INFO=<ID=FA_STATUS,Number=A,Type=String,Description=\"Annotation status per ALT\">"
                        )?;
                    }
                }
                if let Some(w) = &mut vout {
                    writeln!(w, "{line}")?;
                }
            }
            continue;
        }
        if vcf && !seen_chrom {
            bail!("VCF requires #CHROM header before data");
        }
        if line.len() > max_bytes {
            bail!("input record exceeds batch byte budget");
        }
        if !records.is_empty()
            && (records.len() >= opts.batch_size || batch_bytes + line.len() > max_bytes)
        {
            process(
                &engine,
                &pool,
                &mut variants,
                &mut records,
                &mut table,
                &mut vout,
                protocols,
                nastring,
                csv,
                extra_width,
                opts,
                &mut report,
            )?;
            batch_bytes = 0;
        }
        if vcf {
            let (fields, vs) = parse_vcf_record(&line, i + 1, records.len())?;
            if fields.len() != vcf_columns {
                bail!("VCF column count differs from header at line {}", i + 1);
            }
            if let Some(r) = &reference {
                let pos = fields[1].parse::<u64>()? - 1;
                if r.sequence(
                    &crate::model::normalize_chrom(&fields[0]),
                    pos,
                    pos + fields[3].len() as u64,
                )? != fields[3].to_ascii_uppercase().as_bytes()
                {
                    bail!("VCF reference mismatch at line {}", i + 1);
                }
            }
            for mut v in vs {
                if supported_alt(&v.alternate) || v.alternate.is_empty() {
                    if let Some(r) = &reference {
                        r.normalize(&mut v)?;
                    }
                }
                variants.push(v);
            }
            records.push(fields);
        } else {
            let mut v = parse_avinput_record(&line, i + 1)?;
            if let Some(r) = &reference {
                r.normalize(&mut v)?;
            }
            variants.push(v);
            records.push(Vec::new());
        }
        batch_bytes += line.len();
        report.records += 1;
    }
    process(
        &engine,
        &pool,
        &mut variants,
        &mut records,
        &mut table,
        &mut vout,
        protocols,
        nastring,
        csv,
        extra_width,
        opts,
        &mut report,
    )?;
    if vcf && !seen_chrom {
        bail!("missing VCF header");
    }
    table.flush()?;
    if let Some(w) = &mut vout {
        w.flush()?;
    }
    if let Some(w) = vout {
        w.commit()?;
    }
    table.commit()?;
    report.total_seconds = start.elapsed().as_secs_f64();
    if let Some(path) = &opts.report_json {
        let mut w = AtomicOutput::new(path)?;
        serde_json::to_writer_pretty(&mut w, &report)?;
        writeln!(w)?;
        w.commit()?;
    }
    eprintln!(
        "annotated {} records, {} alleles; unsupported {}",
        report.records, report.alleles, report.unsupported
    );
    Ok(())
}
#[allow(clippy::too_many_arguments)]
fn process(
    engine: &AnnotationEngine,
    pool: &rayon::ThreadPool,
    packed: &mut VariantBatch,
    records: &mut Vec<Vec<String>>,
    table: &mut AtomicOutput,
    vout: &mut Option<AtomicOutput>,
    protocols: &[Protocol],
    nastring: &str,
    csv: bool,
    extra_width: usize,
    opts: &ExecutionOptions,
    report: &mut Report,
) -> Result<()> {
    if records.is_empty() {
        return Ok(());
    }
    let variants = std::mem::take(packed).into_variants();
    let mut statuses: Vec<_> = variants
        .iter()
        .map(|v| {
            if v.source_record.is_some() && !supported_alt(&v.alternate) && !v.alternate.is_empty()
            {
                "unsupported_alt"
            } else {
                "ok"
            }
        })
        .collect();
    let unsupported = statuses.iter().filter(|&&s| s != "ok").count();
    if unsupported > 0 && opts.unsupported == Unsupported::Error {
        bail!(
            "unsupported ALT in batch at input line {}",
            variants[statuses.iter().position(|&s| s != "ok").unwrap()].source_line
        );
    }
    let t = Instant::now();
    let mut result = pool.install(|| engine.annotate(&variants, nastring))?;
    for (protocol, range) in protocols.iter().zip(engine.ranges()) {
        if protocol.operation == crate::pipeline::Operation::Cadd {
            for (variant, row) in variants.iter().zip(&result.rows) {
                let state = &row[range.start + 2];
                *report
                    .cadd_counts
                    .entry(protocol.name.clone())
                    .or_default()
                    .entry(state.clone())
                    .or_default() += 1;
                if opts.cadd_require_all && crate::cadd::is_snv_shape(variant) && state != "scored"
                {
                    bail!(
                        "CADD {}: unscored SNV {}:{} {}>{} ({state}), input line {}",
                        protocol.name,
                        variant.chrom,
                        variant.end,
                        variant.reference,
                        variant.alternate,
                        variant.source_line
                    );
                }
            }
        }
    }
    for (i, row) in result.rows.iter().enumerate() {
        if statuses[i] == "ok" {
            for (protocol, range) in protocols.iter().zip(engine.ranges()) {
                if protocol.operation == crate::pipeline::Operation::Gene
                    && row[range.start + 3].split(',').any(|s| s == "unknown")
                {
                    statuses[i] = "coding_unknown";
                }
            }
        }
    }
    report.annotation_seconds += t.elapsed().as_secs_f64();
    let t = Instant::now();
    for (i, row) in result.rows.iter_mut().enumerate() {
        if statuses[i] == "unsupported_alt" {
            for (protocol, range) in protocols.iter().zip(engine.ranges()) {
                if protocol.operation != crate::pipeline::Operation::Cadd {
                    for field in &mut row[range] {
                        *field = nastring.into();
                    }
                }
            }
        }
        row.insert(engine.headers().len(), statuses[i].into());
        row.resize(headers_len(engine, extra_width), nastring.into());
        write_record(table, row, if csv { ',' } else { '\t' })?;
    }
    if let Some(w) = vout {
        let mut cursor = 0;
        for (id, fields) in records.iter_mut().enumerate() {
            let from = cursor;
            while cursor < variants.len() && variants[cursor].source_record == Some(id) {
                cursor += 1;
            }
            let mut additions = Vec::new();
            for (p, range) in protocols.iter().zip(engine.ranges()) {
                let values = result.rows[from..cursor]
                    .iter()
                    .map(|r| {
                        r[range.clone()]
                            .iter()
                            .map(|s| escape(s))
                            .collect::<Vec<_>>()
                            .join("|")
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                additions.push(format!("FA_{}={values}", info_id(&p.name)));
            }
            additions.push(format!("FA_STATUS={}", statuses[from..cursor].join(",")));
            if fields[7] == "." {
                fields[7] = additions.join(";");
            } else {
                fields[7].push(';');
                fields[7].push_str(&additions.join(";"));
            }
            writeln!(w, "{}", fields.join("\t"))?;
        }
    }
    report.write_seconds += t.elapsed().as_secs_f64();
    report.alleles += variants.len() as u64;
    report.unsupported += unsupported as u64;
    report.batches += 1;
    records.clear();
    Ok(())
}
fn headers_len(engine: &AnnotationEngine, extra: usize) -> usize {
    engine.headers().len() + 1 + extra
}
fn info_id(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect()
}
fn escape(s: &str) -> String {
    s.bytes()
        .map(|b| {
            if matches!(b, b'%' | b';' | b'=' | b',' | b'|' | b'"' | b'\\') || b <= 32 || b >= 127 {
                format!("%{b:02X}")
            } else {
                (b as char).to_string()
            }
        })
        .collect()
}
