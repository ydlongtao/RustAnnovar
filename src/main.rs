use anyhow::{Context, Result, bail};
use clap::{Args, Parser, Subcommand, ValueEnum};
use rust_annovar::database::IndexManifest;
use rust_annovar::io::read_avinput;
use rust_annovar::pipeline::{Operation, resolve_protocols};
use sha2::{Digest, Sha256};
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(
    version,
    about = "Fast genomic variant annotation with ANNOVAR-compatible databases"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Convert VCF to ANNOVAR's five-column input format.
    Convert(ConvertArgs),
    /// Run one annotation database.
    Annotate(AnnotateArgs),
    /// Combine multiple gene, region, and filter annotations.
    Table(TableArgs),
    /// Manage and validate local database metadata.
    Db(DbArgs),
    /// Extract reference sequence for AVinput intervals.
    Sequence(SequenceArgs),
    /// Report coding consequences (a gene annotation convenience command).
    CodingChange(AnnotateArgs),
    /// Select rows from a multi-annotation table using a simple predicate.
    Reduce(ReduceArgs),
}

#[derive(Args)]
struct ConvertArgs {
    input: PathBuf,
    #[arg(short, long, default_value = "-")]
    output: PathBuf,
    #[arg(long)]
    include_info: bool,
}

#[derive(Clone, ValueEnum)]
enum OpArg {
    Cadd,
    Gene,
    Region,
    Filter,
}
impl From<OpArg> for Operation {
    fn from(value: OpArg) -> Self {
        match value {
            OpArg::Cadd => Self::Cadd,
            OpArg::Gene => Self::Gene,
            OpArg::Region => Self::Region,
            OpArg::Filter => Self::Filter,
        }
    }
}

#[derive(Args)]
struct AnnotateArgs {
    #[command(flatten)]
    execution: rust_annovar::stream::ExecutionOptions,
    input: PathBuf,
    database: PathBuf,
    #[arg(short, long)]
    output: PathBuf,
    #[arg(long, value_enum, default_value = "filter")]
    operation: OpArg,
    #[arg(long, default_value = "annotation")]
    protocol: String,
    #[arg(long)]
    fasta: Option<PathBuf>,
    #[arg(long)]
    vcf_input: bool,
    #[arg(long, default_value = ".")]
    nastring: String,
}

#[derive(Args)]
struct TableArgs {
    #[command(flatten)]
    execution: rust_annovar::stream::ExecutionOptions,
    input: PathBuf,
    db_dir: PathBuf,
    #[arg(long)]
    build: String,
    #[arg(long, value_delimiter = ',')]
    protocol: Vec<String>,
    #[arg(long, value_delimiter = ',')]
    operation: Vec<String>,
    #[arg(short, long)]
    output: PathBuf,
    #[arg(long)]
    vcf_input: bool,
    #[arg(long)]
    csv: bool,
    #[arg(long)]
    vcf_output: Option<PathBuf>,
    #[arg(long, default_value = ".")]
    nastring: String,
}

#[derive(Args)]
struct DbArgs {
    #[command(subcommand)]
    command: DbCommand,
}
#[derive(Subcommand)]
enum DbCommand {
    /// Import official precomputed CADD SNVs into an ANNOVAR-compatible table.
    ImportCadd {
        input: PathBuf,
        output: PathBuf,
        #[arg(long)]
        build: String,
        #[arg(long)]
        version: String,
    },
    Validate {
        database: PathBuf,
        #[arg(long)]
        full: bool,
    },
    Index {
        database: PathBuf,
        #[arg(long)]
        tmp_dir: Option<PathBuf>,
        #[arg(long)]
        reference: Option<PathBuf>,
        #[arg(long,value_enum,default_value_t=rust_annovar::stream::Normalization::Annovar)]
        normalize: rust_annovar::stream::Normalization,
        #[arg(long, default_value = "generic")]
        kind: String,
        #[arg(short, long)]
        output: Option<PathBuf>,
    },
    Check {
        index: PathBuf,
    },
    List {
        db_dir: PathBuf,
        #[arg(long)]
        build: Option<String>,
    },
    /// Download a public database atomically, with optional SHA-256 verification.
    Download {
        url: String,
        output: PathBuf,
        #[arg(long)]
        sha256: Option<String>,
    },
}

#[derive(Args)]
struct SequenceArgs {
    input: PathBuf,
    fasta: PathBuf,
    #[arg(short, long)]
    output: PathBuf,
}

#[derive(Args)]
struct ReduceArgs {
    input: PathBuf,
    #[arg(short, long)]
    output: PathBuf,
    #[arg(long)]
    column: String,
    #[arg(long)]
    equals: String,
}

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error:#}");
        std::process::exit(2);
    }
}

fn run() -> Result<()> {
    match Cli::parse().command {
        Command::Convert(args) => {
            use std::io::BufRead;
            let mut output = rust_annovar::stream::AtomicOutput::new(&args.output)?;
            for (i, line) in rust_annovar::io::open_reader(&args.input)?
                .lines()
                .enumerate()
            {
                let line = line?;
                if line.starts_with('#') || line.trim().is_empty() {
                    continue;
                }
                let (_, variants) = rust_annovar::io::parse_vcf_record(&line, i + 1, i)?;
                for v in variants {
                    let mut fields = v.avinput_fields().to_vec();
                    if args.include_info {
                        fields.extend(v.extra);
                    }
                    writeln!(output, "{}", fields.join("\t"))?;
                }
            }
            output.commit()?;
        }
        Command::Annotate(args) => run_annotate(args, None)?,
        Command::CodingChange(mut args) => {
            args.operation = OpArg::Gene;
            run_annotate(args, None)?;
        }
        Command::Table(args) => run_table(args)?,
        Command::Db(args) => run_db(args.command)?,
        Command::Sequence(args) => run_sequence(args)?,
        Command::Reduce(args) => run_reduce(args)?,
    }
    Ok(())
}

fn run_annotate(args: AnnotateArgs, _reserved: Option<()>) -> Result<()> {
    let operation = Operation::from(args.operation);
    let protocol = rust_annovar::pipeline::Protocol {
        name: args.protocol,
        operation,
        database: args.database,
        fasta: args.fasta,
    };
    rust_annovar::stream::run(
        &args.input,
        args.vcf_input,
        &[protocol],
        &args.output,
        None,
        false,
        &args.nastring,
        &args.execution,
    )
}

fn run_table(args: TableArgs) -> Result<()> {
    let operations = args
        .operation
        .iter()
        .map(|value| Operation::parse(value))
        .collect::<Result<Vec<_>>>()?;
    let protocols = resolve_protocols(&args.db_dir, &args.build, &args.protocol, &operations)?;
    rust_annovar::stream::run(
        &args.input,
        args.vcf_input,
        &protocols,
        &args.output,
        args.vcf_output.as_deref(),
        args.csv,
        &args.nastring,
        &args.execution,
    )
}

fn run_db(command: DbCommand) -> Result<()> {
    match command {
        DbCommand::ImportCadd {
            input,
            output,
            build,
            version,
        } => {
            let count = rust_annovar::cadd::import(&input, &output, &build, &version)?;
            println!("imported {count} CADD SNVs");
        }
        DbCommand::Validate { database, full } => {
            rust_annovar::disk_index::DiskFilter::open(&database, full)?;
            println!("current");
        }
        DbCommand::Index {
            database,
            tmp_dir,
            reference,
            normalize,
            kind,
            output,
        } => {
            if matches!(kind.as_str(), "filter" | "generic") && output.is_none() {
                println!(
                    "{}",
                    rust_annovar::disk_index::DiskFilter::build_normalized(
                        &database,
                        tmp_dir.as_deref(),
                        if normalize == rust_annovar::stream::Normalization::LeftAlign {
                            Some(rust_annovar::reference::ReferenceGenome::open(
                                reference
                                    .as_deref()
                                    .context("left-align requires --reference")?,
                            )?)
                        } else {
                            None
                        }
                        .as_ref()
                    )?
                    .display()
                );
                return Ok(());
            }
            let output = output.unwrap_or_else(|| database.with_extension("fai.json"));
            IndexManifest::build(&database, &kind)?.write(&output)?;
            println!("{}", output.display());
        }
        DbCommand::Check { index } => {
            let manifest: IndexManifest = serde_json::from_reader(fs::File::open(&index)?)?;
            if !manifest.is_current()? {
                bail!("stale index: {}", index.display());
            }
            println!("current");
        }
        DbCommand::List { db_dir, build } => {
            let prefix = build.map(|v| format!("{v}_"));
            let mut entries = fs::read_dir(db_dir)?
                .filter_map(Result::ok)
                .map(|e| e.path())
                .filter(|path| {
                    path.extension().is_some_and(|ext| ext == "txt")
                        && prefix.as_ref().is_none_or(|p| {
                            path.file_name()
                                .unwrap_or_default()
                                .to_string_lossy()
                                .starts_with(p)
                        })
                })
                .collect::<Vec<_>>();
            entries.sort();
            for path in entries {
                println!("{}", path.display());
            }
        }
        DbCommand::Download {
            url,
            output,
            sha256,
        } => download_database(&url, &output, sha256.as_deref())?,
    }
    Ok(())
}

fn download_database(url: &str, output: &Path, expected_sha256: Option<&str>) -> Result<()> {
    let parent = output.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(parent)?;
    let temp = output.with_extension("download.tmp");
    let result = (|| -> Result<()> {
        let response = ureq::get(url)
            .call()
            .with_context(|| format!("failed to download {url}"))?;
        let mut reader = response.into_reader();
        let mut writer = fs::File::create(&temp)?;
        let mut hasher = Sha256::new();
        let mut buffer = [0u8; 64 * 1024];
        loop {
            let size = reader.read(&mut buffer)?;
            if size == 0 {
                break;
            }
            writer.write_all(&buffer[..size])?;
            hasher.update(&buffer[..size]);
        }
        writer.sync_all()?;
        let actual = format!("{:x}", hasher.finalize());
        if let Some(expected) = expected_sha256
            && !actual.eq_ignore_ascii_case(expected)
        {
            bail!("SHA-256 mismatch: expected {expected}, got {actual}");
        }
        fs::rename(&temp, output)?;
        println!("{actual}  {}", output.display());
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temp);
    }
    result
}

fn run_sequence(args: SequenceArgs) -> Result<()> {
    let genome = read_fasta_genome(&args.fasta)?;
    let variants = read_avinput(&args.input)?;
    let mut output = String::new();
    for (index, variant) in variants.iter().enumerate() {
        let sequence = genome
            .get(&variant.chrom)
            .context("chromosome missing from FASTA")?;
        let end = usize::try_from(variant.end)?;
        let start = usize::try_from(variant.start)?;
        if end > sequence.len() {
            bail!("variant {} exceeds chromosome length", index + 1);
        }
        output.push_str(&format!(
            ">variant_{} {}:{}-{}\n{}\n",
            index + 1,
            variant.chrom,
            variant.start + 1,
            variant.end,
            String::from_utf8_lossy(&sequence[start..end])
        ));
    }
    fs::write(args.output, output)?;
    Ok(())
}

fn read_fasta_genome(path: &Path) -> Result<std::collections::HashMap<String, Vec<u8>>> {
    use std::io::BufRead;
    let mut map = std::collections::HashMap::new();
    let mut id = None::<String>;
    let mut seq = Vec::new();
    for line in rust_annovar::io::open_reader(path)?.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = id.replace(rust_annovar::model::normalize_chrom(
                header
                    .split_whitespace()
                    .next()
                    .context("empty FASTA header")?,
            )) {
                map.insert(previous, std::mem::take(&mut seq));
            }
        } else {
            seq.extend(line.trim().as_bytes().iter().map(u8::to_ascii_uppercase));
        }
    }
    if let Some(id) = id {
        map.insert(id, seq);
    }
    Ok(map)
}

fn run_reduce(args: ReduceArgs) -> Result<()> {
    let content = fs::read_to_string(&args.input)?;
    let mut lines = content.lines();
    let header = lines.next().context("empty table")?;
    let headers: Vec<&str> = header.split('\t').collect();
    let column = headers
        .iter()
        .position(|name| *name == args.column)
        .with_context(|| format!("unknown column {}", args.column))?;
    let mut output = String::from(header);
    output.push('\n');
    for line in lines {
        if line.split('\t').nth(column) == Some(args.equals.as_str()) {
            output.push_str(line);
            output.push('\n');
        }
    }
    fs::write(args.output, output)?;
    Ok(())
}
