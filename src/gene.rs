use crate::model::{Annotation, AnnotationKind, Variant, normalize_chrom};
use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Transcript {
    pub id: String,
    pub gene: String,
    pub chrom: String,
    pub strand: char,
    pub tx_start: u64,
    pub tx_end: u64,
    pub cds_start: u64,
    pub cds_end: u64,
    pub exons: Vec<(u64, u64)>,
}

#[derive(Debug)]
pub struct GeneDatabase {
    by_chrom: HashMap<String, Vec<Transcript>>,
    sequences: HashMap<String, Vec<u8>>,
}

impl GeneDatabase {
    pub fn load(model_path: &Path, fasta_path: Option<&Path>) -> Result<Self> {
        let mut by_chrom: HashMap<String, Vec<Transcript>> = HashMap::new();
        for (line_no, line) in crate::io::open_reader(model_path)?.lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            let offset = usize::from(fields.first().is_some_and(|v| v.parse::<u64>().is_ok()));
            if fields.len() < offset + 10 {
                bail!(
                    "{}:{}: unsupported gene model",
                    model_path.display(),
                    line_no + 1
                );
            }
            let exon_count: usize = fields[offset + 7].parse()?;
            let starts = parse_positions(fields[offset + 8])?;
            let ends = parse_positions(fields[offset + 9])?;
            if starts.len() != exon_count || ends.len() != exon_count {
                bail!("exon count mismatch at line {}", line_no + 1);
            }
            let transcript = Transcript {
                id: fields[offset].to_string(),
                chrom: normalize_chrom(fields[offset + 1]),
                strand: fields[offset + 2]
                    .chars()
                    .next()
                    .context("missing strand")?,
                tx_start: fields[offset + 3].parse()?,
                tx_end: fields[offset + 4].parse()?,
                cds_start: fields[offset + 5].parse()?,
                cds_end: fields[offset + 6].parse()?,
                exons: starts.into_iter().zip(ends).collect(),
                gene: fields
                    .get(offset + 11)
                    .unwrap_or(&fields[offset])
                    .to_string(),
            };
            by_chrom
                .entry(transcript.chrom.clone())
                .or_default()
                .push(transcript);
        }
        for transcripts in by_chrom.values_mut() {
            transcripts.sort_by_key(|tx| tx.tx_start);
        }
        let sequences = fasta_path.map(read_fasta).transpose()?.unwrap_or_default();
        Ok(Self {
            by_chrom,
            sequences,
        })
    }

    pub fn annotate(
        &self,
        variant: &Variant,
        protocol: &str,
        splice: u64,
        flank: u64,
    ) -> Annotation {
        let mut hits = Vec::new();
        if let Some(transcripts) = self.by_chrom.get(&variant.chrom) {
            for transcript in transcripts {
                if transcript.tx_start > variant.end.saturating_add(flank) {
                    break;
                }
                if transcript.tx_end.saturating_add(flank) < variant.start {
                    continue;
                }
                if let Some(hit) = classify(
                    variant,
                    transcript,
                    self.sequences.get(&transcript.id),
                    splice,
                    flank,
                ) {
                    hits.push(hit);
                }
            }
        }
        if hits.is_empty() {
            hits.push(self.nearest_intergenic(variant));
        }
        hits.sort_by_key(|hit| hit.rank);
        let best_rank = hits[0].rank;
        let best: Vec<_> = hits.iter().filter(|hit| hit.rank == best_rank).collect();
        Annotation {
            protocol: protocol.to_string(),
            kind: AnnotationKind::Gene,
            values: vec![
                join_unique(best.iter().map(|hit| hit.function.as_str())),
                join_unique(best.iter().map(|hit| hit.gene.as_str())),
                join_unique(best.iter().map(|hit| hit.detail.as_str())),
                join_unique(best.iter().map(|hit| hit.exonic_function.as_str())),
                join_unique(best.iter().map(|hit| hit.aa_change.as_str())),
            ],
        }
    }

    fn nearest_intergenic(&self, variant: &Variant) -> GeneHit {
        let Some(transcripts) = self.by_chrom.get(&variant.chrom) else {
            return GeneHit::intergenic();
        };
        let mut left: Option<(&str, u64)> = None;
        let mut right: Option<(&str, u64)> = None;
        for transcript in transcripts {
            if transcript.tx_end <= variant.start {
                let distance =
                    variant.start - transcript.tx_end + u64::from(!variant.reference.is_empty());
                if left.is_none_or(|(_, best)| distance < best) {
                    left = Some((&transcript.gene, distance));
                }
            } else if transcript.tx_start >= variant.end {
                let distance = transcript.tx_start - variant.end + 1;
                if right.is_none_or(|(_, best)| distance < best) {
                    right = Some((&transcript.gene, distance));
                }
            }
        }
        let nearest = [left, right].into_iter().flatten().collect::<Vec<_>>();
        if nearest.is_empty() {
            return GeneHit::intergenic();
        }
        GeneHit {
            rank: 8,
            function: "intergenic".into(),
            gene: nearest
                .iter()
                .map(|(gene, _)| *gene)
                .collect::<Vec<_>>()
                .join(","),
            detail: nearest
                .iter()
                .map(|(_, distance)| format!("dist={distance}"))
                .collect::<Vec<_>>()
                .join(";"),
            exonic_function: ".".into(),
            aa_change: ".".into(),
        }
    }
}

#[derive(Debug)]
struct GeneHit {
    rank: u8,
    function: String,
    gene: String,
    detail: String,
    exonic_function: String,
    aa_change: String,
}

impl GeneHit {
    fn intergenic() -> Self {
        Self {
            rank: 8,
            function: "intergenic".into(),
            gene: ".".into(),
            detail: ".".into(),
            exonic_function: ".".into(),
            aa_change: ".".into(),
        }
    }
}

fn classify(
    variant: &Variant,
    tx: &Transcript,
    sequence: Option<&Vec<u8>>,
    splice: u64,
    flank: u64,
) -> Option<GeneHit> {
    let pos = variant.start;
    if variant.end <= tx.tx_start {
        let function = if tx.strand == '+' {
            "upstream"
        } else {
            "downstream"
        };
        let distance = tx.tx_start - variant.end + 1;
        return (distance <= flank).then(|| flank_hit(function, tx, distance));
    }
    if pos >= tx.tx_end {
        let function = if tx.strand == '+' {
            "downstream"
        } else {
            "upstream"
        };
        let distance = pos - tx.tx_end + u64::from(!variant.reference.is_empty());
        return (distance <= flank).then(|| flank_hit(function, tx, distance));
    }
    let exon_index = tx
        .exons
        .iter()
        .position(|(start, end)| variant.overlaps(*start, *end));
    if exon_index.is_none() {
        let near_boundary = tx
            .exons
            .iter()
            .any(|(start, end)| pos.abs_diff(*start) <= splice || pos.abs_diff(*end) <= splice);
        return Some(basic_hit(
            if near_boundary { 1 } else { 5 },
            if near_boundary {
                "splicing"
            } else {
                "intronic"
            },
            tx,
        ));
    }
    let exon_index = exon_index.unwrap();
    if tx.cds_start == tx.cds_end {
        return Some(basic_hit(2, "ncRNA_exonic", tx));
    }
    if variant.end <= tx.cds_start {
        return Some(basic_hit(
            3,
            if tx.strand == '+' { "UTR5" } else { "UTR3" },
            tx,
        ));
    }
    if variant.start >= tx.cds_end {
        return Some(basic_hit(
            3,
            if tx.strand == '+' { "UTR3" } else { "UTR5" },
            tx,
        ));
    }
    let mut hit = basic_hit(0, "exonic", tx);
    hit.detail = ".".into();
    let reference_len = if variant.reference == "0" {
        variant.end.saturating_sub(variant.start) as usize
    } else {
        variant.reference.len()
    };
    let delta = variant.alternate.len() as isize - reference_len as isize;
    if delta != 0 {
        let frame = if delta.unsigned_abs() % 3 == 0 {
            "nonframeshift"
        } else {
            "frameshift"
        };
        let event = if variant.reference.is_empty() {
            "insertion"
        } else if variant.alternate.is_empty() {
            "deletion"
        } else {
            "substitution"
        };
        hit.exonic_function = format!("{frame} {event}");
        hit.aa_change = format!(
            "{}:{}:exon{}:c.?",
            tx.gene,
            tx.id,
            transcript_exon_number(tx, exon_index)
        );
    } else if variant.reference.len() == 1 && variant.alternate.len() == 1 {
        annotate_snv(&mut hit, variant, tx, sequence);
    } else {
        hit.exonic_function = "nonsynonymous block substitution".into();
    }
    Some(hit)
}

fn annotate_snv(hit: &mut GeneHit, variant: &Variant, tx: &Transcript, sequence: Option<&Vec<u8>>) {
    let Some(sequence) = sequence else {
        hit.exonic_function = "unknown".into();
        return;
    };
    let Some(cdna_pos) = genomic_to_cdna(tx, variant.start) else {
        hit.exonic_function = "unknown".into();
        return;
    };
    let Some(cds_start) = genomic_to_cdna(
        tx,
        if tx.strand == '+' {
            tx.cds_start
        } else {
            tx.cds_end - 1
        },
    ) else {
        hit.exonic_function = "unknown".into();
        return;
    };
    let coding_pos = cdna_pos.abs_diff(cds_start);
    let codon_start = cds_start + (coding_pos / 3) * 3;
    if codon_start + 3 > sequence.len() {
        hit.exonic_function = "unknown".into();
        return;
    }
    let old = &sequence[codon_start..codon_start + 3];
    let mut new = old.to_vec();
    let offset = coding_pos % 3;
    let alt = if tx.strand == '+' {
        variant.alternate.as_bytes()[0]
    } else {
        complement(variant.alternate.as_bytes()[0])
    };
    new[offset] = alt;
    let old_aa = translate(old);
    let new_aa = translate(&new);
    hit.exonic_function = match (old_aa, new_aa) {
        (a, b) if a == b => "synonymous SNV",
        (_, b'*') => "stopgain",
        (b'*', _) => "stoploss",
        _ => "nonsynonymous SNV",
    }
    .into();
    let cdna_number = coding_pos + 1;
    let protein_number = coding_pos / 3 + 1;
    let reference = if tx.strand == '+' {
        variant.reference.as_bytes()[0]
    } else {
        complement(variant.reference.as_bytes()[0])
    } as char;
    let alternate = alt as char;
    let exon = tx
        .exons
        .iter()
        .position(|(start, end)| variant.overlaps(*start, *end))
        .map(|index| transcript_exon_number(tx, index))
        .unwrap_or(0);
    hit.aa_change = format!(
        "{}:{}:exon{}:c.{}{}{}:p.{}{}{}",
        tx.gene,
        tx.id,
        exon,
        reference,
        cdna_number,
        alternate,
        aa_name(old_aa),
        protein_number,
        aa_name(new_aa)
    );
}

fn basic_hit(rank: u8, function: &str, tx: &Transcript) -> GeneHit {
    GeneHit {
        rank,
        function: function.into(),
        gene: tx.gene.clone(),
        detail: if matches!(function, "splicing" | "UTR5" | "UTR3") {
            tx.id.clone()
        } else {
            ".".into()
        },
        exonic_function: ".".into(),
        aa_change: ".".into(),
    }
}

fn flank_hit(function: &str, tx: &Transcript, distance: u64) -> GeneHit {
    let mut hit = basic_hit(6, function, tx);
    hit.detail = format!("dist={distance}");
    hit
}

fn genomic_to_cdna(tx: &Transcript, genomic: u64) -> Option<usize> {
    let mut offset = 0usize;
    let iter: Box<dyn Iterator<Item = &(u64, u64)>> = if tx.strand == '+' {
        Box::new(tx.exons.iter())
    } else {
        Box::new(tx.exons.iter().rev())
    };
    for (start, end) in iter {
        if genomic >= *start && genomic < *end {
            return Some(
                offset
                    + if tx.strand == '+' {
                        (genomic - start) as usize
                    } else {
                        (end - 1 - genomic) as usize
                    },
            );
        }
        offset += (end - start) as usize;
    }
    None
}

fn transcript_exon_number(tx: &Transcript, genomic_index: usize) -> usize {
    if tx.strand == '+' {
        genomic_index + 1
    } else {
        tx.exons.len() - genomic_index
    }
}
fn parse_positions(value: &str) -> Result<Vec<u64>> {
    value
        .trim_end_matches(',')
        .split(',')
        .filter(|v| !v.is_empty())
        .map(|v| v.parse().map_err(Into::into))
        .collect()
}

fn read_fasta(path: &Path) -> Result<HashMap<String, Vec<u8>>> {
    let mut result = HashMap::new();
    let mut id = None::<String>;
    let mut sequence = Vec::new();
    for line in crate::io::open_reader(path)?.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = id.replace(
                header
                    .split_whitespace()
                    .next()
                    .context("empty FASTA header")?
                    .to_string(),
            ) {
                result.insert(previous, std::mem::take(&mut sequence));
            }
        } else {
            sequence.extend(line.trim().as_bytes().iter().map(u8::to_ascii_uppercase));
        }
    }
    if let Some(id) = id {
        result.insert(id, sequence);
    }
    Ok(result)
}

fn join_unique<'a>(values: impl Iterator<Item = &'a str>) -> String {
    let mut result: Vec<&str> = Vec::new();
    for value in values {
        if !result.contains(&value) {
            result.push(value);
        }
    }
    result.join(",")
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}
fn translate(codon: &[u8]) -> u8 {
    match codon {
        b"TTT" | b"TTC" => b'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"ATG" => b'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"TAT" | b"TAC" => b'Y',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"TGT" | b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        _ => b'X',
    }
}
fn aa_name(aa: u8) -> &'static str {
    match aa {
        b'A' => "A",
        b'R' => "R",
        b'N' => "N",
        b'D' => "D",
        b'C' => "C",
        b'Q' => "Q",
        b'E' => "E",
        b'G' => "G",
        b'H' => "H",
        b'I' => "I",
        b'L' => "L",
        b'K' => "K",
        b'M' => "M",
        b'F' => "F",
        b'P' => "P",
        b'S' => "S",
        b'T' => "T",
        b'W' => "W",
        b'Y' => "Y",
        b'V' => "V",
        b'*' => "X",
        _ => "X",
    }
}
