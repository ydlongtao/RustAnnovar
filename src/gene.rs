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
struct TranscriptMap {
    prefix: Vec<usize>,
    length: usize,
    cds: Option<(usize, usize)>,
}
impl TranscriptMap {
    fn new(tx: &Transcript) -> Self {
        let mut prefix = Vec::new();
        let mut length = 0;
        for &(s, e) in &tx.exons {
            prefix.push(length);
            length += (e - s) as usize;
        }
        let cds = coding_bounds(tx);
        Self {
            prefix,
            length,
            cds,
        }
    }
    fn position(&self, tx: &Transcript, g: u64) -> Option<usize> {
        let i = tx.exons.partition_point(|(s, _)| *s <= g).checked_sub(1)?;
        let (s, e) = tx.exons[i];
        if g >= e {
            return None;
        }
        let pos = self.prefix[i] + (g - s) as usize;
        Some(if tx.strand == '+' {
            pos
        } else {
            self.length - 1 - pos
        })
    }
}

#[derive(Debug)]
pub struct GeneDatabase {
    maps: HashMap<String, Vec<TranscriptMap>>,
    indexes: HashMap<String, crate::interval::IntervalIndex>,
    by_end: HashMap<String, Vec<usize>>,
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
            if !matches!(transcript.strand, '+' | '-')
                || transcript.tx_start >= transcript.tx_end
                || transcript.cds_start > transcript.cds_end
                || transcript.cds_start < transcript.tx_start
                || transcript.cds_end > transcript.tx_end
                || transcript.exons.is_empty()
                || transcript
                    .exons
                    .iter()
                    .any(|&(s, e)| s >= e || s < transcript.tx_start || e > transcript.tx_end)
                || transcript.exons.windows(2).any(|p| p[0].1 > p[1].0)
            {
                bail!("invalid transcript intervals at line {}", line_no + 1);
            }
            by_chrom
                .entry(transcript.chrom.clone())
                .or_default()
                .push(transcript);
        }
        for transcripts in by_chrom.values_mut() {
            transcripts.sort_by_key(|tx| tx.tx_start);
        }
        let sequences = fasta_path.map(read_fasta).transpose()?.unwrap_or_default();
        let indexes = by_chrom
            .iter()
            .map(|(chrom, txs)| {
                (
                    chrom.clone(),
                    crate::interval::IntervalIndex::new(txs.iter().map(|t| (t.tx_start, t.tx_end))),
                )
            })
            .collect();
        let by_end = by_chrom
            .iter()
            .map(|(chrom, txs)| {
                let mut ids: Vec<_> = (0..txs.len()).collect();
                ids.sort_by_key(|&i| txs[i].tx_end);
                (chrom.clone(), ids)
            })
            .collect();
        let maps = by_chrom
            .iter()
            .map(|(c, ts)| (c.clone(), ts.iter().map(TranscriptMap::new).collect()))
            .collect();
        Ok(Self {
            maps,
            by_chrom,
            sequences,
            indexes,
            by_end,
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
            for index in self.indexes[&variant.chrom].query(
                variant.start.saturating_sub(flank),
                variant.end.saturating_add(flank),
            ) {
                let transcript = &transcripts[index];
                if transcript.tx_start > variant.end.saturating_add(flank) {
                    break;
                }
                if transcript.tx_end.saturating_add(flank) < variant.start {
                    continue;
                }
                if let Some(hit) = classify(
                    variant,
                    transcript,
                    &self.maps[&variant.chrom][index],
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
        let ends = &self.by_end[&variant.chrom];
        let left_stop = ends.partition_point(|&i| transcripts[i].tx_end <= variant.start);
        let right_start = transcripts.partition_point(|t| t.tx_start < variant.end);
        let mut nearest: Vec<(String, u64)> = Vec::new();
        if left_stop > 0 {
            let end = transcripts[ends[left_stop - 1]].tx_end;
            let start = ends[..left_stop].partition_point(|&i| transcripts[i].tx_end < end);
            let names = transcripts[ends[start]].gene.clone();
            nearest.push((
                names,
                variant.start - end + u64::from(!variant.reference.is_empty()),
            ));
        } else {
            nearest.push(("NONE".into(), 0));
        }
        if right_start < transcripts.len() {
            let start = transcripts[right_start].tx_start;
            let _stop = transcripts.partition_point(|t| t.tx_start <= start);
            nearest.push((
                transcripts[right_start].gene.clone(),
                start - variant.end + 1,
            ));
        } else {
            nearest.push(("NONE".into(), 0));
        }
        GeneHit {
            rank: 8,
            function: "intergenic".into(),
            gene: nearest
                .iter()
                .map(|(gene, _)| gene.as_str())
                .collect::<Vec<_>>()
                .join(","),
            detail: nearest
                .iter()
                .map(|(gene, distance)| {
                    if gene == "NONE" {
                        "dist=NONE".into()
                    } else {
                        format!("dist={distance}")
                    }
                })
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
    mapping: &TranscriptMap,
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
        let near_boundary = tx.exons.windows(2).any(|pair| {
            variant.overlaps(pair[0].1, pair[0].1.saturating_add(splice))
                || variant.overlaps(pair[1].0.saturating_sub(splice), pair[1].0)
        });
        let mut hit = basic_hit(
            if tx.cds_start == tx.cds_end {
                2
            } else if near_boundary {
                1
            } else {
                5
            },
            if near_boundary {
                if tx.cds_start == tx.cds_end {
                    "ncRNA_splicing"
                } else {
                    "splicing"
                }
            } else {
                if tx.cds_start == tx.cds_end {
                    "ncRNA_intronic"
                } else {
                    "intronic"
                }
            },
            tx,
        );
        if near_boundary {
            hit.detail =
                noncoding_detail(variant, tx, mapping, true).unwrap_or_else(|| tx.id.clone());
        }
        return Some(hit);
    }
    let exon_index = exon_index.unwrap();
    if tx.cds_start == tx.cds_end {
        return Some(basic_hit(2, "ncRNA_exonic", tx));
    }
    if variant.end <= tx.cds_start || variant.start >= tx.cds_end {
        let before = variant.end <= tx.cds_start;
        let mut hit = basic_hit(
            3,
            if before == (tx.strand == '+') {
                "UTR5"
            } else {
                "UTR3"
            },
            tx,
        );
        hit.detail = noncoding_detail(variant, tx, mapping, false).unwrap_or_else(|| tx.id.clone());
        return Some(hit);
    }
    let mut hit = basic_hit(0, "exonic", tx);
    hit.detail = ".".into();
    let contained = variant.start >= tx.exons[exon_index].0
        && variant.end <= tx.exons[exon_index].1
        && variant.start >= tx.cds_start
        && variant.end <= tx.cds_end
        && !(variant.start == variant.end
            && (variant.start == tx.exons[exon_index].0
                || variant.end == tx.exons[exon_index].1
                || variant.start == tx.cds_start
                || variant.end == tx.cds_end));
    if !contained || variant.reference == "0" {
        hit.exonic_function = "unknown".into();
        hit.detail = "unsupported_coding_boundary".into();
    } else if variant.reference.len() == 1 && variant.alternate.len() == 1 {
        annotate_snv(&mut hit, variant, tx, mapping, sequence);
    } else {
        annotate_small_change(&mut hit, variant, tx, mapping, sequence, exon_index);
    }
    Some(hit)
}

fn annotate_snv(
    hit: &mut GeneHit,
    variant: &Variant,
    tx: &Transcript,
    mapping: &TranscriptMap,
    sequence: Option<&Vec<u8>>,
) {
    let Some(sequence) = sequence else {
        hit.exonic_function = "unknown".into();
        hit.detail = "missing_transcript_sequence".into();
        return;
    };
    let Some(cdna_pos) = mapping.position(tx, variant.start) else {
        hit.exonic_function = "unknown".into();
        return;
    };
    let Some((cds_start, cds_end)) = mapping.cds else {
        hit.exonic_function = "unknown".into();
        hit.detail = "invalid_cds".into();
        return;
    };
    if cds_end > sequence.len() || cds_end <= cds_start || (cds_end - cds_start) % 3 != 0 {
        hit.exonic_function = "unknown".into();
        hit.detail = "invalid_cds".into();
        return;
    }
    let coding_pos = cdna_pos.abs_diff(cds_start);
    let codon_start = cds_start + (coding_pos / 3) * 3;
    if codon_start + 3 > sequence.len() {
        hit.exonic_function = "unknown".into();
        return;
    }
    let expected = if tx.strand == '+' {
        variant.reference.as_bytes()[0]
    } else {
        complement(variant.reference.as_bytes()[0])
    };
    if sequence.get(cdna_pos) != Some(&expected) {
        hit.exonic_function = "unknown".into();
        hit.detail = "reference_mismatch".into();
        return;
    }
    let old = &sequence[codon_start..codon_start + 3];
    if old.iter().any(|b| !b"ACGT".contains(b))
        || !b"ACGT".contains(&variant.alternate.as_bytes()[0])
    {
        hit.exonic_function = "unknown".into();
        hit.detail = "ambiguous_sequence".into();
        return;
    }
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

fn oriented(s: &str, strand: char) -> Vec<u8> {
    if strand == '+' {
        s.as_bytes().to_vec()
    } else {
        s.bytes().rev().map(complement).collect()
    }
}
fn coding_bounds(tx: &Transcript) -> Option<(usize, usize)> {
    let a = genomic_to_cdna(
        tx,
        if tx.strand == '+' {
            tx.cds_start
        } else {
            tx.cds_end.checked_sub(1)?
        },
    )?;
    let b = genomic_to_cdna(
        tx,
        if tx.strand == '+' {
            tx.cds_end.checked_sub(1)?
        } else {
            tx.cds_start
        },
    )? + 1;
    Some((a, b))
}
fn annotate_small_change(
    hit: &mut GeneHit,
    v: &Variant,
    tx: &Transcript,
    mapping: &TranscriptMap,
    seq: Option<&Vec<u8>>,
    exon: usize,
) {
    hit.exonic_function = "unknown".into();
    let Some(seq) = seq else {
        hit.detail = "missing_transcript_sequence".into();
        return;
    };
    let Some((cs, ce)) = mapping.cds else {
        hit.detail = "invalid_cds".into();
        return;
    };
    if ce > seq.len() || ce <= cs || (ce - cs) % 3 != 0 {
        hit.detail = "invalid_cds".into();
        return;
    }
    let g = if tx.strand == '+' {
        v.start
    } else if v.start == v.end {
        v.start - 1
    } else {
        v.end - 1
    };
    let Some(at) = mapping.position(tx, g) else {
        hit.detail = "unsupported_coding_boundary".into();
        return;
    };
    let r = oriented(&v.reference, tx.strand);
    let a = oriented(&v.alternate, tx.strand);
    if at < cs || at + r.len() > ce || seq.get(at..at + r.len()) != Some(r.as_slice()) {
        hit.detail = "reference_mismatch".into();
        return;
    }
    if seq[cs..ce]
        .iter()
        .chain(a.iter())
        .any(|b| !b"ACGT".contains(b))
    {
        hit.detail = "ambiguous_sequence".into();
        return;
    }
    let mut mutant = seq[cs..ce].to_vec();
    let offset = at - cs;
    mutant.splice(offset..offset + r.len(), a.iter().copied());
    let old: Vec<_> = seq[cs..ce].chunks_exact(3).map(translate).collect();
    let new: Vec<_> = mutant.chunks_exact(3).map(translate).collect();
    let delta = a.len() as isize - r.len() as isize;
    hit.exonic_function = if delta == 0 {
        if old == new {
            "synonymous block substitution".into()
        } else if new.iter().position(|&b| b == b'*').unwrap_or(usize::MAX)
            < old.iter().position(|&b| b == b'*').unwrap_or(usize::MAX)
        {
            "stopgain".into()
        } else if old.iter().position(|&b| b == b'*').unwrap_or(usize::MAX)
            < new.iter().position(|&b| b == b'*').unwrap_or(usize::MAX)
        {
            "stoploss".into()
        } else {
            "nonsynonymous block substitution".into()
        }
    } else {
        format!(
            "{} {}",
            if delta.unsigned_abs() % 3 == 0 {
                "nonframeshift"
            } else {
                "frameshift"
            },
            if r.is_empty() {
                "insertion"
            } else if a.is_empty() {
                "deletion"
            } else {
                "substitution"
            }
        )
    };
    let c = if r.is_empty() {
        format!(
            "{}_{}ins{}",
            offset,
            offset + 1,
            String::from_utf8_lossy(&a)
        )
    } else if a.is_empty() {
        format!("{}_{}del", offset + 1, offset + r.len())
    } else {
        format!(
            "{}_{}delins{}",
            offset + 1,
            offset + r.len(),
            String::from_utf8_lossy(&a)
        )
    };
    let first = old
        .iter()
        .zip(&new)
        .position(|(a, b)| a != b)
        .unwrap_or(old.len().min(new.len()));
    let protein = if old == new {
        "=".into()
    } else if delta.unsigned_abs() % 3 != 0 {
        format!(
            "{}{}{}fs",
            aa_name(*old.get(first).unwrap_or(&b'X')),
            first + 1,
            aa_name(*new.get(first).unwrap_or(&b'X'))
        )
    } else {
        let mut end_old = old.len();
        let mut end_new = new.len();
        while end_old > first && end_new > first && old[end_old - 1] == new[end_new - 1] {
            end_old -= 1;
            end_new -= 1;
        }
        let changed = &new[first..end_new];
        let residue = |i: usize| format!("{}{}", aa_name(*old.get(i).unwrap_or(&b'X')), i + 1);
        if end_old == first {
            format!(
                "{}_{}ins{}",
                residue(first.saturating_sub(1)),
                residue(first),
                String::from_utf8_lossy(changed)
            )
        } else {
            let span = if end_old == first + 1 {
                residue(first)
            } else {
                format!("{}_{}", residue(first), residue(end_old - 1))
            };
            if changed.is_empty() {
                format!("{span}del")
            } else if end_old == first + 1 && changed.len() == 1 {
                format!("{span}{}", aa_name(changed[0]))
            } else {
                format!("{span}delins{}", String::from_utf8_lossy(changed))
            }
        }
    };
    hit.aa_change = format!(
        "{}:{}:exon{}:c.{}:p.{}",
        tx.gene,
        tx.id,
        transcript_exon_number(tx, exon),
        c,
        protein
    );
}
fn noncoding_detail(
    v: &Variant,
    tx: &Transcript,
    mapping: &TranscriptMap,
    intron: bool,
) -> Option<String> {
    if v.reference.len() != 1 || v.alternate.len() != 1 || tx.cds_start == tx.cds_end {
        return None;
    }
    let (cs, ce) = mapping.cds?;
    let r = oriented(&v.reference, tx.strand);
    let a = oriented(&v.alternate, tx.strand);
    if !intron {
        let at = mapping.position(tx, v.start)?;
        let c = if at < cs {
            format!("-{}", cs - at)
        } else {
            format!("*{}", at + 1 - ce)
        };
        return Some(format!(
            "{}:c.{}{}>{}",
            tx.id, c, r[0] as char, a[0] as char
        ));
    }
    let (exon, anchor) = tx
        .exons
        .iter()
        .enumerate()
        .flat_map(|(i, (s, e))| [(i, *s), (i, e - 1)])
        .min_by_key(|(_, g)| v.start.abs_diff(*g))?;
    let at = mapping.position(tx, anchor)?;
    if at < cs || at >= ce {
        return None;
    }
    let sign = if (v.start > anchor) == (tx.strand == '+') {
        "+"
    } else {
        "-"
    };
    Some(format!(
        "{}:exon{}:c.{}{}{}{}>{}",
        tx.id,
        transcript_exon_number(tx, exon),
        at - cs + 1,
        sign,
        v.start.abs_diff(anchor),
        r[0] as char,
        a[0] as char
    ))
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
