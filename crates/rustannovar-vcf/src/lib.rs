//! Input layer. The workspace MVP deliberately exposes AVinput only;
//! the existing `rust-annovar` binary retains full VCF handling.
use anyhow::{Context, Result};
use rustannovar_core::Variant;
use std::io::BufRead;

pub struct Record {
    pub variant: Variant,
    pub raw: String,
    /// End of the five core TSV columns when raw bytes can be emitted directly.
    pub tsv_core_end: Option<usize>,
}
pub struct AvinputReader<R> {
    reader: R,
    line: u64,
    row: u64,
    buffer: String,
}
impl<R: BufRead> AvinputReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line: 0,
            row: 0,
            buffer: String::new(),
        }
    }
    pub fn next_record(&mut self) -> Result<Option<Record>> {
        loop {
            self.buffer.clear();
            if self
                .reader
                .read_line(&mut self.buffer)
                .with_context(|| format!("reading input near line {}", self.line + 1))?
                == 0
            {
                return Ok(None);
            }
            self.line += 1;
            let line = self.buffer.trim_end_matches(['\r', '\n']);
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let variant = Variant::parse_avinput(line, self.row)
                .with_context(|| format!("input line {}", self.line))?;
            self.row += 1;
            return Ok(Some(Record {
                variant,
                raw: line.to_owned(),
                tsv_core_end: if line.contains('\t') && !line.contains(['"', '\r']) {
                    Some(
                        line.bytes()
                            .enumerate()
                            .filter(|(_, b)| *b == b'\t')
                            .nth(4)
                            .map_or(line.len(), |(i, _)| i),
                    )
                } else {
                    None
                },
            }));
        }
    }
}
