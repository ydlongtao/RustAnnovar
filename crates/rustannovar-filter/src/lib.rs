//! Coarse ordered chunks; no per-variant key copies or result string formatting.
use anyhow::Result;
use rayon::prelude::*;
use rustannovar_core::{AnnotationResult, Variant};
use rustannovar_db::Database;

pub struct Engine {
    pool: Option<rayon::ThreadPool>,
}
impl Engine {
    pub fn new(threads: usize) -> Result<Self> {
        anyhow::ensure!(threads > 0, "threads must be positive");
        Ok(Self {
            pool: if threads == 1 {
                None
            } else {
                Some(
                    rayon::ThreadPoolBuilder::new()
                        .num_threads(threads)
                        .build()?,
                )
            },
        })
    }
    pub fn annotate<'a, T: Sync>(
        &self,
        db: &'a Database,
        records: &[T],
        key: impl Fn(&T) -> &Variant + Sync,
    ) -> Vec<AnnotationResult<'a>> {
        let annotate = || {
            records
                .par_chunks(4096)
                .map(|chunk| chunk.iter().map(|r| db.lookup(key(r))).collect::<Vec<_>>())
                .collect::<Vec<_>>()
                .into_iter()
                .flatten()
                .collect()
        };
        match &self.pool {
            Some(pool) if records.len() >= 8192 => pool.install(annotate),
            _ => records.iter().map(|r| db.lookup(key(r))).collect(),
        }
    }
}
