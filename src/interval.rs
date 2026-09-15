//! Static interval index. Input order is preserved among query hits.
#[derive(Debug, Default)]
pub struct IntervalIndex {
    entries: Vec<(u64, u64, usize)>,
    max_end: Vec<u64>,
}
impl IntervalIndex {
    pub fn new(intervals: impl Iterator<Item = (u64, u64)>) -> Self {
        let mut entries: Vec<_> = intervals.enumerate().map(|(i, (s, e))| (s, e, i)).collect();
        entries.sort_by_key(|v| (v.0, v.2));
        let mut index = Self {
            max_end: vec![0; entries.len().saturating_mul(4)],
            entries,
        };
        if !index.entries.is_empty() {
            index.build(1, 0, index.entries.len());
        }
        index
    }
    fn build(&mut self, node: usize, lo: usize, hi: usize) -> u64 {
        let value = if hi - lo == 1 {
            self.entries[lo].1
        } else {
            let mid = (lo + hi) / 2;
            self.build(node * 2, lo, mid)
                .max(self.build(node * 2 + 1, mid, hi))
        };
        self.max_end[node] = value;
        value
    }
    /// Inclusive candidate bounds deliberately include insertion boundaries;
    /// the caller applies its exact overlap policy.
    pub fn query(&self, start: u64, end: u64) -> Vec<usize> {
        let mut hits = Vec::new();
        if !self.entries.is_empty() {
            self.visit(1, 0, self.entries.len(), start, end, &mut hits);
        }
        hits.sort_unstable();
        hits
    }
    fn visit(
        &self,
        node: usize,
        lo: usize,
        hi: usize,
        start: u64,
        end: u64,
        hits: &mut Vec<usize>,
    ) {
        if self.max_end[node] < start || self.entries[lo].0 > end {
            return;
        }
        if hi - lo == 1 {
            hits.push(self.entries[lo].2);
            return;
        }
        let mid = (lo + hi) / 2;
        self.visit(node * 2, lo, mid, start, end, hits);
        self.visit(node * 2 + 1, mid, hi, start, end, hits);
    }
}
