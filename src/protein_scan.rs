use std::collections::HashMap;
use std::sync::Arc;
use std::thread;

#[derive(Clone)]
struct Node {
    next: [usize; 26],
    fail: usize,
    out: Vec<usize>,
}

#[derive(Clone)]
pub struct SeedAutomaton {
    nodes: Vec<Node>,
}

impl SeedAutomaton {
    pub fn build(patterns: &[Vec<u8>]) -> Self {
        let mut nodes = vec![Node { next: [usize::MAX; 26], fail: 0, out: Vec::new() }];
        for (id, pattern) in patterns.iter().enumerate() {
            if pattern.is_empty() { continue; }
            let mut state = 0;
            for &b in pattern {
                let c = (b.to_ascii_uppercase().saturating_sub(b'A')) as usize;
                if c >= 26 { continue; }
                let next = nodes[state].next[c];
                if next == usize::MAX {
                    let n = nodes.len();
                    nodes.push(Node { next: [usize::MAX; 26], fail: 0, out: Vec::new() });
                    nodes[state].next[c] = n;
                    state = n;
                } else {
                    state = next;
                }
            }
            nodes[state].out.push(id);
        }

        let mut queue = std::collections::VecDeque::new();
        for c in 0..26 {
            let next = nodes[0].next[c];
            if next == usize::MAX {
                nodes[0].next[c] = 0;
            } else {
                nodes[next].fail = 0;
                queue.push_back(next);
            }
        }

        while let Some(v) = queue.pop_front() {
            let f = nodes[v].fail;
            let inherited = nodes[f].out.clone();
            nodes[v].out.extend(inherited);
            for c in 0..26 {
                let next = nodes[v].next[c];
                if next == usize::MAX {
                    nodes[v].next[c] = nodes[f].next[c];
                } else {
                    nodes[next].fail = nodes[f].next[c];
                    queue.push_back(next);
                }
            }
        }

        Self { nodes }
    }

    fn find(&self, seq: &[u8], mut hit: impl FnMut(usize, usize)) {
        let mut state = 0;
        for (pos, &b) in seq.iter().enumerate() {
            let c = (b.saturating_sub(b'A')) as usize;
            if c >= 26 {
                state = 0;
                continue;
            }
            state = self.nodes[state].next[c];
            for &pattern_id in &self.nodes[state].out {
                hit(pattern_id, pos + 1);
            }
        }
    }
}

#[derive(Clone)]
pub struct PeptideSeed {
    pub row: usize,
    pub offset: usize,
    pub len: usize,
}

pub fn build_seed_automaton(peptides: &[String], max_mismatches: usize) -> (Arc<SeedAutomaton>, Vec<PeptideSeed>) {
    let mut patterns = Vec::new();
    let mut seeds = Vec::new();
    let n_parts = max_mismatches + 1;

    for (row, peptide) in peptides.iter().enumerate() {
        let bytes = peptide.as_bytes();
        if bytes.is_empty() { continue; }
        for part in 0..n_parts.min(bytes.len()) {
            let start = part * bytes.len() / n_parts;
            let end = (part + 1) * bytes.len() / n_parts;
            if end > start {
                patterns.push(bytes[start..end].to_vec());
                seeds.push(PeptideSeed { row, offset: start, len: end - start });
            }
        }
    }

    (Arc::new(SeedAutomaton::build(&patterns)), seeds)
}

pub type ProteinHit = (usize, String, String, Vec<usize>);

pub fn scan_batch(
    records: &[(String, String)],
    automaton: Arc<SeedAutomaton>,
    seeds: Arc<Vec<PeptideSeed>>,
    peptides: Arc<Vec<Vec<u8>>>,
    max_mismatches: usize,
    threads: usize,
) -> Vec<ProteinHit> {
    let workers = threads.max(1).min(records.len().max(1));
    let mut result = Vec::new();

    thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        for chunk in records.chunks(records.len().div_ceil(workers).max(1)) {
            let automaton = Arc::clone(&automaton);
            let seeds = Arc::clone(&seeds);
            let peptides = Arc::clone(&peptides);
            handles.push(scope.spawn(move || {
                let mut local = Vec::new();
                for (header, seq) in chunk {
                    let bytes = seq.as_bytes();
                    let mut candidates: HashMap<usize, Vec<usize>> = HashMap::new();
                    automaton.find(bytes, |seed_id, end| {
                        let seed = &seeds[seed_id];
                        if end >= seed.offset + seed.len {
                            let start = end - seed.offset - seed.len;
                            candidates.entry(seed.row).or_default().push(start);
                        }
                    });

                    let mut accepted_hits = Vec::new();
                    for (row, mut positions) in candidates {
                        positions.sort_unstable();
                        positions.dedup();
                        let peptide = &peptides[row];
                        let mut accepted = Vec::new();
                        for pos in positions {
                            if pos + peptide.len() > bytes.len() { continue; }
                            let mismatches = bytes[pos..pos + peptide.len()]
                                .iter()
                                .zip(peptide)
                                .filter(|(a, b)| a != b)
                                .count();
                            if mismatches <= max_mismatches {
                                accepted.push(pos + 1);
                            }
                        }
                        accepted.sort_unstable();
                        accepted.dedup();
                        if !accepted.is_empty() {
                            accepted_hits.push((row, accepted));
                        }
                    }
                    if !accepted_hits.is_empty() {
                        let (acc, _name, organism, _gene) = crate::parse_uniprot_header(header);
                        let organism = if organism.is_empty() { "unknown".to_string() } else { organism };
                        for (row, positions) in accepted_hits {
                            local.push((row, acc.clone(), organism.clone(), positions));
                        }
                    }
                }
                local
            }));
        }
        for handle in handles {
            result.extend(handle.join().expect("protein scan worker panicked"));
        }
    });

    result
}
