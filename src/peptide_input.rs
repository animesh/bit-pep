/// peptide_input.rs — parse peptide lists for bit-pop protein mode
///
/// Accepts two formats:
///   1. One peptide per line:   ACDEFGHIKLM
///   2. Space-separated per line: ACDE FGHIK LMNPQ
///
/// Lines starting with '#' are comments.
/// Lines starting with '>' trigger FASTA-style sequence accumulation
/// (so a FASTA file of peptides also works).
///
/// Sequences are uppercased and validated against the 20 standard AAs
/// plus ambiguous codes (B, X, Z, U, O).  Sequences shorter than MIN_LEN
/// (default 6) or containing only ambiguous residues are skipped with a warning.

use std::fs;
use std::io::{BufRead, BufReader};

pub const MIN_PEPTIDE_LEN: usize = 6;
pub const MAX_PEPTIDE_LEN: usize = 50; // soft guard; longer sequences accepted but warned

const VALID_AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWYBXZUO*";

#[derive(Debug, Clone)]
pub struct Peptide {
    /// Auto-assigned identifier: "pep_N" or FASTA header if available.
    pub id: String,
    /// Upper-cased, validated sequence (no stop codon).
    pub sequence: String,
}

/// Parse a peptide input file.  Returns all valid peptides found.
pub fn parse_peptide_file(path: &str) -> anyhow::Result<Vec<Peptide>> {
    let file = fs::File::open(path)
        .map_err(|e| anyhow::anyhow!("Cannot open peptide file '{}': {}", path, e))?;
    let reader = BufReader::new(file);
    parse_peptide_reader(reader)
}

/// Parse peptides from a string (useful for tests and stdin piping).
pub fn parse_peptide_str(input: &str) -> anyhow::Result<Vec<Peptide>> {
    parse_peptide_reader(BufReader::new(input.as_bytes()))
}

fn parse_peptide_reader<R: std::io::Read>(reader: BufReader<R>) -> anyhow::Result<Vec<Peptide>> {
    let mut peptides: Vec<Peptide> = Vec::new();
    let mut idx: usize = 0;

    // FASTA mode state
    let mut fasta_mode = false;
    let mut current_id: Option<String> = None;
    let mut current_seq = String::new();

    let flush_fasta = |id: &Option<String>, seq: &str, peptides: &mut Vec<Peptide>, idx: &mut usize| {
        if !seq.is_empty() {
            if let Some(validated) = validate_and_clean(seq) {
                peptides.push(Peptide {
                    id: id.clone().unwrap_or_else(|| format!("pep_{}", idx)),
                    sequence: validated,
                });
                *idx += 1;
            } else {
                eprintln!("WARN: skipping invalid/short FASTA entry '{}'", id.as_deref().unwrap_or("?"));
            }
        }
    };

    for line in reader.lines() {
        let raw = line?;
        let line = raw.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        if line.starts_with('>') {
            // FASTA header encountered
            if fasta_mode {
                flush_fasta(&current_id, &current_seq, &mut peptides, &mut idx);
            }
            fasta_mode = true;
            // Take everything after '>' up to first space as the ID
            let header = &line[1..];
            current_id = Some(header.split_whitespace().next().unwrap_or("pep").to_string());
            current_seq.clear();
            continue;
        }

        if fasta_mode {
            // Accumulate sequence lines (multi-line FASTA)
            current_seq.push_str(line);
            continue;
        }

        // Plain mode: one or more space-separated peptides per line
        for token in line.split_whitespace() {
            let upper = token.to_uppercase();
            match validate_and_clean(&upper) {
                Some(seq) => {
                    peptides.push(Peptide {
                        id: format!("pep_{}", idx),
                        sequence: seq,
                    });
                    idx += 1;
                }
                None => {
                    eprintln!("WARN: skipping '{}' (too short, too long, or invalid chars)", token);
                }
            }
        }
    }

    // Flush last FASTA entry if in FASTA mode
    if fasta_mode {
        flush_fasta(&current_id, &current_seq, &mut peptides, &mut idx);
    }

    if peptides.is_empty() {
        anyhow::bail!("No valid peptides found in input file");
    }

    eprintln!("INFO: loaded {} peptides", peptides.len());
    Ok(peptides)
}

/// Returns Some(cleaned_seq) if valid, None if it should be skipped.
fn validate_and_clean(seq: &str) -> Option<String> {
    // Strip stop codon suffix if present
    let seq = seq.trim_end_matches('*');

    if seq.len() < MIN_PEPTIDE_LEN {
        return None;
    }
    if seq.len() > MAX_PEPTIDE_LEN {
        eprintln!("WARN: peptide length {} > {} — accepted but may be slow: {}…",
                  seq.len(), MAX_PEPTIDE_LEN, &seq[..20]);
    }

    // Validate characters
    let valid = seq.bytes().all(|c| VALID_AA.contains(&c.to_ascii_uppercase()));
    if !valid {
        return None;
    }

    Some(seq.to_uppercase())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_one_per_line() {
        let input = "ACDEFGHIK\nLMNPQRSTV\n";
        let peps = parse_peptide_str(input).unwrap();
        assert_eq!(peps.len(), 2);
        assert_eq!(peps[0].sequence, "ACDEFGHIK");
        assert_eq!(peps[0].id, "pep_0");
    }

    #[test]
    fn parse_space_separated() {
        let input = "ACDEFG HIJKLM NPQRST";
        let peps = parse_peptide_str(input).unwrap();
        assert_eq!(peps.len(), 3);
    }

    #[test]
    fn parse_fasta_mode() {
        let input = ">seq1 some description\nACDEFGHIKL\n>seq2\nMNPQRSTVWY\n";
        let peps = parse_peptide_str(input).unwrap();
        assert_eq!(peps.len(), 2);
        assert_eq!(peps[0].id, "seq1");
        assert_eq!(peps[1].sequence, "MNPQRSTVWY");
    }

    #[test]
    fn skips_too_short() {
        let input = "ACD\nACDEFGHIKL\n"; // "ACD" is 3 < MIN_PEPTIDE_LEN
        let peps = parse_peptide_str(input).unwrap();
        assert_eq!(peps.len(), 1);
    }

    #[test]
    fn strips_stop_codon() {
        let input = "ACDEFGHIK*\n";
        let peps = parse_peptide_str(input).unwrap();
        assert_eq!(peps[0].sequence, "ACDEFGHIK");
    }

    #[test]
    fn comments_ignored() {
        let input = "# this is a comment\nACDEFGHIK\n";
        let peps = parse_peptide_str(input).unwrap();
        assert_eq!(peps.len(), 1);
    }
}
