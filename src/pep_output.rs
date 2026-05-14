/// pep_output.rs — TSV output for peptide-to-proteome mapping results
///
/// Replaces sam.rs.  SAM is genomic-coordinate-centric and maps poorly
/// to proteomics.  We output a simple tab-separated file that downstream
/// tools (R, Python, Pandas) can consume directly.
///
/// Columns:
///   peptide_id      — input ID (e.g. "pep_0" or FASTA header)
///   sequence        — peptide sequence as submitted
///   protein_acc     — UniProt accession (e.g. "P04637")
///   protein_name    — description from FASTA header
///   proteome_id     — UPID the protein came from (e.g. "UP000005640")
///   start           — 1-based start position in protein sequence
///   end             — 1-based end position (inclusive)
///   score           — alignment score [0.0, 1.0]
///   mismatches      — number of residue mismatches
///   status          — "unique"  → maps to exactly 1 protein
///                     "shared"  → maps to >1 protein in same proteome
///                     "xproteome" → maps across different proteomes
///                     "unmapped"  → no hit above threshold

use std::io::{BufWriter, Write};
use std::fs;

pub const HEADER: &str =
    "peptide_id\tsequence\tprotein_acc\tprotein_name\tproteome_id\tstart\tend\tscore\tmismatches\tstatus";

#[derive(Debug, Clone)]
pub struct PepHit {
    pub peptide_id:   String,
    pub sequence:     String,
    pub protein_acc:  String,
    pub protein_name: String,
    pub proteome_id:  String,
    pub start:        usize,   // 1-based
    pub end:          usize,   // 1-based inclusive
    pub score:        f32,
    pub mismatches:   u32,
    pub status:       HitStatus,
}

#[derive(Debug, Clone, PartialEq)]
pub enum HitStatus {
    Unique,
    Shared,
    CrossProteome,
    Unmapped,
}

impl HitStatus {
    pub fn as_str(&self) -> &'static str {
        match self {
            HitStatus::Unique       => "unique",
            HitStatus::Shared       => "shared",
            HitStatus::CrossProteome => "xproteome",
            HitStatus::Unmapped     => "unmapped",
        }
    }
}

/// Assign status based on how many proteins a peptide matched.
pub fn assign_status(hits: &mut Vec<PepHit>) {
    if hits.is_empty() { return; }

    let n_proteins: std::collections::HashSet<&str> =
        hits.iter().map(|h| h.protein_acc.as_str()).collect();
    let n_proteomes: std::collections::HashSet<&str> =
        hits.iter().map(|h| h.proteome_id.as_str()).collect();

    let status = if n_proteomes.len() > 1 {
        HitStatus::CrossProteome
    } else if n_proteins.len() > 1 {
        HitStatus::Shared
    } else {
        HitStatus::Unique
    };

    for h in hits.iter_mut() {
        h.status = status.clone();
    }
}

/// Write all hits to a TSV file.
pub fn write_tsv(hits: &[PepHit], path: &str) -> anyhow::Result<()> {
    let file = fs::File::create(path)
        .map_err(|e| anyhow::anyhow!("Cannot create output file '{}': {}", path, e))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "{}", HEADER)?;
    for h in hits {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{}",
            h.peptide_id, h.sequence, h.protein_acc, h.protein_name,
            h.proteome_id, h.start, h.end, h.score, h.mismatches,
            h.status.as_str()
        )?;
    }
    Ok(())
}

/// Write a brief mapping summary to stderr.
pub fn print_summary(hits: &[PepHit], n_peptides: usize) {
    let mapped: std::collections::HashSet<&str> =
        hits.iter()
            .filter(|h| h.status != HitStatus::Unmapped)
            .map(|h| h.peptide_id.as_str())
            .collect();

    let unique = hits.iter().filter(|h| h.status == HitStatus::Unique).count();
    let shared = hits.iter().filter(|h| h.status == HitStatus::Shared).count();
    let xprot  = hits.iter().filter(|h| h.status == HitStatus::CrossProteome).count();

    eprintln!("─────────────────────────────────────────");
    eprintln!("Peptides submitted : {}", n_peptides);
    eprintln!("Mapped             : {} ({:.1}%)",
              mapped.len(), 100.0 * mapped.len() as f64 / n_peptides as f64);
    eprintln!("  unique           : {}", unique);
    eprintln!("  shared (1 prot.) : {}", shared);
    eprintln!("  cross-proteome   : {}", xprot);
    eprintln!("Unmapped           : {}", n_peptides - mapped.len());
    eprintln!("─────────────────────────────────────────");
}
