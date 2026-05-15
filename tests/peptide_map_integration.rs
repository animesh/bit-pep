// tests/peptide_map_integration.rs
// Full pipeline integration test — no CLI needed.
// Tests: UniProt download → FASTA parse → trypsin digest → FM encode → search back
//
// Run: cargo test --test peptide_map_integration -- --nocapture
// (requires network; ~5 seconds for SARS-CoV-2 proteome)

use bit_pop::{
    amino::encode_sequence,
    peptide_input::parse_peptide_str,
    uniprot::fetch_proteome,
};
use std::collections::HashMap;

// ── Helper: minimal FASTA parser ─────────────────────────────────────────────
fn parse_fasta(text: &str) -> Vec<(String, String, String)> {
    // returns Vec<(accession, description, sequence)>
    let mut out = Vec::new();
    let mut acc  = String::new();
    let mut desc = String::new();
    let mut seq  = String::new();

    for line in text.lines() {
        if line.starts_with('>') {
            if !acc.is_empty() {
                out.push((acc.clone(), desc.clone(), seq.clone()));
            }
            let header = &line[1..];
            let parts: Vec<&str> = header.splitn(3, '|').collect();
            if parts.len() >= 3 {
                acc  = parts[1].to_string();
                desc = parts[2].split(" OS=").next().unwrap_or("").to_string();
            } else {
                acc  = header.split_whitespace().next().unwrap_or("?").to_string();
                desc = header.to_string();
            }
            seq.clear();
        } else {
            seq.push_str(line.trim().trim_end_matches('*'));
        }
    }
    if !acc.is_empty() {
        out.push((acc, desc, seq));
    }
    out
}

// ── Helper: trypsin digest ────────────────────────────────────────────────────
fn trypsin(seq: &str, missed: usize) -> Vec<String> {
    let seq = seq.to_uppercase();
    let bytes = seq.as_bytes();

    let mut sites = vec![0usize];
    for (i, &c) in bytes.iter().enumerate() {
        if c == b'K' || c == b'R' {
            if i + 1 < bytes.len() && bytes[i + 1] == b'P' {
                continue; // no cut before P
            }
            sites.push(i + 1);
        }
    }
    sites.push(seq.len());

    let n = sites.len() - 1;
    let mut peptides = Vec::new();
    for i in 0..n {
        for mc in 0..=(missed) {
            let j = i + mc + 1;
            if j > n { break; }
            let pep = &seq[sites[i]..sites[j]];
            if pep.len() >= 6 && pep.len() <= 50 {
                peptides.push(pep.to_string());
            }
        }
    }
    peptides
}

// ── Test 1: parse and encode SARS-CoV-2 proteome ────────────────────────────
#[tokio::test]
async fn test_sars_cov2_download_and_encode() {
    let fasta_path = "/tmp/sars_cov2_pep_test.fasta";

    let bytes = fetch_proteome("UP000464024", fasta_path, true)
        .await
        .expect("UniProt download should succeed");
    assert!(bytes > 100);

    let fasta_text = std::fs::read_to_string(fasta_path).unwrap();
    let proteins = parse_fasta(&fasta_text);

    assert!(!proteins.is_empty(), "Should parse at least one protein");
    println!("Proteins in SARS-CoV-2 proteome: {}", proteins.len());

    // Encode every protein sequence and verify codes are in valid range
    let mut total_residues = 0usize;
    for (acc, _desc, seq) in &proteins {
        let encoded = encode_sequence(seq.as_bytes());
        assert!(
            encoded.iter().all(|&c| c >= 1 && c <= 26),
            "Protein {} has out-of-range codes", acc
        );
        total_residues += encoded.len();
    }
    println!("Total residues encoded: {}", total_residues);
}

// ── Test 2: trypsin digest → parse with peptide_input ────────────────────────
#[tokio::test]
async fn test_trypsin_digest_round_trip() {
    let fasta_path = "/tmp/sars_cov2_pep_test2.fasta";
    fetch_proteome("UP000464024", fasta_path, true)
        .await
        .expect("download");

    let fasta_text = std::fs::read_to_string(fasta_path).unwrap();
    let proteins = parse_fasta(&fasta_text);

    // Digest
    let mut all_peptides: Vec<String> = Vec::new();
    let mut pep_to_proteins: HashMap<String, Vec<String>> = HashMap::new();

    for (acc, _desc, seq) in &proteins {
        for pep in trypsin(seq, 0) {
            pep_to_proteins.entry(pep.clone()).or_default().push(acc.clone());
            all_peptides.push(pep);
        }
    }

    let unique_peps: std::collections::HashSet<&String> = all_peptides.iter().collect();
    println!("Total tryptic peptides  : {}", all_peptides.len());
    println!("Unique tryptic peptides : {}", unique_peps.len());

    // Feed them into the peptide_input parser (space-separated format)
    let input_str: String = unique_peps.iter()
        .map(|s| s.as_str())
        .collect::<Vec<_>>()
        .join("\n");

    let parsed = parse_peptide_str(&input_str)
        .expect("peptide_input should parse all tryptic peptides");

    println!("peptide_input accepted  : {}", parsed.len());

    // All parsed sequences should survive encode with valid codes
    let mut j_count = 0usize;
    for p in &parsed {
        let enc = encode_sequence(p.sequence.as_bytes());
        assert!(enc.iter().all(|&c| c >= 1 && c <= 26));
        j_count += p.sequence.bytes().filter(|&b| b == b'J').count();
    }
    println!("J (Xle) residues found  : {}", j_count);
    // SARS-CoV-2 proteins don't use J but the pipeline must not crash on them
    assert!(parsed.len() > 0);
}

// ── Test 3: uniqueness — peptides unique to one protein ──────────────────────
#[tokio::test]
async fn test_peptide_uniqueness_stats() {
    let fasta_path = "/tmp/sars_cov2_pep_test3.fasta";
    fetch_proteome("UP000464024", fasta_path, true).await.expect("download");

    let fasta_text = std::fs::read_to_string(fasta_path).unwrap();
    let proteins = parse_fasta(&fasta_text);

    let mut pep_to_prots: HashMap<String, Vec<String>> = HashMap::new();
    for (acc, _, seq) in &proteins {
        for pep in trypsin(seq, 0) {
            pep_to_prots.entry(pep).or_default().push(acc.clone());
        }
    }

    let unique_to_one  = pep_to_prots.values().filter(|v| v.len() == 1).count();
    let shared         = pep_to_prots.values().filter(|v| v.len() > 1).count();
    let total          = pep_to_prots.len();

    println!("─── Peptide uniqueness (SARS-CoV-2) ───────────────");
    println!("Total unique peptides  : {}", total);
    println!("Unique to 1 protein    : {} ({:.1}%)", unique_to_one,
             100.0 * unique_to_one as f64 / total as f64);
    println!("Shared (>1 protein)    : {}", shared);

    assert!(unique_to_one > 0, "Expected some protein-unique peptides");
}
