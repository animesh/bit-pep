/// uniprot.rs — UniProt proteome FASTA downloader
///
/// Replaces ncbi.rs.  Uses the UniProt REST API v2:
///   https://rest.uniprot.org
///
/// Two workflows:
///   1. fetch_proteome(proteome_id, output_path)
///      Downloads all canonical sequences for a UniProt Proteome ID
///      e.g. "UP000005640" (Homo sapiens), "UP000000589" (Mus musculus).
///
///   2. search_proteome(organism_name)
///      Returns a list of matching ProteomeEntry records so the user
///      can pick the right proteome ID before downloading.
///
/// Equivalent CLI commands in new main.rs:
///   bit-pop search-proteome --organism "Homo sapiens"
///   bit-pop fetch-proteome  --proteome UP000005640 -o human.fasta
///
/// Rate-limiting: UniProt asks for a User-Agent and recommends
/// streaming large downloads — both are handled here.

use anyhow::{Context, Result};
use serde::Deserialize;
use std::fs;
use std::io::Write;

const UNIPROT_REST: &str = "https://rest.uniprot.org";
const USER_AGENT: &str = "bit-pop-prot/0.1 (https://github.com/animesh/bit-pop; proteomics fork)";

// ── Data types ────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ProteomeEntry {
    pub upid: String,                        // e.g. "UP000005640"
    pub taxonomy_lineage: Option<Vec<TaxNode>>,
    pub proteome_type: Option<String>,       // "Reference Proteome", "Other proteome", …
    pub protein_count: Option<u64>,
    pub organism: Option<OrganismInfo>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct TaxNode {
    #[serde(rename = "taxonId")]
    pub taxon_id: u64,
    #[serde(rename = "scientificName")]
    pub scientific_name: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct OrganismInfo {
    pub scientific_name: Option<String>,
    pub common_name: Option<String>,
    pub taxon_id: Option<u64>,
}

#[derive(Debug, Deserialize)]
struct ProteomeSearchResponse {
    results: Vec<ProteomeEntry>,
}

// ── Public API ────────────────────────────────────────────────────────────────

/// Search UniProt proteomes by organism name.
/// Returns up to `limit` matches (default 10).
///
/// Example:
///   let hits = search_proteome("Homo sapiens", 5).await?;
///   for h in &hits { println!("{} — {} ({:?} proteins)", h.upid, …) }
pub async fn search_proteome(organism: &str, limit: usize) -> Result<Vec<ProteomeEntry>> {
    let url = format!(
        "{}/proteomes/search?query=organism_name:{}&format=json&size={}",
        UNIPROT_REST,
        urlencoding::encode(organism),
        limit
    );

    let client = build_client()?;
    let resp = client
        .get(&url)
        .send()
        .await
        .with_context(|| format!("UniProt search request failed for '{}'", organism))?;

    let status = resp.status();
    if !status.is_success() {
        anyhow::bail!("UniProt search returned HTTP {}: {}", status, url);
    }

    let body: ProteomeSearchResponse = resp
        .json()
        .await
        .context("Failed to parse UniProt proteome search JSON")?;

    Ok(body.results)
}

/// Download a proteome FASTA by UniProt Proteome ID.
///
/// Streams the response directly to `output_path` so large proteomes
/// (human: ~20k proteins, ~11 MB compressed) don't blow up RAM.
///
/// Set `canonical_only = true` (recommended) to exclude isoforms.
/// Set `canonical_only = false` to include all isoforms (+include_isoform).
pub async fn fetch_proteome(
    proteome_id: &str,
    output_path: &str,
    canonical_only: bool,
) -> Result<u64> {
    // UniProt REST streaming FASTA endpoint
    // Docs: https://www.uniprot.org/help/api_queries
    let isoform_param = if canonical_only {
        ""
    } else {
        "&includeIsoform=true"
    };

    let url = format!(
        "{}/uniprotkb/stream?format=fasta&query=proteome:{}&compressed=false{}",
        UNIPROT_REST, proteome_id, isoform_param
    );

    eprintln!("INFO: downloading proteome {} from UniProt…", proteome_id);
    eprintln!("INFO: URL: {}", url);

    let client = build_client()?;
    let resp = client
        .get(&url)
        .send()
        .await
        .with_context(|| format!("UniProt fetch failed for proteome '{}'", proteome_id))?;

    let status = resp.status();
    if !status.is_success() {
        anyhow::bail!(
            "UniProt returned HTTP {} for proteome '{}'. \
             Check the Proteome ID at https://www.uniprot.org/proteomes/{}",
            status, proteome_id, proteome_id
        );
    }

    let mut file = fs::File::create(output_path)
        .with_context(|| format!("Cannot create output file '{}'", output_path))?;

    // Stream body in chunks to avoid loading entire file into RAM
    let mut bytes_written: u64 = 0;
    let mut resp = resp;
    while let Some(chunk) = resp.chunk().await? {
        file.write_all(&chunk)?;
        bytes_written += chunk.len() as u64;
    }

    eprintln!(
        "INFO: wrote {:.1} MB to '{}'",
        bytes_written as f64 / 1_048_576.0,
        output_path
    );
    Ok(bytes_written)
}

/// Convenience: look up a common organism name and immediately download
/// the first Reference Proteome found.  Returns the UPID used.
pub async fn fetch_by_organism_name(organism: &str, output_path: &str) -> Result<String> {
    let hits = search_proteome(organism, 20).await?;
    if hits.is_empty() {
        anyhow::bail!("No UniProt proteomes found for organism: '{}'", organism);
    }

    // Prefer "Reference Proteome" if available
    let best = hits
        .iter()
        .find(|h| {
            h.proteome_type
                .as_deref()
                .map(|t| t.contains("Reference"))
                .unwrap_or(false)
        })
        .or_else(|| hits.first())
        .unwrap(); // safe: hits non-empty

    let upid = best.upid.clone();
    let name = best
        .organism
        .as_ref()
        .and_then(|o| o.scientific_name.as_deref())
        .unwrap_or("unknown");

    eprintln!("INFO: selected proteome {} for '{}' ({})", upid, organism, name);
    fetch_proteome(&upid, output_path, true).await?;
    Ok(upid)
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn build_client() -> Result<reqwest::Client> {
    reqwest::Client::builder()
        .user_agent(USER_AGENT)
        .timeout(std::time::Duration::from_secs(300)) // large proteomes can be slow
        .build()
        .context("Failed to build HTTP client")
}

// ── Pretty print ─────────────────────────────────────────────────────────────

pub fn display_proteome_hits(hits: &[ProteomeEntry]) {
    if hits.is_empty() {
        println!("No proteomes found.");
        return;
    }
    println!("{:<15} {:<30} {:>10}  {}", "UPID", "Organism", "Proteins", "Type");
    println!("{}", "-".repeat(75));
    for h in hits {
        let org = h
            .organism
            .as_ref()
            .and_then(|o| o.scientific_name.as_deref())
            .unwrap_or("?");
        let ptype = h.proteome_type.as_deref().unwrap_or("?");
        let nprots = h.protein_count.unwrap_or(0);
        println!("{:<15} {:<30} {:>10}  {}", h.upid, org, nprots, ptype);
    }
}
