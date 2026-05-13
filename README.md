# Bit-Pep: Multi-Proteome Peptide Classification

Adapting [bit-pop](https://github.com/mladenpop-oss/bit-pop/) to map Peptide(s) and Proteome(s)

## Plan @claude-code

# bit-pop → pepmap: Peptide Search Across Proteomes

Adapts the bit-pop FM-index + multi-reference classifier for peptide-to-protein mapping.
The FM-index in `fm.rs` is **unchanged** — it operates on `&[u8]`, so encoding amino acids
as bytes (0–20) makes it alphabet-agnostic. The XOR/SW/Myers alignment pipeline collapses
to a single exact FM backward search per peptide (O(m), m = peptide length).

---

## Files to ADD

### `src/aa.rs` — Amino acid alphabet (replaces 2-bit DNA encoding in lib.rs)

```rust
/// Encode a single amino acid character to a byte value 1–20.
/// Unknown/ambiguous residues are mapped to nearest canonical AA.
/// Separator '$' maps to 0.
pub fn encode_aa(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 1,  b'C' => 2,  b'D' => 3,  b'E' => 4,
        b'F' => 5,  b'G' => 6,  b'H' => 7,  b'I' => 8,
        b'K' => 9,  b'L' => 10, b'M' => 11, b'N' => 12,
        b'P' => 13, b'Q' => 14, b'R' => 15, b'S' => 16,
        b'T' => 17, b'V' => 18, b'W' => 19, b'Y' => 20,
        b'U' => 11, // selenocysteine → Met (common in UniProt)
        b'B' => 12, // Asx ambiguity → Asn
        b'Z' => 14, // Glx ambiguity → Gln
        b'X' => 1,  // unknown → Ala
        b'$' => 0,  // separator
        _    => 0,
    }
}

pub fn decode_aa(v: u8) -> char {
    match v {
        1  => 'A', 2  => 'C', 3  => 'D', 4  => 'E',
        5  => 'F', 6  => 'G', 7  => 'H', 8  => 'I',
        9  => 'K', 10 => 'L', 11 => 'M', 12 => 'N',
        13 => 'P', 14 => 'Q', 15 => 'R', 16 => 'S',
        17 => 'T', 18 => 'V', 19 => 'W', 20 => 'Y',
        _  => '$',
    }
}

/// Encode a full protein/peptide string to a Vec<u8> of AA indices.
pub fn encode_sequence(seq: &str) -> Vec<u8> {
    seq.bytes().map(encode_aa).collect()
}

/// Validate that a string looks like a peptide (all standard AA chars).
pub fn is_valid_peptide(s: &str) -> bool {
    !s.is_empty()
        && s.len() >= 5       // discard fragments shorter than 5 aa
        && s.bytes().all(|c| {
            matches!(c.to_ascii_uppercase(),
                b'A'..=b'Z') // broad check; encode_aa handles ambiguous ones
        })
}
```

---

### `src/peptide.rs` — Peptide file parser

```rust
use std::fs;
use crate::aa::is_valid_peptide;

/// Parse a peptide input file.
/// Accepts:
///   - One peptide per line
///   - Multiple peptides on a line separated by whitespace
///   - Lines starting with '#' are comments
///   - Empty lines are skipped
///
/// Returns a deduplicated, uppercased Vec of valid peptide strings.
pub fn parse_peptide_file(path: &str) -> anyhow::Result<Vec<String>> {
    let content = fs::read_to_string(path)
        .map_err(|e| anyhow::anyhow!("Cannot read peptide file {}: {}", path, e))?;

    let mut peptides: Vec<String> = content
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.trim_start().starts_with('#'))
        .flat_map(|line| line.split_whitespace())
        .map(|p| p.to_ascii_uppercase())
        .filter(|p| is_valid_peptide(p))
        .collect();

    // Deduplicate while preserving order
    let mut seen = std::collections::HashSet::new();
    peptides.retain(|p| seen.insert(p.clone()));

    eprintln!("[peptide] Loaded {} unique peptides from {}", peptides.len(), path);
    Ok(peptides)
}
```

---

### `src/uniprot.rs` — UniProt proteome downloader (replaces ncbi.rs)

```rust
use std::fs;
use std::path::Path;
use anyhow::Result;

const UNIPROT_REST: &str = "https://rest.uniprot.org/uniprotkb/search";

/// Download a UniProt reference proteome as FASTA.
///
/// `proteome_id` — UniProt proteome ID, e.g. "UP000005640" (human)
///               — or a taxon name, e.g. "Homo sapiens"
/// `out_path`    — local path to write the FASTA file
///
/// Uses UniProt REST API with cursor-based pagination (500 entries/page).
pub fn download_proteome(proteome_id: &str, out_path: &str) -> Result<()> {
    // Detect whether we got a proteome ID (UP*) or an organism name
    let query = if proteome_id.starts_with("UP") {
        format!("proteome:{}", proteome_id)
    } else {
        format!("proteome:* AND organism_name:{}", proteome_id)
    };

    eprintln!("[uniprot] Downloading proteome: {} → {}", proteome_id, out_path);

    let client = reqwest::blocking::Client::builder()
        .user_agent("pepmap/0.1 (proteomics-tool; contact via GitHub)")
        .build()?;

    let mut all_fasta = String::new();
    let mut cursor: Option<String> = None;
    let mut page = 0usize;

    loop {
        page += 1;
        let mut req = client.get(UNIPROT_REST)
            .query(&[
                ("query",  query.as_str()),
                ("format", "fasta"),
                ("size",   "500"),
            ]);

        if let Some(ref c) = cursor {
            req = req.query(&[("cursor", c.as_str())]);
        }

        let resp = req.send()?;

        // Extract next cursor from Link header (UniProt pagination)
        cursor = resp.headers()
            .get("link")
            .and_then(|v| v.to_str().ok())
            .and_then(|s| parse_next_cursor(s));

        let body = resp.text()?;
        if body.trim().is_empty() { break; }

        all_fasta.push_str(&body);
        eprintln!("[uniprot] Page {} fetched ({} chars total)", page, all_fasta.len());

        if cursor.is_none() { break; }
    }

    fs::write(out_path, &all_fasta)?;
    eprintln!("[uniprot] Written {} bytes to {}", all_fasta.len(), out_path);
    Ok(())
}

/// Parse 'cursor=XYZ' from a Link: <url?cursor=XYZ>; rel="next" header.
fn parse_next_cursor(link_header: &str) -> Option<String> {
    if !link_header.contains("rel=\"next\"") { return None; }
    link_header.split(',')
        .find(|part| part.contains("rel=\"next\""))
        .and_then(|part| {
            let url = part.split('<').nth(1)?.split('>').next()?;
            url.split('&')
               .chain(url.split('?'))
               .find(|seg| seg.starts_with("cursor="))
               .map(|seg| seg["cursor=".len()..].to_owned())
        })
}

/// List available reference proteomes for an organism name.
pub fn search_proteomes(organism: &str) -> Result<Vec<(String, String)>> {
    let client = reqwest::blocking::Client::builder()
        .user_agent("pepmap/0.1")
        .build()?;

    let resp = client.get("https://rest.uniprot.org/proteomes/search")
        .query(&[
            ("query",  organism),
            ("format", "tsv"),
            ("fields", "upid,organism,protein_count,busco"),
            ("size",   "10"),
        ])
        .send()?
        .text()?;

    let results: Vec<(String, String)> = resp.lines()
        .skip(1)  // header row
        .filter_map(|line| {
            let mut cols = line.splitn(2, '\t');
            let upid = cols.next()?.to_owned();
            let org  = cols.next()?.to_owned();
            Some((upid, org))
        })
        .collect();

    Ok(results)
}
```

---

### `src/proteome.rs` — Protein FASTA loader with UniProt header parsing

```rust
use std::fs;
use crate::aa::encode_sequence;

/// A single protein entry from a UniProt FASTA.
#[derive(Debug, Clone)]
pub struct ProteinEntry {
    pub accession:  String,   // e.g. "P12345"
    pub entry_name: String,   // e.g. "GENE_HUMAN"
    pub gene:       String,   // e.g. "GENE"
    pub organism:   String,   // e.g. "Homo sapiens"
    pub description: String,  // full description string
    pub sequence:   String,   // raw AA sequence (uppercase)
    pub offset:     usize,    // byte offset in concatenated index text
}

/// Parsed proteome: all proteins + the concatenated encoded sequence
/// ready for FM-index construction.
pub struct Proteome {
    pub name:       String,
    pub proteins:   Vec<ProteinEntry>,
    /// Concatenated encoded sequence: proteins joined by '$' (0).
    /// This is the text handed to the FM-index builder.
    pub text:       Vec<u8>,
}

impl Proteome {
    pub fn from_fasta(path: &str, proteome_name: &str) -> anyhow::Result<Self> {
        let content = fs::read_to_string(path)
            .map_err(|e| anyhow::anyhow!("Cannot read FASTA {}: {}", path, e))?;

        let mut proteins: Vec<ProteinEntry> = Vec::new();
        let mut text: Vec<u8> = Vec::new();
        let mut current_header: Option<String> = None;
        let mut current_seq: String = String::new();

        for line in content.lines() {
            if line.starts_with('>') {
                if let Some(hdr) = current_header.take() {
                    let offset = text.len();
                    let entry = parse_uniprot_header(&hdr, &current_seq, offset);
                    let encoded = encode_sequence(&current_seq);
                    text.extend_from_slice(&encoded);
                    text.push(0u8); // '$' separator
                    proteins.push(entry);
                    current_seq.clear();
                }
                current_header = Some(line[1..].to_owned());
            } else {
                current_seq.push_str(line.trim());
            }
        }
        // Flush last entry
        if let Some(hdr) = current_header {
            let offset = text.len();
            let entry = parse_uniprot_header(&hdr, &current_seq, offset);
            let encoded = encode_sequence(&current_seq);
            text.extend_from_slice(&encoded);
            text.push(0u8);
            proteins.push(entry);
        }

        eprintln!("[proteome] {} — {} proteins, {} encoded chars",
            proteome_name, proteins.len(), text.len());

        Ok(Proteome {
            name: proteome_name.to_owned(),
            proteins,
            text,
        })
    }

    /// Given a position in the concatenated text, find which protein it belongs to
    /// and the offset within that protein.
    pub fn locate(&self, pos: usize) -> Option<(&ProteinEntry, usize)> {
        // Binary search: find the protein whose offset <= pos < offset+len+1
        let idx = self.proteins.partition_point(|p| p.offset <= pos);
        if idx == 0 { return None; }
        let protein = &self.proteins[idx - 1];
        let local_pos = pos - protein.offset;
        if local_pos < protein.sequence.len() {
            Some((protein, local_pos))
        } else {
            None // hit the separator
        }
    }
}

/// Parse a UniProt FASTA header line (without the leading '>').
///
/// Formats handled:
///   sp|P12345|GENE_HUMAN Description OS=Homo sapiens OX=9606 GN=GENE PE=1 SV=1
///   tr|A0A000|GENE_HUMAN ...
///   Any non-UniProt FASTA (e.g. contaminant db): accession = first token
fn parse_uniprot_header(header: &str, seq: &str, offset: usize) -> ProteinEntry {
    let mut accession  = String::new();
    let mut entry_name = String::new();
    let mut gene       = String::new();
    let mut organism   = String::new();
    let mut description = header.to_owned();

    // Try UniProt pipe format
    let parts: Vec<&str> = header.splitn(3, '|').collect();
    if parts.len() == 3 && (parts[0] == "sp" || parts[0] == "tr") {
        accession  = parts[1].to_owned();
        // third field: "GENE_HUMAN Description OS=..."
        let rest = parts[2];
        let space_pos = rest.find(' ').unwrap_or(rest.len());
        entry_name = rest[..space_pos].to_owned();
        description = rest[space_pos..].trim().to_owned();

        // Extract OS= field
        if let Some(os_start) = description.find("OS=") {
            let os_end = description[os_start..]
                .find(" OX=").map(|i| os_start + i)
                .unwrap_or(description.len());
            organism = description[os_start + 3..os_end].to_owned();
        }

        // Extract GN= field
        if let Some(gn_start) = description.find("GN=") {
            let gn_end = description[gn_start..]
                .find(' ').map(|i| gn_start + i)
                .unwrap_or(description.len());
            gene = description[gn_start + 3..gn_end].to_owned();
        }
    } else {
        // Generic FASTA: first whitespace-delimited token is the accession
        accession = header.split_whitespace().next().unwrap_or("unknown").to_owned();
    }

    ProteinEntry {
        accession,
        entry_name,
        gene,
        organism,
        description,
        sequence: seq.to_ascii_uppercase(),
        offset,
    }
}
```

---

### `src/search.rs` — Peptide FM-index search (exact + fuzzy)

```rust
use crate::aa::encode_sequence;
use crate::proteome::Proteome;
// Reuse the existing FmIndex from fm.rs
use crate::fm::FmIndex;

#[derive(Debug)]
pub struct PeptideHit {
    pub peptide:    String,
    pub accession:  String,
    pub entry_name: String,
    pub gene:       String,
    pub organism:   String,
    pub proteome:   String,
    pub position:   usize,   // 0-based position in protein
    pub match_type: MatchType,
}

#[derive(Debug)]
pub enum MatchType {
    Exact,
    OneMismatch, // future: fuzzy
}

/// Search all peptides against all loaded proteomes.
/// Returns one row per (peptide, hit).
pub fn search_all(
    peptides: &[String],
    proteomes: &[(Proteome, FmIndex)],
    max_mismatches: u8,
) -> Vec<PeptideHit> {
    use rayon::prelude::*;   // reuse existing rayon dependency

    peptides.par_iter().flat_map(|pep| {
        let encoded = encode_sequence(pep);
        let mut hits: Vec<PeptideHit> = Vec::new();

        for (proteome, fm) in proteomes {
            // FM backward search returns a range [lo, hi) of suffix array positions
            if let Some((lo, hi)) = fm.backward_search(&encoded) {
                for sa_pos in lo..hi {
                    let text_pos = fm.sa[sa_pos] as usize;
                    if let Some((protein, local_pos)) = proteome.locate(text_pos) {
                        hits.push(PeptideHit {
                            peptide:    pep.clone(),
                            accession:  protein.accession.clone(),
                            entry_name: protein.entry_name.clone(),
                            gene:       protein.gene.clone(),
                            organism:   protein.organism.clone(),
                            proteome:   proteome.name.clone(),
                            position:   local_pos,
                            match_type: MatchType::Exact,
                        });
                    }
                }
            }
        }
        hits
    }).collect()
}

/// Write results as TSV to stdout or a file.
pub fn write_tsv(hits: &[PeptideHit], out_path: Option<&str>) -> anyhow::Result<()> {
    use std::io::Write;

    let header = "peptide\taccession\tentry\tgene\torganism\tproteome\tposition\tmatch_type\n";
    let mut out: Box<dyn Write> = match out_path {
        Some(p) => Box::new(std::fs::File::create(p)?),
        None    => Box::new(std::io::stdout()),
    };

    out.write_all(header.as_bytes())?;
    for h in hits {
        writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}",
            h.peptide, h.accession, h.entry_name, h.gene,
            h.organism, h.proteome, h.position, h.match_type)?;
    }
    Ok(())
}
```

---

## Files to MODIFY

### `Cargo.toml` — add dependencies

```toml
[dependencies]
# existing deps unchanged...
anyhow   = "1"
reqwest  = { version = "0.11", features = ["blocking"] }
# rayon already present
```

---

### `src/main.rs` — add `peptide-search` subcommand

Add this arm to the existing `match` on subcommands. All existing subcommands are **unchanged**.

```rust
// NEW subcommand: pepmap
("peptide-search", Some(sub)) => {
    let peptide_file = sub.value_of("peptides").unwrap();
    let fasta_paths: Vec<&str> = sub.values_of("proteome")
                                    .unwrap().collect();
    let out_path = sub.value_of("output");
    let threads  = sub.value_of("threads")
                      .and_then(|t| t.parse().ok())
                      .unwrap_or(4usize);

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()?;

    // Optionally download from UniProt
    if let Some(upid) = sub.value_of("uniprot") {
        let local = format!("{}.fasta", upid);
        uniprot::download_proteome(upid, &local)?;
        // fasta_paths would then point to `local`
    }

    // Load all proteomes and build FM-indexes
    let proteomes: Vec<(proteome::Proteome, fm::FmIndex)> = fasta_paths
        .iter()
        .map(|path| {
            let name = std::path::Path::new(path)
                .file_stem().unwrap().to_string_lossy().to_string();
            let prot = proteome::Proteome::from_fasta(path, &name)
                .expect("Failed to load proteome");
            let fm = fm::FmIndex::build(&prot.text); // existing build fn, unchanged
            (prot, fm)
        })
        .collect();

    // Parse peptides
    let peptides = peptide::parse_peptide_file(peptide_file)?;

    // Search
    let hits = search::search_all(&peptides, &proteomes, 0);
    eprintln!("[pepmap] {} hits across {} peptides", hits.len(), peptides.len());

    search::write_tsv(&hits, out_path)?;
}
```

Add the subcommand definition in the clap `App`:

```rust
.subcommand(
    App::new("peptide-search")
        .about("Map peptides to one or more UniProt proteomes")
        .arg(Arg::new("peptides")
            .short('p').long("peptides")
            .value_name("FILE")
            .help("Peptide list: one per line, or space-separated per line")
            .required(true))
        .arg(Arg::new("proteome")
            .short('f').long("proteome")
            .value_name("FASTA")
            .help("Proteome FASTA file(s) (UniProt format)")
            .multiple_occurrences(true)
            .required_unless_present("uniprot"))
        .arg(Arg::new("uniprot")
            .long("uniprot")
            .value_name("PROTEOME_ID")
            .help("UniProt proteome ID to download, e.g. UP000005640"))
        .arg(Arg::new("output")
            .short('o').long("output")
            .value_name("TSV")
            .help("Output TSV file (default: stdout)"))
        .arg(Arg::new("threads")
            .short('t').long("threads")
            .default_value("4"))
)
```

---

## What does NOT change

| File | Status | Reason |
|---|---|---|
| `src/fm.rs` | **Unchanged** | Operates on `&[u8]` — alphabet-agnostic |
| `src/em.rs` | **Unchanged** | EM on abundance vectors — still valid for multi-proteome |
| `src/align.rs` | **Unused** for exact matching; still available for fuzzy | |
| `src/sam.rs` | **Unused** for peptide output | New TSV output in `search.rs` |
| `src/ncbi.rs` | **Unchanged** | DNA classification workflow unaffected |
| `src/fasta.rs` | **Unchanged** | DNA FASTA still needed for `run/build/map` cmds |
| All benchmarks, tests | **Unchanged** | Existing test suite unaffected |

---

## Usage

```bash
# Build
cargo build --release

# Download human proteome and search
./target/release/bit-pop peptide-search \
  --uniprot UP000005640 \
  -p my_peptides.txt \
  -o hits.tsv \
  -t 8

# Search against pre-downloaded FASTAs (human + mouse + contaminants)
./target/release/bit-pop peptide-search \
  -f human.fasta \
  -f mouse.fasta \
  -f contaminants.fasta \
  -p my_peptides.txt \
  -o hits.tsv \
  -t 8
```

### Peptide file formats accepted

```
# Format 1: one per line
PEPTIDER
ACDEFGHIKLM
SAMPLEPEPTIDE

# Format 2: space-separated
PEPTIDER ACDEFGHIKLM
SAMPLEPEPTIDE ANOTHERSEQ

# Format 3: mixed, with comments
# MaxQuant output peptides
PEPTIDER ACDEFGHIKLM
# tryptic fragments
SAMPLEPEPTIDE
```

### Output TSV columns

| Column | Example |
|---|---|
| peptide | PEPTIDER |
| accession | P12345 |
| entry | GENE_HUMAN |
| gene | GENE |
| organism | Homo sapiens |
| proteome | human |
| position | 42 |
| match_type | Exact |

---

## Key design decisions

1. **FM-index reused as-is** — encoding 20 AA as bytes 1–20 fits inside `u8`; libsais SA construction is byte-array agnostic.
2. **No alignment step for exact search** — peptides (from MaxQuant etc.) are exact sequences; FM backward search is sufficient and O(m).
3. **`align.rs` kept available** for future fuzzy matching (e.g. I/L equivalence, deamidation N→D).
4. **EM post-processing** still applies if a peptide maps to multiple proteomes with ambiguity.
5. **Reverse complement logic deleted** only from the new peptide path — DNA workflow is untouched.


## Availability

Source code available under the MIT License.

## Citation

```bibtex
@software{popovic_2026_bitpop,
  author = {Popovi{\'c}, Mladen},
  title = {Bit-Pop: A Proof-of-Concept Tool for Multi-Genome DNA Read Classification},
  year = {2026},
  doi = {10.5281/zenodo.20043593},
  url = {https://github.com/mladenpop-oss/bit-pop}
}
```

Or plain text:

> Popović, M. (2026). Bit-Pop: A Proof-of-Concept Tool for Multi-Genome DNA Read Classification. https://doi.org/10.5281/zenodo.20043593

## License

MIT License
