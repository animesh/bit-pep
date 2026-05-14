# bit-pop → bit-pop-prot: Migration Guide

## New files (drop into src/)

| File               | Replaces     | Purpose                          |
|--------------------|-------------|----------------------------------|
| `amino.rs`         | DNA encoding in `lib.rs` | 5-bit AA encoding, XOR packing |
| `peptide_input.rs` | `fastq.rs`  | Parse peptide lists / FASTA      |
| `uniprot.rs`       | `ncbi.rs`   | UniProt REST API download        |
| `pep_output.rs`    | `sam.rs`    | TSV output (not SAM)             |

---

## Cargo.toml — add/change deps

```toml
[dependencies]
# keep existing: anyhow, rayon, indicatif, clap, serde, libsais-sys, memmap2

# add:
reqwest    = { version = "0.12", features = ["json", "stream"] }
urlencoding = "2.1"
tokio      = { version = "1", features = ["full"] }   # if not already async

# remove or keep unused:
# bio = "..."   # not needed; we roll our own AA encoding
```

---

## src/lib.rs — changes needed

### 1. Replace DNA alphabet constant
```rust
// REMOVE:
pub const ALPHA_SIZE: usize = 5; // ACGT + $
pub fn encode_dna(c: u8) -> u8 { ... }

// ADD:
pub use crate::amino::{ALPHA_SIZE, encode_aa as encode_char, encode_sequence};
```

### 2. Remove reverse-complement logic
```rust
// DELETE the entire rev_comp() function and all call sites.
// Peptides have no strand; this is dead code.
```

### 3. Change default k-mer size
```rust
// DNA default was k=10; for 20-AA alphabet k=5 is roughly equivalent
// discriminative power: 20^5 = 3.2M vs 4^10 = 1M
pub const DEFAULT_K: usize = 5;
```

### 4. Add new modules
```rust
pub mod amino;
pub mod peptide_input;
pub mod uniprot;
pub mod pep_output;
```

---

## src/fm.rs — changes needed

### 1. Alphabet size
```rust
// Change from:
const ALPHA: usize = 5;
// To (pulls from amino.rs):
use crate::amino::ALPHA_SIZE;
const ALPHA: usize = ALPHA_SIZE; // = 22
```

### 2. Sentinel character
The FM-index sentinel '$' should map to code 0.
`encode_aa` already reserves 0 for this — append a `0u8` to the
concatenated sequence before suffix array construction, same as DNA.

### 3. BWT/SA construction
`libsais` works on arbitrary byte alphabets.  The encoded bytes
(range 0–21) are already in the valid range.  **No change needed**
to the SA-IS call itself.

### 4. k-mer lookup
Replace `4usize.pow(k)` with `20usize.pow(k)` in any precomputed
lookup table sizing.  With k=5 that's 3.2M entries — still RAM-safe.

---

## src/align.rs — changes needed

### 1. XOR alignment
The 2-bit DNA XOR trick packed 16 bases per 32-bit word.
5-bit AA encoding packs **6 residues** per 32-bit word (see `amino::pack_u32`).

```rust
// Change chunk size constant:
const CHUNK_SIZE: usize = 6; // was 16 for DNA

// Change XOR mismatch counter to use amino::count_mismatches_u32()
```

### 2. Smith-Waterman scoring matrix
Optionally replace +1/−1 identity scoring with BLOSUM62 for
biologically meaningful approximate matching.

```rust
// Minimal change: keep identity scoring (exact or 1 mismatch)
// Better: load BLOSUM62 as a 20×20 i8 matrix at startup
// This is a separate feature — not required for the port to compile.
```

### 3. Remove quality score handling
Phred-scaled quality penalties are meaningless for peptides.
Delete quality-score branches in SW refinement.

### 4. Myers edit distance
No change needed — Myers operates on raw bytes and works on any alphabet.
This is the preferred aligner for peptides due to short sequence length.

---

## src/fasta.rs — changes needed

Protein FASTA uses the same `>header\nsequence\n` format as genome FASTA.
The main change:

```rust
// REMOVE: reverse-complement sequence generation
// CHANGE: call encode_sequence() from amino.rs instead of encode_dna()
// ADD: parse UniProt FASTA header to extract accession and protein name
//
// UniProt FASTA header format:
//   >sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens ...
//   >tr|A0A023GPI8|... (TrEMBL entry)
//
// Parser:
pub fn parse_uniprot_header(header: &str) -> (String, String) {
    // Returns (accession, protein_name)
    let parts: Vec<&str> = header.splitn(3, '|').collect();
    if parts.len() == 3 {
        let acc  = parts[1].to_string();
        let name = parts[2].split(" OS=").next().unwrap_or(parts[2]).to_string();
        (acc, name)
    } else {
        (header.to_string(), String::new())
    }
}
```

---

## src/main.rs — CLI changes

### Remove commands
- `bit-pop fetch` (NCBI fetch) → replaced by `fetch-proteome`
- `bit-pop search` (NCBI search) → replaced by `search-proteome`

### Rename/change commands
```
bit-pop run genome.fna reads.fastq
→ bit-pop run proteome.fasta peptides.txt [--proteome UP000005640]

bit-pop search --organism "Escherichia coli"
→ bit-pop search-proteome --organism "Homo sapiens"

bit-pop fetch --accession NC_000913.3
→ bit-pop fetch-proteome --proteome UP000005640 -o human.fasta
         OR
   bit-pop fetch-proteome --organism "Homo sapiens" -o human.fasta
           (auto-selects Reference Proteome)
```

### New `run` flags
```
--min-score 0.8          (default higher than DNA; peptides are shorter)
--k 5                    (default k for protein search)
--allow-mismatches N     (per-peptide mismatch budget, default 0 = exact)
--output-format tsv      (only option; SAM removed)
--canonical-only         (exclude isoforms from UniProt download)
```

### Remove flags
```
--read-type short/long   (not applicable)
--spaced-seed            (less relevant; remove or keep as no-op)
--mmap                   (keep — still useful for large proteomes)
-1/-2 (paired-end)       (remove entirely)
```

---

## src/rank.rs — minimal change

Multi-proteome ranking logic is **directly reusable** — it already
assigns reads to the best-matching reference.  For peptides:
- "genome" → "proteome" in variable names (cosmetic)
- Scoring formula stays the same (alignment score 85% + k-mer rarity 15%)
- EM post-processing stays as-is; it models proteome abundance just
  like genome abundance

---

## src/em.rs — no change needed

EM is abundance-model agnostic.  Protein abundance estimation
(instead of genome abundance) is mathematically identical.
The "which proteome does this peptide come from" question is the
same as "which genome does this read come from".

---

## What this port does NOT handle (known limitations)

1. **PTMs** — `PEPTM[ox]IDER` and `PEPTMIDER` are treated as different
   sequences.  No modification-aware matching.  This is a future feature.

2. **Isoforms** — downloading canonical-only (`--canonical-only`, default)
   means peptides unique to isoforms will not map.  Pass `--no-canonical-only`
   to include isoforms (much larger FASTA).

3. **Missed cleavages** — not applicable here; we search the submitted
   sequence as-is, not tryptic peptides derived in-tool.

4. **False discovery rate** — bit-pop has no decoy database / FDR model.
   For proteomics-grade results, use this for candidate generation,
   then validate hits through MaxQuant or MSFragger.
