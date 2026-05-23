# Bit-Pep: Multi-Proteome Peptide Classification

Adapting [bit-pop](https://github.com/mladenpop-oss/bit-pop/) to map Peptide(s) and Proteome(s) as [bit-pep](https://github.com/animesh/bit-pep/) 

## test with random peptides generated per protein from 10% of uniprot peptides ranging from 5 to 35
```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
python generate_random_peptides.py 
    Sampled 57455/574551 proteins (10.0%)
    Generated 57455 peptides from sampled proteins in uniprot_sprot.fasta
    Wrote peptide FASTA: uniprot_sprot_0.1_5_35_42.fasta
    Running: cargo run --bin bit-pep -- run-prot uniprot_sprot.fasta -p uniprot_sprot_0.1_5_35_42.fasta
        Finished `dev` profile [unoptimized + debuginfo] target(s) in 19.19s
        Running `target/debug/bit-pep run-prot uniprot_sprot.fasta -p uniprot_sprot_0.1_5_35_42.fasta`
    bit-pep RunProt  (FM-index parallel peptide->proteome search)
    ===============================================================

    [1/4] Resolving proteome...
    uniprot_sprot.fasta (274.2 MB)
    Loading: [00:00:19 ████████████████████████████████████████ 574627/574627]                                                               574627 proteins indexed in 117.9s  (12 threads)

    [3/4] Loading peptides...
    INFO: loaded 57455 peptides
    57455 peptides

    [4/4] Mapping (12 threads)...
    Mapping: [00:00:00 ████████████████████████████████████████ 57455/57455]                                                               
    ===============================================================
    Output  : uniprot_sprot_0.1_5_35_42.pep.tsv
    Mapping : 8.30s  |  Total: 126.57s

    ─────────────────────────────────────────
    Peptides submitted : 57455
    Mapped             : 57455 (100.0%)
    unique            : 27572
    shared (>1 prot.) : 29883
    cross-proteome    : 0
    Unmapped           : 0
    ─────────────────────────────────────────

    Peptide -> Protein -> Species  (top 20 by total peptides)
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    Organism                                      Unique    Shared     Total  %total  ProtsHit  TotalProt
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    Homo sapiens                                    1100      5217      6317    11.0     12986     20431
    Mus musculus                                     915      4926      5841    10.2     11099     17252
    Escherichia coli (strain K12)                    210      4504      4714     8.2      2700      4531
    Arabidopsis thaliana                            1344      2954      4298     7.5      9277     16418
    Rattus norvegicus                                346      3698      4044     7.0      5274      8226
    Escherichia coli O157:H7                          12      3933      3945     6.9      1489      2047
    Escherichia coli O6:H1 (strain CFT073 /...         6      3730      3736     6.5      1302      1702
    Shigella flexneri                                 10      3575      3585     6.2      1278      1696
    Salmonella typhimurium (strain LT2 / SG...        34      3299      3333     5.8      1336      1831
    Bos taurus                                       312      2951      3263     5.7      3658      6052
    Escherichia coli O9:H4 (strain HS)                 1      3174      3175     5.5       805       847
    Escherichia coli O139:H28 (strain E2437...         3      3169      3172     5.5       800       838
    Shigella sonnei (strain Ss046)                     1      3158      3159     5.5       824       879
    Escherichia coli (strain ATCC 8739 / DS...         2      3145      3147     5.5       788       822
    Escherichia coli O6:K15:H31 (strain 536...         2      3112      3114     5.4       834       894
    Salmonella typhi                                   6      3080      3086     5.4      1077      1349
    Shigella boydii serotype 4 (strain Sb227)          4      3040      3044     5.3       778       829
    Escherichia coli (strain SMS-3-5 / SECEC)          2      3013      3015     5.2       759       791
    Escherichia coli O8 (strain IAI1)                  0      2991      2991     5.2       738       761
    Escherichia coli (strain 55989 / EAEC)             0      2979      2979     5.2       739       762
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    ... and 12424 more organisms

    Unique/Shared/Total = peptide counts; %total = % of all submitted peptides
    ProtsHit = proteins with >=1 mapped peptide; TotalProt = proteins in database
    Wrote metadata CSV: uniprot_sprot_0.1_5_35_42.metadata.csv
    Wrote verified report: uniprot_sprot_0.1_5_35_42.verified.txt

cat uniprot_sprot_0.1_5_35_42.verified.txt
    Total peptides: 57455
    Verified peptides: 57455
    Unique peptides: 27548
    Shared peptides: 29907
    Missing/mismatched peptides: 0

awk -F ',' '{print $NF}' uniprot_sprot_0.1_5_35_42.metadata.csv | sort | uniq -c
  29907 shared
  27548 unique
      1 verified
```

### sweep random seeds from 1 to 100 with peptides ranging from 7-30 per protein from ~1/3 (0.33) of uniprot

```
for seed in $(seq 1 100); do
  python3 generate_random_peptides.py \
    --fasta uniprot_sprot.fasta \
    --min-len 7 \
    --max-len 30 \
    --sample-fraction 0.33 \
    --seed "$seed"
done
```

### peptides from MaxQuant HeLa DIA [results](https://fuzzylife.substack.com/p/proteomics-data-processing-with-maxquant) 

```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
wget "https://zenodo.org/records/14557756/files/proteinGroups.txt"
grep -vE "Peptide|CON__|REV__" proteinGroups.txt| awk -F '\t' '{print $92}'  | sed 's/;/\n/g'  > peptides.txt 
cargo run --bin bit-pep -- run-prot uniprot_sprot.fasta -p peptides.txt -j 12
        Finished `dev` profile [unoptimized + debuginfo] target(s) in 2.91s
        Running `target/debug/bit-pep run-prot uniprot_sprot.fasta -p peptides.txt -j 12`
    Bit-Pop RunProt  (FM-index parallel peptide->proteome search)
    ===============================================================

    [1/4] Resolving proteome...
    uniprot_sprot.fasta (274.2 MB)
    Loading: [00:00:19 ████████████████████████████████████████ 574627/574627]                                                                                                                                                        574627 proteins indexed in 112.6s  (12 threads)

    [3/4] Loading peptides...
    INFO: loaded 58920 peptides
    58920 peptides

    [4/4] Mapping (12 threads)...
    Mapping: [00:00:00 ████████████████████████████████████████ 58920/58920]                                                                                                                                                        
    ===============================================================
    Output  : peptides.pep.tsv
    Mapping : 1.83s  |  Total: 114.59s

    ─────────────────────────────────────────
    Peptides submitted : 58920
    Mapped             : 54626 (92.7%)
    unique            : 24989
    shared (>1 prot.) : 29637
    cross-proteome    : 0
    Unmapped           : 4294
    ─────────────────────────────────────────

    Peptide -> Protein -> Species  (top 20 by total peptides)
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    Organism                                      Unique    Shared     Total  %total  ProtsHit  TotalProt
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    Homo sapiens                                   24871     29560     54431    92.4      6852     20431
    Mus musculus                                      49     23184     23233    39.4      4656     17252
    Rattus norvegicus                                 18     13398     13416    22.8      2333      8226
    Bos taurus                                         4      9333      9337    15.8      1651      6052
    Pongo abelii                                       7      7350      7357    12.5       898      2218
    Gallus gallus                                      3      3594      3597     6.1       634      2314
    Macaca fascicularis                                1      3260      3261     5.5       417      1176
    Sus scrofa                                         0      3158      3158     5.4       439      1462
    Xenopus laevis                                     0      2679      2679     4.5       752      3514
    Oryctolagus cuniculus                              1      2606      2607     4.4       348       979
    Pan troglodytes                                    0      1897      1897     3.2       224       692
    Canis lupus familiaris                             0      1859      1859     3.2       251       857
    Danio rerio                                        1      1698      1699     2.9       502      3369
    Xenopus tropicalis                                 0      1323      1323     2.2       326      1713
    Drosophila melanogaster                            1      1322      1323     2.2       366      3868
    Cricetulus griseus                                 0      1317      1317     2.2       105       248
    Mesocricetus auratus                               0       957       957     1.6       102       277
    Caenorhabditis elegans                             0       873       873     1.5       294      4499
    Equus caballus                                     1       678       679     1.2        73       292
    Arabidopsis thaliana                               1       665       666     1.1       533     16418
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    ... and 2663 more organisms

    Unique/Shared/Total = peptide counts; %total = % of all submitted peptides
    ProtsHit = proteins with >=1 mapped peptide; TotalProt = proteins in database
```

### missing isofom/peptides[LAQPGFPSGGPGGTR,SPIAAARCR,...] in [uniprot_sprot.fasta](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) and deleted protein/peptides[ELQQVTAGEAASIH,QVCQVPASR,...] in [current human proteome](https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&includeIsoform=true&query=%28%28proteome%3AUP000005640%29%29) leading to `Unmapped           : 4294` verify by [including](https://zenodo.org/records/20344308) human proteome from [2024](https://github.com/user-attachments/files/16439166/mqpar.xml.txt)
```
wget "https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&includeIsoform=true&query=%28%28proteome%3AUP000005640%29%29" -O human.fasta
wget https://zenodo.org/records/20344308/files/uniprotkb_proteome_UP000005640_2024_04_18.fasta 
cat uniprotkb_proteome_UP000005640_2024_04_18.fasta >> human.fasta
wget "https://zenodo.org/records/14557756/files/proteinGroups.txt"
grep -vE "Peptide|CON__|REV__" proteinGroups.txt| awk -F '\t' '{print $92}'  | sed 's/;/\n/g'  > peptides.txt 
cargo run --bin bit-pep -- run-prot human.fasta  -p peptides.txt -j 12 
    Finished `dev` profile [unoptimized + debuginfo] target(s) in 3.39s
     Running `target/debug/bit-pep run-prot human.fasta -p peptides.txt -j 12`
    bit-pep RunProt  (FM-index parallel peptide->proteome search)
    ===============================================================

    [1/4] Resolving proteome...
    human.fasta (106.3 MB)
    Loading: [00:00:07 ████████████████████████████████████████ 210230/210230]                                                                                                                                                        210230 proteins indexed in 46.5s  (12 threads)

    [3/4] Loading peptides...
    INFO: loaded 58920 peptides
    58920 peptides

    [4/4] Mapping (12 threads)...
    Mapping: [00:00:00 ████████████████████████████████████████ 58920/58920]                                                                                                                                                        
    ===============================================================
    Output  : peptides.pep.tsv
    Mapping : 4.70s  |  Total: 51.36s

    ─────────────────────────────────────────
    Peptides submitted : 58920
    Mapped             : 58920 (100.0%)
    unique            : 11583
    shared (>1 prot.) : 47337
    cross-proteome    : 0
    Unmapped           : 0
    ─────────────────────────────────────────

    Peptide -> Protein -> Species  (top 20 by total peptides)
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    Organism                                      Unique    Shared     Total  %total  ProtsHit  TotalProt
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    Homo sapiens                                   11583     47337     58920   100.0     35721    210230
    ------------------------------------------  --------  --------  --------  ------  --------  --------

    Unique/Shared/Total = peptide counts; %total = % of all submitted peptides
    ProtsHit = proteins with >=1 mapped peptide; TotalProt = proteins in database
```


### check for overlapping peptides by [Chopping Proteins to Peptides](https://fuzzylife.substack.com/p/chopping-proteins-to-peptides)

```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
wget https://raw.githubusercontent.com/animesh/scripts/d08e6ffb9c9fdc23f37fd12cdd280966d4c9e46e/pepCleave.pl
perl pepCleave.pl uniprot_sprot.fasta 6 7
grep -v "^>" uniprot_sprot.fasta.len6to7.fasta > uniprot_sprot.fasta.len6to7.txt
wc uniprot_sprot.fasta.len6to7.*
  246986978  2411497827 23586051901 uniprot_sprot.fasta.len6to7.fasta
  123493489   123493489   953910896 uniprot_sprot.fasta.len6to7.txt
  370480467  2534991316 24539962797 total
shuf -n 1000000 uniprot_sprot.fasta.len6to7.txt > uniprot_sprot.fasta.len6to7.1M.txt
wc uniprot_sprot.fasta.len6to7.1M.txt
    1000000 1000000 7724587 uniprot_sprot.fasta.len6to7.1M.txt
cargo run --bin bit-pep -- run-prot uniprot_sprot.fasta -p uniprot_sprot.fasta.len6to7.1M.txt -j 12 
        Finished `dev` profile [unoptimized + debuginfo] target(s) in 13.93s
        Running `target/debug/bit-pep run-prot uniprot_sprot.fasta -p uniprot_sprot.fasta.len6to7.1M.txt -j 12`
    bit-pep RunProt  (FM-index parallel peptide->proteome search)
    ===============================================================

    [1/4] Resolving proteome...
    uniprot_sprot.fasta (274.2 MB)
    Loading: [00:00:19 ████████████████████████████████████████ 574627/574627]                                                                                                                                                        574627 proteins indexed in 111.6s  (12 threads)

    [3/4] Loading peptides...
    INFO: loaded 1000000 peptides
    1000000 peptides

    [4/4] Mapping (12 threads)...
    Mapping: [00:00:02 ████████████████████████████████████████ 1000000/1000000]                                                                                                                                                    
    ===============================================================
    Output  : uniprot_sprot.fasta.len6to7.1M.pep.tsv
    Mapping : 28.70s  |  Total: 142.25s

    ─────────────────────────────────────────
    Peptides submitted : 1000000
    Mapped             : 1000000 (100.0%)
    unique            : 555151
    shared (>1 prot.) : 444849
    cross-proteome    : 0
    Unmapped           : 0
    ─────────────────────────────────────────

    Peptide -> Protein -> Species  (top 20 by total peptides)
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    Organism                                      Unique    Shared     Total  %total  ProtsHit  TotalProt
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    Homo sapiens                                   32014    112950    144964    14.5     19744     20431
    Mus musculus                                   24876    106637    131513    13.2     16792     17252
    Arabidopsis thaliana                           39189     54429     93618     9.4     15850     16418
    Rattus norvegicus                               6352     54232     60584     6.1      7987      8226
    Saccharomyces cerevisiae (strain ATCC 2...     16969     25417     42386     4.2      6334      6733
    Caenorhabditis elegans                         15456     20616     36072     3.6      4385      4499
    Bos taurus                                      5282     30584     35866     3.6      5806      6052
    Schizosaccharomyces pombe (strain 972 /...     15849     20016     35865     3.6      4950      5129
    Drosophila melanogaster                        14325     21369     35694     3.6      3759      3868
    Dictyostelium discoideum                       13770     17937     31707     3.2      3924      4163
    Oryza sativa subsp. japonica                    7349     18022     25371     2.5      4083      4197
    Danio rerio                                     8474     16259     24733     2.5      3279      3369
    Xenopus laevis                                  5299     16559     21858     2.2      3421      3514
    Escherichia coli (strain K12)                   3961     16827     20788     2.1      4149      4531
    Bacillus subtilis (strain 168)                  7076     11749     18825     1.9      3865      4191
    Gallus gallus                                   4101     13011     17112     1.7      2221      2314
    Pongo abelii                                     479     14992     15471     1.5      2144      2218
    Mycobacterium tuberculosis (strain ATCC...       787     11260     12047     1.2      2188      2338
    Xenopus tropicalis                              2245      9345     11590     1.2      1667      1713
    Mycobacterium tuberculosis (strain CDC ...        11      9846      9857     1.0      1782      1899
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    ... and 14102 more organisms

    Unique/Shared/Total = peptide counts; %total = % of all submitted peptides
    ProtsHit = proteins with >=1 mapped peptide; TotalProt = proteins in database
```

### confirm position hits by [2 selected peptides](pep2Q9UPZ3.txt) from [Q9UPZ3](https://www.uniprot.org/uniprotkb/Q9UPZ3/entry) [sequences](https://rest.uniprot.org/uniprotkb/Q9UPZ3.fasta)
```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
wget https://raw.githubusercontent.com/animesh/scripts/079de9ce5a6d408d7508fb276f6ad248c58b99d7/pep2protmap.pl
perl pep2protmap.pl uniprot_sprot.fasta pep2Q9UPZ3.txt 1

    Read# 574627 sequences from uniprot_sprot.fasta

    Opening peptide list from pep2Q9UPZ3.txt

    >sp|Q9UPZ3|HPS5_HUMAN   569;
    >sp|Q9UPZ3|HPS5_HUMAN   111;

    Processed 2 Sequences
    Found 2 Matches

cargo run --bin bit-pep -- run-prot uniprot_sprot.fasta  -p pep2Q9UPZ3.txt
        Finished `dev` profile [unoptimized + debuginfo] target(s) in 3.24s
        Running `target/debug/bit-pep run-prot uniprot_sprot.fasta -p pep2Q9UPZ3.txt`
    bit-pep RunProt  (FM-index parallel peptide->proteome search)
    ===============================================================

    [1/4] Resolving proteome...
    uniprot_sprot.fasta (274.2 MB)
    Loading: [00:00:19 ████████████████████████████████████████ 574627/574627]                                                                                                                                                        574627 proteins indexed in 112.6s  (12 threads)

    [3/4] Loading peptides...
    INFO: loaded 2 peptides
    2 peptides

    [4/4] Mapping (12 threads)...
    Mapping: [00:00:00 ████████████████████████████████████████ 2/2]                                                                                                                                                                
    ===============================================================
    Output  : pep2Q9UPZ3.pep.tsv
    Mapping : 0.00s  |  Total: 112.60s

    ─────────────────────────────────────────
    Peptides submitted : 2
    Mapped             : 2 (100.0%)
    unique            : 2
    shared (>1 prot.) : 0
    cross-proteome    : 0
    Unmapped           : 0
    ─────────────────────────────────────────

    Peptide -> Protein -> Species  (top 20 by total peptides)
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    Organism                                      Unique    Shared     Total  %total  ProtsHit  TotalProt
    ------------------------------------------  --------  --------  --------  ------  --------  --------
    Homo sapiens                                       2         0         2   100.0         1     20431
    ------------------------------------------  --------  --------  --------  ------  --------  --------

    Unique/Shared/Total = peptide counts; %total = % of all submitted peptides
    ProtsHit = proteins with >=1 mapped peptide; TotalProt = proteins in database

cat pep2Q9UPZ3.pep.tsv 
    peptide_id      sequence        protein_acc     protein_name    proteome_id     start   end     score   mismatches      status
    pep_0   PELRGDEQSCEEDVSSDTCPK   Q9UPZ3  HPS5_HUMAN BLOC-2 complex member HPS5 [HPS5] Homo sapiens       uniprot_sprot.fasta     570     590     1.0000  0       unique
    pep_1   PEQMYVSSEHK     Q9UPZ3  HPS5_HUMAN BLOC-2 complex member HPS5 [HPS5] Homo sapiens       uniprot_sprot.fasta     112     122     1.0000  0       unique
```

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
