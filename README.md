# Bit-Pep: Fast Multi-Proteome Peptide Mapping

`bit-pep` maps peptide sequences to very large UniProt/proteome FASTA files and reports every matching protein and position. Protein mapping is designed as a streaming, multi-pattern search rather than an FM-index build: the proteome is scanned once, so there is no multi-gigabyte suffix-array construction, no proteome index to retain, and no artificial 2-GB/1.9-billion-residue chunk limit.

The default `run-prot` path is particularly suited to a relatively small peptide list mapped against a very large proteome, which is the common FragPipe/MaxQuant-style use case.

## Build

```bash
cargo build --release
```

Use the release binary for real runs. `cargo run` without `--release` builds an unoptimized debug binary.

## Run

```bash
./target/release/bit-pep run-prot \
    uniprot_trembl.fasta \
    --peptides peptides.txt \
    --output peptides.mapped.tsv \
    --threads 64
```

Short form:

```bash
./target/release/bit-pep run-prot uniprot_trembl.fasta -p peptides.txt -j 64
```

The peptide sequence is read from column 1. Additional input columns are preserved in the output.

A simple one-peptide-per-line file is also accepted:

```text
PEQEHK
PELRGDEQSCEEDVSSDTCPK
```

A TSV with a header is also accepted:

```text
Peptide Sequence	Other information
PEQEHK	foo
PELRGDEQSCEEDVSSDTCPK	bar
```

## Mapping algorithm

`run-prot` uses a compact Aho-Corasick-style multi-pattern automaton implemented in `src/protein_scan.rs`.

For exact mapping (`--fuzzy-mismatches 0`), all submitted peptides are searched simultaneously while each protein is read. This changes the expensive part of the problem from:

```text
build a suffix array for the entire proteome
        +
search the peptides
```

to:

```text
one streaming pass through the proteome
```

For fuzzy mapping, each peptide is split into `m+1` exact seeds when `m` mismatches are allowed. A seed hit creates a candidate start position and the complete peptide is then verified directly against the protein sequence. This uses the pigeonhole property: with at most `m` mismatches, at least one of `m+1` disjoint seeds must match exactly.

Protein sequences and peptide sequences are normalized in the same way:

- case is ignored;
- `I` and `L` are treated as equivalent;
- non-alphabetic characters are removed.

Protein sequences are not reverse-complemented.

## Threads and memory

`-j/--threads` controls the streaming protein scan. The proteome is read in batches and each batch is divided between the requested worker threads. There is no Rayon thread pool involved in `run-prot`.

The mapper does not construct or retain a proteome-wide FM-index. Memory is therefore dominated by the current FASTA batch, the small peptide automaton, and accumulated hit/output metadata rather than by the number of amino acids in the proteome.

`-m/--memory` is retained for command-line compatibility with earlier FM-index versions. It is not used to size the protein-mapping algorithm and normally does not need to be supplied.

The project still contains Rayon because the original DNA/read-mapping functionality elsewhere in the library uses it. It is not part of the `run-prot` execution path.

## Progress and time estimate

`run-prot` reports:

- number of proteins scanned;
- amino acids scanned;
- current protein-throughput;
- an estimated total protein count derived from the FASTA file size and the observed record size;
- an estimated remaining scan time once the first batch has been processed;
- final scan throughput.

The estimate is empirical for the current input and machine. It is not a hard performance prediction.

## Output

The main output is `<peptide-stem>.mapped.tsv` unless `--output` is supplied.

The original input columns are followed by:

- `Protein_IDs` — semicolon-separated matching UniProt accessions;
- `Protein_Match_Count` — number of distinct matching proteins;
- `Match_Positions` — all 1-based peptide positions for each matching protein, formatted as `accession:pos1,pos2,...`.

Unmatched peptide rows are retained.

A second file is written as `<output-stem>.species.tsv` with:

```text
Organism	Unique	Shared	Total	%total	ProtsHit	TotalProt	%Prots
```

`Unique`, `Shared`, and `Total` count submitted peptides associated with each organism. `ProtsHit` is the number of proteins from that organism containing at least one submitted peptide, and `TotalProt` is the number of proteins for that organism in the proteome.

The terminal report also prints the top 100 peptides by number of matching proteins and the top 100 organisms by total mapped peptides.

## Large UniProt databases

The previous implementation constructed an FM-index for the proteome. With the complete UniProt TrEMBL database, this required enormous suffix arrays and forced the proteome into chunks. The current `run-prot` mapper does not do that.

For example, a proteome containing roughly 58 billion amino acids is still processed as a streaming input. Its size affects runtime and I/O, but it does not require a 58-billion-element suffix array in memory.

This is a better fit for peptide-to-protein mapping when the peptide set is small compared with the proteome.

## Performance considerations

The dominant operation is now the sequential read of the proteome plus multi-threaded pattern scanning. Increasing `-j` is therefore not guaranteed to give linear speedup: FASTA decompression, filesystem throughput, memory bandwidth, and CPU/cache behavior can become limiting factors.

For exact peptide mapping, building a full suffix-array/FM-index is unnecessary when the number of query peptides is modest. The multi-pattern automaton searches all queries in one pass and avoids the index-construction cost entirely.

For fuzzy searches, runtime additionally depends on the number of seed candidates that require verification. Short peptides and large mismatch allowances can generate many candidates.

## Validation

The mapper was developed against the previous Perl `pep2protmap` semantics:

- all occurrences are reported;
- positions are 1-based in the new TSV output;
- `I`/`L` equivalence is retained;
- unmatched peptides remain in the result;
- multiple matching proteins are retained rather than truncated.

For a small known example, a peptide list can be compared directly against the original Perl mapper.

## Development notes

The repository also contains the original/general Bit-Pop DNA read-mapping functionality. The protein-specific `run-prot` path is intentionally kept separate from those algorithms because peptide mapping has different performance characteristics.

The important optimization target for `run-prot` is now the streaming scan rather than suffix-array construction. Further optimization should be driven by measured profiling, especially:

- FASTA parsing and allocation;
- memory bandwidth during protein scanning;
- candidate verification for fuzzy searches;
- hit aggregation and output generation.

SIMD verification can be considered for fuzzy matching after those components have been measured.

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

## License

MIT License
