#!/usr/bin/env python3
"""
make_test_peptides.py
─────────────────────
Step 1: Download a UniProt proteome FASTA  (default: SARS-CoV-2, tiny 12 proteins)
Step 2: In-silico trypsin digest (cut after K/R, not before P)
Step 3: Write a peptide list file  → peptides.txt  (one per line)
Step 4: Also write a small sample  → peptides_sample10.txt  (first 10 unique)

Usage:
    python make_test_peptides.py                        # SARS-CoV-2 (fast, 12 proteins)
    python make_test_peptides.py --proteome UP000005640 # Human (slow, ~20k proteins)
    python make_test_peptides.py --fasta my.fasta       # use local FASTA

Requires: requests  (pip install requests)
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("Install requests first:  pip install requests")

# ── CLI ───────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="UniProt download + trypsin digest")
parser.add_argument("--proteome", default="UP000464024",
                    help="UniProt Proteome ID  (default: SARS-CoV-2 UP000464024)")
parser.add_argument("--fasta",    default=None,
                    help="Use local FASTA instead of downloading")
parser.add_argument("--out",      default="peptides.txt",
                    help="Output peptide list  (default: peptides.txt)")
parser.add_argument("--missed",   type=int, default=0,
                    help="Missed cleavages  (default: 0)")
parser.add_argument("--min-len",  type=int, default=6,
                    help="Minimum peptide length  (default: 6)")
parser.add_argument("--max-len",  type=int, default=50,
                    help="Maximum peptide length  (default: 50)")
args = parser.parse_args()

# ── Step 1: FASTA ─────────────────────────────────────────────────────────────
if args.fasta:
    fasta_text = Path(args.fasta).read_text()
    print(f"Read local FASTA: {args.fasta}")
else:
    url = (f"https://rest.uniprot.org/uniprotkb/stream"
           f"?format=fasta&query=proteome:{args.proteome}&compressed=false")
    print(f"Downloading {args.proteome} from UniProt...")
    r = requests.get(url, headers={"User-Agent": "make_test_peptides/1.0"}, timeout=120)
    r.raise_for_status()
    fasta_text = r.text
    fasta_path = f"{args.proteome}.fasta"
    Path(fasta_path).write_text(fasta_text)
    print(f"Saved FASTA → {fasta_path}")

# ── Step 2: Parse FASTA ───────────────────────────────────────────────────────
def parse_fasta(text):
    """Yield (accession, description, sequence) tuples."""
    acc = desc = None
    seq_parts = []
    for line in text.splitlines():
        if line.startswith(">"):
            if acc:
                yield acc, desc, "".join(seq_parts)
            parts = line[1:].split("|")
            if len(parts) >= 3:          # sp|P12345|NAME_HUMAN ...
                acc  = parts[1]
                desc = parts[2].split(" OS=")[0]
            else:
                acc  = line[1:].split()[0]
                desc = line[1:]
            seq_parts = []
        elif line:
            seq_parts.append(line.strip().rstrip("*"))
    if acc:
        yield acc, desc, "".join(seq_parts)

proteins = list(parse_fasta(fasta_text))
print(f"Parsed {len(proteins)} proteins")

# ── Step 3: Trypsin digest ────────────────────────────────────────────────────
def trypsin(seq, missed=0):
    """
    Tryptic digest: cleave after K or R unless followed by P.
    Returns list of peptide strings.
    """
    # Find cleavage sites: after K or R not followed by P
    sites = [0]
    for i, aa in enumerate(seq):
        if aa in ("K", "R"):
            if i + 1 < len(seq) and seq[i + 1] == "P":
                continue        # no cut before P
            sites.append(i + 1)
    sites.append(len(seq))

    # Generate peptides with up to `missed` missed cleavages
    peptides = []
    n = len(sites) - 1
    for i in range(n):
        for mc in range(missed + 1):
            j = i + mc + 1
            if j > n:
                break
            pep = seq[sites[i]:sites[j]]
            peptides.append(pep)
    return peptides

all_peptides = []
protein_map  = defaultdict(list)   # peptide → [accessions]

for acc, desc, seq in proteins:
    for pep in trypsin(seq.upper(), missed=args.missed):
        pep = re.sub(r"[^A-Z]", "", pep)   # strip non-letter chars
        if args.min_len <= len(pep) <= args.max_len:
            all_peptides.append(pep)
            protein_map[pep].append(acc)

total   = len(all_peptides)
unique  = len(set(all_peptides))
unique_to_one = sum(1 for peps in protein_map.values() if len(peps) == 1)

print(f"\nDigest summary:")
print(f"  Total peptides     : {total}")
print(f"  Unique sequences   : {unique}")
print(f"  Unique to 1 protein: {unique_to_one}")
print(f"  Length range       : {args.min_len}–{args.max_len} AA")
print(f"  Missed cleavages   : {args.missed}")

# ── Step 4: Write outputs ─────────────────────────────────────────────────────
unique_peps = sorted(set(all_peptides))

with open(args.out, "w") as f:
    for pep in unique_peps:
        f.write(pep + "\n")
print(f"\nWrote {len(unique_peps)} unique peptides → {args.out}")

# Small sample for quick testing
sample_file = Path(args.out).stem + "_sample10.txt"
with open(sample_file, "w") as f:
    for pep in unique_peps[:10]:
        f.write(pep + "\n")
print(f"Wrote 10-peptide sample  → {sample_file}")

# Space-separated format on one line (tests that input parser too)
spacesep_file = Path(args.out).stem + "_spacesep.txt"
with open(spacesep_file, "w") as f:
    f.write(" ".join(unique_peps[:20]) + "\n")
print(f"Wrote space-separated    → {spacesep_file}")

print("\nDone. Next:")
print(f"  cargo test --test uniprot_integration -- --nocapture")
print(f"  # (after main.rs is wired:)")
print(f"  cargo run -- run {args.proteome}.fasta {args.out}")
