#!/usr/bin/env python3
"""
make_test_peptides.py
─────────────────────
Downloads a UniProt proteome (cached after first run) and produces
a fully-tryptic peptide list for use with bit-pop run-prot.

"Fully tryptic" means every peptide:
  - Is preceded by a cut-site residue (or is position 1 of the protein)
    BUT position-1 peptides are skipped because the FM-index can't
    distinguish them from non-tryptic N-terminal hits.
  - Ends with a cut-site residue (C-terminal protein peptides are skipped
    for the same reason).

This ensures the unique-peptide count from this script matches the
unique count from bit-pop run-prot exactly.

Usage examples:
  python make_test_peptides.py                         # SARS-CoV-2 (fast test)
  python make_test_peptides.py --proteome UP000005640  # Human
  python make_test_peptides.py --cut KR                # Trypsin (default)
  python make_test_peptides.py --cut K                 # Lys-C
  python make_test_peptides.py --cut R                 # Arg-C
  python make_test_peptides.py --missed 1              # allow 1 missed cleavage
  python make_test_peptides.py --fasta my.fasta        # use local FASTA
  python make_test_peptides.py --force                 # re-download even if cached
"""

import argparse, re, sys
from collections import defaultdict
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("Install requests first:  pip install requests")

# ── CLI ───────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="UniProt download + tryptic digest → peptide list for bit-pop",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
parser.add_argument("--proteome", default="UP000464024",
    help="UniProt Proteome ID (default: SARS-CoV-2 UP000464024)")
parser.add_argument("--fasta", default=None,
    help="Use a local FASTA file instead of downloading")
parser.add_argument("--out", default="peptides.txt",
    help="Output peptide list (default: peptides.txt)")
parser.add_argument("--cut", default="KR",
    help=("Residues to cut AFTER (default: KR = trypsin).\n"
          "  KR = trypsin (cuts after K or R, not before P)\n"
          "  K  = Lys-C   (cuts after K only)\n"
          "  R  = Arg-C   (cuts after R only)\n"
          "  Any combination of single-letter codes is accepted."))
parser.add_argument("--missed", type=int, default=0,
    help="Allowed missed cleavages (default: 0)")
parser.add_argument("--min-len", type=int, default=6,
    help="Minimum peptide length (default: 6)")
parser.add_argument("--max-len", type=int, default=50,
    help="Maximum peptide length (default: 50)")
parser.add_argument("--force", action="store_true",
    help="Re-download even if FASTA is already cached")
args = parser.parse_args()

CUT_SITES = args.cut.upper()
print(f"Protease cut sites : {CUT_SITES}")
print(f"  → cuts after    : {' or '.join(list(CUT_SITES))}")
print(f"  → no cut before : P  (proline rule)")
print(f"Missed cleavages   : {args.missed}")
print(f"Length range       : {args.min_len}–{args.max_len} AA")
print()

# ── Step 1: FASTA ─────────────────────────────────────────────────────────────
if args.fasta:
    fasta_text = Path(args.fasta).read_text(encoding="utf-8")
    print(f"Using local FASTA: {args.fasta}")
else:
    fasta_path = f"{args.proteome}.fasta"
    if not args.force and Path(fasta_path).exists() and Path(fasta_path).stat().st_size > 1000:
        print(f"Using cached FASTA: {fasta_path}")
        fasta_text = Path(fasta_path).read_text(encoding="utf-8")
    else:
        url = (f"https://rest.uniprot.org/uniprotkb/stream"
               f"?format=fasta&query=proteome:{args.proteome}&compressed=false")
        print(f"Downloading {args.proteome} from UniProt...")
        r = requests.get(url, headers={"User-Agent": "make_test_peptides/1.0"}, timeout=300)
        r.raise_for_status()
        fasta_text = r.text
        Path(fasta_path).write_text(fasta_text, encoding="utf-8")
        print(f"Saved FASTA -> {fasta_path}")

# ── Step 2: Parse FASTA ───────────────────────────────────────────────────────
def parse_fasta(text):
    acc = desc = None
    seq_parts = []
    for line in text.splitlines():
        if line.startswith(">"):
            if acc:
                yield acc, desc, "".join(seq_parts)
            parts = line[1:].split("|")
            if len(parts) >= 3:
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

# ── Step 3: Digest ────────────────────────────────────────────────────────────
def digest(seq, cut_sites, missed, min_len, max_len):
    """
    Fully-tryptic digest with missed cleavages.

    Returns list of (peptide_str, start_0based) tuples where:
      - The residue at (start-1) is a cut-site residue (N-terminal tryptic)
        — first peptide of protein is EXCLUDED (no preceding K/R)
      - The peptide ends with a cut-site residue (C-terminal tryptic)
        — last peptide of protein is EXCLUDED (no trailing K/R)

    This mirrors the tryptic validation in bit-pop run-prot exactly.
    """
    seq = seq.upper()
    n = len(seq)

    # Find cleavage sites: positions AFTER which we cut
    # i.e., index i such that seq[i] in cut_sites and seq[i+1] != 'P'
    sites = []
    for i, aa in enumerate(seq):
        if aa in cut_sites:
            if i + 1 < n and seq[i + 1] == "P":
                continue   # proline rule: no cut before P
            sites.append(i + 1)   # cut point = start of next peptide

    # Fragment boundaries: [0, site1, site2, ..., n]
    bounds = [0] + sites + [n]
    n_frags = len(bounds) - 1

    peptides = []
    for i in range(n_frags):
        for mc in range(missed + 1):
            j = i + mc + 1
            if j > n_frags:
                break
            pep_start = bounds[i]
            pep_end   = bounds[j]
            pep       = seq[pep_start:pep_end]

            # Enforce fully-tryptic:
            # N-term: must be preceded by a cut-site residue (skip first fragment)
            if pep_start == 0:
                continue   # protein N-terminal peptide — excluded
            n_term_ok = seq[pep_start - 1] in cut_sites

            # C-term: must end with a cut-site residue (skip last fragment)
            if pep_end == n:
                continue   # protein C-terminal peptide — excluded
            c_term_ok = seq[pep_end - 1] in cut_sites

            if not (n_term_ok and c_term_ok):
                continue   # missed cleavage span broke one end

            # Length filter
            if len(pep) < min_len or len(pep) > max_len:
                continue

            # Only standard + IUPAC residues (no digits, spaces, etc.)
            if not re.fullmatch(r"[A-Z]+", pep):
                continue

            peptides.append((pep, pep_start))
    return peptides

all_peptides  = []
pep_to_prots  = defaultdict(set)  # sequence -> set of accessions

for acc, desc, seq in proteins:
    for pep, pos in digest(seq, CUT_SITES, args.missed, args.min_len, args.max_len):
        all_peptides.append(pep)
        pep_to_prots[pep].add(acc)

unique_seqs       = set(all_peptides)
unique_to_one     = sum(1 for s in pep_to_prots.values() if len(s) == 1)
unique_to_multi   = len(unique_seqs) - unique_to_one

print(f"\nDigest summary:")
print(f"  Protease           : {'Trypsin' if CUT_SITES=='KR' else CUT_SITES}")
print(f"  Cut sites          : {CUT_SITES} (not before P)")
print(f"  Missed cleavages   : {args.missed}")
print(f"  Total peptides     : {len(all_peptides)}")
print(f"  Unique sequences   : {len(unique_seqs)}")
print(f"  Unique to 1 protein: {unique_to_one}")
print(f"  Shared (>1 protein): {unique_to_multi}")
print(f"  Length range       : {args.min_len}–{args.max_len} AA")
print()

# ── Step 4: Write outputs ─────────────────────────────────────────────────────
sorted_peps = sorted(unique_seqs)

with open(args.out, "w") as f:
    for pep in sorted_peps:
        f.write(pep + "\n")
print(f"Wrote {len(sorted_peps)} unique peptides -> {args.out}")

stem = Path(args.out).stem
sample_file = f"{stem}_sample10.txt"
with open(sample_file, "w") as f:
    for pep in sorted_peps[:10]:
        f.write(pep + "\n")
print(f"Wrote 10-peptide sample  -> {sample_file}")

spacesep_file = f"{stem}_spacesep.txt"
with open(spacesep_file, "w") as f:
    f.write(" ".join(sorted_peps[:20]) + "\n")
print(f"Wrote space-separated    -> {spacesep_file}")

print(f"""
Done. The unique count here ({unique_to_one} unique-to-1-protein) should match
bit-pop run-prot exactly because both enforce fully-tryptic hits.

Next:
  cargo run --bin bit-pop -- run-prot {args.proteome}.fasta -p {args.out}
""")
