#!/usr/bin/env python3
"""
make_test_peptides.py
─────────────────────
Downloads a UniProt proteome (cached after first run) and produces
a peptide list for use with bit-pep run-prot.

Two modes:

  tryptic  (default)
    Fully-tryptic peptides with proline rule.
    N-terminal and C-terminal protein peptides are excluded so that
    the unique count here matches bit-pep run-prot exactly.

  sliding
    All overlapping subsequences of length --min-len to --max-len.
    Use with: cargo run --bin bit-pep -- run-prot ... --sliding

Usage examples:
  python make_test_peptides.py                           # SARS-CoV-2, tryptic
  python make_test_peptides.py --proteome UP000005640    # Human, tryptic
  python make_test_peptides.py --mode sliding            # all overlapping subseqs
  python make_test_peptides.py --mode sliding --min-len 7 --max-len 25
  python make_test_peptides.py --cut K                   # Lys-C digest
  python make_test_peptides.py --missed 1                # 1 missed cleavage
  python make_test_peptides.py --force                   # re-download FASTA
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
    description="UniProt download + peptide generation for bit-pep",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
parser.add_argument("--proteome", default="UP000464024",
    help="UniProt Proteome ID (default: SARS-CoV-2 UP000464024)")
parser.add_argument("--fasta", default=None,
    help="Use a local FASTA file instead of downloading")
parser.add_argument("--out", default="peptides.txt",
    help="Output peptide list (default: peptides.txt)")
parser.add_argument("--mode", default="tryptic", choices=["tryptic", "sliding"],
    help=("tryptic  = fully-tryptic digest, proline rule, no N/C-term peptides (default)\n"
          "sliding  = all overlapping subsequences of --min-len to --max-len"))
parser.add_argument("--cut", default="KR",
    help="Tryptic mode: residues to cut after (default: KR = trypsin)")
parser.add_argument("--missed", type=int, default=0,
    help="Tryptic mode: allowed missed cleavages (default: 0)")
parser.add_argument("--min-len", type=int, default=6,
    help="Minimum peptide length (default: 6)")
parser.add_argument("--max-len", type=int, default=50,
    help="Maximum peptide length (default: 50)")
parser.add_argument("--force", action="store_true",
    help="Re-download even if FASTA is already cached")
args = parser.parse_args()

CUT_SITES = args.cut.upper()

print(f"Mode               : {args.mode}")
if args.mode == "tryptic":
    print(f"Cut sites          : {CUT_SITES} (not before P)")
    print(f"Missed cleavages   : {args.missed}")
print(f"Length range       : {args.min_len}–{args.max_len} AA")
print()

# ── Step 1: FASTA ─────────────────────────────────────────────────────────────
if args.fasta:
    fasta_text = Path(args.fasta).read_text(encoding="utf-8")
    print(f"Using local FASTA: {args.fasta}")
else:
    fasta_path = f"{args.proteome}.fasta"
    if not args.force and Path(fasta_path).exists() \
            and Path(fasta_path).stat().st_size > 1000:
        print(f"Using cached FASTA: {fasta_path}")
        fasta_text = Path(fasta_path).read_text(encoding="utf-8")
    else:
        url = (f"https://rest.uniprot.org/uniprotkb/stream"
               f"?format=fasta&includeIsoform=true&query=proteome:{args.proteome}&compressed=false")
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
            acc  = parts[1] if len(parts) >= 3 else line[1:].split()[0]
            desc = parts[2].split(" OS=")[0] if len(parts) >= 3 else line[1:]
            seq_parts = []
        elif line:
            seq_parts.append(line.strip().rstrip("*"))
    if acc:
        yield acc, desc, "".join(seq_parts)

proteins = list(parse_fasta(fasta_text))
print(f"Parsed {len(proteins)} proteins\n")

# ── Step 3: Generate peptides ─────────────────────────────────────────────────

def tryptic_digest(seq, cut_sites, missed, min_len, max_len):
    """
    Fully-tryptic digest. Excludes N-terminal and C-terminal protein peptides
    so unique counts match bit-pep run-prot (which applies the same rules).
    """
    seq = seq.upper()
    n   = len(seq)

    # Cleavage sites: index of first residue of next peptide
    sites = []
    for i, aa in enumerate(seq):
        if aa in cut_sites:
            if i + 1 < n and seq[i + 1] == "P":
                continue  # proline rule
            sites.append(i + 1)

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

            if pep_start == 0:  continue   # protein N-terminal → skip
            if pep_end   == n:  continue   # protein C-terminal → skip

            n_term_ok = seq[pep_start - 1] in cut_sites
            c_term_ok = seq[pep_end   - 1] in cut_sites
            proline_n = seq[pep_start] == "P"        # cut before P → invalid
            proline_c = seq[pep_end]   == "P"        # cut before P → invalid

            if not (n_term_ok and c_term_ok and not proline_n and not proline_c):
                continue

            if not (args.min_len <= len(pep) <= args.max_len):
                continue
            if not re.fullmatch(r"[A-Z]+", pep):
                continue

            peptides.append((pep, pep_start))
    return peptides


def sliding_digest(seq, min_len, max_len):
    """
    All overlapping subsequences of length min_len to max_len.
    Use with --sliding in bit-pep run-prot.
    """
    seq = seq.upper()
    n   = len(seq)
    peptides = []
    for start in range(n):
        for length in range(min_len, max_len + 1):
            end = start + length
            if end > n:
                break
            pep = seq[start:end]
            if re.fullmatch(r"[A-Z]+", pep):
                peptides.append((pep, start))
    return peptides


all_peptides = []
pep_to_prots = defaultdict(set)

for acc, desc, seq in proteins:
    if args.mode == "tryptic":
        frags = tryptic_digest(seq, CUT_SITES, args.missed, args.min_len, args.max_len)
    else:
        frags = sliding_digest(seq, args.min_len, args.max_len)
    for pep, pos in frags:
        all_peptides.append(pep)
        pep_to_prots[pep].add(acc)

unique_seqs     = set(all_peptides)
unique_to_one   = sum(1 for s in pep_to_prots.values() if len(s) == 1)
unique_to_multi = len(unique_seqs) - unique_to_one

# ── Step 4: Summary ───────────────────────────────────────────────────────────
print("Digest summary:")
if args.mode == "tryptic":
    print(f"  Protease           : {'Trypsin' if CUT_SITES == 'KR' else CUT_SITES}")
    print(f"  Cut sites          : {CUT_SITES} (not before P)")
    print(f"  Missed cleavages   : {args.missed}")
else:
    print(f"  Window step        : 1 (every position)")
print(f"  Total peptides     : {len(all_peptides)}")
print(f"  Unique sequences   : {len(unique_seqs)}")
print(f"  Unique to 1 protein: {unique_to_one}")
print(f"  Shared (>1 protein): {unique_to_multi}")
print(f"  Length range       : {args.min_len}–{args.max_len} AA")
print()

# ── Step 5: Write outputs ─────────────────────────────────────────────────────
sorted_peps = sorted(unique_seqs)
stem = Path(args.out).stem

with open(args.out, "w") as f:
    for pep in sorted_peps:
        f.write(pep + "\n")
print(f"Wrote {len(sorted_peps)} unique peptides -> {args.out}")

with open(f"{stem}_sample10.txt", "w") as f:
    for pep in sorted_peps[:10]:
        f.write(pep + "\n")
print(f"Wrote 10-peptide sample  -> {stem}_sample10.txt")

with open(f"{stem}_spacesep.txt", "w") as f:
    f.write(" ".join(sorted_peps[:20]) + "\n")
print(f"Wrote space-separated    -> {stem}_spacesep.txt")

run_flag = "--sliding" if args.mode == "sliding" else ""
proteome_arg = args.fasta if args.fasta else f"{args.proteome}.fasta"
print(f"""
The unique count here ({unique_to_one}) should match bit-pep run-prot.

Next:
  cargo run --bin bit-pep -- run-prot {proteome_arg} -p {args.out} {run_flag} -j 12
""")
