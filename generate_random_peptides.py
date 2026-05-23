#!/usr/bin/env python3
import argparse
import csv
import random
import subprocess
from pathlib import Path


def parse_uniprot_fasta(path):
    with open(path, "r", encoding="utf-8") as f:
        header = None
        seq = []
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:].strip()
                seq = []
            else:
                line = line.strip()
                if line:
                    seq.append(line.rstrip("*"))
        if header is not None:
            yield header, "".join(seq)


def parse_protein_accession(header):
    parts = header.split("|")
    if len(parts) >= 3 and parts[0] in {"sp", "tr"}:
        return parts[1]
    return header.split()[0]


def generate_random_peptides(fasta_path, min_len, max_len, sample_fraction, seed):
    if seed is not None:
        random.seed(seed)
    if min_len <= 1:
        raise ValueError("min_len must be > 1")
    if max_len < min_len:
        raise ValueError("max_len must be >= min_len")

    proteins = []
    for header, seq in parse_uniprot_fasta(fasta_path):
        if len(seq) > min_len:
            proteins.append((parse_protein_accession(header), header, seq))

    if not proteins:
        raise RuntimeError(f"No proteins longer than {min_len} aa found in {fasta_path}")

    n_proteins = len(proteins)
    sample_count = max(1, min(n_proteins, int(n_proteins * sample_fraction)))
    if sample_count < n_proteins:
        proteins = random.sample(proteins, sample_count)
    print(f"Sampled {sample_count}/{n_proteins} proteins ({sample_fraction*100:.1f}%)")

    peptides = []
    for idx, (acc, header, seq) in enumerate(proteins):
        max_pep_len = min(max_len, len(seq) - 1)
        if max_pep_len < min_len:
            continue
        pep_len = random.randint(min_len, max_pep_len)
        start = random.randint(0, len(seq) - pep_len)
        peptide = seq[start : start + pep_len].upper()
        if not peptide.isalpha():
            continue
        peptides.append(
            (
                f"pep_{idx}|{acc}|{start+1}-{start+pep_len}",
                peptide,
                acc,
                start + 1,
                start + pep_len,
                header,
            )
        )
    return peptides


def write_peptide_fasta(peptides, path):
    with open(path, "w", encoding="utf-8") as f:
        for pid, seq, *_ in peptides:
            f.write(f">{pid}\n{seq}\n")


def write_metadata(peptides, path, status_list):
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["peptide_id", "sequence", "protein_acc", "start", "end", "source_header", "verified"])
        for (pid, seq, acc, start, end, header), status in zip(peptides, status_list):
            writer.writerow([pid, seq, acc, start, end, header, status])


def read_bitpep_output(tsv_path):
    hits = {}
    with open(tsv_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            hits.setdefault(row["peptide_id"], []).append(row)
    return hits


def verify_hits(peptides, hits):
    statuses = []
    good = 0
    repeated = 0
    missing = 0
    for pid, _, exp_acc, exp_start, exp_end, _ in peptides:
        exp_start = str(exp_start)
        exp_end = str(exp_end)
        rows = hits.get(pid, [])
        if not rows:
            statuses.append("unmapped")
            missing += 1
            continue
        matched = [r for r in rows if r["protein_acc"] == exp_acc and r["start"] == exp_start and r["end"] == exp_end]
        if not matched:
            statuses.append("unmapped")
            missing += 1
            continue
        if len(rows) == 1:
            statuses.append("unique")
        else:
            statuses.append("shared")
            repeated += 1
        good += 1
    return statuses, good, repeated, missing


def write_verified_report(path, total, good, shared, missing):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"Total peptides: {total}\n")
        f.write(f"Verified peptides: {good}\n")
        f.write(f"Unique peptides: {good - shared}\n")
        f.write(f"Shared peptides: {shared}\n")
        f.write(f"Missing/mismatched peptides: {missing}\n")


def run_bitpep(fasta, peptides_file, extra_args):
    cmd = ["cargo", "run", "--bin", "bit-pep", "--", "run-prot", fasta, "-p", peptides_file]
    if extra_args:
        cmd.extend(extra_args)
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description="Generate random peptides from a UniProt FASTA.")
    parser.add_argument("--fasta", default="uniprot_sprot.fasta", help="Proteome FASTA file")
    parser.add_argument("--min-len", type=int, default=5, help="Minimum peptide length (>1)")
    parser.add_argument("--max-len", type=int, default=35, help="Maximum peptide length")
    parser.add_argument("--sample-fraction", type=float, default=0.1, help="Protein sampling fraction")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verify", action="store_true", help="Exit nonzero if any peptide is missing or mismatched")
    parser.add_argument("--no-run", action="store_true", help="Do not run bit-pep mapping")
    parser.add_argument("--bitpep-args", nargs=argparse.REMAINDER, help="Extra args for bit-pep mapping")
    args = parser.parse_args()

    stem = Path(args.fasta).stem
    base_name = f"{stem}_{args.sample_fraction}_{args.min_len}_{args.max_len}_{args.seed}"
    peptides_path = f"{base_name}.fasta"
    metadata_path = f"{base_name}.metadata.csv"
    verified_path = f"{base_name}.verified.txt"

    peptides = generate_random_peptides(
        args.fasta,
        args.min_len,
        args.max_len,
        args.sample_fraction,
        args.seed,
    )
    print(f"Generated {len(peptides)} peptides from sampled proteins in {args.fasta}")
    write_peptide_fasta(peptides, peptides_path)
    print(f"Wrote peptide FASTA: {peptides_path}")

    if args.no_run:
        write_metadata(peptides, metadata_path, ["pending"] * len(peptides))
        print(f"Wrote metadata CSV: {metadata_path}")
        return

    run_bitpep(args.fasta, peptides_path, args.bitpep_args)
    hits = read_bitpep_output(Path(peptides_path).with_suffix('.pep.tsv'))
    statuses, good, shared, missing = verify_hits(peptides, hits)
    write_metadata(peptides, metadata_path, statuses)
    write_verified_report(verified_path, len(peptides), good, shared, missing)
    print(f"Wrote metadata CSV: {metadata_path}")
    print(f"Wrote verified report: {verified_path}")
    if missing and args.verify:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
