#!/usr/bin/env python3
"""
Translate a coding DNA/RNA sequence into an amino acid sequence.

Usage examples:
  # from a FASTA file (find first ATG/AUG, stop at stop codon)
  python ex01_translate.py --in data/gene.fa

  # inline sequence
  python ex01_translate.py --seq ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG

  # translate entire sequence from frame 2 (no start/stop search)
  python ex01_translate.py --in data/gene.fa --no-orf --frame 2

  # write result to file
  python ex01_translate.py --in data/gene.fa --out outputs/protein.faa
"""

from pathlib import Path
import argparse

# --- Genetic code (RNA codons -> 1-letter amino acids) ---
# Stops are represented as '*'
GENETIC_CODE = {
    # U**
    "UUU":"F","UUC":"F","UUA":"L","UUG":"L",
    "UCU":"S","UCC":"S","UCA":"S","UCG":"S",
    "UAU":"Y","UAC":"Y","UAA":"*","UAG":"*",
    "UGU":"C","UGC":"C","UGA":"*","UGG":"W",
    # C**
    "CUU":"L","CUC":"L","CUA":"L","CUG":"L",
    "CCU":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAU":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGU":"R","CGC":"R","CGA":"R","CGG":"R",
    # A**
    "AUU":"I","AUC":"I","AUA":"I","AUG":"M",
    "ACU":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAU":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGU":"S","AGC":"S","AGA":"R","AGG":"R",
    # G**
    "GUU":"V","GUC":"V","GUA":"V","GUG":"V",
    "GCU":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAU":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGU":"G","GGC":"G","GGA":"G","GGG":"G",
}
START_CODON = "AUG"
STOP_CODONS = {"UAA", "UAG", "UGA"}


def read_sequence_from_file(path: str) -> str:
    """Reads a raw/FASTA file and returns one concatenated sequence (uppercased)."""
    seq_parts = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    return "".join(seq_parts).upper()


def sanitize_to_rna(seq: str) -> str:
    """Uppercase, remove spaces, convert DNA (T) to RNA (U)."""
    seq = seq.upper().replace(" ", "")
    return seq.replace("T", "U")


def find_orf_bounds(rna: str) -> tuple[int, int]:
    """
    Return (start_idx, end_idx) for the first ORF:
    - start at first AUG
    - end at the first in-frame stop (UAA/UAG/UGA), end is exclusive of stop.
    If no valid ORF is found, returns (0, len(rna)) so we can at least translate.
    """
    start = rna.find(START_CODON)
    if start == -1:
        return 0, len(rna)

    for i in range(start, len(rna) - 2, 3):
        codon = rna[i:i+3]
        if codon in STOP_CODONS:
            return start, i  # stop codon not included
    return start, len(rna)


def translate_rna(rna: str) -> str:
    """Translate RNA (length multiple of 3 preferred). Unknown/partial codons -> 'X'."""
    protein = []
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        aa = GENETIC_CODE.get(codon, "X")
        protein.append(aa)
    return "".join(protein)


def main():
    ap = argparse.ArgumentParser(description="Translate coding DNA/RNA to amino acids.")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--in", dest="inp", help="Input file (FASTA or raw sequence)")
    src.add_argument("--seq", help="Inline sequence (DNA or RNA)")

    ap.add_argument("--no-orf", action="store_true",
                    help="Translate full sequence from a chosen frame (default: find first AUG and stop at stop codon).")
    ap.add_argument("--frame", type=int, default=1, choices=[1, 2, 3],
                    help="Reading frame to use when --no-orf (1/2/3).")
    ap.add_argument("--out", help="Write protein to file (optional)")
    args = ap.parse_args()

    # 1) get sequence
    if args.inp:
        seq = read_sequence_from_file(args.inp)
    else:
        seq = args.seq.strip().upper()

    # 2) normalize to RNA
    rna = sanitize_to_rna(seq)

    # 3) pick region to translate
    if args.no_orf:
        offset = args.frame - 1
        rna_region = rna[offset:]
    else:
        start, end = find_orf_bounds(rna)
        rna_region = rna[start:end]

    # trim trailing incomplete codon
    trim_len = len(rna_region) - (len(rna_region) % 3)
    rna_region = rna_region[:trim_len]

    # 4) translate
    protein = translate_rna(rna_region)

    # 5) output
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w") as fh:
            # simple FASTA-style header if input file name is known
            header = args.inp if args.inp else "inline_sequence"
            fh.write(f">{header}\n")
            # wrap to 60 chars per line for readability
            for i in range(0, len(protein), 60):
                fh.write(protein[i:i+60] + "\n")
        print(f"Wrote: {args.out}")
    else:
        print(protein)


if __name__ == "__main__":
    main()
