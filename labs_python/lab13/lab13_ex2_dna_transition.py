"""
Lab 13 - Exercise 2
DNA sequence (expected 50 letters) -> transition count matrix + probability matrix -> JSON file.

Stdlib only.
"""

import argparse
import json
import sys
from typing import Dict, List, Tuple


DNA_ALPHABET = ["A", "C", "G", "T"]
DNA_SET = set(DNA_ALPHABET)


def read_text_file(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


def clean_dna(seq: str) -> str:
    # remove whitespace/newlines, uppercase
    s = "".join(ch for ch in seq.upper() if not ch.isspace())
    # validate
    bad = sorted(set(ch for ch in s if ch not in DNA_SET))
    if bad:
        raise ValueError(f"DNA sequence contains invalid characters: {bad}. Allowed: A,C,G,T")
    return s


def build_dna_transition(seq: str) -> Tuple[List[List[int]], List[List[float]]]:
    n = len(DNA_ALPHABET)
    idx: Dict[str, int] = {ch: i for i, ch in enumerate(DNA_ALPHABET)}

    counts: List[List[int]] = [[0 for _ in range(n)] for _ in range(n)]

    for i in range(len(seq) - 1):
        a = idx[seq[i]]
        b = idx[seq[i + 1]]
        counts[a][b] += 1

    probs: List[List[float]] = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        row_sum = sum(counts[i])
        if row_sum == 0:
            continue
        for j in range(n):
            probs[i][j] = counts[i][j] / row_sum

    return counts, probs


def main() -> None:
    parser = argparse.ArgumentParser(description="Exercise 2: DNA transition matrix -> JSON")
    parser.add_argument("--seq", type=str, help="DNA sequence directly (A,C,G,T only).")
    parser.add_argument("--file", type=str, help="Path to a file containing the DNA sequence.")
    parser.add_argument("--out", type=str, default="dna_transition.json", help="Output JSON path.")
    parser.add_argument("--expected_len", type=int, default=50, help="Expected DNA length (default: 50).")
    args = parser.parse_args()

    if not args.seq and not args.file:
        print("Error: provide --seq or --file", file=sys.stderr)
        sys.exit(2)
    if args.seq and args.file:
        print("Error: provide only one of --seq or --file", file=sys.stderr)
        sys.exit(2)

    raw = args.seq if args.seq else read_text_file(args.file)
    seq = clean_dna(raw)

    if len(seq) != args.expected_len:
        print(f"Warning: sequence length is {len(seq)}, expected {args.expected_len}. Continuing...")

    counts, probs = build_dna_transition(seq)

    payload = {
        "type": "dna_transition_matrix",
        "alphabet": DNA_ALPHABET,
        "sequence_length": len(seq),
        "sequence": seq,
        "counts": counts,
        "probabilities": probs,
    }

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print(f"Saved DNA transition matrix to: {args.out}")


if __name__ == "__main__":
    main()
