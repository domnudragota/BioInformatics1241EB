#!/usr/bin/env python3
"""
We build a position-specific model from known motif instances (length L), then scan a
long sequence S with a sliding window of length L.

Matrices:
  1) Count matrix:         counts[b][pos]
  2) Relative frequencies: freq[b][pos] = counts[b][pos] / N
  3) Weight matrix:        w[b][pos] = freq[b][pos] / bg[b]
  4) Log-likelihoods:      ll[b][pos] = ln(w[b][pos])

Scoring:
  score(window) = sum_{pos=1..L} ll[ window[pos] ][pos]

Notes:
- With pseudocount=0 (the default), some frequencies are 0, so log-likelihood can be -inf.
  This makes the scan very strict (any "forbidden" base at a fully conserved position
  kills the window).
- If you want smoother behavior, run with --pseudocount 1.

This lab's known motifs (10 sequences, length 9) are taken from the provided image.
"""

from __future__ import annotations

import argparse
import math
from typing import Dict, List, Tuple


ALPHABET = "ACGT"

# Known exon-intron boundary motif instances (length 9), from the lab screenshot.
DEFAULT_MOTIFS = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT",
]

DEFAULT_SEQUENCE_S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"


Matrix = Dict[str, List[float]]


def validate_motifs(motifs: List[str]) -> int:
    if not motifs:
        raise ValueError("No motif sequences provided.")

    length = len(motifs[0])
    for m in motifs:
        if len(m) != length:
            raise ValueError("All motif sequences must have the same length.")
        for ch in m:
            if ch not in ALPHABET:
                raise ValueError(f"Invalid base '{ch}' in motif '{m}'. Expected only A/C/G/T.")
    return length


def build_count_matrix(motifs: List[str], motif_len: int) -> Matrix:
    counts: Matrix = {b: [0.0] * motif_len for b in ALPHABET}

    for m in motifs:
        for i, ch in enumerate(m):
            counts[ch][i] += 1.0

    return counts


def add_pseudocounts(counts: Matrix, pseudocount: float) -> Matrix:
    if pseudocount == 0:
        return {b: col[:] for b, col in counts.items()}

    adjusted: Matrix = {b: [0.0] * len(counts[b]) for b in ALPHABET}
    for b in ALPHABET:
        for i, val in enumerate(counts[b]):
            adjusted[b][i] = val + pseudocount
    return adjusted


def build_frequency_matrix(counts: Matrix) -> Matrix:
    motif_len = len(next(iter(counts.values())))
    freqs: Matrix = {b: [0.0] * motif_len for b in ALPHABET}

    for i in range(motif_len):
        column_total = 0.0
        for b in ALPHABET:
            column_total += counts[b][i]

        if column_total == 0:
            # This should never happen if motifs are valid.
            for b in ALPHABET:
                freqs[b][i] = 0.0
        else:
            for b in ALPHABET:
                freqs[b][i] = counts[b][i] / column_total

    return freqs


def build_weight_matrix(freqs: Matrix, background: Dict[str, float]) -> Matrix:
    motif_len = len(next(iter(freqs.values())))
    weights: Matrix = {b: [0.0] * motif_len for b in ALPHABET}

    for b in ALPHABET:
        bg = background[b]
        if bg <= 0:
            raise ValueError("Background probabilities must be > 0.")
        for i in range(motif_len):
            weights[b][i] = freqs[b][i] / bg

    return weights


def build_log_likelihood_matrix(weights: Matrix) -> Matrix:
    motif_len = len(next(iter(weights.values())))
    ll: Matrix = {b: [0.0] * motif_len for b in ALPHABET}

    for b in ALPHABET:
        for i in range(motif_len):
            w = weights[b][i]
            if w == 0.0:
                ll[b][i] = float("-inf")
            else:
                ll[b][i] = math.log(w)  # natural log
    return ll


def score_window(window: str, ll: Matrix) -> float:
    score = 0.0
    for pos, ch in enumerate(window):
        if ch not in ll:
            return float("-inf")
        value = ll[ch][pos]
        if value == float("-inf"):
            return float("-inf")
        score += value
    return score


def scan_sequence(sequence: str, ll: Matrix, motif_len: int) -> List[Tuple[int, str, float]]:
    results: List[Tuple[int, str, float]] = []
    for start in range(0, len(sequence) - motif_len + 1):
        window = sequence[start : start + motif_len]
        score = score_window(window, ll)
        results.append((start, window, score))
    return results


def format_number(x: float, digits: int = 3) -> str:
    if x == float("-inf"):
        return "-inf"
    return f"{x:.{digits}f}"


def print_matrix(title: str, matrix: Matrix, digits: int = 3, ints: bool = False) -> None:
    motif_len = len(next(iter(matrix.values())))
    print()
    print(title)
    print("     " + " ".join(f"{i+1:>7}" for i in range(motif_len)))

    for b in ALPHABET:
        row = []
        for i in range(motif_len):
            val = matrix[b][i]
            if ints:
                row.append(f"{int(val):>7}")
            else:
                row.append(f"{format_number(val, digits):>7}")
        print(f"{b:>3} " + " ".join(row))


def main() -> None:
    parser = argparse.ArgumentParser(description="Lab12 motif finding with log-likelihood scanning.")
    parser.add_argument(
        "--pseudocount",
        type=float,
        default=0.0,
        help="Add this pseudocount to each base count before normalizing (default: 0).",
    )
    parser.add_argument(
        "--sequence",
        type=str,
        default=DEFAULT_SEQUENCE_S,
        help="Sequence S to scan (default: the one given in the lab).",
    )
    args = parser.parse_args()

    motifs = DEFAULT_MOTIFS[:]
    motif_len = validate_motifs(motifs)

    # Null/background model (random DNA): P(A)=P(C)=P(G)=P(T)=0.25
    background = {b: 0.25 for b in ALPHABET}

    # 1) Count matrix
    counts = build_count_matrix(motifs, motif_len)

    # Optional: add pseudocounts before normalizing
    counts_pc = add_pseudocounts(counts, args.pseudocount)

    # 2) Relative frequencies matrix
    freqs = build_frequency_matrix(counts_pc)

    # 3) Weight matrix (likelihood ratio vs background)
    weights = build_weight_matrix(freqs, background)

    # 4) Log-likelihoods matrix
    ll = build_log_likelihood_matrix(weights)

    # Display matrices
    print(f"Known motifs: {len(motifs)} sequences, motif length L={motif_len}")
    print(f"Using pseudocount = {args.pseudocount}")
    print_matrix("1) Count Matrix", counts, ints=True)
    print_matrix("2) Relative Frequencies Matrix", freqs, digits=3)
    print_matrix("3) Weight Matrix (freq / background)", weights, digits=3)
    print_matrix("4) Log-Likelihoods Matrix ln(weight)", ll, digits=3)

    # 5) Scan S with a sliding window of length L
    sequence = args.sequence.strip().upper()
    results = scan_sequence(sequence, ll, motif_len)

    # Compute a simple reference threshold: the minimum score among the training motifs
    training_scores = [score_window(m, ll) for m in motifs]
    threshold = min(training_scores)

    print()
    print("5) Sliding window scores")
    print(f"Sequence length: {len(sequence)}")
    print(f"Windows scored:  {len(results)}")
    print(f"Reference threshold (min training motif score): {format_number(threshold, 3)}")
    print()
    print(f"{'start(0b)':>9} {'start(1b)':>9} {'window':>{motif_len}} {'score':>12}")
    print("-" * (9 + 1 + 9 + 1 + motif_len + 1 + 12))

    for start, window, score in results:
        print(f"{start:>9} {start+1:>9} {window:>{motif_len}} {format_number(score, 3):>12}")

    # Show the best hits
    finite_results = [r for r in results if r[2] != float("-inf")]
    finite_results.sort(key=lambda x: x[2], reverse=True)

    print()
    if not finite_results:
        print("No finite scores (everything became -inf). Try --pseudocount 1.")
        return

    best_start, best_window, best_score = finite_results[0]
    print("Top hit:")
    print(f"  start (0-based) = {best_start}")
    print(f"  start (1-based) = {best_start + 1}")
    print(f"  window          = {best_window}")
    print(f"  score           = {format_number(best_score, 3)}")

    # Signals: any hit above threshold?
    hits = [r for r in finite_results if r[2] >= threshold]
    print()
    if hits:
        print("Signals found (windows with score >= threshold):")
        for start, window, score in hits:
            print(f"  start={start:>2} (1-based {start+1:>2})  window={window}  score={format_number(score,3)}")
        print()
        print("Conclusion: YES, S likely contains an exon-intron border motif-like pattern.")
    else:
        print("No window exceeded the threshold.")
        print("Conclusion: NO strong signal found for an exon-intron border in S.")


if __name__ == "__main__":
    main()
