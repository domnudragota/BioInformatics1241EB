"""
Lab: CpG Islands using 1st-order Markov models (transition matrices)

Given:
- S1: CpG+ (ISLAND) training sequence
- S2: CpG- (NON-ISLAND) training sequence
- S : unknown sequence to classify

Steps:
1) Count transition frequencies from S1 and S2
2) Convert counts to transition probabilities (row-normalized)
3) Build log-likelihood / log-odds matrix:
      beta[x][y] = log2( P_plus(y|x) / P_minus(y|x) )
4) Score S by summing beta over adjacent pairs
5) Decide CpG+ if score > 0 else CpG-

Notes:
- With pseudocount=0.0, probabilities can be 0, producing +/-inf in beta and scores.
- If you want to avoid infinities, set PSEUDOCOUNT to a small value like 1e-6 or 1.0.
"""

from __future__ import annotations

import math
from typing import Dict, Tuple, List

ALPHABET = "ACGT"

S1 = "ATCGATTCGATATCATACACGTAT"          # CpG+ (ISLAND)
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"         # CpG- (NON-ISLAND)
S  = "CAGGTTGGAAACGTAA"                   # unknown (classify)

# Change this if you want smoothing to avoid 0 probabilities and +/-inf logs.
# Typical choices:
#   0.0   -> raw (may create inf / -inf)
#   1e-6  -> tiny smoothing (keeps close to raw, avoids inf)
#   1.0   -> Laplace smoothing
PSEUDOCOUNT = 1.0


Matrix = Dict[str, Dict[str, float]]


def check_seq(seq: str, alphabet: str = ALPHABET) -> None:
    bad = sorted(set(seq) - set(alphabet))
    if bad:
        raise ValueError(f"Sequence contains invalid symbols: {bad}")


def transition_counts(seq: str, alphabet: str = ALPHABET) -> Matrix:
    """
    Counts transitions X->Y for all adjacent pairs in the sequence.
    Returns counts[x][y].
    """
    counts: Matrix = {a: {b: 0.0 for b in alphabet} for a in alphabet}
    for x, y in zip(seq, seq[1:]):
        counts[x][y] += 1.0
    return counts


def transition_probs(counts: Matrix, pseudocount: float = 0.0, alphabet: str = ALPHABET) -> Matrix:
    """
    Row-normalize counts into probabilities P(next|current).
    Optional pseudocount avoids exact zeros.
    """
    probs: Matrix = {a: {b: 0.0 for b in alphabet} for a in alphabet}

    for a in alphabet:
        row_sum = sum(counts[a][b] for b in alphabet)
        denom = row_sum + pseudocount * len(alphabet)

        for b in alphabet:
            if denom > 0:
                probs[a][b] = (counts[a][b] + pseudocount) / denom
            else:
                probs[a][b] = 0.0

    return probs


def log2_ratio_matrix(p_plus: Matrix, p_minus: Matrix, alphabet: str = ALPHABET) -> Matrix:
    """
    beta[x][y] = log2( P_plus(y|x) / P_minus(y|x) )

    Zero-handling:
    - if pp==0 and pm>0  => -inf
    - if pm==0 and pp>0  => +inf
    - if both 0          => 0 (uninformative)
    """
    beta: Matrix = {a: {b: 0.0 for b in alphabet} for a in alphabet}

    for a in alphabet:
        for b in alphabet:
            pp = p_plus[a][b]
            pm = p_minus[a][b]

            if pp == 0.0 and pm == 0.0:
                beta[a][b] = 0.0
            elif pm == 0.0:
                beta[a][b] = float("inf")
            elif pp == 0.0:
                beta[a][b] = float("-inf")
            else:
                beta[a][b] = math.log2(pp / pm)

    return beta


def score_sequence(seq: str, beta: Matrix, show_steps: bool = True) -> float:
    """
    Score(seq) = sum_{i=1..n-1} beta[seq[i]][seq[i+1]]
    """
    total = 0.0
    steps: List[Tuple[int, str, float]] = []

    for i, (x, y) in enumerate(zip(seq, seq[1:]), start=1):
        v = beta[x][y]
        total += v
        steps.append((i, x + y, v))

    if show_steps:
        print("Per-transition contributions:")
        for i, pair, v in steps:
            if math.isinf(v):
                s = "inf" if v > 0 else "-inf"
            else:
                s = f"{v:.6f}"
            print(f"{i:>2}. {pair}: {s}")
        print()

    return total


def fmt(v: float, decimals: int = 3) -> str:
    if math.isinf(v):
        return "inf" if v > 0 else "-inf"
    return f"{v:.{decimals}f}"


def print_matrix(m: Matrix, title: str, decimals: int = 3) -> None:
    print(title)
    header = "     " + "  ".join(f"{b:>10}" for b in ALPHABET)
    print(header)
    for a in ALPHABET:
        row = "  ".join(f"{fmt(m[a][b], decimals):>10}" for b in ALPHABET)
        print(f"{a:>2} | {row}")
    print()


def main() -> None:
    check_seq(S1)
    check_seq(S2)
    check_seq(S)

    # 1) Counts
    c_plus = transition_counts(S1)
    c_minus = transition_counts(S2)

    # 2) Probabilities
    p_plus = transition_probs(c_plus, pseudocount=PSEUDOCOUNT)
    p_minus = transition_probs(c_minus, pseudocount=PSEUDOCOUNT)

    # 3) Log-likelihood / log-odds matrix
    beta = log2_ratio_matrix(p_plus, p_minus)

    # Print results
    print(f"PSEUDOCOUNT = {PSEUDOCOUNT}\n")
    print_matrix(c_plus, "Counts (CpG+ / ISLAND) from S1", decimals=0)
    print_matrix(c_minus, "Counts (CpG- / NON-ISLAND) from S2", decimals=0)

    print_matrix(p_plus, "Transition Probabilities M1: CpG+ (ISLAND)", decimals=6)
    print_matrix(p_minus, "Transition Probabilities M2: CpG- (NON-ISLAND)", decimals=6)

    print_matrix(beta, "Log-likelihood / Log-odds matrix β = log2(P+/P-)", decimals=6)

    # 4) Score and classify S
    score = score_sequence(S, beta, show_steps=True)

    print(f"Total log-odds score for S = {score}")
    if score > 0:
        print("Decision: CpG+ (ISLAND)  [score > 0]")
    else:
        print("Decision: CpG- (NON-ISLAND)  [score <= 0]")


if __name__ == "__main__":
    main()
