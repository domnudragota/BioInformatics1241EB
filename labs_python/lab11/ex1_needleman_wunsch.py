#!/usr/bin/env python3
"""
Lab11 - Exercise 1
Needleman–Wunsch global alignment (DNA) - minimal dependencies (stdlib only).
"""

import argparse


def read_first_fasta_sequence(path: str) -> str:
    """Reads the first FASTA record from a file and returns it as a single uppercase string."""
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # If we already collected a sequence and we meet a new header -> stop (first record only)
                if seq_parts:
                    break
                continue
            seq_parts.append(line)

    seq = "".join(seq_parts).upper()
    # keep only DNA letters (basic cleanup)
    cleaned = []
    for ch in seq:
        if ch in ("A", "C", "G", "T", "N"):
            cleaned.append(ch)
    return "".join(cleaned)


def needleman_wunsch(s1: str, s2: str, match: int, mismatch: int, gap: int):
    """
    Returns:
      score_matrix, traceback_matrix, aligned_s1, aligned_s2, final_score
    traceback_matrix contains:
      'D' = diagonal, 'U' = up, 'L' = left
    """
    s1 = s1.upper().replace(" ", "").replace("\n", "")
    s2 = s2.upper().replace(" ", "").replace("\n", "")

    m = len(s1)
    n = len(s2)

    # DP score matrix
    score = []
    for _ in range(m + 1):
        score.append([0] * (n + 1))

    # Traceback matrix (stores best direction to come from)
    tb = []
    for _ in range(m + 1):
        tb.append([None] * (n + 1))

    # Initialize first column (align s1 with gaps)
    score[0][0] = 0
    tb[0][0] = None

    for i in range(1, m + 1):
        score[i][0] = score[i - 1][0] + gap
        tb[i][0] = "U"  # came from up

    # Initialize first row (align s2 with gaps)
    for j in range(1, n + 1):
        score[0][j] = score[0][j - 1] + gap
        tb[0][j] = "L"  # came from left

    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            a = s1[i - 1]
            b = s2[j - 1]

            if a == b:
                diag = score[i - 1][j - 1] + match
            else:
                diag = score[i - 1][j - 1] + mismatch

            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap

            best = max(diag, up, left)
            score[i][j] = best

            # Tie-breaking: prefer diagonal, then up, then left (simple + stable)
            if best == diag:
                tb[i][j] = "D"
            elif best == up:
                tb[i][j] = "U"
            else:
                tb[i][j] = "L"

    # Traceback to build alignment
    i = m
    j = n
    aligned1 = []
    aligned2 = []

    while i > 0 or j > 0:
        direction = tb[i][j]

        if direction == "D":
            aligned1.append(s1[i - 1])
            aligned2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif direction == "U":
            aligned1.append(s1[i - 1])
            aligned2.append("-")
            i -= 1
        else:  # "L"
            aligned1.append("-")
            aligned2.append(s2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()

    aln1 = "".join(aligned1)
    aln2 = "".join(aligned2)
    final_score = score[m][n]
    return score, tb, aln1, aln2, final_score


def make_match_line(aln1: str, aln2: str) -> str:
    """
    Builds a middle line for display:
      '|' match
      '.' mismatch
      ' ' gap
    """
    out = []
    for a, b in zip(aln1, aln2):
        if a == "-" or b == "-":
            out.append(" ")
        elif a == b:
            out.append("|")
        else:
            out.append(".")
    return "".join(out)


def alignment_identity(aln1: str, aln2: str):
    matches = 0
    length = 0
    for a, b in zip(aln1, aln2):
        if a == "-" or b == "-":
            length += 1
            continue
        length += 1
        if a == b:
            matches += 1
    pct = (matches / length * 100.0) if length > 0 else 0.0
    return matches, length, pct


def print_alignment(aln1: str, aln2: str, block: int = 60):
    mid = make_match_line(aln1, aln2)

    for start in range(0, len(aln1), block):
        a1 = aln1[start : start + block]
        mm = mid[start : start + block]
        a2 = aln2[start : start + block]
        print(a1)
        print(mm)
        print(a2)
        print()


def print_score_matrix(score, s1: str, s2: str):
    """Optional debug: prints the DP matrix with headers."""
    s1 = " " + s1  # for row header alignment
    s2 = " " + s2  # for col header alignment

    # header row
    header = ["    "]
    for ch in s2:
        header.append(f"{ch:>4}")
    print("".join(header))

    for i in range(len(score)):
        row = [f"{s1[i]:>4}"]
        for j in range(len(score[0])):
            row.append(f"{score[i][j]:>4}")
        print("".join(row))


def main():
    parser = argparse.ArgumentParser(description="Needleman–Wunsch global alignment (DNA) - stdlib only.")
    parser.add_argument("--s1", type=str, default="ACCGTGAAGCCAATAC", help="Sequence 1 (default: lab statement)")
    parser.add_argument("--s2", type=str, default="AGCGTGCAGCCAATAC", help="Sequence 2 (default: lab statement)")
    parser.add_argument("--f1", type=str, default=None, help="Optional FASTA file for S1 (uses first record)")
    parser.add_argument("--f2", type=str, default=None, help="Optional FASTA file for S2 (uses first record)")
    parser.add_argument("--match", type=int, default=1, help="Match score (default: 1)")
    parser.add_argument("--mismatch", type=int, default=-1, help="Mismatch score (default: -1)")
    parser.add_argument("--gap", type=int, default=-1, help="Gap penalty (default: -1)")
    parser.add_argument("--block", type=int, default=60, help="Alignment print width (default: 60)")
    parser.add_argument("--matrix", action="store_true", help="Print DP score matrix (debug)")
    args = parser.parse_args()

    s1 = args.s1
    s2 = args.s2

    if args.f1 is not None:
        s1 = read_first_fasta_sequence(args.f1)
    if args.f2 is not None:
        s2 = read_first_fasta_sequence(args.f2)

    score, tb, aln1, aln2, final_score = needleman_wunsch(
        s1=s1,
        s2=s2,
        match=args.match,
        mismatch=args.mismatch,
        gap=args.gap,
    )

    matches, length, pct = alignment_identity(aln1, aln2)

    print("Needleman–Wunsch (global) alignment")
    print(f"S1 length: {len(s1)} | S2 length: {len(s2)}")
    print(f"Scoring: match={args.match}, mismatch={args.mismatch}, gap={args.gap}")
    print(f"Final score: {final_score}")
    print(f"Identity: {matches}/{length} ({pct:.2f}%)\n")

    print_alignment(aln1, aln2, block=args.block)

    if args.matrix:
        print("DP score matrix:\n")
        print_score_matrix(score, s1, s2)


if __name__ == "__main__":
    main()
