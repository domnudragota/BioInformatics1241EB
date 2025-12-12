#!/usr/bin/env python3
"""
Lab11 - Exercise 2 (extended for console visualization + tables)
Pairwise alignment of multiple HIV-1 genomes (multi-FASTA) using a piecewise strategy,
and show alignment + a proportional "rectangular bar" similarity visualization in console.

Also produces:
- summary.csv (pair list)
- summary_table.txt (pretty console-style table)
- identity_matrix.csv (NxN)
- score_matrix.csv (NxN)

Stdlib only.
"""

import argparse
import bisect
import csv
import os
from itertools import combinations


# -----------------------------
# FASTA parsing
# -----------------------------

IUPAC_DNA = set("ACGTRYSWKMBDHVN")


def parse_fasta(path: str):
    """
    Parses a multi-FASTA file.
    Returns: list of dicts: {"id": str, "header": str, "seq": str}
    """
    records = []
    header = None
    seq_parts = []

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_parts).upper()
                    records.append(
                        {
                            "id": header.split()[0],
                            "header": header,
                            "seq": seq,
                        }
                    )
                header = line[1:].strip()
                seq_parts = []
            else:
                cleaned = []
                for ch in line.upper():
                    if ch in IUPAC_DNA:
                        cleaned.append(ch)
                if cleaned:
                    seq_parts.append("".join(cleaned))

    if header is not None:
        seq = "".join(seq_parts).upper()
        records.append(
            {
                "id": header.split()[0],
                "header": header,
                "seq": seq,
            }
        )

    return records


# -----------------------------
# Needleman–Wunsch (small pieces)
# -----------------------------

def nw_align(s1: str, s2: str, match: int, mismatch: int, gap: int):
    """
    Global alignment via Needleman–Wunsch.
    Intended for small pieces (e.g., <= 600-1000 bases).

    Returns: (aligned1, aligned2)
    """
    m = len(s1)
    n = len(s2)

    score = [[0] * (n + 1) for _ in range(m + 1)]
    tb = [[0] * (n + 1) for _ in range(m + 1)]  # 1=diag, 2=up, 3=left

    for i in range(1, m + 1):
        score[i][0] = score[i - 1][0] + gap
        tb[i][0] = 2

    for j in range(1, n + 1):
        score[0][j] = score[0][j - 1] + gap
        tb[0][j] = 3

    for i in range(1, m + 1):
        a = s1[i - 1]
        row = score[i]
        for j in range(1, n + 1):
            b = s2[j - 1]

            diag = score[i - 1][j - 1] + (match if a == b else mismatch)
            up = score[i - 1][j] + gap
            left = row[j - 1] + gap

            best = diag
            direction = 1  # prefer diagonal on ties
            if up > best:
                best = up
                direction = 2
            if left > best:
                best = left
                direction = 3

            row[j] = best
            tb[i][j] = direction

    i, j = m, n
    a1 = []
    a2 = []
    while i > 0 or j > 0:
        d = tb[i][j]
        if d == 1:
            a1.append(s1[i - 1])
            a2.append(s2[j - 1])
            i -= 1
            j -= 1
        elif d == 2:
            a1.append(s1[i - 1])
            a2.append("-")
            i -= 1
        else:
            a1.append("-")
            a2.append(s2[j - 1])
            j -= 1

    a1.reverse()
    a2.reverse()
    return "".join(a1), "".join(a2)


# -----------------------------
# Anchors (unique k-mers) + chaining
# -----------------------------

def build_unique_kmer_pos(seq: str, k: int):
    """
    Returns dict {kmer: position_in_seq} for k-mers that appear exactly once in seq.
    Skips kmers containing ambiguous letters (anything not A/C/G/T).
    """
    counts = {}
    L = len(seq)

    for i in range(L - k + 1):
        kmer = seq[i : i + k]
        if any(ch not in "ACGT" for ch in kmer):
            continue
        counts[kmer] = counts.get(kmer, 0) + 1

    pos = {}
    for i in range(L - k + 1):
        kmer = seq[i : i + k]
        if any(ch not in "ACGT" for ch in kmer):
            continue
        if counts.get(kmer, 0) == 1:
            pos[kmer] = i

    return pos


def unique_kmer_counts(seq: str, k: int):
    counts = {}
    L = len(seq)
    for i in range(L - k + 1):
        kmer = seq[i : i + k]
        if any(ch not in "ACGT" for ch in kmer):
            continue
        counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def find_anchors_unique(seq1: str, seq2: str, k: int):
    """
    Finds anchors as exact matches of unique k-mers in BOTH sequences.
    Returns list of (pos1, pos2).
    """
    pos1 = build_unique_kmer_pos(seq1, k)
    c2 = unique_kmer_counts(seq2, k)

    anchors = []
    L2 = len(seq2)
    for j in range(L2 - k + 1):
        kmer = seq2[j : j + k]
        if any(ch not in "ACGT" for ch in kmer):
            continue
        if c2.get(kmer, 0) == 1:
            p1 = pos1.get(kmer)
            if p1 is not None:
                anchors.append((p1, j))
    return anchors


def lis_chain(anchors):
    """Longest increasing subsequence on pos2 after sorting by pos1."""
    if not anchors:
        return []

    anchors = sorted(anchors, key=lambda x: (x[0], x[1]))

    tails = []
    tails_idx = []
    prev = [-1] * len(anchors)

    for i, (_, p2) in enumerate(anchors):
        pos = bisect.bisect_left(tails, p2)
        if pos == len(tails):
            tails.append(p2)
            tails_idx.append(i)
        else:
            tails[pos] = p2
            tails_idx[pos] = i

        prev[i] = tails_idx[pos - 1] if pos > 0 else -1

    k = tails_idx[-1]
    chain = []
    while k != -1:
        chain.append(anchors[k])
        k = prev[k]
    chain.reverse()
    return chain


def thin_chain(chain, k: int, anchor_step: int):
    """Reduce anchor density to avoid too many tiny segments."""
    out = []
    last1 = -10**9
    last2 = -10**9
    min_step = max(k, anchor_step)

    for p1, p2 in chain:
        if p1 >= last1 + min_step and p2 >= last2 + min_step:
            out.append((p1, p2))
            last1, last2 = p1, p2

    return out


# -----------------------------
# Piecewise alignment
# -----------------------------

def count_non_gaps(aln: str) -> int:
    return sum(1 for c in aln if c != "-")


def align_by_windows(s1: str, s2: str, match: int, mismatch: int, gap: int, max_piece: int, overlap: int = 80):
    """
    Fallback step-by-step alignment:
    Align windows of size max_piece, advance by consumed chars (minus overlap).
    """
    i = 0
    j = 0
    out1 = []
    out2 = []

    max_iters = 200000
    it = 0

    while (i < len(s1) or j < len(s2)) and it < max_iters:
        it += 1
        p1 = s1[i : i + max_piece]
        p2 = s2[j : j + max_piece]
        if not p1 and not p2:
            break

        a1, a2 = nw_align(p1, p2, match, mismatch, gap)
        out1.append(a1)
        out2.append(a2)

        di = count_non_gaps(a1)
        dj = count_non_gaps(a2)
        if di == 0 and dj == 0:
            break

        step_i = max(di - overlap, 1) if i + di < len(s1) else di
        step_j = max(dj - overlap, 1) if j + dj < len(s2) else dj
        i = min(i + step_i, len(s1))
        j = min(j + step_j, len(s2))

    return "".join(out1), "".join(out2)


def align_piecewise(
    s1: str,
    s2: str,
    match: int,
    mismatch: int,
    gap: int,
    k: int,
    min_k: int,
    max_piece: int,
    anchor_step: int,
    depth: int = 0,
    max_depth: int = 30,
):
    """
    Strategy:
    - if small enough -> NW
    - else find anchor chain -> recursively align gaps between anchors -> stitch
    - if no anchors -> decrease k; if still none -> window fallback
    """
    if len(s1) <= max_piece and len(s2) <= max_piece:
        return nw_align(s1, s2, match, mismatch, gap)

    if depth >= max_depth:
        return align_by_windows(s1, s2, match, mismatch, gap, max_piece=max_piece, overlap=80)

    anchors = find_anchors_unique(s1, s2, k)
    chain = thin_chain(lis_chain(anchors), k, anchor_step)

    if not chain:
        if k > min_k:
            return align_piecewise(
                s1, s2, match, mismatch, gap,
                k=k - 2, min_k=min_k,
                max_piece=max_piece, anchor_step=anchor_step,
                depth=depth + 1, max_depth=max_depth,
            )
        return align_by_windows(s1, s2, match, mismatch, gap, max_piece=max_piece, overlap=80)

    out1 = []
    out2 = []
    prev1 = 0
    prev2 = 0

    for p1, p2 in chain:
        seg1 = s1[prev1:p1]
        seg2 = s2[prev2:p2]

        a1, a2 = align_piecewise(
            seg1, seg2, match, mismatch, gap,
            k=k, min_k=min_k,
            max_piece=max_piece, anchor_step=anchor_step,
            depth=depth + 1, max_depth=max_depth,
        )
        out1.append(a1)
        out2.append(a2)

        anchor_seq = s1[p1:p1 + k]
        out1.append(anchor_seq)
        out2.append(anchor_seq)

        prev1 = p1 + k
        prev2 = p2 + k

    a1, a2 = align_piecewise(
        s1[prev1:], s2[prev2:], match, mismatch, gap,
        k=k, min_k=min_k,
        max_piece=max_piece, anchor_step=anchor_step,
        depth=depth + 1, max_depth=max_depth,
    )
    out1.append(a1)
    out2.append(a2)

    return "".join(out1), "".join(out2)


# -----------------------------
# Console visualization helpers
# -----------------------------

def make_match_line(aln1: str, aln2: str) -> str:
    """'|' match, '.' mismatch, ' ' gap"""
    out = []
    for a, b in zip(aln1, aln2):
        if a == "-" or b == "-":
            out.append(" ")
        elif a == b:
            out.append("|")
        else:
            out.append(".")
    return "".join(out)


def alignment_score(aln1: str, aln2: str, match: int, mismatch: int, gap: int) -> int:
    s = 0
    for a, b in zip(aln1, aln2):
        if a == "-" or b == "-":
            s += gap
        elif a == b:
            s += match
        else:
            s += mismatch
    return s


def identity_excluding_gaps(aln1: str, aln2: str):
    matches = 0
    compared = 0
    for a, b in zip(aln1, aln2):
        if a == "-" or b == "-":
            continue
        compared += 1
        if a == b:
            matches += 1
    pct = (matches / compared * 100.0) if compared else 0.0
    return matches, compared, pct


def similarity_bar(aln1: str, aln2: str, width: int = 100, win: int = 60, threshold: float = 0.85):
    """
    Proportional "rectangular bar" similarity view:
    draw '█' where windowed identity >= threshold.
    """
    L = len(aln1)
    if L == 0:
        return ""

    bar = []
    for b in range(width):
        start = (b * L) // width
        end = ((b + 1) * L) // width
        if end <= start:
            end = min(start + 1, L)

        mid = (start + end) // 2
        w_start = max(0, mid - win // 2)
        w_end = min(L, w_start + win)

        matches = 0
        compared = 0
        for i in range(w_start, w_end):
            a = aln1[i]
            c = aln2[i]
            if a == "-" or c == "-":
                continue
            compared += 1
            if a == c:
                matches += 1

        ident = (matches / compared) if compared else 0.0
        bar.append("█" if ident >= threshold else " ")

    return "".join(bar)


def print_alignment_blocks(aln1: str, aln2: str, block: int = 80, head_tail_blocks: int = 3):
    """
    Prints alignment in blocks.
    If huge, print head+tail blocks.
    """
    mid = make_match_line(aln1, aln2)
    total_blocks = (len(aln1) + block - 1) // block

    def print_block(bi: int):
        start = bi * block
        end = min((bi + 1) * block, len(aln1))
        print(aln1[start:end])
        print(mid[start:end])
        print(aln2[start:end])
        print()

    if total_blocks <= head_tail_blocks * 2:
        for bi in range(total_blocks):
            print_block(bi)
        return

    for bi in range(head_tail_blocks):
        print_block(bi)

    omitted = total_blocks - 2 * head_tail_blocks
    print(f"... ({omitted} blocks omitted) ...\n")

    for bi in range(total_blocks - head_tail_blocks, total_blocks):
        print_block(bi)


# -----------------------------
# File output + tables
# -----------------------------

def safe_filename(name: str) -> str:
    return "".join(ch if (ch.isalnum() or ch in "._-") else "_" for ch in name)


def write_alignment_txt(path: str, id1: str, id2: str, aln1: str, aln2: str, score: int, ident_pct: float, block: int = 80):
    mid = make_match_line(aln1, aln2)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"{id1} vs {id2}\n")
        f.write(f"alignment_length={len(aln1)} score={score} identity_excl_gaps={ident_pct:.2f}%\n\n")

        for start in range(0, len(aln1), block):
            f.write(aln1[start : start + block] + "\n")
            f.write(mid[start : start + block] + "\n")
            f.write(aln2[start : start + block] + "\n\n")


def write_matrix_csv(path: str, ids, matrix, fmt=str):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([""] + ids)
        for i, rid in enumerate(ids):
            row = [rid]
            for j in range(len(ids)):
                val = matrix[i][j]
                row.append("" if val is None else fmt(val))
            w.writerow(row)


def format_console_table(headers, rows):
    col_w = [len(h) for h in headers]
    for r in rows:
        for i, cell in enumerate(r):
            col_w[i] = max(col_w[i], len(cell))

    def line(sep="-"):
        return sep.join(sep * (w + 2) for w in col_w)

    out = []
    out.append(" | ".join(headers[i].ljust(col_w[i]) for i in range(len(headers))))
    out.append(line(sep="-"))
    for r in rows:
        out.append(" | ".join(r[i].ljust(col_w[i]) for i in range(len(headers))))
    return "\n".join(out)


# -----------------------------
# Main
# -----------------------------

def main():
    parser = argparse.ArgumentParser(description="Lab11 Ex2: piecewise pairwise alignment + console bars + tables (stdlib only).")
    parser.add_argument("-i", "--input", required=True, help="Path to multi-FASTA with ~10 HIV-1 genomes")
    parser.add_argument("-o", "--out", default="labs_python/lab11/out/ex2_hiv1_pairwise", help="Output directory")

    parser.add_argument("--match", type=int, default=1)
    parser.add_argument("--mismatch", type=int, default=-1)
    parser.add_argument("--gap", type=int, default=-1)

    parser.add_argument("--k", type=int, default=20, help="Anchor k-mer size (default: 20)")
    parser.add_argument("--min-k", type=int, default=12, help="Minimum k if anchors are missing (default: 12)")
    parser.add_argument("--max-piece", type=int, default=600, help="Max piece length to run NW directly (default: 600)")
    parser.add_argument("--anchor-step", type=int, default=250, help="Minimum spacing between anchors (default: 250)")
    parser.add_argument("--block", type=int, default=80, help="Alignment block width for printing/writing (default: 80)")

    # Console visualization controls
    parser.add_argument("--print-bar", action="store_true", help="Print similarity bar in console for each pair")
    parser.add_argument("--bar-width", type=int, default=100, help="Similarity bar width in characters (default: 100)")
    parser.add_argument("--bar-window", type=int, default=60, help="Similarity window size used per bar bin (default: 60)")
    parser.add_argument("--bar-threshold", type=float, default=0.85, help="Identity threshold for bar blocks (default: 0.85)")

    parser.add_argument("--print-align", action="store_true", help="Print alignment blocks in console for each pair")
    parser.add_argument("--head-tail", type=int, default=3, help="If huge, print N head blocks + N tail blocks (default: 3)")

    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    records = parse_fasta(args.input)
    if len(records) < 2:
        raise SystemExit("Need at least 2 FASTA records.")

    print(f"Loaded {len(records)} sequences from: {args.input}")
    for r in records:
        print(f"  - {r['id']}: {len(r['seq'])} bp")

    # prepare matrices + pretty table rows
    ids = [r["id"] for r in records]
    id_to_idx = {rid: i for i, rid in enumerate(ids)}
    n = len(ids)

    score_mat = [[None] * n for _ in range(n)]
    ident_mat = [[None] * n for _ in range(n)]
    for i in range(n):
        score_mat[i][i] = 0
        ident_mat[i][i] = 100.0

    table_rows = []

    summary_path = os.path.join(args.out, "summary.csv")
    with open(summary_path, "w", newline="", encoding="utf-8") as csvf:
        writer = csv.writer(csvf)
        writer.writerow(
            [
                "seq1",
                "seq2",
                "len1",
                "len2",
                "alignment_len",
                "matches_no_gaps",
                "compared_no_gaps",
                "identity_no_gaps_pct",
                "score",
            ]
        )

        pair_count = 0
        total_pairs = (len(records) * (len(records) - 1)) // 2

        for a, b in combinations(records, 2):
            pair_count += 1
            id1 = a["id"]
            id2 = b["id"]
            s1 = a["seq"]
            s2 = b["seq"]

            print(f"\n[{pair_count:02d}/{total_pairs:02d}] {id1}  VS  {id2}")

            aln1, aln2 = align_piecewise(
                s1,
                s2,
                match=args.match,
                mismatch=args.mismatch,
                gap=args.gap,
                k=args.k,
                min_k=args.min_k,
                max_piece=args.max_piece,
                anchor_step=args.anchor_step,
            )

            score = alignment_score(aln1, aln2, args.match, args.mismatch, args.gap)
            matches, compared, ident_pct = identity_excluding_gaps(aln1, aln2)

            # store in matrices
            i = id_to_idx[id1]
            j = id_to_idx[id2]
            score_mat[i][j] = score
            score_mat[j][i] = score
            ident_mat[i][j] = ident_pct
            ident_mat[j][i] = ident_pct

            # row for pretty console table
            table_rows.append([
                str(pair_count),
                id1,
                id2,
                f"{ident_pct:.2f}",
                str(score),
                str(len(aln1)),
            ])

            # Console bar visualization
            if args.print_bar:
                bar = similarity_bar(
                    aln1, aln2,
                    width=args.bar_width,
                    win=args.bar_window,
                    threshold=args.bar_threshold,
                )
                print(f"score={score}  identity={ident_pct:.2f}%  aln_len={len(aln1)}")
                print(f"[0]{bar}[end]")
                print("    " + "█ = similar region (windowed identity >= threshold)")

            # Console alignment (head+tail blocks)
            if args.print_align:
                print()
                print_alignment_blocks(aln1, aln2, block=args.block, head_tail_blocks=args.head_tail)

            # Write per-pair alignment file
            out_name = safe_filename(f"{id1}__VS__{id2}.txt")
            out_path = os.path.join(args.out, out_name)
            write_alignment_txt(out_path, id1, id2, aln1, aln2, score, ident_pct, block=args.block)

            # Summary CSV
            writer.writerow(
                [
                    id1,
                    id2,
                    len(s1),
                    len(s2),
                    len(aln1),
                    matches,
                    compared,
                    f"{ident_pct:.4f}",
                    score,
                ]
            )

    # ---- pretty table + matrices ----
    headers = ["#", "seq1", "seq2", "similarity_%", "nw_score", "aln_len"]
    table_text = format_console_table(headers, table_rows)

    print("\n=== Pairwise Similarity Table (identity excluding gaps %) ===")
    print(table_text)

    summary_table_path = os.path.join(args.out, "summary_table.txt")
    with open(summary_table_path, "w", encoding="utf-8") as f:
        f.write(table_text + "\n")

    identity_matrix_path = os.path.join(args.out, "identity_matrix.csv")
    score_matrix_path = os.path.join(args.out, "score_matrix.csv")
    write_matrix_csv(identity_matrix_path, ids, ident_mat, fmt=lambda x: f"{x:.2f}")
    write_matrix_csv(score_matrix_path, ids, score_mat, fmt=lambda x: str(int(x)))

    print(f"\nSaved summary table: {summary_table_path}")
    print(f"Saved identity matrix: {identity_matrix_path}")
    print(f"Saved score matrix: {score_matrix_path}")

    print(f"\nDone.")
    print(f"Alignments written to: {args.out}")
    print(f"Summary CSV: {summary_path}")


if __name__ == "__main__":
    main()
