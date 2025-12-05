#!/usr/bin/env python3
"""
Lab 10 - DNA promoter patterns using C+G% and Kappa Index of Coincidence

Usage examples
--------------
# 1) Just test the sequence from the lab statement
python lab10.py

# 2) Use one or more promoter sequences from a FASTA file
python lab10.py data/promoters.fasta --window 30 --step 1 --prefix promoters
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple

import matplotlib.pyplot as plt

try:
    from Bio import SeqIO  # Biopython
except ImportError:
    SeqIO = None  # Script still works for the built-in test sequence


# ---------------------------------------------------------------------------
# Constants from the lab statement
# ---------------------------------------------------------------------------

TEST_SEQUENCE = (
    "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
)

DEFAULT_WINDOW = 30
DEFAULT_STEP = 1


# ---------------------------------------------------------------------------
# Core computations
# ---------------------------------------------------------------------------

def gc_percent(sequence: str) -> float:
    """
    Compute C+G% for a whole sequence.

    Returns a value in [0, 100], rounded to 2 decimals.
    For TEST_SEQUENCE this should be 29.27.
    """
    seq = sequence.upper()
    length = len(seq)
    if length == 0:
        return 0.0

    c_count = seq.count("C")
    g_count = seq.count("G")
    cg = c_count + g_count

    return round(cg / length * 100.0, 2)


def kappa_ic(sequence: str) -> float:
    """
    Compute Kappa Index of Coincidence (KIC) for a single DNA sequence.

    Implementation follows the algorithm described in the promoter
    analysis papers (PromKappa):contentReference[oaicite:1]{index=1}:

    1. For every shift d from 1 to N-1:
       - compare sequence[i] with sequence[i + d] for all i
       - compute the percentage of matches for that shift
    2. Average those percentages over all shifts.

    The result is rounded to 2 decimals.
    """
    seq = sequence.upper()
    L = len(seq)
    if L < 2:
        return 0.0

    # Number of different shifts
    N = L - 1
    total = 0.0

    for d in range(1, N + 1):
        matches = 0
        window_len = L - d  # number of pairs compared at this shift

        for i in range(window_len):
            if seq[i] == seq[i + d]:
                matches += 1

        total += (matches / window_len) * 100.0

    ic_value = round(total / N, 2)
    return ic_value


def sliding_windows(sequence: str, window_size: int, step: int = 1) -> Iterable[str]:
    """
    Generate all sliding windows of fixed size along the sequence.
    """
    seq = sequence.upper()
    for start in range(0, len(seq) - window_size + 1, step):
        yield seq[start:start + window_size]


@dataclass
class Pattern:
    gc_values: List[float]
    kappa_values: List[float]

    @property
    def center_of_weight(self) -> Tuple[float, float]:
        """
        Center of weight of the pattern: mean GC% and mean Kappa IC
        over all sliding windows.
        """
        if not self.gc_values:
            return float("nan"), float("nan")
        cx = sum(self.gc_values) / len(self.gc_values)
        cy = sum(self.kappa_values) / len(self.kappa_values)
        return cx, cy


def build_pattern(
    sequence: str,
    window_size: int,
    step: int = 1,
    normalize_gc: bool = False,
) -> Pattern:
    """
    For a given sequence, generate the DNA pattern points.

    For each sliding window:
      - compute GC% in that window
      - compute Kappa IC in that window

    If normalize_gc=True, GC% is rescaled relative to the total GC%
    of the whole sequence (as described in the DNA patterns papers:contentReference[oaicite:2]{index=2}).
    For the lab, using the raw GC% is usually enough, so default=False.
    """
    seq = sequence.upper()
    gc_tot = gc_percent(seq)  # used only if normalize_gc is True

    gc_values: List[float] = []
    kappa_values: List[float] = []

    for window in sliding_windows(seq, window_size, step):
        raw_gc = gc_percent(window)

        if normalize_gc and gc_tot > 0:
            gc_value = raw_gc / gc_tot * 100.0
        else:
            gc_value = raw_gc

        kappa_value = kappa_ic(window)

        gc_values.append(gc_value)
        kappa_values.append(kappa_value)

    return Pattern(gc_values=gc_values, kappa_values=kappa_values)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_pattern(pattern: Pattern, output_path: Path, title: str) -> None:
    """
    Scatter plot of (GC%, Kappa IC) for all sliding windows.
    """
    fig, ax = plt.subplots(figsize=(6, 6))

    ax.scatter(pattern.gc_values, pattern.kappa_values, s=10, alpha=0.7)

    # mark center of weight
    cx, cy = pattern.center_of_weight
    ax.scatter([cx], [cy], s=60, marker="x")
    ax.text(cx, cy, "  center", fontsize=8, va="center")

    ax.set_xlabel("C+G% (per window)")
    ax.set_ylabel("Kappa IC")
    ax.set_title(title)
    ax.grid(True, linestyle="--", alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_centers(
    centers: List[Tuple[float, float]],
    output_path: Path,
    title: str,
) -> None:
    """
    Plot the center of weight for each pattern on a second chart.
    """
    if not centers:
        return

    xs = [c[0] for c in centers]
    ys = [c[1] for c in centers]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(xs, ys, s=40)

    # label points with 1, 2, 3, ...
    for idx, (x, y) in enumerate(centers, start=1):
        ax.text(x, y, f" {idx}", fontsize=8, va="center")

    ax.set_xlabel("C+G% (center of weight)")
    ax.set_ylabel("Kappa IC (center of weight)")
    ax.set_title(title)
    ax.grid(True, linestyle="--", alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_sequences_from_fasta(path: Path) -> List[str]:
    """
    Load promoter sequences from a FASTA file (all records).
    """
    if SeqIO is None:
        raise RuntimeError("Biopython is not installed, can't read FASTA files.")

    sequences: List[str] = []
    for record in SeqIO.parse(str(path), "fasta"):
        sequences.append(str(record.seq).upper())

    return sequences


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Lab 10: DNA promoter patterns using C+G% and Kappa IC",
    )
    parser.add_argument(
        "fasta",
        nargs="?",
        help="Optional FASTA file with one or more promoter sequences. "
             "If omitted, the built-in test sequence from the lab statement is used.",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=DEFAULT_WINDOW,
        help=f"Sliding window size (default: {DEFAULT_WINDOW})",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=DEFAULT_STEP,
        help=f"Sliding window step (default: {DEFAULT_STEP})",
    )
    parser.add_argument(
        "--prefix",
        default="lab10_pattern",
        help="Prefix for output PNG files (default: lab10_pattern)",
    )
    parser.add_argument(
        "--normalize-gc",
        action="store_true",
        help="Normalize per-window GC%% relative to total GC%% (PromKappa style).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.fasta:
        sequences = load_sequences_from_fasta(Path(args.fasta))
        print(f"Loaded {len(sequences)} promoter sequence(s) from {args.fasta}")
    else:
        sequences = [TEST_SEQUENCE]
        print("No FASTA provided, using test sequence from lab statement.")
        print(f"Test sequence length: {len(TEST_SEQUENCE)} bp")

    centers: List[Tuple[float, float]] = []

    for idx, seq in enumerate(sequences, start=1):
        print(f"\n=== Sequence {idx} ===")
        cg_tot = gc_percent(seq)
        kappa_tot = kappa_ic(seq)
        print(f"Total C+G%: {cg_tot:.2f}")
        print(f"Total Kappa IC: {kappa_tot:.2f}")

        pattern = build_pattern(
            seq,
            window_size=args.window,
            step=args.step,
            normalize_gc=args.normalize_gc,
        )

        cx, cy = pattern.center_of_weight
        centers.append((cx, cy))
        print(f"Sliding windows: {len(pattern.gc_values)}")
        print(f"Center of weight: GC%={cx:.2f}, Kappa={cy:.2f}")

        out_png = Path(f"{args.prefix}_seq{idx}.png")
        plot_pattern(pattern, out_png, title=f"Promoter pattern #{idx}")
        print(f"Pattern plot saved to: {out_png}")

    # Second chart: centers of all patterns
    centers_png = Path(f"{args.prefix}_centers.png")
    plot_centers(centers, centers_png, title="Centers of promoter patterns")
    print(f"\nCenters plot saved to: {centers_png}")


if __name__ == "__main__":
    main()
