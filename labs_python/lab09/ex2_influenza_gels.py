#!/usr/bin/env python3
"""
Lab 09 – Influenza genomes + restriction enzymes + "difference" gel

Usage:
    python labs_python/lab09/ex2_influenza_gels.py labs_python/lab09/data/influenza
"""

import sys
import os
import glob
from typing import List, Tuple, Dict, Set

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# ---------- Helpers from previous exercise ----------

def read_fasta_sequence(path: str) -> str:
    """Read the first sequence from a FASTA file and return it as a string."""
    with open(path, "r") as f:
        lines = f.readlines()

    seq_lines = [line.strip() for line in lines if not line.startswith(">")]
    seq = "".join(seq_lines).upper()
    return seq


# cut_index = index inside recognition site where the cut happens (0-based)
ENZYMES = {
    "EcoRI":  {"site": "GAATTC", "cut_index": 1},  # G|AATTC
    "BamHI":  {"site": "GGATCC", "cut_index": 1},  # G|GATCC
    "HindIII": {"site": "AAGCTT", "cut_index": 1}, # A|AGCTT
    "TaqI":   {"site": "TCGA",   "cut_index": 1},  # T|CGA
    "HaeIII": {"site": "GGCC",   "cut_index": 2},  # GG|CC
}


def find_cut_positions(seq: str, recog: str, cut_index: int) -> List[int]:
    """Find all cut positions for one enzyme. Returns list of 0-based cut indices."""
    seq = seq.upper()
    recog = recog.upper()
    cuts = []
    start = 0

    while True:
        i = seq.find(recog, start)
        if i == -1:
            break
        cuts.append(i + cut_index)
        start = i + 1  # continue after this position

    return cuts


def digest_sequence(seq: str, recog: str, cut_index: int) -> Tuple[List[int], List[int]]:
    """
    Return:
      - internal cut positions (0-based) – not including 0 and len(seq)
      - list of fragment lengths
    """
    n = len(seq)
    cut_positions = [0] + sorted(find_cut_positions(seq, recog, cut_index)) + [n]

    fragments_lengths: List[int] = []
    internal_cuts = cut_positions[1:-1]

    for i in range(len(cut_positions) - 1):
        start = cut_positions[i]
        end = cut_positions[i + 1]
        fragments_lengths.append(end - start)

    return internal_cuts, fragments_lengths


# ---------- Simple gel plotting (reused) ----------

def plot_gel(enzymes_fragments: List[Tuple[str, List[int]]],
             title: str = "Simulated DNA gel") -> None:
    """
    Draw a simple fake gel:
      - one lane per enzyme
      - smaller fragments migrate further down
      - y-axis ~ fragment size (bp)
    """
    # Collect all fragment lengths to scale the band positions
    all_lengths = [L for _, frags in enzymes_fragments for L in frags]
    if not all_lengths:
        print("Nothing to plot (no fragment lengths).")
        return

    max_len = max(all_lengths)
    min_len = min(all_lengths)
    span = max_len - min_len if max_len != min_len else 1

    fig, ax = plt.subplots(figsize=(6, 6))
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    lane_count = len(enzymes_fragments)
    x_positions = [i * 2 for i in range(lane_count)]

    ax.set_xlim(-1, (lane_count - 1) * 2 + 1)
    ax.set_ylim(0, 1)
    ax.invert_yaxis()  # top = wells

    # Draw lanes
    for lane_index, (name, fragment_lengths) in enumerate(enzymes_fragments):
        x_center = x_positions[lane_index]

        # Well at the top
        ax.add_line(Line2D([x_center - 0.4, x_center + 0.4],
                           [0.02, 0.02],
                           linewidth=8,
                           color="grey"))

        # Bands
        for length in fragment_lengths:
            y = (max_len - length) / span * 0.9 + 0.05
            ax.add_line(Line2D([x_center - 0.3, x_center + 0.3],
                               [y, y],
                               linewidth=4,
                               color="white"))

    # X labels = enzyme names
    ax.set_xticks(x_positions)
    ax.set_xticklabels(
        [name for name, _ in enzymes_fragments],
        color="white",
        fontsize=9
    )

    # Y axis ~ fragment size
    step = max_len // 5 or 1
    tick_lengths = list(range(step, max_len + 1, step))
    if len(tick_lengths) > 6:
        tick_lengths = [tick_lengths[i]
                        for i in range(0, len(tick_lengths),
                                       max(1, len(tick_lengths) // 5))]

    y_positions = [(max_len - L) / span * 0.9 + 0.05 for L in tick_lengths]
    ax.set_yticks(y_positions)
    ax.set_yticklabels([str(L) for L in tick_lengths], color="white", fontsize=8)
    ax.set_ylabel("Fragment size (bp)", color="white")

    # Cosmetics
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="x", colors="white")
    ax.tick_params(axis="y", colors="white")
    ax.set_title(title, color="white", fontsize=10)

    plt.tight_layout()
    plt.savefig("influenza_difference_gel.png", dpi=300)
    plt.show()


# ---------- Main logic for this exercise ----------

def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python labs_python/lab09/ex2_influenza_gels.py path/to/influenza_folder")
        sys.exit(1)

    folder = sys.argv[1]

    # Find all FASTA files (*.fa, *.fasta)
    fasta_files = sorted(
        glob.glob(os.path.join(folder, "*.fa")) +
        glob.glob(os.path.join(folder, "*.fasta")) +
        glob.glob(os.path.join(folder, "*.fna"))
    )

    if not fasta_files:
        print(f"No FASTA files found in folder: {folder}")
        sys.exit(1)

    print(f"Found {len(fasta_files)} influenza genomes:")
    for path in fasta_files:
        print("  -", os.path.basename(path))

    # per_enzyme_per_genome[enzyme_name] = list of sets of fragment lengths, one set per genome
    per_enzyme_per_genome: Dict[str, List[Set[int]]] = {name: [] for name in ENZYMES.keys()}

    # ----- Digest each genome with each enzyme -----
    for idx, path in enumerate(fasta_files):
        seq = read_fasta_sequence(path)
        print(f"\nGenome {idx + 1}: {os.path.basename(path)} (length = {len(seq)} bp)")

        for enzyme_name, info in ENZYMES.items():
            cuts, frag_lengths = digest_sequence(seq, info["site"], info["cut_index"])

            # store as set to avoid counting identical length twice in the same genome
            per_enzyme_per_genome[enzyme_name].append(set(frag_lengths))

            print(
                f"  {enzyme_name}: cuts = {len(cuts)}, "
                f"fragments = {len(frag_lengths)}"
            )

    # ----- Compute "difference bands" per enzyme -----
    difference_lanes: List[Tuple[str, List[int]]] = []

    print("\n=== DIFFERENCE ANALYSIS (bands NOT common to all genomes) ===")

    num_genomes = len(fasta_files)

    for enzyme_name, frag_sets_for_genomes in per_enzyme_per_genome.items():
        # all fragment lengths that ever appear for this enzyme
        all_lengths: Set[int] = set().union(*frag_sets_for_genomes)

        # fragment lengths that are present in EVERY genome for this enzyme
        common_lengths: Set[int] = set()
        for L in all_lengths:
            count = sum(1 for s in frag_sets_for_genomes if L in s)
            if count == num_genomes:
                common_lengths.add(L)

        # remove common bands -> keep only differences
        diff_lengths = sorted(all_lengths - common_lengths)

        if diff_lengths:
            difference_lanes.append((enzyme_name, diff_lengths))
            print(f"\n{enzyme_name}:")
            print("  Common bands (in ALL 10 genomes):", sorted(common_lengths) if common_lengths else "none")
            print("  Unique / differing bands (shown on gel):", diff_lengths)
        else:
            print(f"\n{enzyme_name}: all fragment sizes are common to all genomes.")
            print("  -> This enzyme will not have a lane in the difference gel.")

    if not difference_lanes:
        print("\nNo differences found: all enzymes give identical patterns across the 10 genomes.")
        return

    # ----- Plot the "difference gel" -----
    plot_gel(
        difference_lanes,
        title="Influenza genomes – electrophoresis bands differing between strains"
    )


if __name__ == "__main__":
    main()
