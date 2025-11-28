"""
Lab 09 – Restriction enzymes + gel simulation

Usage:
    python labs_python/lab09/ex1_restriction_digest.py labs_python/lab09/data/sequence.fasta
"""

import sys
from typing import List, Tuple

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# ---------- STEP 1: Read DNA sequence from FASTA ----------

def read_fasta_sequence(path: str) -> str:
    """Read the first sequence from a FASTA file and return it as a string."""
    with open(path, "r") as f:
        lines = f.readlines()

    # Skip header lines starting with '>'
    seq_lines = [line.strip() for line in lines if not line.startswith(">")]
    seq = "".join(seq_lines).upper()

    return seq


# ---------- STEP 2: Define restriction enzymes ----------

# cut_index = where we cut inside the recognition site (0-based, from the left)
ENZYMES = {
    "EcoRI":  {"site": "GAATTC", "cut_index": 1},  # G|AATTC
    "BamHI":  {"site": "GGATCC", "cut_index": 1},  # G|GATCC
    "HindIII": {"site": "AAGCTT", "cut_index": 1},  # A|AGCTT
    "TaqI":   {"site": "TCGA",   "cut_index": 1},  # T|CGA
    "HaeIII": {"site": "GGCC",   "cut_index": 2},  # GG|CC
}


# ---------- STEP 3: Find cuts and fragments ----------

def find_cut_positions(seq: str, recog: str, cut_index: int) -> List[int]:
    """
    Find all cut positions for one enzyme.
    Returns a list of positions (0-based) where the DNA is cut.
    """
    seq = seq.upper()
    recog = recog.upper()
    cuts = []
    start = 0

    while True:
        i = seq.find(recog, start)
        if i == -1:
            break
        cuts.append(i + cut_index)
        start = i + 1  # continue search after this position

    return cuts


def digest_sequence(seq: str, recog: str, cut_index: int) -> Tuple[List[int], List[Tuple[int, int, int]]]:
    """
    Return:
      - list of cut positions (0-based, without 0 and len(seq))
      - list of fragments as (start, end, length), with 0 <= start < end <= len(seq)
    """
    n = len(seq)
    cut_positions = [0] + sorted(find_cut_positions(seq, recog, cut_index)) + [n]

    fragments = []
    for i in range(len(cut_positions) - 1):
        start = cut_positions[i]
        end = cut_positions[i + 1]
        fragments.append((start, end, end - start))

    # internal cuts only (no 0 and n)
    internal_cuts = cut_positions[1:-1]
    return internal_cuts, fragments


def pretty_print_digest(enzyme_name: str, cuts: List[int], fragments: List[Tuple[int, int, int]]) -> List[int]:
    """
    Print info for one enzyme.
    Returns the list of fragment lengths.
    """
    print(f"\n==== {enzyme_name} ====")
    print(f"Number of cuts: {len(cuts)}")

    # positions as 1-based for humans
    cuts_1based = [c + 1 for c in cuts]
    print("Cleavage positions (1-based):", cuts_1based if cuts_1based else "No cuts")

    print("Fragments (1-based inclusive ranges):")
    fragment_lengths = []
    for start, end, length in fragments:
        # start: 0-based, end: 0-based exclusive -> convert to 1-based inclusive
        start_1 = start + 1
        end_1 = end
        fragment_lengths.append(length)
        print(f"  {start_1:5d}-{end_1:5d}   length = {length}")

    return fragment_lengths


# ---------- STEP 4: Simple gel simulation with matplotlib ----------

def plot_gel(enzymes_fragments: List[Tuple[str, List[int]]]) -> None:
    """
    Draw a simple fake gel:
      - one lane per enzyme
      - smaller fragments migrate further down
      - y-axis shows approximate fragment size (bp)
    """
    # Collect all fragment lengths to scale the positions
    all_lengths = [L for _, frags in enzymes_fragments for L in frags]
    max_len = max(all_lengths)
    min_len = min(all_lengths)
    span = max_len - min_len if max_len != min_len else 1

    fig, ax = plt.subplots(figsize=(6, 6))
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")

    # Layout
    lane_count = len(enzymes_fragments)
    x_positions = [i * 2 for i in range(lane_count)]

    ax.set_xlim(-1, (lane_count - 1) * 2 + 1)
    ax.set_ylim(0, 1)
    ax.invert_yaxis()  # top = wells, bottom = migrated fragments

    # Draw lanes
    for lane_index, (name, fragment_lengths) in enumerate(enzymes_fragments):
        x_center = x_positions[lane_index]

        # "Well" at the top
        ax.add_line(Line2D([x_center - 0.4, x_center + 0.4],
                           [0.02, 0.02],
                           linewidth=8,
                           color="grey"))

        # Bands
        for length in fragment_lengths:
            # map length to y-position: big fragments stay near the top, small go down
            y = (max_len - length) / span * 0.9 + 0.05  # between 0.05 and 0.95
            ax.add_line(Line2D([x_center - 0.3, x_center + 0.3],
                               [y, y],
                               linewidth=4,
                               color="white"))

    # ----- X axis: enzyme names under each lane -----
    ax.set_xticks(x_positions)
    ax.set_xticklabels(
        [name for name, _ in enzymes_fragments],
        color="white",
        fontsize=9
    )

    # ----- Y axis: rough fragment sizes (bp) -----
    # choose ~5 ticks between min and max
    step = max_len // 5 or 1
    tick_lengths = list(range(step, max_len + 1, step))
    # keep a few evenly spaced values
    if len(tick_lengths) > 6:
        tick_lengths = [tick_lengths[i] for i in range(0, len(tick_lengths), len(tick_lengths) // 5)]

    y_positions = [(max_len - L) / span * 0.9 + 0.05 for L in tick_lengths]

    ax.set_yticks(y_positions)
    ax.set_yticklabels([str(L) for L in tick_lengths], color="white", fontsize=8)
    ax.set_ylabel("Fragment size (bp)", color="white")

    # Cosmetics: remove frame, set tick colors, title
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.tick_params(axis="x", colors="white")
    ax.tick_params(axis="y", colors="white")

    ax.set_title(
        "Simulated DNA gel (shorter fragments run lower)",
        color="white",
        fontsize=10
    )

    plt.tight_layout()
    plt.savefig("gel_simulation_lab09.png", dpi=300)
    plt.show()



# ---------- STEP 5: Main ----------

def main():
    if len(sys.argv) != 2:
        print("Usage: python labs_python/lab09/ex1_restriction_digest.py path/to/sequence.fasta")
        sys.exit(1)

    fasta_path = sys.argv[1]
    dna_seq = read_fasta_sequence(fasta_path)

    print(f"Sequence length: {len(dna_seq)} bases")

    enzymes_fragments = []

    for name, info in ENZYMES.items():
        cuts, fragments = digest_sequence(dna_seq, info["site"], info["cut_index"])
        fragment_lengths = pretty_print_digest(name, cuts, fragments)
        enzymes_fragments.append((name, fragment_lengths))

    # draw the gel
    plot_gel(enzymes_fragments)


if __name__ == "__main__":
    main()
