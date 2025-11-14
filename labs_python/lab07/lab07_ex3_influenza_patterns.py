"""
Lab 07 - Exercise: Influenza strains
------------------------------------

For each Influenza genome (FASTA file) in data/influenza:

- detect all motifs with length 6–10 bases,
- count how many times each motif appears,
- select the top 20 most frequent motifs,
- plot a bar chart with those 20 motifs (one chart per strain).
"""

from pathlib import Path
import matplotlib.pyplot as plt

# --- configuration ---

MIN_LEN = 6          # minimum motif length
MAX_LEN = 10         # maximum motif length
TOP_K = 20           # how many motifs to show in the chart

# folder with all influenza FASTA files
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
INFLUENZA_DIR = DATA_DIR / "influenza"


# --- helpers ---

def read_fasta(path: Path) -> str:
    """
    Read a FASTA file and return the concatenated DNA sequence.
    Header lines (starting with '>') are ignored.
    """
    parts = []
    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            parts.append(line.upper())
    return "".join(parts)


def count_motifs(sequence: str, min_len: int, max_len: int) -> dict[str, int]:
    """
    Count all motifs of length min_len..max_len in the sequence.

    Returns a dictionary: motif -> count.
    (We don't store positions here, just total counts.)
    """
    counts: dict[str, int] = {}
    n = len(sequence)

    for L in range(min_len, max_len + 1):
        if n < L:
            continue

        for i in range(0, n - L + 1):
            motif = sequence[i:i + L]
            counts[motif] = counts.get(motif, 0) + 1

    return counts


def plot_top_motifs(motif_counts: dict[str, int],
                    title: str,
                    out_path: Path,
                    top_k: int = TOP_K):
    """
    Plot a bar chart with the top_k motifs by count.
    Save figure to out_path as PNG.
    """
    # sort motifs by count (descending), then by motif string
    sorted_items = sorted(
        motif_counts.items(),
        key=lambda kv: (-kv[1], kv[0])
    )

    top_items = sorted_items[:top_k]

    motifs = [m for m, _ in top_items]
    counts = [c for _, c in top_items]

    # create figure
    plt.figure(figsize=(12, 6))
    x = range(len(motifs))
    plt.bar(x, counts)
    plt.xticks(x, motifs, rotation=45, ha="right")
    plt.xlabel("Motif (6–10 bp)")
    plt.ylabel("Count in genome")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def find_fasta_files(directory: Path) -> list[Path]:
    """
    Return a sorted list of FASTA-like files in the directory.
    We look for several common extensions.
    """
    patterns = ("*.fasta", "*.fa", "*.fna", "*.ffn")
    files: set[Path] = set()
    for pat in patterns:
        files.update(directory.glob(pat))
    return sorted(files)


# --- main ---

def main():
    if not INFLUENZA_DIR.exists():
        raise FileNotFoundError(
            f"Folder with influenza genomes not found: {INFLUENZA_DIR}\n"
            f"Create it and put your 10 FASTA files inside."
        )

    fasta_files = find_fasta_files(INFLUENZA_DIR)
    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA files (*.fasta, *.fa, *.fna, *.ffn) found in {INFLUENZA_DIR}"
        )

    print(f"Found {len(fasta_files)} influenza genomes:")
    for fp in fasta_files:
        print("  -", fp.name)

    # folder for plots
    plots_dir = INFLUENZA_DIR / "plots"
    plots_dir.mkdir(exist_ok=True)

    for fasta_path in fasta_files:
        print(f"\nProcessing {fasta_path.name} ...")

        seq = read_fasta(fasta_path)
        print(f"  sequence length: {len(seq)} bases")

        motif_counts = count_motifs(seq, MIN_LEN, MAX_LEN)
        print(f"  distinct motifs (6–10 bp): {len(motif_counts)}")

        # save top 20 as a chart
        strain_name = fasta_path.stem
        out_png = plots_dir / f"{strain_name}_top{TOP_K}_motifs.png"

        plot_top_motifs(motif_counts, title=strain_name, out_path=out_png)

        print(f"  chart saved to: {out_png}")

    print("\nDone. One PNG chart per strain was created in:", plots_dir)


if __name__ == "__main__":
    main()
