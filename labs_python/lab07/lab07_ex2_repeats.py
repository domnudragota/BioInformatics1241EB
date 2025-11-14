"""
Lab 07 - Exercise 2
-------------------

Detect repetitions of length 6–10 bases in the DNA sequence
we downloaded from NCBI in Exercise 1.

A "repetition" here means: a substring that appears at least
twice in the same sequence.
"""

from pathlib import Path

# configuration
FASTA_FILENAME = "ncbi_sequence.fasta"  # same file used in ex1
MIN_REPEAT_LEN = 6
MAX_REPEAT_LEN = 10
MIN_OCCURRENCES = 2  # how many times a motif must appear to be considered a repeat


def read_fasta(path: Path) -> str:
    """
    Read a FASTA file and return the concatenated DNA sequence as one string.
    Header lines (starting with '>') are ignored.
    """
    sequence_parts = []

    with path.open("r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue  # skip empty lines

            if line.startswith(">"):
                continue  # skip header

            sequence_parts.append(line.upper())

    return "".join(sequence_parts)


def find_repeats(sequence: str,
                 min_len: int,
                 max_len: int,
                 min_occurrences: int):
    """
    Find repeated substrings between min_len and max_len (inclusive).

    Returns a dictionary:
        { length_L: { motif: [pos1, pos2, ...], ... }, ... }

    where each motif appears at least min_occurrences times.
    """
    result = {}

    seq_len = len(sequence)

    for L in range(min_len, max_len + 1):
        motifs = {}  # motif -> list of start positions

        # slide a window of length L across the sequence
        for i in range(0, seq_len - L + 1):
            fragment = sequence[i:i + L]

            if fragment not in motifs:
                motifs[fragment] = [i]
            else:
                motifs[fragment].append(i)

        # keep only those motifs that occur at least min_occurrences times
        filtered = {
            motif: positions
            for motif, positions in motifs.items()
            if len(positions) >= min_occurrences
        }

        result[L] = filtered

    return result


def save_repeats_to_file(repeats, sequence_length: int, output_path: Path):
    """
    Save a human-readable report of the repeats to a text file.
    """
    with output_path.open("w") as f:
        f.write(f"DNA sequence length: {sequence_length} bases\n")
        f.write(
            f"Motifs of length {MIN_REPEAT_LEN}–{MAX_REPEAT_LEN} "
            f"appearing at least {MIN_OCCURRENCES} times.\n\n"
        )

        for L in sorted(repeats.keys()):
            motifs = repeats[L]
            f.write(f"=== Motifs of length {L} ===\n")
            if not motifs:
                f.write("No repetitions found for this length.\n\n")
                continue

            # sort motifs by descending number of occurrences
            sorted_items = sorted(
                motifs.items(),
                key=lambda item: len(item[1]),
                reverse=True,
            )

            for motif, positions in sorted_items:
                f.write(
                    f"Motif: {motif} | "
                    f"occurrences: {len(positions)} | "
                    f"positions (0-based): {positions}\n"
                )
            f.write("\n")


def main():
    base_dir = Path(__file__).parent
    data_dir = base_dir / "data"
    fasta_path = data_dir / FASTA_FILENAME

    if not fasta_path.exists():
        raise FileNotFoundError(
            f"FASTA file not found at {fasta_path}. "
            f"Either create it or change FASTA_FILENAME."
        )

    # 1) read the sequence
    sequence = read_fasta(fasta_path)
    print(f"Sequence length: {len(sequence)} bases")

    # 2) find repeats
    print(
        f"Searching for repeats of length {MIN_REPEAT_LEN}–{MAX_REPEAT_LEN} "
        f"appearing at least {MIN_OCCURRENCES} times..."
    )

    repeats = find_repeats(
        sequence,
        MIN_REPEAT_LEN,
        MAX_REPEAT_LEN,
        MIN_OCCURRENCES,
    )

    # 3) show a short summary in terminal
    for L in range(MIN_REPEAT_LEN, MAX_REPEAT_LEN + 1):
        count_motifs = len(repeats[L])
        print(f"Length {L}: {count_motifs} distinct repeated motifs found.")

    # 4) save detailed report
    output_path = data_dir / "repeats_6_10bp.txt"
    save_repeats_to_file(repeats, len(sequence), output_path)
    print(f"Detailed repeats report saved to: {output_path}")


if __name__ == "__main__":
    main()
