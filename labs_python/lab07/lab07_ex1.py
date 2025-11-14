"""
Lab 07 - Exercise 1
-------------------

1. Read a DNA sequence from an NCBI FASTA file (1000–3000 bp).
2. Generate 2000 random reads of length 100–150 from this sequence.
3. Rebuild (assemble) the original sequence with a simple greedy algorithm.
"""

from pathlib import Path
import random


# -------- Configuration (change if needed) --------

FASTA_FILENAME = "ncbi_sequence.fasta"  # file in labs_python/lab07/data/
N_READS = 2000
MIN_READ_LEN = 100
MAX_READ_LEN = 150
MIN_OVERLAP = 30  # minimum overlap length when merging


# -------- Helper functions --------

def read_fasta(path: Path) -> str:
    """
    Read a FASTA file and return the concatenated DNA sequence as one string.
    Header lines (starting with '>') are ignored.
    """
    sequence_parts = []

    with path.open("r") as f:
        for line in f:
            line = line.strip()

            # skip empty lines
            if not line:
                continue

            # skip header
            if line.startswith(">"):
                continue

            # add sequence (uppercased)
            sequence_parts.append(line.upper())

    sequence = "".join(sequence_parts)
    return sequence


def generate_reads(sequence: str,
                   n_reads: int,
                   min_len: int,
                   max_len: int) -> list[str]:
    """
    Generate n_reads random substrings ("reads") from the given sequence.
    Each read has a random length between min_len and max_len.
    Sampling is done with replacement (a position can be used multiple times).
    """
    reads = []
    seq_len = len(sequence)

    if seq_len < min_len:
        raise ValueError("Sequence is shorter than the minimum read length.")

    for _ in range(n_reads):
        read_len = random.randint(min_len, max_len)

        # ensure we don't go out of bounds
        max_start = seq_len - read_len
        start = random.randint(0, max_start)

        read = sequence[start:start + read_len]
        reads.append(read)

    return reads


def overlap(a: str, b: str, min_length: int) -> int:
    """
    Return the length of the longest suffix of 'a' that matches
    a prefix of 'b', with length at least min_length.
    This version uses str.find and is much faster than the naive O(L^2) one.
    """
    start = 0
    while True:
        # look for the first chunk of b (prefix of length min_length)
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0  # no possible overlap

        # check if the rest matches
        if b.startswith(a[start:]):
            return len(a) - start

        # try again from next position
        start += 1


def assemble_by_successor(reads: list[str], min_overlap: int) -> str:
    """
    Faster O(n^2) assembly:

    1. For each read i, find the best successor j (max overlap(i, j)).
    2. Build a chain starting from a read with no incoming edges.
    3. Append any leftover reads separated by Ns.

    Not perfect, but good enough for the lab and much faster than the
    full n^3 greedy merge.
    """
    n = len(reads)
    if n == 0:
        return ""

    # For each read: (best_successor_index, overlap_length)
    successors: list[tuple[int | None, int]] = [(None, 0)] * n
    in_degree = [0] * n

    print("Precomputing best overlaps for each read...")

    for i in range(n):
        best_j = None
        best_olen = 0

        for j in range(n):
            if i == j:
                continue

            olen = overlap(reads[i], reads[j], min_overlap)
            if olen > best_olen:
                best_olen = olen
                best_j = j

        successors[i] = (best_j, best_olen)
        if best_j is not None and best_olen >= min_overlap:
            in_degree[best_j] += 1

        if (i + 1) % 100 == 0 or i == n - 1:
            print(f"  processed {i + 1}/{n} reads")

    # choose start: a read that is not a best_successor of anyone
    start = 0
    for i in range(n):
        if in_degree[i] == 0:
            start = i
            break

    print(f"\nStarting assembly chain from read index {start}.")

    used = set()
    current = start
    assembled = reads[current]
    used.add(current)

    while True:
        next_idx, olen = successors[current]
        if next_idx is None or olen < min_overlap or next_idx in used:
            break

        assembled += reads[next_idx][olen:]
        used.add(next_idx)
        current = next_idx

    # attach leftover reads (if any) so we don't lose them completely
    # separated by a run of 'N's (unknown bases)
    for i in range(n):
        if i not in used:
            assembled += "N" * min_overlap + reads[i]

    return assembled


# -------- Main script --------

def main():
    base_dir = Path(__file__).parent
    data_dir = base_dir / "data"
    fasta_path = data_dir / FASTA_FILENAME

    if not fasta_path.exists():
        raise FileNotFoundError(
            f"FASTA file not found at {fasta_path}. "
            f"Either create it or change FASTA_FILENAME."
        )

    # 1) Read original sequence
    sequence = read_fasta(fasta_path)
    print(f"Original sequence length: {len(sequence)} bases")

    if not (1000 <= len(sequence) <= 3000):
        print("WARNING: Sequence length is not in [1000, 3000] as requested.")

    # 2) Generate reads
    random.seed(42)  # for reproducible results
    reads = generate_reads(sequence, N_READS, MIN_READ_LEN, MAX_READ_LEN)
    print(f"Generated {len(reads)} reads.")
    print(f"Example read (first 60 bases): {reads[0][:60]}")

    # Optional: save reads to a file for inspection
    data_dir.mkdir(exist_ok=True)
    reads_file = data_dir / "reads.txt"
    with reads_file.open("w") as f:
        for r in reads:
            f.write(r + "\n")
    print(f"Reads saved to: {reads_file}")

    # 3) Assemble using faster algorithm
    print("\nStarting assembly (successor-based)...")
    assembled = assemble_by_successor(reads, MIN_OVERLAP)
    print("\nAssembly finished.")
    print(f"Assembled length: {len(assembled)} bases")

    # Optional: save assembled sequence
    assembled_file = data_dir / "assembled_sequence.txt"
    with assembled_file.open("w") as f:
        f.write(assembled)
    print(f"Assembled sequence saved to: {assembled_file}")

    # 4) Check how good the assembly is
    contains_original = sequence in assembled
    print(f"\nDoes the assembled sequence contain the original exactly? {contains_original}")

    # quick comparison of start/end
    print("\nOriginal (first 100 bases):")
    print(sequence[:100])
    print("\nAssembled (first 100 bases):")
    print(assembled[:100])
    print("\nOriginal (last 100 bases):")
    print(sequence[-100:])
    print("\nAssembled (last 100 bases):")
    print(assembled[-100:])


if __name__ == "__main__":
    main()
