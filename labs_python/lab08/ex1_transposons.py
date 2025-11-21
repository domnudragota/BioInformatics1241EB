from pathlib import Path
import random
import argparse
from typing import List, Tuple

NUCLEOTIDES = "ACGT"


def random_dna(length: int) -> str:
    """Generate a random DNA string of given length."""
    return "".join(random.choice(NUCLEOTIDES) for _ in range(length))


def generate_te_families(num_families: int) -> List[str]:
    """Create num_families random TE motifs (8–15 bp each)."""
    families = []
    for _ in range(num_families):
        motif_len = random.randint(8, 15)
        families.append(random_dna(motif_len))
    return families


def simulate_sequence(
    target_min: int = 200,
    target_max: int = 400,
    min_families: int = 3,
    max_families: int = 4,
) -> Tuple[str, List[str], List[int]]:
    """
    Build one artificial DNA sequence:
    - total length in [target_min, target_max]
    - 3–4 TE families
    - each family repeated 2–3 times
    """
    num_families = random.randint(min_families, max_families)

    while True:
        # pick TE motifs and how many copies of each
        te_families = generate_te_families(num_families)
        copies = [random.randint(2, 3) for _ in range(num_families)]
        total_te_len = sum(len(m) * c for m, c in zip(te_families, copies))

        target_len = random.randint(target_min, target_max)

        # require at least 50 bp of non-TE background so it doesn't become “all repeats”
        if total_te_len <= target_len - 50:
            break

    background_len = target_len - total_te_len
    seq = random_dna(background_len)

    # insert each TE copy at a random position inside the background
    for motif, n_copies in zip(te_families, copies):
        for _ in range(n_copies):
            pos = random.randint(0, len(seq))
            seq = seq[:pos] + motif + seq[pos:]

    assert target_min <= len(seq) <= target_max
    return seq, te_families, copies


def find_all_positions(seq: str, motif: str) -> List[int]:
    """Find all occurrences of motif in seq (simple substring search)."""
    positions = []
    start = 0
    while True:
        idx = seq.find(motif, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def write_fasta(seq: str, path: Path, header: str = "artificial_TE_sequence") -> None:
    """Write sequence to FASTA, wrapped at 60 bp/line."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + "\n")


def write_annotations(
    seq: str,
    te_families: List[str],
    copies: List[int],
    path: Path,
) -> None:
    """
    Save a small TSV with TE positions:
    family_id, motif, copy_index, start, end
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("family_id\tmotif\tcopy_index\tstart\tend\n")
        for i, motif in enumerate(te_families):
            family_id = f"TE{i + 1}"
            positions = find_all_positions(seq, motif)
            # keep only as many as we intended to insert (2–3 per family)
            for copy_index, start in enumerate(positions[:copies[i]], start=1):
                end = start + len(motif)
                f.write(f"{family_id}\t{motif}\t{copy_index}\t{start}\t{end}\n")


def main():
    here = Path(__file__).resolve().parent
    default_out = here / "data" / "ex1_artificial_te.fasta"
    default_ann = here / "data" / "ex1_artificial_te_annotations.tsv"

    parser = argparse.ArgumentParser(
        description="Exercise 1: simulate DNA with transposable elements."
    )
    parser.add_argument(
        "-o",
        "--output-fasta",
        type=Path,
        default=default_out,
        help=f"Output FASTA path (default: {default_out})",
    )
    parser.add_argument(
        "--output-ann",
        type=Path,
        default=default_ann,
        help=f"Output TSV with TE positions (default: {default_ann})",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )
    parser.add_argument(
        "--min-len",
        type=int,
        default=200,
        help="Minimum sequence length (default: 200)",
    )
    parser.add_argument(
        "--max-len",
        type=int,
        default=400,
        help="Maximum sequence length (default: 400)",
    )
    args = parser.parse_args()

    random.seed(args.seed)

    seq, te_families, copies = simulate_sequence(
        target_min=args.min_len,
        target_max=args.max_len,
    )

    header = (
        f"artificial_TE_sequence len={len(seq)} "
        f"families={len(te_families)} "
        f"copies={copies}"
    )
    write_fasta(seq, args.output_fasta, header=header)
    write_annotations(seq, te_families, copies, args.output_ann)

    print("Generated artificial DNA sequence with transposable elements")
    print(f"  Length: {len(seq)} bp")
    print("  TE families (3–4 simulated):")
    for i, (motif, n_copies) in enumerate(zip(te_families, copies), start=1):
        print(f"    TE{i}: motif={motif} (len={len(motif)}), copies={n_copies}")
    print(f"FASTA saved to      : {args.output_fasta}")
    print(f"Annotations saved to: {args.output_ann}")


if __name__ == "__main__":
    main()
