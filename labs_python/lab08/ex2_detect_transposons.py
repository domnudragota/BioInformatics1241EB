from pathlib import Path
import argparse
from typing import Dict, List, Set

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def read_fasta(path: Path) -> str:
    """Read a (single-sequence) FASTA file and return the sequence as a string."""
    seq_lines: List[str] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # header line, skip
                continue
            seq_lines.append(line)
    return "".join(seq_lines).upper()


def load_motifs_from_annotations(path: Path) -> Dict[str, Set[str]]:
    """
    Load TE motifs from the annotations TSV created in Exercise 1.

    Expected header (from ex1):
        family_id  motif  copy_index  start  end

    We only care about family_id and motif, and we deduplicate motifs per family.
    """
    motifs_by_family: Dict[str, Set[str]] = {}

    with path.open() as f:
        header = f.readline().strip().split("\t")
        if not header:
            raise ValueError(f"Empty or invalid annotations file: {path}")

        try:
            family_idx = header.index("family_id")
            motif_idx = header.index("motif")
        except ValueError:
            raise ValueError(
                f"Annotations file {path} must have 'family_id' and 'motif' columns"
            )

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= max(family_idx, motif_idx):
                continue

            family_id = parts[family_idx]
            motif = parts[motif_idx].upper()

            if family_id not in motifs_by_family:
                motifs_by_family[family_id] = set()
            motifs_by_family[family_id].add(motif)

    return motifs_by_family


def find_all_positions(seq: str, motif: str) -> List[int]:
    """
    Find all (possibly overlapping) occurrences of motif in seq.
    Returns 0-based starting indices.
    """
    positions: List[int] = []
    start = 0
    while True:
        idx = seq.find(motif, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1  # move one base forward to allow overlapping hits if any
    return positions


def detect_transposons(
    seq: str, motifs_by_family: Dict[str, Set[str]]
) -> List[dict]:
    """
    Detect positions of each TE motif in the sequence.

    Returns a list of dictionaries with:
        - family_id
        - motif
        - start_0
        - end_0
        - start_1
        - end_1
    """
    detections: List[dict] = []

    for family_id, motifs in motifs_by_family.items():
        for motif in motifs:
            hits = find_all_positions(seq, motif)
            for start0 in hits:
                end0 = start0 + len(motif)
                detections.append(
                    {
                        "family_id": family_id,
                        "motif": motif,
                        "start_0": start0,         # 0-based
                        "end_0": end0,             # 0-based exclusive
                        "start_1": start0 + 1,     # 1-based (more “bio” style)
                        "end_1": end0,             # 1-based inclusive
                    }
                )

    # sort by genomic position
    detections.sort(key=lambda d: d["start_0"])
    return detections


def write_detections(detections: List[dict], path: Path) -> None:
    """
    Save detections to TSV.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(
            "family_id\tmotif\tstart_0_based\tend_0_based_exclusive\t"
            "start_1_based\tend_1_based\n"
        )
        for d in detections:
            f.write(
                f"{d['family_id']}\t{d['motif']}\t"
                f"{d['start_0']}\t{d['end_0']}\t"
                f"{d['start_1']}\t{d['end_1']}\n"
            )


def plot_transposon_map(seq_len: int, detections: List[dict], out_path: Path) -> None:
    """
    Plot a simple map of the artificial gene with transposons as boxes.
    """
    if not detections:
        return

    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 2))

    # Base "gene" line
    ax.hlines(y=0.5, xmin=0, xmax=seq_len, linewidth=2)

    # Draw each TE as a rectangle on top of the line.
    for d in detections:
        start = d["start_0"]
        length = d["end_0"] - d["start_0"]
        rect = Rectangle((start, 0.25), width=length, height=0.5)
        ax.add_patch(rect)
        ax.text(
            start + length / 2,
            0.5,
            d["family_id"],
            ha="center",
            va="center",
            fontsize=8,
        )

    ax.set_ylim(0, 1)
    ax.set_xlim(-5, seq_len + 5)
    ax.set_xlabel("Position along gene (bp)")
    ax.set_yticks([])
    ax.set_title("Map of transposable elements in artificial gene")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]


def plot_example_insertion(
    seq: str,
    detection: dict,
    flank_size: int,
    out_path: Path,
) -> None:
    """
    Plot a conceptual diagram of one TE insertion:

    - New gene (with TE): left flank + TE + right flank
    - Original gene (without TE): left flank + right flank ("gap filled")
    - Show candidate direct repeats (small k-mers flanking the TE)
    - Show candidate inverted repeats inside the TE
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)

    start = detection["start_0"]
    end = detection["end_0"]
    family_id = detection["family_id"]
    motif = detection["motif"]

    # Window around the insertion site
    left_start = max(0, start - flank_size)
    right_end = min(len(seq), end + flank_size)

    left_flank = seq[left_start:start]
    te_seq = seq[start:end]
    right_flank = seq[end:right_end]

    # Direct repeat candidates: short segments immediately flanking TE
    dr_size = min(4, len(left_flank), len(right_flank))
    if dr_size > 0:
        left_dr = left_flank[-dr_size:]
        right_dr = right_flank[:dr_size]
    else:
        left_dr = ""
        right_dr = ""

    # Inverted repeat candidates: short segments at TE ends
    ir_size = min(4, len(te_seq) // 2)
    if ir_size > 0:
        ir_left = te_seq[:ir_size]
        ir_right = te_seq[-ir_size:]
        ir_right_rc = revcomp(ir_right)
    else:
        ir_left = ""
        ir_right = ""
        ir_right_rc = ""

    # For the schematic we use a simple relative coordinate system
    total_span = len(left_flank) + len(te_seq) + len(right_flank)

    fig, ax = plt.subplots(figsize=(10, 3))

    # New gene (with TE)
    ax.text(
        -1,
        1.0,
        "New gene (with TE)",
        ha="right",
        va="center",
        fontsize=9,
    )
    # Left flank
    ax.add_patch(Rectangle((0, 0.9), len(left_flank), 0.15))
    # TE
    ax.add_patch(Rectangle((len(left_flank), 0.9), len(te_seq), 0.15))
    # Right flank
    ax.add_patch(
        Rectangle(
            (len(left_flank) + len(te_seq), 0.9),
            len(right_flank),
            0.15,
        )
    )
    ax.text(
        len(left_flank) + len(te_seq) / 2,
        1.07,
        f"{family_id} (TE)",
        ha="center",
        va="bottom",
        fontsize=8,
    )

    # Original gene (without TE) – TE removed, flanks joined
    ax.text(
        -1,
        0.6,
        "Original gene",
        ha="right",
        va="center",
        fontsize=9,
    )
    ax.add_patch(Rectangle((0, 0.55), len(left_flank), 0.15))
    ax.add_patch(Rectangle((len(left_flank), 0.55), len(right_flank), 0.15))
    ax.text(
        len(left_flank),
        0.75,
        "gap filled",
        ha="center",
        va="bottom",
        fontsize=8,
    )

    ax.set_ylim(0.4, 1.2)
    ax.set_xlim(0, total_span)
    ax.set_yticks([])
    ax.set_xlabel("Relative position (bp)")
    ax.set_title("Example of transposon insertion vs original gene")

    # Info text: flanks, direct repeat candidates, inverted repeat candidates
    info_lines: List[str] = [
        f"TE motif: {motif}",
        f"Left flank (len {len(left_flank)}): {left_flank}",
        f"Right flank (len {len(right_flank)}): {right_flank}",
    ]

    if dr_size > 0:
        info_lines.append(
            f"Direct repeat candidates ({dr_size} bp): "
            f"left='{left_dr}', right='{right_dr}'"
        )

    if ir_size > 0:
        info_lines.append(
            f"Inverted repeat candidates in TE ({ir_size} bp): "
            f"IR-L='{ir_left}', IR-R='{ir_right}', IR-R (revcomp)='{ir_right_rc}'"
        )

    fig.text(
        0.01,
        0.02,
        "\n".join(info_lines),
        fontsize=7,
        va="bottom",
        ha="left",
    )

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    here = Path(__file__).resolve().parent

    default_fasta = here / "data" / "ex1_artificial_te.fasta"
    default_ann = here / "data" / "ex1_artificial_te_annotations.tsv"
    default_out = here / "data" / "ex2_detected_transposons.tsv"
    default_map = here / "data" / "ex2_te_map.png"
    default_example = here / "data" / "ex2_example_insertion.png"

    parser = argparse.ArgumentParser(
        description=(
            "Exercise 2: detect positions of transposable elements "
            "inside the artificial DNA sequence and generate simple diagrams."
        )
    )
    parser.add_argument(
        "-i",
        "--input-fasta",
        type=Path,
        default=default_fasta,
        help=f"FASTA file with artificial DNA (default: {default_fasta})",
    )
    parser.add_argument(
        "-a",
        "--annotations",
        type=Path,
        default=default_ann,
        help=(
            "Annotations TSV from Exercise 1, from which TE motifs are read "
            f"(default: {default_ann})"
        ),
    )
    parser.add_argument(
        "-o",
        "--output-tsv",
        type=Path,
        default=default_out,
        help=f"Output TSV with detected TE positions (default: {default_out})",
    )
    parser.add_argument(
        "--map-png",
        type=Path,
        default=default_map,
        help=f"Output PNG with gene + TE map (default: {default_map})",
    )
    parser.add_argument(
        "--example-png",
        type=Path,
        default=default_example,
        help=(
            "Output PNG with example insertion (original vs new gene, "
            f"direct/inverted repeats) (default: {default_example})"
        ),
    )
    parser.add_argument(
        "--flank-size",
        type=int,
        default=15,
        help="Number of bases to show on each side of the TE in the example diagram.",
    )
    args = parser.parse_args()

    if not args.input_fasta.exists():
        raise SystemExit(f"ERROR: FASTA file not found: {args.input_fasta}")
    if not args.annotations.exists():
        raise SystemExit(f"ERROR: annotations file not found: {args.annotations}")

    seq = read_fasta(args.input_fasta)
    motifs_by_family = load_motifs_from_annotations(args.annotations)
    detections = detect_transposons(seq, motifs_by_family)
    write_detections(detections, args.output_tsv)

    print(f"Sequence length: {len(seq)} bp")
    print("Loaded TE motifs per family:")
    for family_id, motifs in motifs_by_family.items():
        print(f"  {family_id}: {len(motifs)} motif(s)")

    print(f"Total detected TE occurrences: {len(detections)}")
    if detections:
        first = detections[0]
        print(
            f"First hit: {first['family_id']} ({first['motif']}) "
            f"at {first['start_1']}–{first['end_1']} (1-based)"
        )

    # Diagrams
    if detections:
        print(f"Generating TE map image: {args.map_png}")
        plot_transposon_map(len(seq), detections, args.map_png)

        print(f"Generating example insertion diagram: {args.example_png}")
        plot_example_insertion(
            seq,
            detections[0],
            flank_size=args.flank_size,
            out_path=args.example_png,
        )
    else:
        print("No detections -> no diagrams generated.")

    print(f"Detections saved to: {args.output_tsv}")


if __name__ == "__main__":
    main()
