from pathlib import Path
import argparse
from typing import Dict, List, Tuple


def revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]


def read_fasta_multi(path: Path) -> List[Tuple[str, str]]:
    """
    Read a FASTA file that may contain multiple sequences (contigs).
    Returns a list of (header, sequence).
    """
    records: List[Tuple[str, str]] = []
    header: str | None = None
    seq_lines: List[str] = []

    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # store previous record
                if header is not None:
                    records.append((header, "".join(seq_lines).upper()))
                    seq_lines = []
                # keep only first token as contig id (simpler)
                header = line[1:].split()[0]
            else:
                seq_lines.append(line)

    if header is not None:
        records.append((header, "".join(seq_lines).upper()))

    return records


def build_kmer_positions(
    seq: str,
    min_ir_len: int = 5,
    max_ir_len: int = 6,
) -> Dict[int, Dict[str, List[int]]]:
    """
    Index all k-mers of length in [min_ir_len, max_ir_len] along the sequence.

    Returns: {k: {kmer: [positions...]}}
    """
    kmer_positions: Dict[int, Dict[str, List[int]]] = {}
    n = len(seq)

    for k in range(min_ir_len, max_ir_len + 1):
        inner: Dict[str, List[int]] = {}
        if n < k:
            kmer_positions[k] = inner
            continue

        for i in range(0, n - k + 1):
            kmer = seq[i : i + k]
            if "N" in kmer:  # skip ambiguous bases
                continue
            if kmer not in inner:
                inner[kmer] = []
            inner[kmer].append(i)

        kmer_positions[k] = inner

    return kmer_positions


def find_te_candidates(
    seq: str,
    min_ir_len: int = 5,
    max_ir_len: int = 6,
    min_te_len: int = 100,
    max_te_len: int = 5000,
) -> List[dict]:
    """
    Detect candidate TEs as segments bounded by terminal inverted repeats (TIRs)
    of length 5–6 bases.

    A TE is defined here as:
        TE = [start, end) where seq[start:start+k] and seq[j:j+k] are inverted
        repeats (reverse complements), with:
            - k in [min_ir_len, max_ir_len]
            - min_te_len <= end - start <= max_te_len
    """
    kmer_positions = build_kmer_positions(seq, min_ir_len, max_ir_len)
    candidates: List[dict] = []
    seen: set[tuple[int, int, int]] = set()  # (start, end, ir_len)

    for k, inner in kmer_positions.items():
        if not inner:
            continue

        # precompute rc -> positions as needed
        for motif, pos_list in inner.items():
            rc = revcomp(motif)
            if rc not in inner:
                continue

            rc_positions = inner[rc]

            # Two-pointer scan so we don't do O(n^2) blindly
            j_start = 0
            for i_pos in pos_list:
                # we only want rc_positions strictly after i_pos
                while j_start < len(rc_positions) and rc_positions[j_start] <= i_pos:
                    j_start += 1

                j_idx = j_start
                while j_idx < len(rc_positions):
                    j_pos = rc_positions[j_idx]
                    te_start = i_pos
                    te_end = j_pos + k
                    te_len = te_end - te_start

                    if te_len < min_te_len:
                        j_idx += 1
                        continue
                    if te_len > max_te_len:
                        # rc_positions is sorted; further j_pos will only be larger
                        break

                    key = (te_start, te_end, k)
                    if key not in seen:
                        seen.add(key)
                        candidates.append(
                            {
                                "start": te_start,
                                "end": te_end,
                                "length": te_len,
                                "ir_len": k,
                                "left_ir_seq": motif,
                                "right_ir_seq": rc,
                            }
                        )

                    j_idx += 1

    # sort by genomic coordinates
    candidates.sort(key=lambda d: (d["start"], d["end"]))
    return candidates


def classify_embedding_overlaps(tes: List[dict]) -> None:
    """
    Given a list of TE dicts with 'start' and 'end', annotate them with:
      - 'id'          : TE index (1..N per contig)
      - 'embedded_in' : set of IDs that contain this TE fully
      - 'contains'    : set of IDs fully contained in this TE
      - 'overlaps'    : set of IDs that partially overlap
      - 'relation'    : "simple" / "embedded" / "container" / "overlap" / combos
    """
    # assign IDs
    for idx, te in enumerate(tes):
        te["id"] = idx + 1
        te["embedded_in"] = set()
        te["contains"] = set()
        te["overlaps"] = set()

    # already sorted by start, then end
    n = len(tes)
    for i in range(n):
        a = tes[i]
        for j in range(i + 1, n):
            b = tes[j]
            if b["start"] >= a["end"]:
                # beyond the region where 'a' can overlap anyone
                break

            # We know: b.start < a.end  => they overlap in some way on x-axis
            if b["end"] <= a["end"]:
                # b fully inside a -> embedding
                b["embedded_in"].add(a["id"])
                a["contains"].add(b["id"])
            elif b["start"] <= a["start"] and b["end"] >= a["end"]:
                # a fully inside b (less common due to sorting, but keep it)
                a["embedded_in"].add(b["id"])
                b["contains"].add(a["id"])
            else:
                # partial overlap (neither fully contains the other)
                a["overlaps"].add(b["id"])
                b["overlaps"].add(a["id"])

    # summarize relation label
    for te in tes:
        flags: List[str] = []
        if te["embedded_in"]:
            flags.append("embedded")
        if te["contains"]:
            flags.append("container")
        if te["overlaps"]:
            flags.append("overlap")
        te["relation"] = ",".join(flags) if flags else "simple"


def list_fasta_files(input_path: Path) -> List[Path]:
    """
    If input_path is a file -> [file].
    If it's a directory -> all *.fna, *.fa, *.fasta inside (non-recursive).
    """
    if input_path.is_file():
        return [input_path]

    if not input_path.is_dir():
        raise SystemExit(f"ERROR: input path is neither file nor directory: {input_path}")

    files: List[Path] = []
    for pattern in ("*.fna", "*.fa", "*.fasta"):
        files.extend(sorted(input_path.glob(pattern)))
    if not files:
        raise SystemExit(f"ERROR: no FASTA files found in directory: {input_path}")
    return files


def process_genome_file(
    fasta_path: Path,
    min_ir_len: int,
    max_ir_len: int,
    min_te_len: int,
    max_te_len: int,
) -> List[dict]:
    """
    Process one FASTA genome file (possibly multi-contig).
    Returns a list of TE records with metadata.
    """
    genome_name = fasta_path.name
    print(f"\nProcessing genome file: {genome_name}")

    records = read_fasta_multi(fasta_path)
    print(f"  Contigs found: {len(records)}")

    all_rows: List[dict] = []

    for contig_idx, (contig_id, seq) in enumerate(records, start=1):
        print(f"    Contig {contig_idx}: {contig_id}, length={len(seq)} bp")

        tes = find_te_candidates(
            seq,
            min_ir_len=min_ir_len,
            max_ir_len=max_ir_len,
            min_te_len=min_te_len,
            max_te_len=max_te_len,
        )

        print(f"      TE candidates detected: {len(tes)}")

        if not tes:
            continue

        classify_embedding_overlaps(tes)

        # convert to output-friendly rows
        for te in tes:
            row = {
                "genome_file": genome_name,
                "contig_id": contig_id,
                "contig_index": contig_idx,
                "te_id": te["id"],
                "start_0_based": te["start"],
                "end_0_based_exclusive": te["end"],
                "start_1_based": te["start"] + 1,
                "end_1_based": te["end"],
                "length": te["length"],
                "ir_len": te["ir_len"],
                "left_ir_seq": te["left_ir_seq"],
                "right_ir_seq": te["right_ir_seq"],
                "relation": te["relation"],  # simple / embedded / container / overlap / combo
            }
            all_rows.append(row)

    return all_rows


def write_results(rows: List[dict], out_path: Path) -> None:
    """Write all TE rows to a TSV file."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        # header
        f.write(
            "genome_file\tcontig_id\tcontig_index\tte_id\t"
            "start_0_based\tend_0_based_exclusive\tstart_1_based\tend_1_based\t"
            "length\tir_len\tleft_ir_seq\tright_ir_seq\trelation\n"
        )
        for r in rows:
            f.write(
                f"{r['genome_file']}\t{r['contig_id']}\t{r['contig_index']}\t{r['te_id']}\t"
                f"{r['start_0_based']}\t{r['end_0_based_exclusive']}\t"
                f"{r['start_1_based']}\t{r['end_1_based']}\t"
                f"{r['length']}\t{r['ir_len']}\t"
                f"{r['left_ir_seq']}\t{r['right_ir_seq']}\t"
                f"{r['relation']}\n"
            )


def main():
    here = Path(__file__).resolve().parent
    default_input = here / "data"
    default_out = here / "data" / "ex3_bacteria_tes.tsv"

    parser = argparse.ArgumentParser(
        description=(
            "Exercise 3: detect transposable elements in bacterial genomes "
            "based on inverted repeats of length 5–6 bp."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=default_input,
        help=(
            "Input FASTA/.fna file OR directory containing multiple genomes "
            f"(default: {default_input})"
        ),
    )
    parser.add_argument(
        "-o",
        "--output-tsv",
        type=Path,
        default=default_out,
        help=f"Output TSV with detected TEs (default: {default_out})",
    )
    parser.add_argument(
        "--min-ir-len",
        type=int,
        default=5,
        help="Minimum inverted repeat length (default: 5)",
    )
    parser.add_argument(
        "--max-ir-len",
        type=int,
        default=6,
        help="Maximum inverted repeat length (default: 6)",
    )
    parser.add_argument(
        "--min-te-len",
        type=int,
        default=100,
        help="Minimum TE length (bp) between IRs (default: 100)",
    )
    parser.add_argument(
        "--max-te-len",
        type=int,
        default=5000,
        help="Maximum TE length (bp) between IRs (default: 5000)",
    )
    args = parser.parse_args()

    fasta_files = list_fasta_files(args.input)

    all_rows: List[dict] = []
    for fasta_path in fasta_files:
        rows = process_genome_file(
            fasta_path,
            min_ir_len=args.min_ir_len,
            max_ir_len=args.max_ir_len,
            min_te_len=args.min_te_len,
            max_te_len=args.max_te_len,
        )
        all_rows.extend(rows)

    print(f"\nTotal TE records across all genomes: {len(all_rows)}")
    write_results(all_rows, args.output_tsv)
    print(f"Results written to: {args.output_tsv}")


if __name__ == "__main__":
    main()
