
"""
Gel electrophoresis toy simulation (lab06)

Steps implemented:
1) Read an arbitrary DNA sequence from a FASTA (1000–3000 nt recommended).
2) Take 10 random samples with random lengths in [100, 3000] but clamped to the source length.
3) Store samples in an array and write them to out/samples.fasta.
4) Render a simple gel image to out/gel.png:
   - Small fragments migrate farther (appear lower).
   - Large fragments migrate less (appear higher).
No external internet access needed; you provide the FASTA file.
"""

import argparse
import math
import os
import random
from pathlib import Path
from typing import List, Tuple

from PIL import Image, ImageDraw, ImageFont  # pillow

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "out"
DATA_DIR = ROOT / "data"

# ----------------------------- FASTA utils ----------------------------- #
def read_fasta_sequence(fp: Path) -> str:
    if not fp.exists():
        raise FileNotFoundError(f"Input FASTA not found: {fp}")
    seq = []
    with fp.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append("".join(c for c in line.upper() if c in "ACGTN"))
    s = "".join(seq)
    if not s:
        raise ValueError("No DNA letters (ACGTN) found in the FASTA.")
    return s


def write_fasta(records: List[Tuple[str, str]], fp: Path) -> None:
    fp.parent.mkdir(parents=True, exist_ok=True)
    with fp.open("w") as f:
        for hdr, s in records:
            f.write(f">{hdr}\n")
            # wrap at 80 columns
            for i in range(0, len(s), 80):
                f.write(s[i : i + 80] + "\n")


# ----------------------------- Sampling ----------------------------- #
def take_random_samples(
    sequence: str,
    n_samples: int = 10,
    min_len: int = 100,
    max_len: int = 3000,
    rng: random.Random = random,
) -> List[Tuple[int, int, str]]:
    """Return list of (start, length, fragment_seq)."""
    L = len(sequence)
    upper = min(max_len, L)
    lower = min_len if L >= min_len else L  # clamp if sequence is short

    samples = []
    for _ in range(n_samples):
        frag_len = rng.randint(lower, upper)
        start_max = L - frag_len
        if start_max < 0:
            # Shouldn't happen after clamping, but guard anyway
            frag_len = L
            start_max = 0
        start = rng.randint(0, start_max)
        frag = sequence[start : start + frag_len]
        samples.append((start, frag_len, frag))
    return samples


# ----------------------------- Gel rendering ----------------------------- #
def _log_migration_y(size: int, smin: int, smax: int, top: int, bottom: int) -> int:
    """
    Map fragment size to vertical position using a log relationship:
      larger size -> smaller migration -> nearer to 'top'
      smaller size -> larger migration -> nearer to 'bottom'
    """
    # avoid division by zero
    if smax == smin:
        return (top + bottom) // 2
    num = math.log10(smax) - math.log10(size)
    den = math.log10(smax) - math.log10(smin)
    frac = max(0.0, min(1.0, num / den))  # 0..1
    return int(top + frac * (bottom - top))


def draw_gel(
    sizes: List[int],
    out_png: Path,
    ladder_sizes: List[int] = (3000, 2000, 1500, 1000, 700, 500, 300, 200, 100),
) -> None:
    """
    Create a grayscale gel-like image with one lane of sample bands and a simple 'ladder'.
    We avoid matplotlib so there are no chart color/style constraints.
    """
    W, H = 600, 900                       # canvas
    GEL_X, GEL_W = 240, 120               # sample lane position and width
    LADDER_X, LADDER_W = 120, 60          # ladder lane
    WELL_H = 16
    TOP_MARGIN, BOT_MARGIN = 80, 60

    img = Image.new("L", (W, H), color=0)  # black background
    d = ImageDraw.Draw(img)

    # lanes (light grey rectangles)
    d.rectangle([GEL_X, TOP_MARGIN - WELL_H, GEL_X + GEL_W, H - BOT_MARGIN], fill=32)
    d.rectangle([LADDER_X, TOP_MARGIN - WELL_H, LADDER_X + LADDER_W, H - BOT_MARGIN], fill=32)

    # wells (brighter)
    d.rectangle([GEL_X + 8, TOP_MARGIN - WELL_H, GEL_X + GEL_W - 8, TOP_MARGIN], fill=128)
    d.rectangle([LADDER_X + 8, TOP_MARGIN - WELL_H, LADDER_X + LADDER_W - 8, TOP_MARGIN], fill=128)

    # scales
    smin = max(1, min(sizes))
    smax = max(sizes)

    # draw ladder ticks and labels
    try:
        font = ImageFont.load_default()
    except Exception:
        font = None

    for lad in ladder_sizes:
        y = _log_migration_y(lad, min(ladder_sizes), max(ladder_sizes), TOP_MARGIN, H - BOT_MARGIN)
        # tick band
        d.rectangle([LADDER_X + 10, y - 2, LADDER_X + LADDER_W - 10, y + 2], fill=200)
        # label
        label = f"{lad} bp"
        d.text((20, y - 6), label, fill=200, font=font)

    # draw sample bands (same intensity/width)
    for sz in sizes:
        y = _log_migration_y(sz, smin, smax, TOP_MARGIN, H - BOT_MARGIN)
        d.rectangle([GEL_X + 10, y - 3, GEL_X + GEL_W - 10, y + 3], fill=240)

    # light frame
    d.rectangle([10, 10, W - 10, H - 10], outline=180, width=1)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    img.save(out_png)


# ----------------------------- Main CLI ----------------------------- #
def main():
    ap = argparse.ArgumentParser(description="Gel electrophoresis toy simulation")
    ap.add_argument("--input", "-i", type=str, default=str(DATA_DIR / "sequence.fasta"),
                    help="Path to input FASTA downloaded from NCBI")
    ap.add_argument("--samples", "-n", type=int, default=10, help="Number of random fragments")
    ap.add_argument("--minlen", type=int, default=100, help="Minimum fragment length")
    ap.add_argument("--maxlen", type=int, default=3000, help="Maximum fragment length")
    ap.add_argument("--seed", type=int, default=42, help="PRNG seed for reproducibility")
    args = ap.parse_args()

    rng = random.Random(args.seed)

    # 1) Read sequence
    seq = read_fasta_sequence(Path(args.input))
    if not (1000 <= len(seq) <= 3000):
        # Soft enforcement: if the user picked a longer sequence, take a 1000–3000 nt window.
        if len(seq) > 3000:
            win = rng.randint(1000, 3000)
            start = rng.randint(0, len(seq) - win)
            seq = seq[start : start + win]
        else:
            raise ValueError(
                f"Sequence length is {len(seq)} nt (<1000). Please pick an NCBI sequence between 1000 and 3000 nt."
            )

    # 2–3) Random samples
    samples = take_random_samples(
        sequence=seq,
        n_samples=args.samples,
        min_len=args.minlen,
        max_len=args.maxlen,
        rng=rng,
    )

    # Prepare outputs
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fasta_out = OUT_DIR / "samples.fasta"
    tsv_out = OUT_DIR / "sizes.tsv"
    png_out = OUT_DIR / "gel.png"

    # write samples.fasta
    records = []
    sizes = []
    for idx, (start, frag_len, frag) in enumerate(samples, 1):
        hdr = f"sample_{idx}|start={start}|len={frag_len}"
        records.append((hdr, frag))
        sizes.append(frag_len)
    write_fasta(records, fasta_out)

    # write sizes.tsv for reference
    with tsv_out.open("w") as f:
        f.write("sample\tstart\tlength\n")
        for idx, (start, frag_len, _) in enumerate(samples, 1):
            f.write(f"{idx}\t{start}\t{frag_len}\n")

    # 4) Draw gel
    draw_gel(sizes=sizes, out_png=png_out)

    print(f"[ok] Input length: {len(seq)} nt")
    print(f"[ok] Wrote {len(samples)} samples to: {fasta_out}")
    print(f"[ok] Sizes table: {tsv_out}")
    print(f"[ok] Gel image: {png_out}")


if __name__ == "__main__":
    main()
