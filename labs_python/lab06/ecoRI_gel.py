#!/usr/bin/env python3
"""
EcoRI digest + gel simulation (lab06, task 2) — robust reader

- Input: either
    A) a folder with FASTA files  -> one lane per file, OR
    B) a single multi-FASTA file  -> one lane per record.
- Accepts: .fa/.fna/.fasta (any case) and .gz variants transparently.
- Enzyme: EcoRI (GAATTC), cleavage G^AATTC (linear molecules).
- Output:
  * out/ecoRI_influenza.png : multi-lane gel (ladder + one lane per sample/record)
  * out/ecoRI_fragments.tsv : fragment sizes per lane

Dependency: Pillow (installed for task 1).
"""

from pathlib import Path
from typing import List, Tuple, Dict
import argparse
import math
import gzip

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "out"


# ---------------------- FASTA utils ---------------------- #
def _open_text_maybe_gzip(fp: Path):
    """Open text; if gzipped, read as gzip text."""
    # Trust extension first; also sniff header as fallback
    lower = fp.name.lower()
    if lower.endswith(".gz"):
        return gzip.open(fp, "rt", encoding="utf-8", errors="ignore")
    # sniff magic number
    with fp.open("rb") as bf:
        start = bf.read(2)
    if start == b"\x1f\x8b":  # gzip magic
        return gzip.open(fp, "rt", encoding="utf-8", errors="ignore")
    return fp.open("r", encoding="utf-8", errors="ignore")


def read_fasta_records(fp: Path) -> List[Tuple[str, str]]:
    """Return list of (header, sequence). Keeps ACGTN; converts U->T. Skips empty records."""
    records: List[Tuple[str, str]] = []
    if not fp.exists() or not fp.is_file():
        return records
    header = None
    seq = []
    with _open_text_maybe_gzip(fp) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None and seq:
                    s = _clean("".join(seq))
                    if s:
                        records.append((header, s))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if header is not None and seq:
            s = _clean("".join(seq))
            if s:
                records.append((header, s))
    return records


def _clean(s: str) -> str:
    s = s.upper().replace("U", "T")  # handle RNA sources
    # Map wide IUPAC alphabet to DNA/N; keep A/C/G/T/N, replace others with N
    out = []
    for c in s:
        if c in "ACGTN":
            out.append(c)
        elif c in "RYSWKMBDHV":  # ambiguous -> N
            out.append("N")
        # ignore any other symbols/spaces/digits
    return "".join(out)


def short_label(hdr: str, idx: int) -> str:
    """Create a compact lane label from FASTA header."""
    if not hdr:
        return f"rec{idx:02d}"
    token = hdr.split()[0]
    if len(token) > 14:
        token = token[:14]
    return token or f"rec{idx:02d}"


# -------------------- Restriction digest ------------------ #
def ecoRI_fragments_linear(seq: str) -> List[int]:
    """
    EcoRI recognizes GAATTC and cuts G^AATTC (after G, index+1).
    For a linear molecule, fragments are lengths between consecutive cut positions.
    """
    motif = "GAATTC"
    cut_shift = 1
    cuts = []
    i = 0
    while True:
        j = seq.find(motif, i)
        if j == -1:
            break
        cuts.append(j + cut_shift)
        i = j + 1  # allow overlapping/adjacent sites
    pos = [0] + cuts + [len(seq)]
    pos = sorted(set(p for p in pos if 0 <= p <= len(seq)))
    frags = [b - a for a, b in zip(pos, pos[1:]) if b - a > 0]
    return frags if frags else [len(seq)]


# ---------------------- Gel drawing ----------------------- #
def _log_migration_y(size: int, smin: int, smax: int, top: int, bottom: int) -> int:
    if smax == smin:
        return (top + bottom) // 2
    import math as _m
    num = _m.log10(smax) - _m.log10(max(size, 1))
    den = _m.log10(smax) - _m.log10(max(smin, 1))
    frac = max(0.0, min(1.0, num / den))
    return int(top + frac * (bottom - top))


def draw_multi_lane_gel(
    lanes: Dict[str, List[int]],
    out_png: Path,
    ladder_sizes: List[int] = (10000, 7000, 5000, 3000, 2000, 1500, 1000, 700, 500, 300, 200, 100),
):
    labels = list(lanes.keys())
    sizes_all = [s for lst in lanes.values() for s in lst] or [100, 3000]
    smin = min(min(ladder_sizes), max(1, min(sizes_all)))
    smax = max(max(ladder_sizes), max(sizes_all))

    # layout
    n = len(labels)
    lane_w = 110
    gap = 40
    left_pad = 50
    right_pad = 40
    top = 90
    bottom_pad = 100
    ladder_w = 90
    wells_h = 18

    W = left_pad + ladder_w + gap + n * lane_w + (n - 1) * gap + right_pad
    H = 900

    img = Image.new("L", (W, H), color=0)
    d = ImageDraw.Draw(img)
    try:
        font = ImageFont.load_default()
    except Exception:
        font = None

    # frame
    d.rectangle([10, 10, W - 10, H - 10], outline=180, width=1)

    # ladder background
    ladder_x = left_pad
    d.rectangle([ladder_x, top - wells_h, ladder_x + ladder_w, H - bottom_pad], fill=32)
    d.rectangle([ladder_x + 6, top - wells_h, ladder_x + ladder_w - 6, top], fill=128)

    # ladder ticks + labels
    for lad in ladder_sizes:
        y = _log_migration_y(lad, smin, smax, top, H - bottom_pad)
        d.rectangle([ladder_x + 12, y - 2, ladder_x + ladder_w - 12, y + 2], fill=200)
        d.text((15, y - 6), f"{lad} bp", fill=200, font=font)

    # lanes
    cur_x = ladder_x + ladder_w + gap
    for label in labels:
        d.rectangle([cur_x, top - wells_h, cur_x + lane_w, H - bottom_pad], fill=32)
        d.rectangle([cur_x + 8, top - wells_h, cur_x + lane_w - 8, top], fill=128)
        for sz in lanes[label]:
            y = _log_migration_y(sz, smin, smax, top, H - bottom_pad)
            d.rectangle([cur_x + 10, y - 3, cur_x + lane_w - 10, y + 3], fill=240)
        # label under lane
        txt = label[:14]
        tw = d.textlength(txt, font=font) if hasattr(d, "textlength") else len(txt) * 6
        d.text((cur_x + (lane_w - tw) / 2, H - bottom_pad + 25), txt, fill=200, font=font)
        cur_x += lane_w + gap

    out_png.parent.mkdir(parents=True, exist_ok=True)
    img.save(out_png)


# -------------------------- CLI --------------------------- #
def main():
    ap = argparse.ArgumentParser(description="EcoRI digest + gel simulation for FASTA input (dir or multi-FASTA file)")
    ap.add_argument("--input", "-i", type=str,
                    default=str(ROOT / "data" / "influenza"),
                    help="Folder with FASTAs (one lane per file) OR a single multi-FASTA file.")
    ap.add_argument("--out_png", type=str, default=str(OUT_DIR / "ecoRI_influenza.png"),
                    help="Output gel image.")
    ap.add_argument("--out_tsv", type=str, default=str(OUT_DIR / "ecoRI_fragments.tsv"),
                    help="Output TSV with fragment sizes.")
    args = ap.parse_args()

    in_path = Path(args.input)
    lanes: Dict[str, List[int]] = {}

    def add_lane(label: str, seqs: List[str]):
        frags_all: List[int] = []
        for seq in seqs:
            if not seq:
                continue
            frags_all.extend(ecoRI_fragments_linear(seq))
        if frags_all:
            lanes[label] = sorted(frags_all)

    if in_path.is_dir():
        # collect likely FASTA files (case-insensitive, with or without .gz)
        patterns = ["*.fa", "*.fna", "*.fasta", "*.FA", "*.FNA", "*.FASTA",
                    "*.fa.gz", "*.fna.gz", "*.fasta.gz", "*.FA.GZ", "*.FNA.GZ", "*.FASTA.GZ"]
        files = []
        for pat in patterns:
            files.extend(in_path.glob(pat))
        files = sorted(set(files))
        if not files:
            print(f"[warn] No FASTA-like files in: {in_path}")
        for fp in files:
            recs = read_fasta_records(fp)
            seqs = [s for _, s in recs]
            print(f"[info] {fp.name}: {len(seqs)} record(s)")
            add_lane(fp.stem, seqs)
    elif in_path.is_file():
        recs = read_fasta_records(in_path)
        if not recs:
            print(f"[warn] No FASTA records in file: {in_path}")
        for i, (hdr, seq) in enumerate(recs, 1):
            label = short_label(hdr, i)
            print(f"[info] record {i}: {label}, length={len(seq)}")
            add_lane(label, [seq])
    else:
        raise FileNotFoundError(f"Input not found: {in_path}")

    if not lanes:
        raise SystemExit("No lanes produced; check your input (are they FASTA? gzipped?).")

    # write sizes table
    tsv = Path(args.out_tsv)
    tsv.parent.mkdir(parents=True, exist_ok=True)
    with tsv.open("w") as f:
        f.write("lane\tfragment_index\tlength_bp\n")
        for label in lanes:
            for k, L in enumerate(lanes[label], 1):
                f.write(f"{label}\t{k}\t{L}\n")

    # draw gel
    draw_multi_lane_gel(lanes, Path(args.out_png))

    print(f"[ok] Lanes: {len(lanes)}")
    print(f"[ok] TSV:   {tsv}")
    print(f"[ok] Gel:   {args.out_png}")


if __name__ == "__main__":
    main()
