#!/usr/bin/env python3
"""
Exercise 1 (Python, Lab 04): Random read sampling + greedy overlap assembly

Steps:
1) Load an input DNA sequence (FASTA, 1,000–3,000 nt). Default: ../data/reference.fasta
2) Sample 2,000 reads of length 100–150 uniformly at random (no sequencing errors).
3) Store reads in a Python list and write reads.fasta for inspection.
4) Rebuild the original sequence using a simple greedy overlap-layout strategy:
   - Start from a seed read.
   - Extend contig to the right/left by picking the unused read with the largest
     suffix/prefix overlap above a threshold k (default k=30).
   - Continue until no further extension is possible.
5) Write assembled contig to outputs/assembled.fasta and a small report.
6) Write algorithm_limitations.txt explaining the main difficulties (repeats, low complexity, etc.).
"""
import argparse
import os
import random
from collections import Counter, defaultdict

def read_fasta_single(path: str) -> str:
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line: 
                continue
            if line.startswith(">"):
                continue
            seq.append(line.upper())
    dna = "".join(seq)
    # Filter to A/C/G/T only (soft sanitize)
    dna = "".join([c for c in dna if c in "ACGT"])
    return dna

def write_fasta(path: str, header: str, seq: str, width: int = 70) -> None:
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i+width] + "\n")

def sample_reads(seq: str, n_reads: int = 2000, min_len: int = 100, max_len: int = 150, seed: int = 42):
    rng = random.Random(seed)
    reads = []
    for _ in range(n_reads):
        L = rng.randint(min_len, max_len)
        start = rng.randint(0, len(seq) - L)
        read = seq[start:start+L]
        reads.append(read)
    return reads

def write_reads_fasta(path: str, reads):
    with open(path, "w") as f:
        for i, r in enumerate(reads, start=1):
            f.write(f">read_{i}\n{r}\n")

def build_prefix_index(reads, k):
    """Map k-length prefixes to list of read indices for quick candidate lookup."""
    index = defaultdict(list)
    for i, r in enumerate(reads):
        if len(r) >= k:
            index[r[:k]].append(i)
    return index

def build_suffix_index(reads, k):
    """Map k-length suffixes to list of read indices for quick candidate lookup."""
    index = defaultdict(list)
    for i, r in enumerate(reads):
        if len(r) >= k:
            index[r[-k:]].append(i)
    return index

def overlap_len(a: str, b: str, min_len: int) -> int:
    """
    Longest suffix of 'a' that matches a prefix of 'b' (length >= min_len).
    Returns the overlap length, or 0 if none >= min_len.
    """
    start = len(a) - min_len
    if start < 0:
        return 0
    # Search decreasing start positions to allow longer overlaps first
    for i in range(start, -1, -1):
        # Quick check first char to prune
        if a[i] == b[0] and a[i:] == b[:len(a)-i]:
            return len(a) - i
    return 0

def greedy_extend_right(contig: str, reads, used, prefix_index, k: int) -> tuple[str, int]:
    """
    Try to extend contig to the right using maximum overlap with an unused read.
    Returns (new_contig, used_read_index or -1).
    """
    if len(contig) < k:
        return contig, -1
    suffix = contig[-k:]
    best_idx = -1
    best_ov = 0
    # Candidate reads must begin with contig's last k-mer
    for idx in prefix_index.get(suffix, []):
        if used[idx]:
            continue
        r = reads[idx]
        ov = overlap_len(contig, r, k)
        if ov > best_ov:
            best_ov = ov
            best_idx = idx
    if best_idx == -1:
        return contig, -1
    # Merge
    r = reads[best_idx]
    new_contig = contig + r[best_ov:]
    return new_contig, best_idx

def greedy_extend_left(contig: str, reads, used, suffix_index, k: int) -> tuple[str, int]:
    """
    Try to extend contig to the left using maximum overlap with an unused read.
    Returns (new_contig, used_read_index or -1).
    """
    if len(contig) < k:
        return contig, -1
    prefix = contig[:k]
    best_idx = -1
    best_ov = 0
    for idx in suffix_index.get(prefix, []):
        if used[idx]:
            continue
        r = reads[idx]
        # Overlap where r's suffix matches contig's prefix
        ov = overlap_len(r, contig, k)
        if ov > best_ov:
            best_ov = ov
            best_idx = idx
    if best_idx == -1:
        return contig, -1
    r = reads[best_idx]
    new_contig = r + contig[best_ov:]
    return new_contig, best_idx

def assemble_greedy(reads, k: int = 30):
    """
    Very simple greedy contig assembler:
    - Choose the longest read as seed to improve stability.
    - Extend right, then left (alternating until no progress).
    - Mark reads as used as we consume them.
    - If we get stuck before covering the whole reference, we could relax k,
      but here we keep it simple.
    """
    if not reads:
        return ""
    # Optional: collapse exact duplicates to reduce work (but keep counts)
    counts = Counter(reads)
    unique_reads = list(counts.keys())
    # Keep mapping from unique read to its count, but we treat "used" at unique level.
    used = [False] * len(unique_reads)

    # Seed: longest read
    seed_idx = max(range(len(unique_reads)), key=lambda i: len(unique_reads[i]))
    contig = unique_reads[seed_idx]
    used[seed_idx] = True

    prefix_index = build_prefix_index(unique_reads, k)
    suffix_index = build_suffix_index(unique_reads, k)

    progress = True
    while progress:
        progress = False
        # Try to extend RIGHT as much as possible
        while True:
            contig2, idx = greedy_extend_right(contig, unique_reads, used, prefix_index, k)
            if idx == -1:
                break
            contig = contig2
            used[idx] = True
            progress = True
        # Try to extend LEFT as much as possible
        while True:
            contig2, idx = greedy_extend_left(contig, unique_reads, used, suffix_index, k)
            if idx == -1:
                break
            contig = contig2
            used[idx] = True
            progress = True
    return contig

def best_shift_identity(a: str, b: str) -> tuple[int, float]:
    """
    Crude identity estimate by sliding the shorter across the longer and counting exact matches.
    Returns (best_offset, best_identity_fraction).
    """
    if len(a) < len(b):
        short, long_ = a, b
        sign = -1  # indicates a inside b
    else:
        short, long_ = b, a
        sign = 1   # indicates b inside a

    best_id = 0.0
    best_off = 0
    for off in range(-len(short), len(long_)):
        matches = 0
        total = 0
        for i in range(len(short)):
            j = i + off
            if 0 <= j < len(long_):
                total += 1
                if short[i] == long_[j]:
                    matches += 1
        if total > 0:
            ident = matches / total
            if ident > best_id:
                best_id = ident
                best_off = off * sign
    return best_off, best_id

def main():
    parser = argparse.ArgumentParser(description="Exercise 1: random sampling + greedy assembly")
    parser.add_argument("--ref", type=str, default="../data/reference.fasta",
                        help="Path to FASTA (single sequence, 1k–3k nt). Default: ../data/reference.fasta")
    parser.add_argument("--reads", type=int, default=2000, help="Number of reads to sample (default: 2000)")
    parser.add_argument("--minlen", type=int, default=100, help="Minimum read length (default: 100)")
    parser.add_argument("--maxlen", type=int, default=150, help="Maximum read length (default: 150)")
    parser.add_argument("--k", type=int, default=30, help="Minimum overlap length (default: 30)")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed for reproducibility")
    parser.add_argument("--outdir", type=str, default="../outputs", help="Outputs folder")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    ref_seq = read_fasta_single(args.ref)
    if not (1000 <= len(ref_seq) <= 3000):
        raise SystemExit(f"Reference length must be 1000–3000 nt. Got {len(ref_seq)} nt from {args.ref}.")

    # 1) Sample reads
    reads = sample_reads(ref_seq, n_reads=args.reads, min_len=args.minlen, max_len=args.maxlen, seed=args.seed)

    # 2) Save reads (also kept in a Python list as requested)
    reads_fa = os.path.join(args.outdir, "reads.fasta")
    write_reads_fasta(reads_fa, reads)

    # 3) Assemble
    contig = assemble_greedy(reads, k=args.k)

    # 4) Write contig
    assembled_fa = os.path.join(args.outdir, "assembled.fasta")
    write_fasta(assembled_fa, "greedy_assembled_contig", contig)

    # 5) Simple evaluation
    # Check if ref is substring of contig or vice versa; also crude identity
    ref_in_contig = ref_seq in contig
    contig_in_ref = contig in ref_seq
    offset, ident = best_shift_identity(ref_seq, contig)

    report_path = os.path.join(args.outdir, "assembly_report.txt")
    with open(report_path, "w") as f:
        f.write("=== Assembly Report (Exercise 1) ===\n")
        f.write(f"Reference length: {len(ref_seq)}\n")
        f.write(f"Reads: {args.reads} (len {args.minlen}–{args.maxlen})\n")
        f.write(f"Greedy min-overlap k: {args.k}\n")
        f.write(f"Assembled contig length: {len(contig)}\n")
        f.write(f"Reference ⊂ Contig? {'YES' if ref_in_contig else 'NO'}\n")
        f.write(f"Contig ⊂ Reference? {'YES' if contig_in_ref else 'NO'}\n")
        f.write(f"Best shift (ref vs contig): {offset}\n")
        f.write(f"Crude identity at best shift: {ident:.4f}\n")

    # 6) Algorithm limitations text (kept short and to the point for the lab)
    limitations_txt = os.path.join(args.outdir, "algorithm_limitations.txt")
    with open(limitations_txt, "w") as f:
        f.write(
"""Main problem when reconstructing from short random samples (reads)

1) Repeats longer than the overlap (k) or longer than the read length:
   - If the reference has repeated segments (e.g., ...ACGTT...ACGTT...) longer than k,
     those parts are indistinguishable from one another using only local overlaps.
   - The greedy assembler may either collapse multiple copies into one (under-assembly)
     or join the wrong copies (mis-assembly). This is the dominant difficulty.

2) Low-complexity / homopolymer regions (e.g., AAAAA..., ATATAT..., or GC-rich microsatellites):
   - Many reads look almost identical, so overlaps are not unique.
   - The graph has ambiguous branches; greedy choices can pick the wrong path.

3) Branching due to near-identical repeats (paralogs, duplicated elements):
   - Even a single mismatch is enough to break a unique overlap if k is large,
     but if k is small we risk false overlaps. Tuning k is a trade-off.

4) Coverage gaps or uneven sampling (less likely here but relevant in practice):
   - If some region is under-sampled, we cannot bridge it reliably, so contigs
     stop early or remain disconnected.

5) Orientation and errors (not simulated here):
   - Real reads can be forward or reverse-complement and contain sequencing errors.
     Without handling orientation and error-tolerant overlaps, assembly quality degrades.

In short: repeats and low-complexity sequences create ambiguous overlaps that a simple greedy
overlap-based algorithm cannot resolve uniquely. Longer reads (or larger k), paired-end links,
or de Bruijn/Eulerian methods with coverage-aware heuristics are typically used to mitigate this.
"""
        )

    print("Done.")
    print(f"Wrote reads:      {reads_fa}")
    print(f"Wrote contig:     {assembled_fa}")
    print(f"Wrote report:     {report_path}")
    print(f"Wrote limitations:{limitations_txt}")

if __name__ == "__main__":
    main()
