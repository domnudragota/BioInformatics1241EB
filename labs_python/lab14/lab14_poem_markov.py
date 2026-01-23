#!/usr/bin/env python3
"""
Lab: Poem author-like detection with word-transition Markov models (LLR trend)

Inputs:
  - Eminescu poem text file (plain text)
  - Nichita Stanescu poem text file (plain text)
  - Test text file (plain text)

Method:
  1) Tokenize into words
  2) Replace each unique word with a single character (Unicode private-use chars)
  3) Build transition counts and probabilities P(next|current) for each poem
  4) Log-likelihood / log-odds score:
       beta(a->b) = log2( P_E(b|a) / P_S(b|a) )
     with Laplace smoothing (pseudocount alpha) to avoid zeros
  5) Sliding window scan of the test text; plot score trend:
       positive => more like Eminescu
       negative => more like Stanescu

Run example:
  python3 lab14_poem_markov.py \
    --eminescu eminescu.txt \
    --stanescu stanescu.txt \
    --test test.txt \
    --window 40 \
    --step 5 \
    --alpha 1.0 \
    --out scores.png

"""

from __future__ import annotations

import argparse
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt


# -----------------------------
# Tokenization / normalization
# -----------------------------
WORD_RE = re.compile(r"[A-Za-zĂÂÎȘŞȚŢăâîșşțţ]+", re.UNICODE)

def normalize_ro(s: str) -> str:
    # Normalize common Romanian cedilla variants to comma-below variants (optional)
    return (s.replace("ş", "ș")
             .replace("Ş", "Ș")
             .replace("ţ", "ț")
             .replace("Ţ", "Ț"))

def tokenize_words(text: str) -> List[str]:
    text = normalize_ro(text).lower()
    return WORD_RE.findall(text)


def read_text_file(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


# -----------------------------
# Word -> single-character mapping
# -----------------------------
PUA_START = 0xE000  # Unicode Private Use Area start

@dataclass
class VocabMap:
    word_to_char: Dict[str, str]
    char_to_word: Dict[str, str]
    unk_char: str

def build_vocab_map(tokens_a: List[str], tokens_b: List[str]) -> VocabMap:
    """
    Build a shared vocab from BOTH poems.
    Any test word not in vocab becomes <UNK>.
    """
    vocab = sorted(set(tokens_a) | set(tokens_b))
    # Reserve <UNK>
    word_to_char: Dict[str, str] = {}
    char_to_word: Dict[str, str] = {}

    unk_char = chr(PUA_START)  # first PUA char reserved for UNK
    word_to_char["<UNK>"] = unk_char
    char_to_word[unk_char] = "<UNK>"

    for i, w in enumerate(vocab, start=1):
        ch = chr(PUA_START + i)
        word_to_char[w] = ch
        char_to_word[ch] = w

    return VocabMap(word_to_char=word_to_char, char_to_word=char_to_word, unk_char=unk_char)

def encode_tokens(tokens: List[str], vocab: VocabMap) -> str:
    out = []
    for w in tokens:
        out.append(vocab.word_to_char.get(w, vocab.unk_char))
    return "".join(out)

def decode_chars(chars: str, vocab: VocabMap) -> List[str]:
    return [vocab.char_to_word.get(ch, "<UNK>") for ch in chars]


# -----------------------------
# Transition counts + scoring
# -----------------------------
Counts = Dict[str, Dict[str, int]]
Totals = Dict[str, int]

def transition_counts(seq: str) -> Tuple[Counts, Totals]:
    """
    Counts transitions x->y over adjacent symbols in seq.
    Sparse dict-of-dicts.
    """
    counts: Counts = defaultdict(lambda: defaultdict(int))
    totals: Totals = defaultdict(int)

    for x, y in zip(seq, seq[1:]):
        counts[x][y] += 1
        totals[x] += 1

    return counts, totals

def prob_next(counts: Counts, totals: Totals, prev: str, nxt: str, V: int, alpha: float) -> float:
    """
    Laplace smoothing over next-symbol vocabulary size V:
      P(nxt|prev) = (c(prev->nxt) + alpha) / (total(prev) + alpha*V)
    If prev never appears (total=0), this becomes uniform 1/V.
    """
    c = counts.get(prev, {}).get(nxt, 0)
    t = totals.get(prev, 0)
    return (c + alpha) / (t + alpha * V)

def beta_llr(prev: str, nxt: str,
             cE: Counts, tE: Totals,
             cS: Counts, tS: Totals,
             V: int, alpha: float) -> float:
    """
    beta = log2( P_E(nxt|prev) / P_S(nxt|prev) )
    """
    pE = prob_next(cE, tE, prev, nxt, V, alpha)
    pS = prob_next(cS, tS, prev, nxt, V, alpha)
    return math.log2(pE / pS)

def score_window(seq: str, start: int, window_len: int,
                 cE: Counts, tE: Totals,
                 cS: Counts, tS: Totals,
                 V: int, alpha: float) -> float:
    """
    Sum beta over transitions inside the window [start, start+window_len).
    window_len is in *tokens* (characters = words).
    """
    end = start + window_len
    if end > len(seq):
        return float("nan")
    if window_len < 2:
        return 0.0

    total = 0.0
    for x, y in zip(seq[start:end], seq[start+1:end]):
        total += beta_llr(x, y, cE, tE, cS, tS, V, alpha)
    return total


# -----------------------------
# CLI / main
# -----------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--eminescu", required=True, help="Path to Eminescu poem (.txt, plain text)")
    ap.add_argument("--stanescu", required=True, help="Path to Nichita Stanescu poem (.txt, plain text)")
    ap.add_argument("--test", required=True, help="Path to test/mixed text (.txt, plain text)")
    ap.add_argument("--window", type=int, default=40, help="Sliding window size (words). Default: 40")
    ap.add_argument("--step", type=int, default=5, help="Step size (words). Default: 5")
    ap.add_argument("--alpha", type=float, default=1.0, help="Laplace pseudocount. Default: 1.0")
    ap.add_argument("--out", default="scores.png", help="Output chart filename. Default: scores.png")
    ap.add_argument("--topk", type=int, default=5, help="Print top-k strongest windows per author. Default: 5")
    args = ap.parse_args()

    # Read and tokenize
    em_text = read_text_file(args.eminescu)
    st_text = read_text_file(args.stanescu)
    te_text = read_text_file(args.test)

    em_tokens = tokenize_words(em_text)
    st_tokens = tokenize_words(st_text)
    te_tokens = tokenize_words(te_text)

    if len(em_tokens) < 2 or len(st_tokens) < 2 or len(te_tokens) < 2:
        raise SystemExit("Error: each input must contain at least 2 word tokens after tokenization.")

    # Build shared vocab (from training poems only), map to single charactershook
    vocab = build_vocab_map(em_tokens, st_tokens)
    V = len(vocab.word_to_char)  # includes <UNK>

    em_seq = encode_tokens(em_tokens, vocab)
    st_seq = encode_tokens(st_tokens, vocab)
    te_seq = encode_tokens(te_tokens, vocab)

    # Transition counts
    cE, tE = transition_counts(em_seq)
    cS, tS = transition_counts(st_seq)

    # Sliding scan
    w = args.window
    step = args.step
    alpha = args.alpha

    positions: List[int] = []
    scores: List[float] = []

    for start in range(0, len(te_seq) - w + 1, step):
        sc = score_window(te_seq, start, w, cE, tE, cS, tS, V, alpha)
        positions.append(start + w // 2)  # midpoint (word index)
        scores.append(sc)

    # Plot
    plt.figure()
    plt.plot(positions, scores)
    plt.axhline(0.0)
    plt.title("Sliding-window log-odds: + Eminescu, - Stanescu")
    plt.xlabel("Word index (window midpoint)")
    plt.ylabel("Log-odds score (sum of log2 ratios)")
    plt.tight_layout()
    plt.savefig(args.out, dpi=200)

    # Print summary + strongest windows
    print(f"Vocab size (incl <UNK>): {V}")
    print(f"Test length (words): {len(te_tokens)}")
    print(f"Window={w} Step={step} Alpha={alpha}")
    print(f"Chart saved to: {args.out}")
    print()

    # Collect windows for ranking
    window_records: List[Tuple[float, int]] = []  # (score, start)
    for start in range(0, len(te_seq) - w + 1, step):
        sc = score_window(te_seq, start, w, cE, tE, cS, tS, V, alpha)
        window_records.append((sc, start))

    # Top Eminescu-like (highest scores) and Stanescu-like (lowest scores)
    topk = args.topk
    window_records_sorted = sorted(window_records, key=lambda x: x[0])
    most_stanescu = window_records_sorted[:topk]
    most_eminescu = window_records_sorted[-topk:][::-1]

    def snippet(start: int) -> str:
        words = te_tokens[start:start+w]
        return " ".join(words[:min(20, len(words))]) + (" ..." if len(words) > 20 else "")

    print(f"Top {topk} most Stanescu-like windows (lowest scores):")
    for sc, start in most_stanescu:
        print(f"  score={sc:.3f}  start={start}  snippet: {snippet(start)}")
    print()

    print(f"Top {topk} most Eminescu-like windows (highest scores):")
    for sc, start in most_eminescu:
        print(f"  score={sc:.3f}  start={start}  snippet: {snippet(start)}")
    print()


if __name__ == "__main__":
    main()
