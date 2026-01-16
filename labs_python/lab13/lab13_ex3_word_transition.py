"""
Lab 13 - Exercise 3
English text (~300 letters) -> word transition probabilities -> JSON file.
Each unique word is mapped to a unique symbol (S0, S1, ...).

Stdlib only.
"""

import argparse
import json
import re
import sys
from typing import Dict, List


def read_text_file(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


def normalize_text_to_words(text: str) -> List[str]:
    # Lowercase, replace non-letters with spaces, collapse spaces
    lowered = text.lower()
    cleaned = re.sub(r"[^a-z]+", " ", lowered)
    cleaned = cleaned.strip()
    if not cleaned:
        return []
    return cleaned.split()


def count_letters(text: str) -> int:
    return sum(1 for ch in text if ch.isalpha())


def build_symbol_mapping(words: List[str]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    next_id = 0
    for w in words:
        if w not in mapping:
            mapping[w] = f"S{next_id}"
            next_id += 1
    return mapping


def build_word_transition_probs(words: List[str], word_to_sym: Dict[str, str]) -> Dict[str, Dict[str, float]]:
    # counts[from_sym][to_sym] = int
    counts: Dict[str, Dict[str, int]] = {}

    for i in range(len(words) - 1):
        a = word_to_sym[words[i]]
        b = word_to_sym[words[i + 1]]
        if a not in counts:
            counts[a] = {}
        counts[a][b] = counts[a].get(b, 0) + 1

    # normalize to probabilities
    probs: Dict[str, Dict[str, float]] = {}
    for a, row in counts.items():
        total = sum(row.values())
        probs[a] = {}
        for b, c in row.items():
            probs[a][b] = c / total if total else 0.0

    return probs


def main() -> None:
    parser = argparse.ArgumentParser(description="Exercise 3: word transition probabilities -> JSON (with symbols)")
    parser.add_argument("--text", type=str, help="Text directly (quote it in shell).")
    parser.add_argument("--file", type=str, help="Path to a text file.")
    parser.add_argument("--out", type=str, default="word_transitions.json", help="Output JSON path.")
    parser.add_argument("--expected_letters", type=int, default=300, help="Expected letter count (default: 300).")
    args = parser.parse_args()

    if not args.text and not args.file:
        print("Error: provide --text or --file", file=sys.stderr)
        sys.exit(2)
    if args.text and args.file:
        print("Error: provide only one of --text or --file", file=sys.stderr)
        sys.exit(2)

    raw = args.text if args.text else read_text_file(args.file)

    letters = count_letters(raw)
    if letters != args.expected_letters:
        print(f"Warning: text has {letters} letters, expected {args.expected_letters}. Continuing...")

    words = normalize_text_to_words(raw)
    if len(words) < 2:
        print("Error: need at least 2 words after normalization.", file=sys.stderr)
        sys.exit(2)

    word_to_sym = build_symbol_mapping(words)
    sym_to_word = {sym: word for word, sym in word_to_sym.items()}

    probs = build_word_transition_probs(words, word_to_sym)

    payload = {
        "type": "word_transition_probabilities",
        "letters_count": letters,
        "word_count": len(words),
        "unique_words": len(word_to_sym),
        "word_to_symbol": word_to_sym,
        "symbol_to_word": sym_to_word,
        "transition_probabilities": probs,
    }

    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print(f"Saved word transition probabilities to: {args.out}")


if __name__ == "__main__":
    main()
