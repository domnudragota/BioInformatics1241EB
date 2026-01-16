"""
Exercise 4 (English text)
Load a word transition model JSON (from ex3) and generate a new word sequence
by sampling transitions between symbols, then decoding symbols back to words.

Stdlib only.
"""

import argparse
import json
import random
import sys
from typing import Dict, List, Optional


def load_json(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def sample_from_sparse_dist(dist: Dict[str, float], rng: random.Random) -> str:
    """
    Sample a key from a sparse distribution dict: {state: prob}.
    Defensive normalization and safe fallback.
    """
    if not dist:
        raise ValueError("Empty distribution")

    total = sum(dist.values())
    if total <= 0.0:
        return rng.choice(list(dist.keys()))

    r = rng.random()
    cum = 0.0
    last_key = next(iter(dist.keys()))
    for k, p in dist.items():
        cum += p / total
        if r <= cum:
            return k
        last_key = k
    return last_key


def choose_start_symbol(
    transitions: Dict[str, Dict[str, float]],
    symbol_to_word: Dict[str, str],
    start_word: Optional[str],
    start_symbol: Optional[str],
    rng: random.Random,
) -> str:
    if start_symbol:
        if start_symbol not in symbol_to_word:
            raise ValueError("start_symbol not found in model")
        return start_symbol

    if start_word:
        # reverse lookup
        for sym, w in symbol_to_word.items():
            if w == start_word.lower():
                return sym
        raise ValueError("start_word not found in model")

    # default: pick a symbol that has outgoing transitions
    candidates = list(transitions.keys())
    if not candidates:
        # fallback: any symbol
        return rng.choice(list(symbol_to_word.keys()))
    return rng.choice(candidates)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate English-like text using a word transition JSON model.")
    parser.add_argument("--model", default="words.json", help="Path to word model JSON (ex3 output).")
    parser.add_argument("--words", type=int, default=None, help="Number of words to generate (default: model word_count).")
    parser.add_argument("--letters", type=int, default=None, help="Target number of letters (stops when reached).")
    parser.add_argument("--start_word", type=str, default=None, help="Optional start word (must exist in training vocab).")
    parser.add_argument("--start_symbol", type=str, default=None, help="Optional start symbol (e.g., S0).")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")
    parser.add_argument("--out", type=str, default=None, help="Optional output file (otherwise prints).")
    args = parser.parse_args()

    model = load_json(args.model)

    if model.get("type") != "word_transition_probabilities":
        print("Error: this does not look like a word transition model JSON.", file=sys.stderr)
        sys.exit(2)

    transitions: Dict[str, Dict[str, float]] = model["transition_probabilities"]
    symbol_to_word: Dict[str, str] = model["symbol_to_word"]

    if args.words is None and args.letters is None:
        args.words = int(model.get("word_count", 100))

    rng = random.Random(args.seed)

    try:
        current = choose_start_symbol(
            transitions=transitions,
            symbol_to_word=symbol_to_word,
            start_word=args.start_word,
            start_symbol=args.start_symbol,
            rng=rng,
        )
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(2)

    out_syms: List[str] = [current]

    def current_letters_count(words: List[str]) -> int:
        return sum(sum(1 for c in w if c.isalpha()) for w in words)

    # Generate by letters OR by word count
    while True:
        if args.letters is not None:
            decoded = [symbol_to_word.get(s, "") for s in out_syms]
            if current_letters_count(decoded) >= args.letters:
                break
        else:
            if len(out_syms) >= args.words:
                break

        row = transitions.get(current)
        if not row:
            # dead end: restart from a random valid state
            current = rng.choice(list(transitions.keys())) if transitions else rng.choice(list(symbol_to_word.keys()))
            out_syms.append(current)
            continue

        nxt = sample_from_sparse_dist(row, rng)
        out_syms.append(nxt)
        current = nxt

    out_words = [symbol_to_word.get(s, "") for s in out_syms]
    result = " ".join(w for w in out_words if w)

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(result + "\n")
        print(f"Saved generated text to: {args.out}")
    else:
        print(result)


if __name__ == "__main__":
    main()
