#!/usr/bin/env python3
"""
Exercise 4 (DNA)
Load a DNA transition model JSON (from ex2) and generate a new DNA sequence
by sampling from transition probabilities.

Stdlib only.
"""

import argparse
import json
import random
import sys
from typing import Dict, List


def load_json(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def sample_next_from_row(states: List[str], probs_row: List[float], rng: random.Random) -> str:
    """
    Sample next state from a dense probability row aligned with `states`.
    Handles rows that sum to ~0 by falling back to uniform sampling.
    """
    total = sum(probs_row)
    if total <= 0.0:
        return rng.choice(states)

    r = rng.random()
    cum = 0.0
    last_state = states[-1]
    for s, p in zip(states, probs_row):
        cum += p / total  # normalize defensively
        if r <= cum:
            return s
        last_state = s
    return last_state


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate DNA using a transition model JSON.")
    parser.add_argument("--model", default="dna.json", help="Path to DNA model JSON (ex2 output).")
    parser.add_argument("--length", type=int, default=None, help="Output length (default: model sequence_length).")
    parser.add_argument("--start", type=str, default=None, help="Starting nucleotide (A/C/G/T).")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")
    parser.add_argument("--out", type=str, default=None, help="Optional output file (otherwise prints to stdout).")
    args = parser.parse_args()

    model = load_json(args.model)

    if model.get("type") != "dna_transition_matrix":
        print("Error: this does not look like a DNA transition model JSON.", file=sys.stderr)
        sys.exit(2)

    alphabet: List[str] = model["alphabet"]
    probs: List[List[float]] = model["probabilities"]

    if args.length is None:
        args.length = int(model.get("sequence_length", 50))
    if args.length <= 0:
        print("Error: --length must be > 0", file=sys.stderr)
        sys.exit(2)

    rng = random.Random(args.seed)

    # pick start
    if args.start:
        start = args.start.upper().strip()
        if start not in alphabet:
            print(f"Error: --start must be one of {alphabet}", file=sys.stderr)
            sys.exit(2)
        current = start
    else:
        current = rng.choice(alphabet)

    index: Dict[str, int] = {ch: i for i, ch in enumerate(alphabet)}

    out_chars = [current]
    while len(out_chars) < args.length:
        row = probs[index[current]]
        nxt = sample_next_from_row(alphabet, row, rng)
        out_chars.append(nxt)
        current = nxt

    result = "".join(out_chars)

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(result + "\n")
        print(f"Saved generated DNA to: {args.out}")
    else:
        print(result)


if __name__ == "__main__":
    main()
