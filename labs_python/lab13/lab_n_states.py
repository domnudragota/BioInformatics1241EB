#!/usr/bin/env python3
"""
Lab - n-state prediction (5 discrete steps)

Given:
  - A: n x n square matrix
  - x0: vector of length n

Compute:
  x1 = A * x0
  x2 = A * x1
  ...
  x5 = A * x4

Input is read from console (space or comma-separated values supported).
"""

from typing import List


def parse_number_list(line: str) -> List[float]:
    """
    Parses a line containing numbers separated by spaces and/or commas.
    Example accepted inputs:
      "1 2 3"
      "1,2,3"
      "1, 2, 3"
      "1  2,3"
    """
    cleaned = line.replace(",", " ")
    parts = cleaned.split()

    values: List[float] = []
    for p in parts:
        values.append(float(p))
    return values


def read_int(prompt: str) -> int:
    while True:
        raw = input(prompt).strip()
        try:
            value = int(raw)
            if value <= 0:
                print("n must be a positive integer.")
                continue
            return value
        except ValueError:
            print("Please enter a valid integer.")


def read_matrix(n: int) -> List[List[float]]:
    print(f"\nEnter the {n}x{n} matrix A (one row per line).")
    print("Numbers can be space-separated or comma-separated.\n")

    matrix: List[List[float]] = []
    for i in range(n):
        while True:
            row_raw = input(f"A[{i}] row ({n} values): ").strip()
            try:
                row = parse_number_list(row_raw)
            except ValueError:
                print("Invalid row: please use only numbers.")
                continue

            if len(row) != n:
                print(f"Row must have exactly {n} numbers (you entered {len(row)}).")
                continue

            matrix.append(row)
            break

    return matrix


def read_vector(n: int) -> List[float]:
    print(f"\nEnter the initial vector x0 ({n} values).")
    print("Numbers can be space-separated or comma-separated.\n")

    while True:
        raw = input("x0: ").strip()
        try:
            vec = parse_number_list(raw)
        except ValueError:
            print("Invalid vector: please use only numbers.")
            continue

        if len(vec) != n:
            print(f"Vector must have exactly {n} numbers (you entered {len(vec)}).")
            continue

        return vec


def mat_vec_mul(matrix: List[List[float]], vector: List[float]) -> List[float]:
    """
    Computes y = matrix * vector
    where matrix is n x n and vector is length n.
    """
    n = len(matrix)
    result: List[float] = [0.0] * n

    for i in range(n):
        s = 0.0
        for j in range(n):
            s += matrix[i][j] * vector[j]
        result[i] = s

    return result


def predict_steps(matrix: List[List[float]], x0: List[float], steps: int = 5) -> List[List[float]]:
    """
    Returns a list of vectors:
      [x0, x1, x2, ..., x_steps]
    """
    all_states: List[List[float]] = [x0[:]]  # copy x0
    current = x0[:]

    for _ in range(steps):
        current = mat_vec_mul(matrix, current)
        all_states.append(current[:])

    return all_states


def format_vector(vec: List[float], precision: int = 6) -> str:
    fmt = f"{{:.{precision}f}}"
    return "[" + ", ".join(fmt.format(v) for v in vec) + "]"


def main() -> None:
    print("=== n-state prediction (5 discrete steps) ===")

    n = read_int("Enter n (number of states): ")
    A = read_matrix(n)
    x0 = read_vector(n)

    steps = 5
    states = predict_steps(A, x0, steps=steps)

    print("\n--- Results ---")
    for t, xt in enumerate(states):
        print(f"x{t} = {format_vector(xt)}")

    print("\nDone.")


if __name__ == "__main__":
    main()
