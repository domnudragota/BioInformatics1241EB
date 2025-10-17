# labs/lab03 — minimal, like lab02

.PHONY: ex1 ex2 fmt clippy clean

# ex1: interactive / args for Tm of a single sequence
ex1:
	cargo run --bin ex01_tm

# ex2: sliding window (k=9) from FASTA (.fa/.fasta or .gz) + PNG chart
ex2:
	cargo run --bin ex02_tm_window

ex3:
	cargo run --bin ex03_threshold


fmt:
	cargo fmt --all

clippy:
	cargo clippy --all-targets --all-features -D warnings

clean:
	cargo clean
