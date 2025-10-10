use std::collections::HashMap;

const S: &str = "ATTGTCCCAATCTGTTG";
const ALPHABET: [char; 4] = ['A', 'C', 'G', 'T'];

fn main() {
    println!("[lab02] ex01_kmer_freq");
    println!("Sequence S (len={}): {}", S.len(), S);
    println!();

    // Di-nucleotides (k = 2)
    run_for_k(2);

    println!();

    // Tri-nucleotides (k = 3)
    run_for_k(3);
}

/// Orchestrates: generate all k-mers, count overlapping occurrences, print table.
fn run_for_k(k: usize) {
    let kmers = generate_kmers(k);
    let (counts, total_windows) = count_overlapping(S, k);

    println!("k = {} ({} total windows)", k, total_windows);
    for kmer in kmers {
        let c = counts.get(kmer.as_str()).copied().unwrap_or(0);
        let freq = if total_windows > 0 {
            (c as f64) / (total_windows as f64)
        } else {
            0.0
        };
        println!("{}: count = {:>2}, rel_freq = {:.4}", kmer, c, freq);
    }

}

/// Generate every k-mer over {A,C,G,T} in lexicographic order (AA, AC, ..., TT for k=2).
fn generate_kmers(k: usize) -> Vec<String> {
    let mut out = Vec::new();
    let mut buf = String::with_capacity(k);
    product(&mut out, &mut buf, k);
    out
}

fn product(out: &mut Vec<String>, buf: &mut String, k_left: usize) {
    if k_left == 0 {
        out.push(buf.clone());
        return;
    }
    for ch in ALPHABET {
        buf.push(ch);
        product(out, buf, k_left - 1);
        buf.pop();
    }
}

/// Count overlapping k-mer occurrences in `seq` (slide by 1).
/// Relative frequency later is count / (len(seq) - k + 1).
fn count_overlapping(seq: &str, k: usize) -> (HashMap<String, usize>, usize) {
    let n = seq.len();
    if n < k {
        return (HashMap::new(), 0);
    }

    let total = n - k + 1;
    let mut map: HashMap<String, usize> = HashMap::new();

    for i in 0..=n - k {
        // Safe since S is ASCII A/C/G/T
        let kmer = &seq[i..i + k];
        *map.entry(kmer.to_string()).or_insert(0) += 1;
    }
    (map, total)
}
