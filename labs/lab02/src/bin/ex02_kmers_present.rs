use std::collections::HashSet;
use std::env;

/// Default sequence; override by passing a CLI arg:
///   make ex2 args="ABAA"  (if you wire a generic run target)
/// or simply:
///   cargo run --bin ex02_kmers_present ABAA
const DEFAULT_S: &str = "ATTGTCCCAATCTGTTG";

fn main() {
    let args: Vec<String> = env::args().collect();
    let s = if args.len() > 1 { &args[1] } else { DEFAULT_S };

    println!("[lab02] ex02_kmers_present");
    println!("Sequence S (len={}): {}", s.len(), s);

    report_unique_kmers(s, 2);
    report_unique_kmers(s, 3);
}

/// Scan S once and list unique k-mers (k=2 or 3) in discovery order.
/// NOTE: assumes ASCII nucleotides; slicing uses byte indices.
fn report_unique_kmers(s: &str, k: usize) {
    let n = s.len();
    if n < k {
        println!("\nk = {} → 0 distinct (sequence too short)", k);
        return;
    }

    let mut seen: HashSet<&str> = HashSet::new();
    let mut uniques: Vec<&str> = Vec::new();

    for i in 0..=n - k {
        let kmer = &s[i..i + k];
        if seen.insert(kmer) {
            uniques.push(kmer);
        }
    }

    println!("\nk = {} → {} distinct", k, uniques.len());
    for u in uniques {
        println!("{u}");
    }
}
