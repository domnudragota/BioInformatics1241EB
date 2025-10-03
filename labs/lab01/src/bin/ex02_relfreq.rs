fn main() {
    let s = "ATTTCGCCGATA";

    // Alphabet in first-appearance order
    let mut alphabet = String::new();
    for ch in s.chars() {
        if !alphabet.contains(ch) {
            alphabet.push(ch);
        }
    }

    // Counts per symbol aligned with `alphabet`
    let mut counts = vec![0usize; alphabet.len()];
    let mut total = 0usize;
    for ch in s.chars() {
        if let Some(pos) = alphabet.find(ch) {
            counts[pos] += 1;
            total += 1;
        }
    }

    println!("Sequence: {}", s);
    println!("Alphabet: {}", alphabet);
    println!("Total symbols: {}", total);
    for (i, ch) in alphabet.chars().enumerate() {
        let rel = counts[i] as f64 / total as f64;
        println!("{} -> count = {}, rel = {:.4}", ch, counts[i], rel);
    }
}
