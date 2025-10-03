fn main() {
    let s = "ATTTCG*CCGATA";
    let mut text = String::new();

    // 1st exercise
    for ch in s.chars() {
        if !text.contains(ch) {
            text.push(ch);
        }
    }
    // 2nd exercise
    let mut counts = vec![0usize; text.len()];
    let mut total = 0usize;
    for ch in s.chars() {
        if let Some(pos) = text.find(ch) {
            counts[pos] += 1;
            total += 1;
        }
    }

    println!("Alphabet: {}", text);
    println!("Size: {}", text.len());
    println!("Total symbols: {}", total);
    for (i, ch) in text.chars().enumerate() {
        let rel = counts[i] as f64 / total as f64;
        println!("{} = count = {}, rel = {:.4}", ch, counts[i], rel);
    }
}
