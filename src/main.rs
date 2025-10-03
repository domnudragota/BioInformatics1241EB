fn main() {
    let s = "ATTTCGCCGATA";
    let mut text = String::new();

    for ch in s.chars() {
        if !text.contains(ch) {
            text.push(ch);
        }
    }

    println!("Alphabet: {}", text);
    println!("Size: {}", text.len());
}
