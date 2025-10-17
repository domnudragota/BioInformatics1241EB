use std::env;
use std::io::{self, Write};

/// Compute Tm using two formulas:
/// 1) Wallace rule (good for short primers): Tm = 4*(G+C) + 2*(A+T)
/// 2) Salt-adjusted: Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) – 600/length
fn main() {
    println!("[lab03] ex01_melting_temp");

    // Parse CLI args if provided: <DNA> [Na+ in molar]
    let mut args = env::args().skip(1);
    let mut seq = args.next();
    let mut na_in_m = args
        .next()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.05); // default 50 mM

    // If no DNA provided, prompt the user (keeps Makefile super simple)
    if seq.is_none() {
        seq = Some(prompt("Enter DNA sequence (A,C,G,T only): "));
        let maybe_na = prompt("Optional [Na+] in mol/L (press Enter for 0.05): ");
        if !maybe_na.trim().is_empty() {
            na_in_m = maybe_na.trim().parse::<f64>().unwrap_or(0.05);
        }
    }

    let seq = seq.unwrap_or_default().trim().to_string();
    if seq.is_empty() {
        eprintln!("Error: empty DNA sequence.");
        std::process::exit(1);
    }

    // Normalize to uppercase and validate characters
    let seq = seq.to_uppercase();
    if let Some(bad) = seq.chars().find(|&c| !matches!(c, 'A' | 'C' | 'G' | 'T')) {
        eprintln!("Error: invalid character '{bad}'. Allowed: A,C,G,T.");
        std::process::exit(1);
    }

    let (a, c, g, t) = count_bases(&seq);
    let len = (a + c + g + t) as usize;
    if len == 0 {
        eprintln!("Error: sequence length is zero after validation.");
        std::process::exit(1);
    }

    // Wallace rule
    let tm_wallace = 4.0 * (g + c) as f64 + 2.0 * (a + t) as f64;

    // Salt-adjusted formula
    let gc_percent = 100.0 * (g + c) as f64 / (len as f64);
    let na = if na_in_m > 0.0 { na_in_m } else { 0.05 }; // guard log10
    let tm_salt = 81.5 + 16.6 * na.log10() + 0.41 * gc_percent - 600.0 / (len as f64);

    // Output
    println!();
    println!("Sequence length: {len}");
    println!("A:{a}  C:{c}  G:{g}  T:{t}  GC%: {:.2}", gc_percent);
    println!("-----------------------------------------------");
    println!("Tm (Wallace rule):           {:.2} °C", tm_wallace);
    println!("Tm (salt-adjusted, [Na+]={:.3} M): {:.2} °C", na, tm_salt);

    println!();
    println!("Tip: you can also pass args directly:");
    println!("  cargo run --bin ex01_tm -- ACTG 0.1");
}

fn prompt(msg: &str) -> String {
    print!("{msg}");
    let _ = io::stdout().flush();
    let mut s = String::new();
    io::stdin().read_line(&mut s).ok();
    s
}

fn count_bases(seq: &str) -> (u32, u32, u32, u32) {
    let mut a = 0u32;
    let mut c = 0u32;
    let mut g = 0u32;
    let mut t = 0u32;
    for ch in seq.chars() {
        match ch {
            'A' => a += 1,
            'C' => c += 1,
            'G' => g += 1,
            'T' => t += 1,
            _ => {} // already validated
        }
    }
    (a, c, g, t)
}
