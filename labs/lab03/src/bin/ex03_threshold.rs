use anyhow::{bail, Context, Result};
use flate2::read::GzDecoder;
use plotters::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

/// Tm formulas (same as ex02)
#[derive(Clone, Copy)]
enum TmMethod {
    Wallace,           // Tm = 4*(G+C) + 2*(A+T)
    SaltAdjusted(f64), // Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) – 600/len
}

fn main() -> Result<()> {
    println!("[lab03] ex03_threshold (single-pass, k=9, cutoff bars)");

    // CLI:
    // ex03_threshold <input.fa[.gz]>
    //   [--method wallace|salt] [--na 0.05]
    //   [--out tm_cutoff.png] [--width 1200] [--binw 10000]
    //   [--th v1[,v2,...]] [--progress]
    let mut args = std::env::args().skip(1);
    let input_path = args
        .next()
        .unwrap_or_else(|| prompt("Path to FASTA (.fa/.fasta/.gz): ").trim().to_string());

    let mut na_m: f64 = 0.05;
    let mut method = TmMethod::Wallace;
    let mut out_path: String = "tm_cutoff.png".to_string();
    let mut width_target: usize = 1200;      // final plotted bins
    let mut windows_per_bin: usize = 10_000; // streaming bin size (pre-merge)
    let mut thresholds: Option<Vec<f64>> = None;
    let mut progress = false;

    while let Some(flag) = args.next() {
        match flag.as_str() {
            "--method" => {
                let v = args.next().unwrap_or_else(|| "wallace".to_string());
                method = if v.eq_ignore_ascii_case("salt") || v.eq_ignore_ascii_case("salt-adjusted")
                {
                    TmMethod::SaltAdjusted(na_m)
                } else {
                    TmMethod::Wallace
                };
            }
            "--na" => {
                na_m = args
                    .next()
                    .and_then(|s| s.parse::<f64>().ok())
                    .unwrap_or(0.05);
                if let TmMethod::SaltAdjusted(_) = method {
                    method = TmMethod::SaltAdjusted(na_m);
                }
            }
            "--out" => {
                out_path = args.next().unwrap_or_else(|| "tm_cutoff.png".to_string());
            }
            "--width" => {
                width_target = args
                    .next()
                    .and_then(|s| s.parse::<usize>().ok())
                    .unwrap_or(1200)
                    .max(32);
            }
            "--binw" => {
                windows_per_bin = args
                    .next()
                    .and_then(|s| s.parse::<usize>().ok())
                    .unwrap_or(10_000)
                    .max(256);
            }
            "--th" | "--threshold" => {
                let raw = args.next().unwrap_or_default();
                let parsed: Vec<f64> = raw
                    .split(',')
                    .filter_map(|t| t.trim().parse::<f64>().ok())
                    .collect();
                if !parsed.is_empty() {
                    thresholds = Some(parsed);
                }
            }
            "--progress" => progress = true,
            _ => {}
        }
    }

    let k: usize = 9;

    // Single-pass stream (shared logic with ex02)
    let stats = compute_single_pass(&input_path, k, method, windows_per_bin, progress)?;

    if stats.total_valid < k {
        bail!(
            "Sequence has only {} valid A/C/G/T bases (< k={}).",
            stats.total_valid,
            k
        );
    }

    // Downsample to desired width
    let (bins_final, y_min, y_max) = downsample_bins(&stats.bins, width_target);

    // Build X,Y for plotting (avg per bin)
    let xs: Vec<i32> = (0..bins_final.len()).map(|i| i as i32).collect();
    let ys: Vec<f64> = bins_final
        .iter()
        .map(|b| if b.count > 0 { b.sum / (b.count as f64) } else { 0.0 })
        .collect();

    // Default cutoff if none provided: mean Tm
    let ths: Vec<f64> = match thresholds {
        Some(v) => v,
        None => {
            let mean = if ys.is_empty() { 0.0 } else { ys.iter().sum::<f64>() / (ys.len() as f64) };
            vec![mean]
        }
    };

    plot_with_thresholds(
        &xs,
        &ys,
        (y_min, y_max),
        &ths,
        &out_path,
        "Tm (k=9) — cutoff bars",
        "Position (bin)",
        "Tm (°C)",
    )?;

    println!(
        "Done. Valid bases: {} | windows: {} | bins: {} → {} | Tm range: {:.2}..{:.2} °C",
        stats.total_valid,
        stats.n_windows,
        stats.bins.len(),
        xs.len(),
        y_min,
        y_max
    );
    println!("Thresholds: {:?}", ths);
    println!("Saved chart → {out_path}");
    Ok(())
}

/// ==== streaming & plotting infra (optimized, borrow-safe) ====

#[derive(Clone, Copy)]
struct BinAgg {
    sum: f64,
    min: f64,
    max: f64,
    count: u32,
}
impl BinAgg {
    fn new() -> Self {
        Self { sum: 0.0, min: f64::INFINITY, max: f64::NEG_INFINITY, count: 0 }
    }
    fn push(&mut self, v: f64) {
        self.sum += v;
        if v < self.min { self.min = v; }
        if v > self.max { self.max = v; }
        self.count += 1;
    }
}

struct StreamStats {
    bins: Vec<BinAgg>,
    total_valid: usize,
    n_windows: usize,
}

fn compute_single_pass<P: AsRef<Path>>(
    path: P,
    k: usize,
    method: TmMethod,
    windows_per_bin: usize,
    progress: bool,
) -> Result<StreamStats> {
    let mut br = open_reader(path)?;
    let mut bins: Vec<BinAgg> = Vec::with_capacity(4096);
    let mut current = BinAgg::new();

    // Precomputed constants
    let wallace_base = (2 * k) as f64; // 2*k + 2*GC
    let (salt_c0, salt_c1) = match method {
        TmMethod::SaltAdjusted(na) => (
            81.5 + 16.6 * na.log10() - 600.0 / (k as f64),
            41.0 / (k as f64),
        ),
        _ => (0.0, 0.0),
    };

    // GC ring (k=9 -> fixed-size array for speed)
    let mut ring_gc = [0u8; 9];
    let mut ring_filled = 0usize;
    let mut ring_idx = 0usize;
    let mut gc: i32 = 0;

    #[inline]
    fn gc_flag(b: u8) -> Option<u8> {
        match b {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(1),
            b'T' | b't' => Some(0),
            _ => None,
        }
    }

    // FASTA header skipping (chunk-friendly)
    let mut at_line_start = true;
    let mut in_header = false;

    let mut total_valid: usize = 0;
    let mut n_windows: usize = 0;

    let mut processed_bytes: usize = 0;
    let mut next_print: usize = 64 * 1024 * 1024; // 64 MiB

    loop {
        let consumed;
        {
            let buf = br.fill_buf()?;
            if buf.is_empty() { break; }

            processed_bytes += buf.len();
            if progress && processed_bytes >= next_print {
                eprintln!(" .. processed ~{} MiB", processed_bytes / (1024 * 1024));
                next_print += 64 * 1024 * 1024;
            }

            let mut i = 0;
            while i < buf.len() {
                let b = buf[i];

                // Header lines start with '>' at start of a line and go until newline
                if in_header {
                    if b == b'\n' { in_header = false; at_line_start = true; }
                    i += 1;
                    continue;
                }
                if at_line_start && b == b'>' {
                    in_header = true;
                    i += 1;
                    continue;
                }
                if b == b'\n' {
                    at_line_start = true;
                    i += 1;
                    continue;
                } else if b == b'\r' {
                    i += 1;
                    continue;
                } else if at_line_start {
                    at_line_start = false;
                }

                if let Some(f) = gc_flag(b) {
                    total_valid += 1;

                    if ring_filled < k {
                        ring_gc[ring_filled] = f;
                        gc += f as i32;
                        ring_filled += 1;

                        if ring_filled == k {
                            let tm = tm_from_gc(gc, wallace_base, salt_c0, salt_c1, method);
                            current.push(tm);
                            n_windows += 1;

                            if current.count as usize >= windows_per_bin {
                                bins.push(current);
                                current = BinAgg::new();
                            }
                        }
                    } else {
                        let f_old = ring_gc[ring_idx];
                        gc -= f_old as i32;
                        ring_gc[ring_idx] = f;
                        gc += f as i32;
                        ring_idx = (ring_idx + 1) % k;

                        let tm = tm_from_gc(gc, wallace_base, salt_c0, salt_c1, method);
                        current.push(tm);
                        n_windows += 1;

                        if current.count as usize >= windows_per_bin {
                            bins.push(current);
                            current = BinAgg::new();
                        }
                    }
                }

                i += 1;
            }

            consumed = buf.len(); // record length while borrow is active
        } // buf borrow ends here

        br.consume(consumed); // safe to consume now
    }

    if current.count > 0 {
        bins.push(current);
    }

    Ok(StreamStats { bins, total_valid, n_windows })
}

fn downsample_bins(bins: &[BinAgg], target: usize) -> (Vec<BinAgg>, f64, f64) {
    if bins.is_empty() {
        return (vec![], 0.0, 1.0);
    }
    if bins.len() <= target {
        let mut ymin = f64::INFINITY;
        let mut ymax = f64::NEG_INFINITY;
        for b in bins {
            if b.count > 0 {
                let v = b.sum / (b.count as f64);
                if v < ymin { ymin = v; }
                if v > ymax { ymax = v; }
            }
        }
        return (bins.to_vec(), ymin, ymax);
    }

    let group = (bins.len() + target - 1) / target; // ceil
    let mut out = Vec::with_capacity(target);
    let mut ymin = f64::INFINITY;
    let mut ymax = f64::NEG_INFINITY;

    let mut i = 0;
    while i < bins.len() {
        let mut agg = BinAgg::new();
        let end = (i + group).min(bins.len());
        for j in i..end {
            let b = bins[j];
            if b.count > 0 {
                agg.sum += b.sum;
                agg.count += b.count;
                if b.min < agg.min { agg.min = b.min; }
                if b.max > agg.max { agg.max = b.max; }
            }
        }
        if agg.count > 0 {
            let v = agg.sum / (agg.count as f64);
            if v < ymin { ymin = v; }
            if v > ymax { ymax = v; }
        }
        out.push(agg);
        i = end;
    }
    (out, ymin, ymax)
}

#[inline]
fn tm_from_gc(gc: i32, wallace_base: f64, salt_c0: f64, salt_c1: f64, method: TmMethod) -> f64 {
    match method {
        TmMethod::Wallace => wallace_base + 2.0 * (gc as f64), // 2*k + 2*GC
        TmMethod::SaltAdjusted(_) => salt_c0 + salt_c1 * (gc as f64),
    }
}

fn open_reader<P: AsRef<Path>>(path: P) -> Result<BufReader<Box<dyn Read>>> {
    let p = path.as_ref();
    let f = File::open(p).with_context(|| format!("open {}", p.display()))?;
    let boxed: Box<dyn Read> = match p.extension().and_then(|e| e.to_str()) {
        Some(ext) if ext.eq_ignore_ascii_case("gz") || ext.eq_ignore_ascii_case("gzip") => {
            Box::new(GzDecoder::new(f))
        }
        _ => Box::new(f),
    };
    Ok(BufReader::with_capacity(1 << 20, boxed)) // 1 MiB buffer
}

fn plot_with_thresholds(
    xs: &[i32],
    ys: &[f64],
    y_range: (f64, f64),
    thresholds: &[f64],
    out: &str,
    title: &str,
    xlab: &str,
    ylab: &str,
) -> Result<()> {
    let w: u32 = 1200;
    let h: u32 = 420;

    let root = BitMapBackend::new(out, (w, h)).into_drawing_area();
    root.fill(&WHITE)?;

    let (ymin, ymax) = y_range;
    let pad = ((ymax - ymin) * 0.05).max(1.0);
    let y_low = (ymin - pad).floor();
    let y_high = (ymax + pad).ceil();

    let x_max = (xs.len().saturating_sub(1)) as i32;

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 24))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0..x_max, y_low..y_high)?;

    chart.configure_mesh().x_desc(xlab).y_desc(ylab).draw()?;

    // Main series
    chart.draw_series(LineSeries::new(
        xs.iter().zip(ys.iter()).map(|(x, y)| (*x, *y)),
        &BLUE,
    ))?;

    // Horizontal threshold bars
    for (i, &t) in thresholds.iter().enumerate() {
        let color = match i % 5 {
            0 => &RED,
            1 => &GREEN,
            2 => &BLACK,
            3 => &MAGENTA,
            _ => &CYAN,
        };
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(0, t), (x_max, t)],
            ShapeStyle::from(color).stroke_width(3),
        )))?;
    }

    root.present()?;
    Ok(())
}

fn prompt(msg: &str) -> String {
    print!("{msg}");
    let _ = io::stdout().flush();
    let mut s = String::new();
    io::stdin().read_line(&mut s).ok();
    s
}
