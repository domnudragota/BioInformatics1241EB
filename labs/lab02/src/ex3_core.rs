// labs/lab02/src/ex3_core.rs

use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[cfg(feature = "gui")]
use slint::{Image, Rgba8Pixel, SharedPixelBuffer};

// Limits / layout
pub const MAX_WINDOWS_FULL: usize = 2_000_000; // above this, stream + downsample
// Supersampling factor for sharper charts without resizing the UI.
// 1 = current quality, 2 = ~8 MB RGBA buffer, 3 = ~18 MB (use only if you have RAM/CPU headroom).
pub const CHART_SS: u32 = 2;
pub const CHART_BASE_W: u32 = 1200;
pub const CHART_BASE_H: u32 = 420;
// Actual render size (higher than what we display -> smoother)
pub const CHART_W: u32 = CHART_BASE_W * CHART_SS;
pub const CHART_H: u32 = CHART_BASE_H * CHART_SS;

#[inline]
pub fn num_windows(n_bases: usize, win: usize, step: usize) -> usize {
    if n_bases < win { 0 } else { 1 + (n_bases - win) / step }
}

// ---------- FASTA reader (.gz supported) ----------
pub fn parse_fasta_bytes<P: AsRef<Path>>(path: P) -> Result<Vec<u8>, String> {
    let p = path.as_ref();
    let file = File::open(p).map_err(|e| format!("open error: {e}"))?;
    let gz = matches!(p.extension().and_then(|e| e.to_str()), Some("gz"));

    let mut reader: Box<dyn BufRead> = if gz {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    // Parse FASTA: skip headers ('>'), keep only A/C/G/T (uppercase)
    let mut seq: Vec<u8> = Vec::with_capacity(1 << 20);
    let mut line = Vec::<u8>::with_capacity(1 << 12);
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line).map_err(|e| format!("read error: {e}"))?;
        if n == 0 { break; }
        if !line.is_empty() && line[0] == b'>' { continue; } // header
        for &b in &line {
            let b = b.to_ascii_uppercase();
            match b { b'A' | b'C' | b'G' | b'T' => seq.push(b), _ => {} }
        }
    }
    if seq.is_empty() { return Err("No A/C/G/T symbols found in file".into()); }
    Ok(seq)
}

// ---------- Sliding-window (rolling) ----------
#[derive(Debug, Clone)]
pub struct Series { pub a: Vec<f32>, pub c: Vec<f32>, pub g: Vec<f32>, pub t: Vec<f32> }

#[inline]
fn idx(b: u8) -> usize {
    match b { b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => unreachable!(), }
}

pub fn sliding_freqs(seq: &[u8], win: usize, step: usize) -> Result<Series, String> {
    if win == 0 || step == 0 { return Err("Window and step must be > 0".into()); }
    if seq.len() < win { return Err(format!("Sequence too short for window={win}")); }

    let n = seq.len();
    let mut counts = [0usize; 4];
    for &b in &seq[..win] { counts[idx(b)] += 1; }

    let cap = 1 + (n - win) / step;
    let mut a = Vec::with_capacity(cap);
    let mut c = Vec::with_capacity(cap);
    let mut g = Vec::with_capacity(cap);
    let mut t = Vec::with_capacity(cap);

    let wlen = win as f32;
    a.push(counts[0] as f32 / wlen);
    c.push(counts[1] as f32 / wlen);
    g.push(counts[2] as f32 / wlen);
    t.push(counts[3] as f32 / wlen);

    let mut i = step;
    while i + win <= n {
        for &b in &seq[i - step .. i]     { counts[idx(b)] -= 1; }
        for &b in &seq[i + win - step .. i + win] { counts[idx(b)] += 1; }

        a.push(counts[0] as f32 / wlen);
        c.push(counts[1] as f32 / wlen);
        g.push(counts[2] as f32 / wlen);
        t.push(counts[3] as f32 / wlen);

        i += step;
    }
    Ok(Series { a, c, g, t })
}

// Streamed version: compute bucket-averages; memory O(width).
pub fn sliding_freqs_downsampled(seq: &[u8], win: usize, step: usize, width: usize) -> Result<Series, String> {
    if win == 0 || step == 0 { return Err("Window and step must be > 0".into()); }
    if seq.len() < win { return Err(format!("Sequence too short for window={win}")); }
    if width == 0 { return Err("Width must be > 0".into()); }

    let total = num_windows(seq.len(), win, step);
    let bucket = (total + width - 1) / width; // ceil

    let mut counts = [0usize; 4];
    for &b in &seq[..win] { counts[idx(b)] += 1; }

    let mut a = Vec::with_capacity(width);
    let mut c = Vec::with_capacity(width);
    let mut g = Vec::with_capacity(width);
    let mut t = Vec::with_capacity(width);

    let wlen = win as f32;
    let mut acc = [0f32; 4];
    let mut in_bucket = 0usize;

    let flush = |acc: [f32;4], cnt: usize, a:&mut Vec<f32>, c:&mut Vec<f32>, g:&mut Vec<f32>, t:&mut Vec<f32>| {
        if cnt == 0 { return; }
        let f = 1.0 / cnt as f32;
        a.push(acc[0]*f); c.push(acc[1]*f); g.push(acc[2]*f); t.push(acc[3]*f);
    };

    acc[0] += counts[0] as f32 / wlen;
    acc[1] += counts[1] as f32 / wlen;
    acc[2] += counts[2] as f32 / wlen;
    acc[3] += counts[3] as f32 / wlen;
    in_bucket += 1;

    let mut i = step;
    while i + win <= seq.len() {
        for &b in &seq[i - step .. i]          { counts[idx(b)] -= 1; }
        for &b in &seq[i + win - step .. i + win] { counts[idx(b)] += 1; }

        acc[0] += counts[0] as f32 / wlen;
        acc[1] += counts[1] as f32 / wlen;
        acc[2] += counts[2] as f32 / wlen;
        acc[3] += counts[3] as f32 / wlen;
        in_bucket += 1;

        if in_bucket == bucket {
            flush(acc, in_bucket, &mut a, &mut c, &mut g, &mut t);
            acc = [0.0; 4];
            in_bucket = 0;
        }
        i += step;
    }
    if in_bucket > 0 { flush(acc, in_bucket, &mut a, &mut c, &mut g, &mut t); }
    Ok(Series { a, c, g, t })
}

// ---------- Chart rendering (no extra crates) ----------
#[inline]
fn clamp01(x: f32) -> f32 { if x < 0.0 { 0.0 } else if x > 1.0 { 1.0 } else { x } }

fn downsample_to(v: &[f32], w: usize) -> Vec<f32> {
    if v.is_empty() || w == 0 { return Vec::new(); }
    if v.len() <= w { return v.to_vec(); }
    let bucket = (v.len() + w - 1) / w; // ceil
    let mut out = Vec::with_capacity(w);
    let mut i = 0usize;
    while i < v.len() && out.len() < w {
        let end = (i + bucket).min(v.len());
        let mut s = 0.0; let mut c = 0usize;
        for &x in &v[i..end] { s += x; c += 1; }
        out.push(clamp01(s / c as f32));
        i = end;
    }
    out
}

fn draw_line(pix: &mut [u8], stride: usize, w: i32, h: i32,
             mut x0: i32, mut y0: i32, x1: i32, y1: i32, col: [u8;4]) {
    let dx = (x1 - x0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let dy = -(y1 - y0).abs();
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    loop {
        if (0..w).contains(&x0) && (0..h).contains(&y0) {
            let idx = (y0 as usize) * stride + (x0 as usize) * 4;
            pix[idx]   = col[0]; pix[idx+1] = col[1]; pix[idx+2] = col[2]; pix[idx+3] = col[3];
        }
        if x0 == x1 && y0 == y1 { break; }
        let e2 = 2 * err;
        if e2 >= dy { err += dy; x0 += sx; }
        if e2 <= dx { err += dx; y0 += sy; }
    }
}

#[cfg(feature = "gui")]
pub fn render_chart_image(series: &Series, width: u32, height: u32) -> Image {
    let w = width as i32; let h = height as i32;
    let margin = 8i32;
    let plot_w = (w - 2*margin).max(1);
    let plot_h = (h - 2*margin).max(1);

    // downsample each to <= width
    let a = downsample_to(&series.a, width as usize);
    let c = downsample_to(&series.c, width as usize);
    let g = downsample_to(&series.g, width as usize);
    let t = downsample_to(&series.t, width as usize);

    let mut buf = SharedPixelBuffer::<Rgba8Pixel>::new(width, height);
    let stride = (width as usize) * 4;
    let pixels = buf.make_mut_bytes();

    // background
    for y in 0..(height as usize) {
        let row = &mut pixels[y*stride..(y+1)*stride];
        for x in (0..row.len()).step_by(4) { row[x]=0x12; row[x+1]=0x12; row[x+2]=0x12; row[x+3]=0xFF; }
    }

    // grid (0.0, 0.5, 1.0)
    let grid_col = [0x33,0x33,0x33,0xFF];
    for gv in [0.0f32, 0.25, 0.5, 0.75, 1.0] {
        let y = margin + (plot_h as f32 - 1.0 - gv * (plot_h as f32 - 1.0)).round() as i32;
        draw_line(pixels, stride, w, h, margin, y, w - margin - 1, y, grid_col);
    }

    // colors for A/C/G/T (A changed to purple for better contrast with T/red)
    let col_a = [0x8E, 0x44, 0xAD, 0xFF]; // A = purple
    let col_c = [0x26, 0x8B, 0xD2, 0xFF]; // C = blue
    let col_g = [0x85, 0x99, 0x00, 0xFF]; // G = green
    let col_t = [0xDC, 0x32, 0x2F, 0xFF]; // T = red


    let mut plot = |v: &[f32], col: [u8;4]| {
        if v.len() < 2 { return; }
        let last = (v.len()-1) as i32;
        for i in 0..last {
            let x0 = margin + (i * (plot_w - 1)) / last;
            let x1 = margin + ((i+1) * (plot_w - 1)) / last;
            let y0 = margin + (plot_h - 1) - ((v[i as usize] * (plot_h - 1) as f32).round() as i32);
            let y1 = margin + (plot_h - 1) - ((v[(i+1) as usize] * (plot_h - 1) as f32).round() as i32);
            draw_line(pixels, stride, w, h, x0, y0, x1, y1, col);
        }
    };
    plot(&a, col_a); plot(&c, col_c); plot(&g, col_g); plot(&t, col_t);

    Image::from_rgba8(buf)
}
