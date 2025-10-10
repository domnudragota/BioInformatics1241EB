// labs/lab02/src/bin/ex3_fasta_gui.rs

#[cfg(not(feature = "gui"))]
fn main() {
    eprintln!("[lab02] ex3_fasta_gui: build with `--features gui`");
    eprintln!("Example: make ex3   # (uses --features gui --release)");
}

#[cfg(feature = "gui")]
mod gui_app {
    use flate2::read::GzDecoder;
    use rfd::FileDialog;
    use slint::{
        invoke_from_event_loop, ComponentHandle, Image, Rgba8Pixel, SharedPixelBuffer, SharedString,
    };
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::Path;

    // -------------------- Limits / layout ------------------------------------

    // Limit how much we keep/display in memory
    const MAX_WINDOWS_FULL: usize = 2_000_000; // above this, stream + downsample
    const SPARK_MAX_POINTS: usize = 800; // cap sparkline width

    // Chart image size (in pixels)
    const CHART_W: u32 = 860;
    const CHART_H: u32 = 220;

    #[inline]
    fn num_windows(n_bases: usize, win: usize, step: usize) -> usize {
        if n_bases < win {
            0
        } else {
            1 + (n_bases - win) / step
        }
    }

    // -------------------- FASTA reader (.gz supported) -----------------------

    fn parse_fasta_bytes<P: AsRef<Path>>(path: P) -> Result<Vec<u8>, String> {
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
            let n = reader
                .read_until(b'\n', &mut line)
                .map_err(|e| format!("read error: {e}"))?;
            if n == 0 {
                break;
            }
            if !line.is_empty() && line[0] == b'>' {
                continue; // header line
            }
            for &b in &line {
                let b = b.to_ascii_uppercase();
                match b {
                    b'A' | b'C' | b'G' | b'T' => seq.push(b),
                    _ => {}
                }
            }
        }
        if seq.is_empty() {
            return Err("No A/C/G/T symbols found in file".into());
        }
        Ok(seq)
    }

    // -------------------- Sliding-window (rolling) ---------------------------

    #[derive(Debug, Clone)]
    struct Series {
        a: Vec<f32>,
        c: Vec<f32>,
        g: Vec<f32>,
        t: Vec<f32>,
    }

    #[inline]
    fn idx(b: u8) -> usize {
        match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => unreachable!(),
        }
    }

    fn sliding_freqs(seq: &[u8], win: usize, step: usize) -> Result<Series, String> {
        if win == 0 || step == 0 {
            return Err("Window and step must be > 0".into());
        }
        if seq.len() < win {
            return Err(format!("Sequence too short for window={win}"));
        }

        let n = seq.len();
        let mut counts = [0usize; 4];

        // init first window
        for &b in &seq[..win] {
            counts[idx(b)] += 1;
        }

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

        // slide by `step`
        let mut i = step;
        while i + win <= n {
            // remove outgoing
            for &b in &seq[i - step..i] {
                counts[idx(b)] -= 1;
            }
            // add incoming
            for &b in &seq[i + win - step..i + win] {
                counts[idx(b)] += 1;
            }

            a.push(counts[0] as f32 / wlen);
            c.push(counts[1] as f32 / wlen);
            g.push(counts[2] as f32 / wlen);
            t.push(counts[3] as f32 / wlen);

            i += step;
        }

        Ok(Series { a, c, g, t })
    }

    // Streamed version: compute bucket-averages; memory O(width).
    fn sliding_freqs_downsampled(
        seq: &[u8],
        win: usize,
        step: usize,
        width: usize,
    ) -> Result<Series, String> {
        if win == 0 || step == 0 {
            return Err("Window and step must be > 0".into());
        }
        if seq.len() < win {
            return Err(format!("Sequence too short for window={win}"));
        }
        if width == 0 {
            return Err("Width must be > 0".into());
        }

        let total = num_windows(seq.len(), win, step);
        let bucket = (total + width - 1) / width; // ceil

        // rolling counts for current window
        let mut counts = [0usize; 4];
        for &b in &seq[..win] {
            counts[idx(b)] += 1;
        }

        let mut a = Vec::with_capacity(width);
        let mut c = Vec::with_capacity(width);
        let mut g = Vec::with_capacity(width);
        let mut t = Vec::with_capacity(width);

        let wlen = win as f32;

        let mut acc = [0f32; 4];
        let mut in_bucket = 0usize;

        // helper to flush a bucket average
        let mut flush = |acc: [f32; 4],
                         cnt: usize,
                         a: &mut Vec<f32>,
                         c: &mut Vec<f32>,
                         g: &mut Vec<f32>,
                         t: &mut Vec<f32>| {
            if cnt == 0 {
                return;
            }
            let f = 1.0 / cnt as f32;
            a.push(acc[0] * f);
            c.push(acc[1] * f);
            g.push(acc[2] * f);
            t.push(acc[3] * f);
        };

        // first window
        acc[0] += counts[0] as f32 / wlen;
        acc[1] += counts[1] as f32 / wlen;
        acc[2] += counts[2] as f32 / wlen;
        acc[3] += counts[3] as f32 / wlen;
        in_bucket += 1;

        // slide through the rest
        let mut i = step;
        while i + win <= seq.len() {
            // remove outgoing step bytes
            for &b in &seq[i - step..i] {
                counts[idx(b)] -= 1;
            }
            // add incoming step bytes
            for &b in &seq[i + win - step..i + win] {
                counts[idx(b)] += 1;
            }

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

        // tail
        if in_bucket > 0 {
            flush(acc, in_bucket, &mut a, &mut c, &mut g, &mut t);
        }

        Ok(Series { a, c, g, t })
    }

    // -------------------- Sparklines (downsampled) ---------------------------

    fn sparkline_downsample(v: &[f32], max_points: usize) -> String {
        const BLOCKS: &[char] = &['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];
        if v.is_empty() || max_points == 0 {
            return String::new();
        }

        let len = v.len();
        if len <= max_points {
            let mut out = String::with_capacity(len);
            for &x in v {
                let xx = x.clamp(0.0, 1.0);
                let idx = (xx * (BLOCKS.len() as f32 - 1.0)).round() as usize;
                out.push(BLOCKS[idx]);
            }
            return out;
        }

        let bucket = (len + max_points - 1) / max_points; // ceil
        let mut out = String::with_capacity(max_points);
        let mut i = 0usize;
        while i < len {
            let end = (i + bucket).min(len);
            let mut sum = 0.0f32;
            let mut cnt = 0usize;
            for &x in &v[i..end] {
                sum += x;
                cnt += 1;
            }
            let avg = (sum / cnt as f32).clamp(0.0, 1.0);
            let idx = (avg * (BLOCKS.len() as f32 - 1.0)).round() as usize;
            out.push(BLOCKS[idx]);
            i = end;
        }
        out
    }

    // -------------------- Chart rendering (no extra crates) ------------------

    #[inline]
    fn clamp01(x: f32) -> f32 {
        if x < 0.0 {
            0.0
        } else if x > 1.0 {
            1.0
        } else {
            x
        }
    }

    fn downsample_to(v: &[f32], w: usize) -> Vec<f32> {
        if v.is_empty() || w == 0 {
            return Vec::new();
        }
        if v.len() <= w {
            return v.to_vec();
        }
        let bucket = (v.len() + w - 1) / w; // ceil
        let mut out = Vec::with_capacity(w);
        let mut i = 0usize;
        while i < v.len() && out.len() < w {
            let end = (i + bucket).min(v.len());
            let mut s = 0.0;
            let mut c = 0usize;
            for &x in &v[i..end] {
                s += x;
                c += 1;
            }
            out.push(clamp01(s / c as f32));
            i = end;
        }
        out
    }

    fn draw_line(
        pix: &mut [u8],
        stride: usize,
        w: i32,
        h: i32,
        mut x0: i32,
        mut y0: i32,
        mut x1: i32,
        mut y1: i32,
        col: [u8; 4],
    ) {
        let mut dx = (x1 - x0).abs();
        let sx = if x0 < x1 { 1 } else { -1 };
        let mut dy = -(y1 - y0).abs();
        let sy = if y0 < y1 { 1 } else { -1 };
        let mut err = dx + dy;
        loop {
            if (0..w).contains(&x0) && (0..h).contains(&y0) {
                let idx = (y0 as usize) * stride + (x0 as usize) * 4;
                pix[idx] = col[0];
                pix[idx + 1] = col[1];
                pix[idx + 2] = col[2];
                pix[idx + 3] = col[3];
            }
            if x0 == x1 && y0 == y1 {
                break;
            }
            let e2 = 2 * err;
            if e2 >= dy {
                err += dy;
                x0 += sx;
            }
            if e2 <= dx {
                err += dx;
                y0 += sy;
            }
        }
    }

    fn render_chart_image(series: &Series, width: u32, height: u32) -> Image {
        let w = width as i32;
        let h = height as i32;
        let margin = 8i32;
        let plot_w = (w - 2 * margin).max(1);
        let plot_h = (h - 2 * margin).max(1);

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
            let row = &mut pixels[y * stride..(y + 1) * stride];
            for x in (0..row.len()).step_by(4) {
                row[x] = 0x12;
                row[x + 1] = 0x12;
                row[x + 2] = 0x12;
                row[x + 3] = 0xFF;
            }
        }

        // grid (0.0, 0.5, 1.0)
        let grid_col = [0x33, 0x33, 0x33, 0xFF];
        for gv in [0.0f32, 0.5, 1.0] {
            let y = margin + (plot_h as f32 - 1.0 - gv * (plot_h as f32 - 1.0)).round() as i32;
            draw_line(pixels, stride, w, h, margin, y, w - margin - 1, y, grid_col);
        }

        // colors for A/C/G/T
        let col_a = [0xE6, 0x4A, 0x19, 0xFF]; // orange
        let col_c = [0x26, 0x8B, 0xD2, 0xFF]; // blue
        let col_g = [0x85, 0x99, 0x00, 0xFF]; // green
        let col_t = [0xDC, 0x32, 0x2F, 0xFF]; // red

        let mut plot = |v: &[f32], col: [u8; 4]| {
            if v.len() < 2 {
                return;
            }
            let last = (v.len() - 1) as i32;
            for i in 0..last {
                let x0 = margin + (i * (plot_w - 1)) / last;
                let x1 = margin + ((i + 1) * (plot_w - 1)) / last;
                let y0 = margin + (plot_h - 1) - ((v[i as usize] * (plot_h - 1) as f32).round() as i32);
                let y1 = margin
                    + (plot_h - 1)
                    - ((v[(i + 1) as usize] * (plot_h - 1) as f32).round() as i32);
                draw_line(pixels, stride, w, h, x0, y0, x1, y1, col);
            }
        };

        plot(&a, col_a);
        plot(&c, col_c);
        plot(&g, col_g);
        plot(&t, col_t);

        Image::from_rgba8(buf)
    }

    // -------------------- Slint UI -------------------------------------------

    slint::slint! {
        import { Button, SpinBox } from "std-widgets.slint";

        export component AppWindow inherits Window {
            width: 900px;
            height: 700px;
            title: "Lab02 - FASTA Sliding Window Frequencies";

            in-out property <int> window_size: 5;
            in-out property <int> step_size: 1;

            in-out property <string> file_path: "";
            in-out property <string> status: "Pick a FASTA and press Compute";

            in-out property <string> vec_a: "";
            in-out property <string> vec_c: "";
            in-out property <string> vec_g: "";
            in-out property <string> vec_t: "";

            in-out property <string> spark_a: "";
            in-out property <string> spark_c: "";
            in-out property <string> spark_g: "";
            in-out property <string> spark_t: "";

            in-out property <image> chart_img; // <-- chart

            callback open_file();
            callback compute();

            VerticalLayout {
                padding: 12px;
                spacing: 10px;

                // Controls row
                HorizontalLayout {
                    spacing: 8px;

                    Text { text: "Window:"; }
                    SpinBox { value: root.window_size; minimum: 1; maximum: 10000; }
                    Text { text: "Step:"; }
                    SpinBox { value: root.step_size; minimum: 1; maximum: 10000; }

                    Rectangle { width: 8px; height: 1px; }

                    Button { text: "Open FASTA…"; clicked => { root.open_file(); } }
                    Button { text: "Compute"; clicked => { root.compute(); } }

                    Text { text: root.file_path; horizontal-stretch: 1; }
                }

                Text { text: root.status; }

                Rectangle { background: #1a1a1a; border-width: 1px; border-color: #303030; width: 100%; height: 1px; }

                // Vectors (raw values, brief preview)
                VerticalLayout {
                    spacing: 6px;

                    HorizontalLayout { spacing: 6px; Text { text: "A:"; } Text { text: root.vec_a; wrap: word-wrap; } }
                    HorizontalLayout { spacing: 6px; Text { text: "C:"; } Text { text: root.vec_c; wrap: word-wrap; } }
                    HorizontalLayout { spacing: 6px; Text { text: "G:"; } Text { text: root.vec_g; wrap: word-wrap; } }
                    HorizontalLayout { spacing: 6px; Text { text: "T:"; } Text { text: root.vec_t; wrap: word-wrap; } }
                }

                Rectangle { background: #1a1a1a; border-width: 1px; border-color: #303030; width: 100%; height: 1px; }

                // Chart image
                Text { text: "Chart:"; }
                Image { source: root.chart_img; width: 100%; height: 240px; }

                Rectangle { background: #1a1a1a; border-width: 1px; border-color: #303030; width: 100%; height: 1px; }

                // Optional: keep the compact sparklines too
                VerticalLayout {
                    spacing: 4px;
                    Text { text: "Sparklines (0..1):"; }
                    Text { text: "A  " + root.spark_a; font-family: "monospace"; }
                    Text { text: "C  " + root.spark_c; font-family: "monospace"; }
                    Text { text: "G  " + root.spark_g; font-family: "monospace"; }
                    Text { text: "T  " + root.spark_t; font-family: "monospace"; }
                }

                Rectangle { horizontal-stretch: 1; vertical-stretch: 1; background: transparent; }
            }
        }
    }

    pub fn main() -> Result<(), Box<dyn std::error::Error>> {
        let app = AppWindow::new()?;

        // Store the chosen file path (if any) on the Rust side
        let chosen_path = std::rc::Rc::new(std::cell::RefCell::new(None::<String>));
        {
            let app_weak = app.as_weak();
            let chosen_path = chosen_path.clone();
            app.on_open_file(move || {
                if let Some(path) = FileDialog::new()
                    .add_filter("FASTA", &["fa", "fasta", "fna", "gz", "txt"])
                    .pick_file()
                {
                    let p = path.display().to_string();
                    if let Some(app) = app_weak.upgrade() {
                        app.set_file_path(SharedString::from(p.clone()));
                        app.set_status(SharedString::from("File selected. Set window/step and press Compute."));
                    }
                    *chosen_path.borrow_mut() = Some(p);
                }
            });
        }

        {
            let app_weak = app.as_weak();
            let chosen_path = chosen_path.clone();
            app.on_compute(move || {
                let (win, step) = {
                    if let Some(app) = app_weak.upgrade() {
                        (app.get_window_size() as usize, app.get_step_size() as usize)
                    } else {
                        return;
                    }
                };

                let path_opt = chosen_path.borrow().clone();
                let app_for_ui = app_weak.clone();

                std::thread::spawn(move || {
                    // Load sequence
                    let seq_res = if let Some(p) = path_opt {
                        parse_fasta_bytes(p)
                    } else {
                        Ok(b"ATTGTCCCAATCTGTTG".to_vec())
                    };

                    // Compute (choose full vs downsampled)
                    let outcome = match seq_res {
                        Err(e) => Err(e),
                        Ok(seq) => {
                            let total = num_windows(seq.len(), win, step);
                            if total == 0 {
                                Err("Sequence shorter than window".to_string())
                            } else if total > MAX_WINDOWS_FULL {
                                sliding_freqs_downsampled(&seq, win, step, SPARK_MAX_POINTS * 2)
                            } else {
                                sliding_freqs(&seq, win, step)
                            }
                        }
                    };

                    let _ = invoke_from_event_loop(move || {
                        if let Some(app) = app_for_ui.upgrade() {
                            match outcome {
                                Err(e) => app.set_status(SharedString::from(format!("Error: {e}"))),
                                Ok(series) => {
                                    // Previews to avoid giant strings
                                    let preview = |v: &[f32], max: usize| -> String {
                                        if v.is_empty() {
                                            return String::new();
                                        }
                                        let mut s = v
                                            .iter()
                                            .take(max)
                                            .map(|x| format!("{:.3}", x))
                                            .collect::<Vec<_>>()
                                            .join(", ");
                                        if v.len() > max {
                                            s.push_str(&format!(" … (+{} more)", v.len() - max));
                                        }
                                        s
                                    };

                                    // Set numeric previews
                                    app.set_vec_a(SharedString::from(preview(&series.a, 64)));
                                    app.set_vec_c(SharedString::from(preview(&series.c, 64)));
                                    app.set_vec_g(SharedString::from(preview(&series.g, 64)));
                                    app.set_vec_t(SharedString::from(preview(&series.t, 64)));

                                    // Set sparklines
                                    app.set_spark_a(SharedString::from(sparkline_downsample(&series.a, SPARK_MAX_POINTS)));
                                    app.set_spark_c(SharedString::from(sparkline_downsample(&series.c, SPARK_MAX_POINTS)));
                                    app.set_spark_g(SharedString::from(sparkline_downsample(&series.g, SPARK_MAX_POINTS)));
                                    app.set_spark_t(SharedString::from(sparkline_downsample(&series.t, SPARK_MAX_POINTS)));

                                    // Render and set chart image
                                    let img = render_chart_image(&series, CHART_W, CHART_H);
                                    app.set_chart_img(img);

                                    app.set_status(SharedString::from(format!(
                                        "Computed {} windows (win={}, step={})",
                                        series.a.len(),
                                        win,
                                        step
                                    )));
                                }
                            }
                        }
                    });
                });
            });
        }

        app.run()?;
        Ok(())
    }
}

#[cfg(feature = "gui")]
fn main() {
    if let Err(e) = gui_app::main() {
        eprintln!("Error: {e}");
    }
}
