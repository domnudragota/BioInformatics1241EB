mod parser;

use parser::{parse_fasta, Stats, RecordMeta};
use slint::{ModelRc, SharedString, VecModel};
use std::path::PathBuf;
use std::sync::{
    atomic::{AtomicBool, Ordering},
    Arc,
};
use std::thread;

slint::include_modules!();

fn human_mb(bytes: u64) -> String {
    format!("{:.1} MB", (bytes as f64) / 1_000_000.0)
}

fn build_frequency_rows(stats: &Stats) -> Vec<Row> {
    let total = stats.total_len as f64;
    if total == 0.0 {
        return Vec::new();
    }

    let mut rows: Vec<(char, u64)> = Vec::new();
    for c in b'A'..=b'Z' {
        let cnt = stats.counts[c as usize];
        if cnt > 0 { rows.push((c as char, cnt)); }
    }
    for i in 0u16..=255 {
        let b = i as u8;
        if (b'A'..=b'Z').contains(&b) { continue; }
        let cnt = stats.counts[b as usize];
        if cnt > 0 {
            let ch = b as char;
            if ch.is_ascii_graphic() { rows.push((ch, cnt)); }
        }
    }

    rows.into_iter()
        .map(|(ch, cnt)| {
            let pct = 100.0 * (cnt as f64) / total;
            Row {
                symbol: SharedString::from(ch.to_string()),
                count: cnt as i32,
                pct: SharedString::from(format!("{:.2}%", pct)),
            }
        })
        .collect()
}

fn build_record_rows(stats: &Stats) -> Vec<Rec> {
    stats
        .records
        .iter()
        .enumerate()
        .map(|(i, r)| Rec {
            idx: (i + 1) as i32,
            header: SharedString::from(r.header.clone()),
            len: r.len as i32,
        })
        .collect()
}

fn main() -> Result<(), slint::PlatformError> {
    let app = AppWindow::new()?;

    let cancel_flag = Arc::new(AtomicBool::new(false));

    // Cancel
    {
        let cancel_flag = cancel_flag.clone();
        let weak = app.as_weak();
        app.on_cancel_clicked(move || {
            cancel_flag.store(true, Ordering::SeqCst);
            if let Some(win) = weak.upgrade() {
                win.set_status(SharedString::from("Cancel requested…"));
            }
        });
    }

    // Open
    {
        let weak = app.as_weak();
        let cancel_flag = cancel_flag.clone();

        app.on_open_file_clicked(move || {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("FASTA (plain)", &["fa", "fasta", "fna", "faa", "fas"])
                .add_filter("GZip", &["gz"])
                .pick_file()
            {
                let pbuf: PathBuf = path.into();
                let display = pbuf.display().to_string();

                if let Some(win) = weak.upgrade() {
                    win.set_file_path(SharedString::from(display.clone()));
                    win.set_status(SharedString::from("Opening file…"));
                    win.set_progress(0.0);
                    win.set_parsing(true);
                    win.set_header_preview(SharedString::from(""));
                    win.set_record_count(0);
                    win.set_total_len(0);
                    win.set_frequencies(ModelRc::new(VecModel::from(Vec::<Row>::new())));
                    win.set_records(ModelRc::new(VecModel::from(Vec::<Rec>::new())));
                }

                cancel_flag.store(false, Ordering::SeqCst);

                let weak_ui = weak.clone();
                let cancel_in_thread = cancel_flag.clone();
                thread::spawn(move || {
                    let progress_cb = |done: u64, total: u64| {
                        let weak_ui = weak_ui.clone();
                        let mb_done = (done as f64) / 1_000_000.0;
                        let status_txt = if total > 0 {
                            let mb_total = (total as f64) / 1_000_000.0;
                            format!("Parsing… {:.1}/{:.1} MB", mb_done, mb_total)
                        } else {
                            // gz (unknown decompressed total)
                            format!("Parsing… {:.1} MB (compressed source)", mb_done)
                        };

                        let frac = if total > 0 { (done as f32) / (total as f32) } else { 0.0 };

                        let _ = slint::invoke_from_event_loop(move || {
                            if let Some(win) = weak_ui.upgrade() {
                                win.set_progress(frac.min(1.0));
                                win.set_status(SharedString::from(status_txt));
                            }
                        });
                    };

                    let result = parse_fasta(&pbuf, &cancel_in_thread, progress_cb);

                    let _ = slint::invoke_from_event_loop(move || {
                        if let Some(win) = weak_ui.upgrade() {
                            match result {
                                Ok(stats) => {
                                    let freq_rows = build_frequency_rows(&stats);
                                    let rec_rows  = build_record_rows(&stats);

                                    win.set_header_preview(SharedString::from(stats.header_preview));
                                    win.set_record_count(stats.record_count as i32);
                                    win.set_total_len(stats.total_len as i32);
                                    win.set_status(SharedString::from(format!(
                                        "Done: {} residues across {} record(s) ({})",
                                        stats.total_len,
                                        stats.record_count,
                                        human_mb(stats.file_size)
                                    )));
                                    win.set_progress(1.0);
                                    win.set_frequencies(ModelRc::new(VecModel::from(freq_rows)));
                                    win.set_records(ModelRc::new(VecModel::from(rec_rows)));
                                }
                                Err(msg) => {
                                    let cancelled = msg == "Canceled";
                                    win.set_status(SharedString::from(if cancelled { "Canceled." } else { &msg }));
                                    if !cancelled { win.set_progress(0.0); }
                                }
                            }
                            win.set_parsing(false);
                        }
                    });
                });
            }
        });
    }

    app.run()
}
