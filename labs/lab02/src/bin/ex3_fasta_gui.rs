// labs/lab02/src/bin/ex3_fasta_gui.rs

#[cfg(not(feature = "gui"))]
fn main() {
    eprintln!("[lab02] ex3_fasta_gui: build with `--features gui`");
    eprintln!("Example: make ex3   # (uses --features gui --release)");
}

#[cfg(feature = "gui")]
#[path = "../ex3_core.rs"]
mod ex3_core;

#[cfg(feature = "gui")]
mod gui_app {
    use super::ex3_core::*;
    use rfd::FileDialog;
    use slint::{invoke_from_event_loop, ComponentHandle, Image, SharedString};
    use std::cell::RefCell;
    use std::rc::Rc;

    slint::slint! {
        import { Button, SpinBox } from "std-widgets.slint";

        export component AppWindow inherits Window {
            width: 1200px;
            height: 800px;
            title: "Lab02 - FASTA Sliding Window Frequencies";

            in-out property <int> window_size: 5;
            in-out property <int> step_size: 1;

            in-out property <string> file_path: "";
            in-out property <string> status: "Pick a FASTA and press Compute";

            // vector previews
            in-out property <string> vec_a: "";
            in-out property <string> vec_c: "";
            in-out property <string> vec_g: "";
            in-out property <string> vec_t: "";

            // chart image
            in-out property <image> chart_img;

            callback open_file();
            callback compute();

            VerticalLayout {
                padding: 12px;
                spacing: 10px;

                // Controls row
                HorizontalLayout {
                    spacing: 8px;

                    Text { text: "Window:"; }
                    // Two-way binding so changes propagate back to root
                    SpinBox { value <=> root.window_size; minimum: 1; maximum: 10000; }

                    Text { text: "Step:"; }
                    SpinBox { value <=> root.step_size; minimum: 1; maximum: 10000; }

                    Rectangle { width: 8px; height: 1px; }

                    Button { text: "Open FASTA…"; clicked => { root.open_file(); } }
                    Button { text: "Compute"; clicked => { root.compute(); } }

                    Text { text: root.file_path; horizontal-stretch: 1; }
                }

                Text { text: root.status; }

                Rectangle { background: #1a1a1a; border-width: 1px; border-color: #303030; width: 100%; height: 1px; }

                // Vectors (raw values preview)
                Text { text: "Vectors (first 64 values):"; }
                Rectangle {
                    background: #101010;
                    border-width: 1px; border-color: #303030;
                    width: 100%;

                    VerticalLayout {
                        padding: 8px;
                        spacing: 6px;

                        HorizontalLayout { spacing: 6px;
                            Text { text: "A:"; color: #E6E6E6; }
                            Text { text: root.vec_a; wrap: word-wrap; color: #E6E6E6; }
                        }
                        HorizontalLayout { spacing: 6px;
                            Text { text: "C:"; color: #E6E6E6; }
                            Text { text: root.vec_c; wrap: word-wrap; color: #E6E6E6; }
                        }
                        HorizontalLayout { spacing: 6px;
                            Text { text: "G:"; color: #E6E6E6; }
                            Text { text: root.vec_g; wrap: word-wrap; color: #E6E6E6; }
                        }
                        HorizontalLayout { spacing: 6px;
                            Text { text: "T:"; color: #E6E6E6; }
                            Text { text: root.vec_t; wrap: word-wrap; color: #E6E6E6; }
                        }
                    }
                }

                Rectangle { background: #1a1a1a; border-width: 1px; border-color: #303030; width: 100%; height: 1px; }

                // Chart header + legend
                HorizontalLayout {
                    spacing: 12px;

                    Text { text: "Chart: Relative nucleotide frequencies"; }

                    // legend entries
                    HorizontalLayout { spacing: 6px;
                        Rectangle { width: 12px; height: 12px; background: #8E44AD; }  // A (purple)

                        Text { text: "A"; }
                    }
                    HorizontalLayout { spacing: 6px;
                        Rectangle { width: 12px; height: 12px; background: #268BD2; }  // C
                        Text { text: "C"; }
                    }
                    HorizontalLayout { spacing: 6px;
                        Rectangle { width: 12px; height: 12px; background: #859900; }  // G
                        Text { text: "G"; }
                    }
                    HorizontalLayout { spacing: 6px;
                        Rectangle { width: 12px; height: 12px; background: #DC322F; }  // T
                        Text { text: "T"; }
                    }

                    Rectangle { width: 8px; height: 1px; }
                    Text { text: "y: 0 → 1"; }
                }


                HorizontalLayout {
    spacing: 8px;

    // Y-axis labels
    VerticalLayout {
        width: 50px;
        spacing: 0px;

        Text { text: "1.00"; horizontal-alignment: right; }
        Rectangle { vertical-stretch: 1; background: transparent; }

        Text { text: "0.75"; horizontal-alignment: right; }
        Rectangle { vertical-stretch: 1; background: transparent; }

        Text { text: "0.50"; horizontal-alignment: right; }
        Rectangle { vertical-stretch: 1; background: transparent; }

        Text { text: "0.25"; horizontal-alignment: right; }
        Rectangle { vertical-stretch: 1; background: transparent; }

        Text { text: "0.00"; horizontal-alignment: right; }
    }

    // Chart image
    Image { source: root.chart_img; width: 100%; height: 420px; }
}


                Rectangle { horizontal-stretch: 1; vertical-stretch: 1; background: transparent; }
            }
        }
    }

    pub fn main() -> Result<(), Box<dyn std::error::Error>> {
        let app = AppWindow::new()?;

        // keep chosen path on Rust side
        let chosen_path: Rc<RefCell<Option<String>>> = Rc::new(RefCell::new(None));

        // --- Open file ---
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

        // --- Compute ---
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
                                sliding_freqs_downsampled(&seq, win, step, 2 * CHART_W as usize)
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
                                    // Render and set chart image
                                    let img: Image = render_chart_image(&series, CHART_W, CHART_H);
                                    app.set_chart_img(img);

                                    // Short previews to avoid huge strings
                                    let preview = |v: &[f32], max: usize| -> String {
                                        if v.is_empty() { return String::new(); }
                                        let mut s = v.iter().take(max).map(|x| format!("{:.3}", x))
                                            .collect::<Vec<_>>().join(", ");
                                        if v.len() > max { s.push_str(&format!(" … (+{} more)", v.len() - max)); }
                                        s
                                    };

                                    app.set_vec_a(SharedString::from(preview(&series.a, 64)));
                                    app.set_vec_c(SharedString::from(preview(&series.c, 64)));
                                    app.set_vec_g(SharedString::from(preview(&series.g, 64)));
                                    app.set_vec_t(SharedString::from(preview(&series.t, 64)));

                                    app.set_status(SharedString::from(format!(
                                        "Computed {} windows (win={}, step={})",
                                        series.a.len(), win, step
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
