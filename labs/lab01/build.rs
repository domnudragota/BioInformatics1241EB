

#[cfg(feature = "gui")]
fn main() {
    slint_build::compile("ui/app.slint").unwrap();
}

#[cfg(not(feature = "gui"))]
fn main() {}
