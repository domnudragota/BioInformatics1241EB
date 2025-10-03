use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};

use flate2::read::GzDecoder;

#[derive(Debug, Clone)]
pub struct RecordMeta {
    pub header: String,
    pub len: u64,
}

#[derive(Debug)]
pub struct Stats {
    pub counts: [u64; 256],
    pub total_len: u64,
    pub record_count: u64,
    pub header_preview: String,
    pub file_size: u64,          // on-disk size (compressed if .gz)
    pub records: Vec<RecordMeta> // per-record info
}

fn is_gzip_file<P: AsRef<Path>>(path: P, file: &mut File) -> bool {
    let ext_is_gz = path
        .as_ref()
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);

    let mut magic = [0u8; 2];
    let mut is_gz_magic = false;
    if file.read_exact(&mut magic).is_ok() {
        is_gz_magic = magic == [0x1f, 0x8b];
    }
    let _ = file.seek(SeekFrom::Start(0));
    ext_is_gz || is_gz_magic
}

pub fn parse_fasta<P: AsRef<Path>>(
    path: P,
    cancel_flag: &AtomicBool,
    mut on_progress: impl FnMut(u64, u64),
) -> Result<Stats, String> {
    let mut file = File::open(&path).map_err(|e| format!("Open error: {e}"))?;
    let on_disk_size = file
        .metadata()
        .map_err(|e| format!("Metadata error: {e}"))?
        .len();

    let gz = is_gzip_file(&path, &mut file);

    enum Source { Plain(File), Gzip(GzDecoder<File>) }
    let source = if gz { Source::Gzip(GzDecoder::new(file)) } else { Source::Plain(file) };

    let mut reader: Box<dyn BufRead> = match source {
        Source::Plain(f) => Box::new(BufReader::with_capacity(4 * 1024 * 1024, f)),
        Source::Gzip(gz) => Box::new(BufReader::with_capacity(4 * 1024 * 1024, gz)),
    };

    let mut counts = [0u64; 256];
    let mut total_len: u64 = 0;
    let mut records: Vec<RecordMeta> = Vec::new();

    let mut header_preview = String::new();
    let mut current_header: Option<String> = None;
    let mut current_len: u64 = 0;

    let mut bytes_read_decompressed: u64 = 0;
    let mut line = String::new();

    loop {
        if cancel_flag.load(Ordering::Relaxed) {
            return Err("Canceled".to_string());
        }

        line.clear();
        let n = reader.read_line(&mut line).map_err(|e| format!("Read error: {e}"))?;
        if n == 0 { break; }
        bytes_read_decompressed = bytes_read_decompressed.saturating_add(n as u64);

        if line.starts_with('>') {
            // push previous record if we had one
            if let Some(hdr) = current_header.take() {
                records.push(RecordMeta { header: hdr, len: current_len });
                current_len = 0;
            }

            // start new record
            let mut hdr = line[1..].trim().to_string();
            if hdr.len() > 120 { hdr.truncate(120); hdr.push('…'); }
            if header_preview.is_empty() {
                header_preview = hdr.clone();
            }
            current_header = Some(hdr);
        } else {
            // sequence content line
            for &b in line.as_bytes() {
                match b {
                    b'\n' | b'\r' | b' ' | b'\t' => continue,
                    _ => {
                        let bb = if (b'a'..=b'z').contains(&b) { b - 32 } else { b }; // ASCII uppercase
                        counts[bb as usize] = counts[bb as usize].saturating_add(1);
                        total_len = total_len.saturating_add(1);
                        current_len = current_len.saturating_add(1);
                    }
                }
            }
        }

        // progress (unknown total for .gz)
        if gz { on_progress(bytes_read_decompressed, 0); }
        else { on_progress(bytes_read_decompressed, on_disk_size); }
    }

    // push last record (if any)
    match (current_header, current_len) {
        (Some(hdr), len) => records.push(RecordMeta { header: hdr, len }),
        (None, len) if len > 0 => {
            // no header but data → synthesize one
            let hdr = "(no header)".to_string();
            if header_preview.is_empty() { header_preview = hdr.clone(); }
            records.push(RecordMeta { header: hdr, len });
        }
        _ => {}
    }

    let record_count = records.len() as u64;

    Ok(Stats {
        counts,
        total_len,
        record_count,
        header_preview,
        file_size: on_disk_size,
        records,
    })
}
