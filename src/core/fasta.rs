use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
#[cfg(unix)]
use std::os::unix::fs::FileExt;
#[cfg(windows)]
use std::os::windows::fs::FileExt;
use std::path::Path;

struct FaiEntry {
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

enum Backend {
    Pread(File),
    Mmap(memmap2::Mmap),
}

/// Positioned read that does not disturb the file cursor.
///
/// `read_at` is `pread(2)`, which only exists on Unix; Windows spells the same
/// operation `seek_read`. Both leave the offset unchanged, so the two backends
/// behave identically.
#[cfg(unix)]
fn pread(file: &File, buf: &mut [u8], offset: u64) -> std::io::Result<usize> {
    file.read_at(buf, offset)
}

#[cfg(windows)]
fn pread(file: &File, buf: &mut [u8], offset: u64) -> std::io::Result<usize> {
    file.seek_read(buf, offset)
}

pub struct FastaReader {
    backend: Backend,
    index: HashMap<String, FaiEntry>,
    seq_names: HashSet<String>,
    /// Sequence names in `.fai` (i.e. FASTA) order, for deterministic output.
    order: Vec<String>,
}

unsafe impl Send for FastaReader {}
unsafe impl Sync for FastaReader {}

impl FastaReader {
    /// Open with pread backend (low RSS, good for random whole-genome access).
    pub fn open<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let path = path.as_ref();
        Self::ensure_fai(path)?;
        let (index, order) = Self::parse_fai(&path.with_extension(
            format!("{}.fai", path.extension().unwrap_or_default().to_string_lossy()),
        ))?;
        let seq_names: HashSet<String> = index.keys().cloned().collect();
        let file = File::open(path)?;
        Ok(Self { backend: Backend::Pread(file), index, seq_names, order })
    }

    /// Open with mmap backend (fast sequential access, RSS ~ touched pages).
    pub fn open_mmap<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let path = path.as_ref();
        Self::ensure_fai(path)?;
        let (index, order) = Self::parse_fai(&path.with_extension(
            format!("{}.fai", path.extension().unwrap_or_default().to_string_lossy()),
        ))?;
        let seq_names: HashSet<String> = index.keys().cloned().collect();
        let file = File::open(path)?;
        let mmap = unsafe { memmap2::MmapOptions::new().map(&file)? };
        Ok(Self { backend: Backend::Mmap(mmap), index, seq_names, order })
    }

    fn ensure_fai(path: &Path) -> std::io::Result<()> {
        let fai_path = path.with_extension(
            format!("{}.fai", path.extension().unwrap_or_default().to_string_lossy()),
        );
        if !fai_path.exists() {
            Self::build_fai(path, &fai_path)?;
        }
        Ok(())
    }

    fn build_fai(fasta_path: &Path, fai_path: &Path) -> std::io::Result<()> {
        let file = File::open(fasta_path)?;
        let mut reader = BufReader::new(file);
        let mut out = std::io::BufWriter::new(File::create(fai_path)?);

        let mut name = String::new();
        let mut seq_len: u64 = 0;
        let mut seq_offset: u64 = 0;
        let mut line_bases: u64 = 0;
        let mut line_width: u64 = 0;
        let mut byte_offset: u64 = 0;
        let mut in_seq = false;

        let mut buf = String::new();
        loop {
            buf.clear();
            let n = reader.read_line(&mut buf)?;
            if n == 0 {
                if in_seq && !name.is_empty() {
                    writeln!(out, "{}\t{}\t{}\t{}\t{}", name, seq_len, seq_offset, line_bases, line_width)?;
                }
                break;
            }

            if buf.starts_with('>') {
                if in_seq && !name.is_empty() {
                    writeln!(out, "{}\t{}\t{}\t{}\t{}", name, seq_len, seq_offset, line_bases, line_width)?;
                }
                name = buf[1..].trim_end().split_whitespace().next().unwrap_or("").to_string();
                seq_len = 0;
                seq_offset = byte_offset + n as u64;
                line_bases = 0;
                line_width = 0;
                in_seq = true;
            } else if in_seq {
                let bases = buf.trim_end_matches(&['\n', '\r'][..]).len() as u64;
                if line_bases == 0 {
                    line_bases = bases;
                    line_width = n as u64;
                }
                seq_len += bases;
            }

            byte_offset += n as u64;
        }
        Ok(())
    }

    /// Parse a `.fai` index.
    ///
    /// Returns the index map plus the sequence names in file order, so that
    /// callers can emit deterministic (FASTA-ordered) output.
    fn parse_fai(fai_path: &Path) -> std::io::Result<(HashMap<String, FaiEntry>, Vec<String>)> {
        let file = File::open(fai_path)?;
        let reader = BufReader::new(file);
        let mut index = HashMap::new();
        let mut order = Vec::new();
        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 5 {
                continue;
            }
            let name = fields[0].to_string();
            let length: u64 = fields[1].parse().map_err(|_|
                std::io::Error::new(std::io::ErrorKind::InvalidData, "bad fai length"))?;
            let offset: u64 = fields[2].parse().map_err(|_|
                std::io::Error::new(std::io::ErrorKind::InvalidData, "bad fai offset"))?;
            let line_bases: u64 = fields[3].parse().map_err(|_|
                std::io::Error::new(std::io::ErrorKind::InvalidData, "bad fai line_bases"))?;
            let line_width: u64 = fields[4].parse().map_err(|_|
                std::io::Error::new(std::io::ErrorKind::InvalidData, "bad fai line_width"))?;
            order.push(name.clone());
            index.insert(name, FaiEntry { length, offset, line_bases, line_width });
        }
        Ok((index, order))
    }

    pub fn resolve_name(&self, chrom: &str) -> Option<String> {
        if self.seq_names.contains(chrom) {
            return Some(chrom.to_string());
        }
        if chrom.starts_with("chr") {
            let alt = &chrom[3..];
            if self.seq_names.contains(alt) {
                return Some(alt.to_string());
            }
        } else {
            let alt = format!("chr{}", chrom);
            if self.seq_names.contains(&alt) {
                return Some(alt);
            }
        }
        None
    }

    #[inline]
    fn byte_offset(entry: &FaiEntry, pos: u64) -> u64 {
        entry.offset + pos / entry.line_bases * entry.line_width + pos % entry.line_bases
    }

    fn extract_seq(buf: &[u8], seq_len: usize) -> Option<String> {
        let mut result = Vec::with_capacity(seq_len);
        for &b in buf {
            if b != b'\n' && b != b'\r' {
                result.push(b.to_ascii_uppercase());
            }
        }
        if result.len() == seq_len {
            Some(unsafe { String::from_utf8_unchecked(result) })
        } else {
            None
        }
    }

    fn read_range(&self, entry: &FaiEntry, start: u64, end: u64) -> Option<String> {
        if start >= end || end > entry.length {
            return None;
        }
        let seq_len = (end - start) as usize;
        let first_byte = Self::byte_offset(entry, start);
        let last_byte = Self::byte_offset(entry, end - 1);
        let raw_len = (last_byte - first_byte + 1) as usize;

        match &self.backend {
            Backend::Pread(file) => {
                let mut buf = vec![0u8; raw_len];
                if pread(file, &mut buf, first_byte).ok()? < raw_len {
                    return None;
                }
                Self::extract_seq(&buf, seq_len)
            }
            Backend::Mmap(mmap) => {
                let start_idx = first_byte as usize;
                let end_idx = start_idx + raw_len;
                if end_idx > mmap.len() {
                    return None;
                }
                Self::extract_seq(&mmap[start_idx..end_idx], seq_len)
            }
        }
    }

    pub fn fetch(&self, chrom: &str, start: u64, end: u64) -> Option<String> {
        let name = self.resolve_name(chrom)?;
        let entry = self.index.get(&name)?;
        self.read_range(entry, start, end)
    }

    /// Batch fetch sorted by file offset for sequential I/O (pread backend).
    /// Results are returned in the original input order.
    pub fn batch_fetch(&self, requests: &[(String, u64, u64)]) -> Vec<Option<String>> {
        let mut results: Vec<Option<String>> = vec![None; requests.len()];

        struct FetchJob {
            idx: usize,
            resolved_name: String,
            start: u64,
            end: u64,
            file_offset: u64,
        }

        let mut jobs: Vec<FetchJob> = Vec::with_capacity(requests.len());

        for (i, (chrom, start, end)) in requests.iter().enumerate() {
            if let Some(name) = self.resolve_name(chrom) {
                if let Some(entry) = self.index.get(&name) {
                    if *end <= entry.length && *start < *end {
                        let fo = Self::byte_offset(entry, *start);
                        jobs.push(FetchJob {
                            idx: i,
                            resolved_name: name,
                            start: *start,
                            end: *end,
                            file_offset: fo,
                        });
                    }
                }
            }
        }

        jobs.sort_unstable_by_key(|j| j.file_offset);

        for job in &jobs {
            if let Some(entry) = self.index.get(&job.resolved_name) {
                results[job.idx] = self.read_range(entry, job.start, job.end);
            }
        }

        results
    }

    /// Sequence names in FASTA (`.fai`) order.
    pub fn references(&self) -> Vec<String> {
        self.order.clone()
    }

    pub fn lengths(&self) -> Vec<usize> {
        self.references()
            .iter()
            .map(|name| {
                self.index.get(name).map(|e| e.length as usize).unwrap_or(0)
            })
            .collect()
    }
}
