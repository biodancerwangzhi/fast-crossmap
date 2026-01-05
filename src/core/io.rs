//! High-performance I/O abstraction layer
//!
//! Provides optimized file reading with configurable buffer sizes
//! and optional memory mapping for large files.

use memmap2::Mmap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

/// Default buffer size for BufReader (128KB)
pub const DEFAULT_BUFFER_SIZE: usize = 128 * 1024;

/// Large buffer size for high-throughput I/O (1MB)
pub const LARGE_BUFFER_SIZE: usize = 1024 * 1024;

/// Threshold for using memory mapping (100MB)
pub const MMAP_THRESHOLD: u64 = 100 * 1024 * 1024;

/// I/O strategy selection
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IoStrategy {
    /// Use buffered reading with configurable buffer size
    Buffered(usize),
    /// Use memory mapping for the entire file
    MemoryMapped,
    /// Automatically select based on file size
    Auto,
}

impl Default for IoStrategy {
    fn default() -> Self {
        IoStrategy::Auto
    }
}

/// A smart reader that automatically selects the optimal I/O strategy
pub enum SmartReader {
    /// Buffered reader for smaller files or streaming
    Buffered(BufReader<File>),
    /// Memory-mapped reader for large files
    Mapped(MappedReader),
}

/// Memory-mapped file reader
pub struct MappedReader {
    #[allow(dead_code)]
    mmap: Mmap,
    position: usize,
}

impl MappedReader {
    /// Create a new memory-mapped reader
    pub fn new(file: &File) -> io::Result<Self> {
        // SAFETY: We assume the file won't be modified while mapped
        let mmap = unsafe { Mmap::map(file)? };
        Ok(Self { mmap, position: 0 })
    }

    /// Get the entire file content as a byte slice
    pub fn as_bytes(&self) -> &[u8] {
        &self.mmap
    }

    /// Get remaining bytes from current position
    pub fn remaining(&self) -> &[u8] {
        &self.mmap[self.position..]
    }

    /// Get file size
    pub fn len(&self) -> usize {
        self.mmap.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.mmap.is_empty()
    }
}

impl Read for MappedReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let remaining = &self.mmap[self.position..];
        let to_read = std::cmp::min(buf.len(), remaining.len());
        buf[..to_read].copy_from_slice(&remaining[..to_read]);
        self.position += to_read;
        Ok(to_read)
    }
}

impl BufRead for MappedReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        Ok(&self.mmap[self.position..])
    }

    fn consume(&mut self, amt: usize) {
        self.position = std::cmp::min(self.position + amt, self.mmap.len());
    }
}

impl SmartReader {
    /// Open a file with the specified I/O strategy
    pub fn open<P: AsRef<Path>>(path: P, strategy: IoStrategy) -> io::Result<Self> {
        let file = File::open(path.as_ref())?;
        let metadata = file.metadata()?;
        let file_size = metadata.len();

        match strategy {
            IoStrategy::Buffered(buf_size) => {
                Ok(SmartReader::Buffered(BufReader::with_capacity(buf_size, file)))
            }
            IoStrategy::MemoryMapped => {
                Ok(SmartReader::Mapped(MappedReader::new(&file)?))
            }
            IoStrategy::Auto => {
                if file_size >= MMAP_THRESHOLD {
                    // Use memory mapping for large files
                    Ok(SmartReader::Mapped(MappedReader::new(&file)?))
                } else {
                    // Use buffered reading for smaller files
                    let buf_size = if file_size > 10 * 1024 * 1024 {
                        LARGE_BUFFER_SIZE
                    } else {
                        DEFAULT_BUFFER_SIZE
                    };
                    Ok(SmartReader::Buffered(BufReader::with_capacity(buf_size, file)))
                }
            }
        }
    }

    /// Open with default auto strategy
    pub fn open_auto<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::open(path, IoStrategy::Auto)
    }

    /// Open with specified buffer size
    pub fn open_buffered<P: AsRef<Path>>(path: P, buffer_size: usize) -> io::Result<Self> {
        Self::open(path, IoStrategy::Buffered(buffer_size))
    }

    /// Check if using memory mapping
    pub fn is_mapped(&self) -> bool {
        matches!(self, SmartReader::Mapped(_))
    }
}

impl Read for SmartReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            SmartReader::Buffered(reader) => reader.read(buf),
            SmartReader::Mapped(reader) => reader.read(buf),
        }
    }
}

impl BufRead for SmartReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            SmartReader::Buffered(reader) => reader.fill_buf(),
            SmartReader::Mapped(reader) => reader.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            SmartReader::Buffered(reader) => reader.consume(amt),
            SmartReader::Mapped(reader) => reader.consume(amt),
        }
    }
}

/// Create a buffered reader with optimal buffer size
pub fn create_buf_reader<P: AsRef<Path>>(path: P) -> io::Result<BufReader<File>> {
    let file = File::open(path)?;
    Ok(BufReader::with_capacity(DEFAULT_BUFFER_SIZE, file))
}

/// Create a buffered reader with custom buffer size
pub fn create_buf_reader_with_capacity<P: AsRef<Path>>(
    path: P,
    capacity: usize,
) -> io::Result<BufReader<File>> {
    let file = File::open(path)?;
    Ok(BufReader::with_capacity(capacity, file))
}

/// Line iterator that reuses a buffer to avoid allocations
pub struct LineIterator<R: BufRead> {
    reader: R,
    buffer: String,
}

impl<R: BufRead> LineIterator<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            buffer: String::with_capacity(1024),
        }
    }

    /// Read the next line into the internal buffer
    /// Returns None at EOF, Some(Ok(&str)) on success, Some(Err) on error
    pub fn next_line(&mut self) -> Option<io::Result<&str>> {
        self.buffer.clear();
        match self.reader.read_line(&mut self.buffer) {
            Ok(0) => None, // EOF
            Ok(_) => {
                // Remove trailing newline
                if self.buffer.ends_with('\n') {
                    self.buffer.pop();
                    if self.buffer.ends_with('\r') {
                        self.buffer.pop();
                    }
                }
                Some(Ok(&self.buffer))
            }
            Err(e) => Some(Err(e)),
        }
    }
}

/// Byte line iterator for zero-copy parsing
pub struct ByteLineIterator<R: BufRead> {
    reader: R,
    buffer: Vec<u8>,
}

impl<R: BufRead> ByteLineIterator<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            buffer: Vec::with_capacity(4096),
        }
    }

    /// Read the next line as bytes
    pub fn next_line(&mut self) -> Option<io::Result<&[u8]>> {
        self.buffer.clear();
        match self.reader.read_until(b'\n', &mut self.buffer) {
            Ok(0) => None, // EOF
            Ok(_) => {
                // Remove trailing newline
                if self.buffer.last() == Some(&b'\n') {
                    self.buffer.pop();
                    if self.buffer.last() == Some(&b'\r') {
                        self.buffer.pop();
                    }
                }
                Some(Ok(&self.buffer))
            }
            Err(e) => Some(Err(e)),
        }
    }

    /// Get the internal buffer for reuse
    pub fn buffer(&self) -> &[u8] {
        &self.buffer
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_default_buffer_size() {
        assert_eq!(DEFAULT_BUFFER_SIZE, 128 * 1024);
    }

    #[test]
    fn test_large_buffer_size() {
        assert_eq!(LARGE_BUFFER_SIZE, 1024 * 1024);
    }

    #[test]
    fn test_io_strategy_default() {
        assert_eq!(IoStrategy::default(), IoStrategy::Auto);
    }

    #[test]
    fn test_smart_reader_buffered() -> io::Result<()> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, "line1\nline2\nline3")?;
        
        let reader = SmartReader::open(temp.path(), IoStrategy::Buffered(1024))?;
        assert!(!reader.is_mapped());
        Ok(())
    }

    #[test]
    fn test_smart_reader_auto_small_file() -> io::Result<()> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, "small file content")?;
        
        let reader = SmartReader::open_auto(temp.path())?;
        // Small file should use buffered reading
        assert!(!reader.is_mapped());
        Ok(())
    }

    #[test]
    fn test_line_iterator() -> io::Result<()> {
        let mut temp = NamedTempFile::new()?;
        writeln!(temp, "line1")?;
        writeln!(temp, "line2")?;
        writeln!(temp, "line3")?;
        temp.flush()?;

        let file = File::open(temp.path())?;
        let reader = BufReader::new(file);
        let mut iter = LineIterator::new(reader);

        assert_eq!(iter.next_line().unwrap()?, "line1");
        assert_eq!(iter.next_line().unwrap()?, "line2");
        assert_eq!(iter.next_line().unwrap()?, "line3");
        assert!(iter.next_line().is_none());
        Ok(())
    }

    #[test]
    fn test_byte_line_iterator() -> io::Result<()> {
        let mut temp = NamedTempFile::new()?;
        temp.write_all(b"line1\nline2\nline3\n")?;
        temp.flush()?;

        let file = File::open(temp.path())?;
        let reader = BufReader::new(file);
        let mut iter = ByteLineIterator::new(reader);

        assert_eq!(iter.next_line().unwrap()?, b"line1");
        assert_eq!(iter.next_line().unwrap()?, b"line2");
        assert_eq!(iter.next_line().unwrap()?, b"line3");
        assert!(iter.next_line().is_none());
        Ok(())
    }

    #[test]
    fn test_mapped_reader_len() -> io::Result<()> {
        let mut temp = NamedTempFile::new()?;
        temp.write_all(b"test content")?;
        temp.flush()?;

        let file = File::open(temp.path())?;
        let reader = MappedReader::new(&file)?;
        
        assert_eq!(reader.len(), 12);
        assert!(!reader.is_empty());
        assert_eq!(reader.as_bytes(), b"test content");
        Ok(())
    }
}
