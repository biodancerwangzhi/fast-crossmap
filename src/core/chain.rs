//! Chain file parsing
//!
//! Parses UCSC chain format files used for coordinate liftover.
//!
//! # Chain File Format
//!
//! ```text
//! chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
//! size dt dq
//! size dt dq
//! size
//! ```
//!
//! - Header line starts with "chain"
//! - Data lines contain: size (alignment block), dt (target gap), dq (query/source gap)
//! - Last data line has only size (no gaps)

use crate::core::Strand;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Error type for chain file parsing
/// 
/// Provides detailed error information including line numbers and
/// descriptive messages for debugging chain file issues.
#[derive(Debug, Clone)]
pub struct ChainParseError {
    /// Human-readable error message
    pub message: String,
    /// Line number where the error occurred (1-based)
    pub line_number: Option<usize>,
    /// The kind of error that occurred
    pub kind: ChainParseErrorKind,
    /// The problematic content (if available)
    pub content: Option<String>,
}

/// Specific kinds of chain parsing errors
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChainParseErrorKind {
    /// Invalid chain header format
    InvalidHeader,
    /// Invalid data line format
    InvalidDataLine,
    /// Invalid strand character (must be '+' or '-')
    InvalidStrand,
    /// Failed to parse a numeric value
    InvalidNumber,
    /// Unexpected end of file
    UnexpectedEof,
    /// I/O error during reading
    IoError,
    /// File not found
    FileNotFound,
    /// Unsupported compression format
    UnsupportedCompression,
    /// Coordinate validation error (e.g., start > end)
    InvalidCoordinates,
}

impl std::fmt::Display for ChainParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.line_number {
            Some(line) => write!(f, "Line {}: {}", line, self.message),
            None => write!(f, "{}", self.message),
        }
    }
}

impl std::error::Error for ChainParseError {}

impl ChainParseError {
    /// Create a new error without line number
    pub fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            line_number: None,
            kind: ChainParseErrorKind::InvalidHeader,
            content: None,
        }
    }

    /// Create a new error with line number
    pub fn with_line(message: impl Into<String>, line_number: usize) -> Self {
        Self {
            message: message.into(),
            line_number: Some(line_number),
            kind: ChainParseErrorKind::InvalidHeader,
            content: None,
        }
    }
    
    /// Create an error with full context
    pub fn with_context(
        message: impl Into<String>,
        line_number: usize,
        kind: ChainParseErrorKind,
        content: Option<String>,
    ) -> Self {
        Self {
            message: message.into(),
            line_number: Some(line_number),
            kind,
            content,
        }
    }
    
    /// Create an invalid header error
    pub fn invalid_header(message: impl Into<String>, line_number: usize, content: &str) -> Self {
        Self::with_context(
            message,
            line_number,
            ChainParseErrorKind::InvalidHeader,
            Some(content.chars().take(100).collect()),
        )
    }
    
    /// Create an invalid data line error
    pub fn invalid_data_line(message: impl Into<String>, line_number: usize, content: &str) -> Self {
        Self::with_context(
            message,
            line_number,
            ChainParseErrorKind::InvalidDataLine,
            Some(content.chars().take(100).collect()),
        )
    }
    
    /// Create an invalid strand error
    pub fn invalid_strand(strand: char, line_number: usize) -> Self {
        Self::with_context(
            format!("Invalid strand character '{}', expected '+' or '-'", strand),
            line_number,
            ChainParseErrorKind::InvalidStrand,
            None,
        )
    }
    
    /// Create an invalid number error
    pub fn invalid_number(field: &str, value: &str, line_number: usize) -> Self {
        Self::with_context(
            format!("Invalid {} value '{}': expected a non-negative integer", field, value),
            line_number,
            ChainParseErrorKind::InvalidNumber,
            None,
        )
    }
    
    /// Create a file not found error
    pub fn file_not_found(path: &Path) -> Self {
        Self {
            message: format!("Chain file not found: {}", path.display()),
            line_number: None,
            kind: ChainParseErrorKind::FileNotFound,
            content: None,
        }
    }
    
    /// Create an invalid coordinates error
    pub fn invalid_coordinates(message: impl Into<String>, line_number: usize) -> Self {
        Self::with_context(
            message,
            line_number,
            ChainParseErrorKind::InvalidCoordinates,
            None,
        )
    }
    
    /// Check if this is a specific kind of error
    pub fn is_kind(&self, kind: ChainParseErrorKind) -> bool {
        self.kind == kind
    }
}

impl From<std::io::Error> for ChainParseError {
    fn from(e: std::io::Error) -> Self {
        Self {
            message: format!("IO error: {}", e),
            line_number: None,
            kind: ChainParseErrorKind::IoError,
            content: None,
        }
    }
}

/// Parsed chain header information
/// 
/// Note on naming convention:
/// In UCSC chain format, "target" (t) is the reference genome and "query" (q) is the aligned genome.
/// For liftover purposes (e.g., GRCh37_to_GRCh38.chain):
/// - "target" in chain = source coordinates (GRCh37) = what user queries
/// - "query" in chain = target coordinates (GRCh38) = mapping result
/// 
/// We use "source" to mean "user input coordinates" and "target" to mean "mapping output coordinates".
/// So we swap the UCSC naming: our source = UCSC target, our target = UCSC query.
#[derive(Debug, Clone)]
pub struct ChainHeader {
    pub score: u64,
    /// Source chromosome (UCSC "target") - the coordinate system user queries
    pub source_name: String,
    pub source_size: u64,
    pub source_strand: Strand,
    pub source_start: u64,
    pub source_end: u64,
    /// Target chromosome (UCSC "query") - the coordinate system we map to
    pub target_name: String,
    pub target_size: u64,
    pub target_strand: Strand,
    pub target_start: u64,
    pub target_end: u64,
    pub chain_id: String,
}


impl ChainHeader {
    /// Parse a chain header line
    ///
    /// Format: chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
    /// 
    /// Note: UCSC "target" (t) = our "source" (user input coordinates)
    ///       UCSC "query" (q) = our "target" (mapping output coordinates)
    pub fn parse(line: &str, line_number: usize) -> Result<Self, ChainParseError> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        
        if fields.len() < 12 {
            return Err(ChainParseError::invalid_header(
                format!("Expected 12+ fields, got {}", fields.len()),
                line_number,
                line,
            ));
        }
        
        if fields[0] != "chain" {
            return Err(ChainParseError::invalid_header(
                format!("Expected 'chain' keyword, got '{}'", fields[0]),
                line_number,
                line,
            ));
        }
        
        let score = fields[1].parse::<u64>().map_err(|_| {
            ChainParseError::invalid_number("score", fields[1], line_number)
        })?;
        
        // UCSC "target" = our "source" (user input coordinates)
        let source_name = fields[2].to_string();
        let source_size = fields[3].parse::<u64>().map_err(|_| {
            ChainParseError::invalid_number("source size", fields[3], line_number)
        })?;
        
        let source_strand_char = fields[4].chars().next().unwrap_or('?');
        let source_strand = Strand::from_char(source_strand_char).ok_or_else(|| {
            ChainParseError::invalid_strand(source_strand_char, line_number)
        })?;
        
        let source_start = fields[5].parse::<u64>().map_err(|_| {
            ChainParseError::invalid_number("source start", fields[5], line_number)
        })?;
        
        let source_end = fields[6].parse::<u64>().map_err(|_| {
            ChainParseError::invalid_number("source end", fields[6], line_number)
        })?;
        
        // Validate source coordinates
        if source_start > source_end {
            return Err(ChainParseError::invalid_coordinates(
                format!("Source start ({}) > source end ({})", source_start, source_end),
                line_number,
            ));
        }
        if source_end > source_size {
            return Err(ChainParseError::invalid_coordinates(
                format!("Source end ({}) > source size ({})", source_end, source_size),
                line_number,
            ));
        }
        
        // UCSC "query" = our "target" (mapping output coordinates)
        let target_name = fields[7].to_string();
        let target_size = fields[8].parse::<u64>().map_err(|_| {
            ChainParseError::invalid_number("target size", fields[8], line_number)
        })?;
        
        let target_strand_char = fields[9].chars().next().unwrap_or('?');
        let target_strand = Strand::from_char(target_strand_char).ok_or_else(|| {
            ChainParseError::invalid_strand(target_strand_char, line_number)
        })?;
        
        let target_start = fields[10].parse::<u64>().map_err(|_| {
            ChainParseError::invalid_number("target start", fields[10], line_number)
        })?;
        
        let target_end = fields[11].parse::<u64>().map_err(|_| {
            ChainParseError::invalid_number("target end", fields[11], line_number)
        })?;
        
        // Validate target coordinates
        if target_start > target_end {
            return Err(ChainParseError::invalid_coordinates(
                format!("Target start ({}) > target end ({})", target_start, target_end),
                line_number,
            ));
        }
        if target_end > target_size {
            return Err(ChainParseError::invalid_coordinates(
                format!("Target end ({}) > target size ({})", target_end, target_size),
                line_number,
            ));
        }
        
        // Chain ID is optional (field 12)
        let chain_id = fields.get(12).map(|s| s.to_string()).unwrap_or_default();
        
        Ok(Self {
            score,
            target_name,
            target_size,
            target_strand,
            target_start,
            target_end,
            source_name,
            source_size,
            source_strand,
            source_start,
            source_end,
            chain_id,
        })
    }
}


/// A single alignment block from a chain file
///
/// Represents a contiguous aligned region between source and target genomes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChainBlock {
    /// Source chromosome name
    pub source_chrom: String,
    /// Source start position (0-based)
    pub source_start: u64,
    /// Source end position (exclusive)
    pub source_end: u64,
    /// Target chromosome name
    pub target_chrom: String,
    /// Target start position (0-based, already flipped for negative strand)
    pub target_start: u64,
    /// Target end position (exclusive, already flipped for negative strand)
    pub target_end: u64,
    /// Target strand direction
    pub target_strand: Strand,
}

/// Data line in a chain file (size, dt, dq)
/// 
/// Note: UCSC "dt" = gap in UCSC target = gap in our source
///       UCSC "dq" = gap in UCSC query = gap in our target
#[derive(Debug, Clone, Copy)]
struct DataLine {
    /// Alignment block size (ungapped)
    size: u64,
    /// Gap in source sequence (UCSC "dt" = gap in UCSC target)
    source_gap: u64,
    /// Gap in target sequence (UCSC "dq" = gap in UCSC query)
    target_gap: u64,
}

impl DataLine {
    /// Parse a data line (middle line with 3 fields or last line with 1 field)
    fn parse(line: &str, line_number: usize) -> Result<Self, ChainParseError> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        
        match fields.len() {
            1 => {
                // Last line: only size
                let size = fields[0].parse::<u64>().map_err(|_| {
                    ChainParseError::invalid_number("block size", fields[0], line_number)
                })?;
                if size == 0 {
                    return Err(ChainParseError::invalid_data_line(
                        "Block size must be greater than 0",
                        line_number,
                        line,
                    ));
                }
                Ok(Self {
                    size,
                    target_gap: 0,
                    source_gap: 0,
                })
            }
            3 => {
                // Middle line: size dt dq
                // dt = gap in UCSC target = gap in our source
                // dq = gap in UCSC query = gap in our target
                let size = fields[0].parse::<u64>().map_err(|_| {
                    ChainParseError::invalid_number("block size", fields[0], line_number)
                })?;
                if size == 0 {
                    return Err(ChainParseError::invalid_data_line(
                        "Block size must be greater than 0",
                        line_number,
                        line,
                    ));
                }
                let source_gap = fields[1].parse::<u64>().map_err(|_| {
                    ChainParseError::invalid_number("source gap (dt)", fields[1], line_number)
                })?;
                let target_gap = fields[2].parse::<u64>().map_err(|_| {
                    ChainParseError::invalid_number("target gap (dq)", fields[2], line_number)
                })?;
                Ok(Self {
                    size,
                    source_gap,
                    target_gap,
                })
            }
            _ => Err(ChainParseError::invalid_data_line(
                format!("Expected 1 or 3 fields, got {}", fields.len()),
                line_number,
                line,
            )),
        }
    }
}


/// Result of parsing a chain file
#[derive(Debug, Clone)]
pub struct ChainFile {
    /// All alignment blocks
    pub blocks: Vec<ChainBlock>,
    /// Target chromosome sizes
    pub target_chrom_sizes: HashMap<String, u64>,
    /// Source chromosome sizes
    pub source_chrom_sizes: HashMap<String, u64>,
}

impl ChainFile {
    /// Create a new empty ChainFile
    pub fn new() -> Self {
        Self {
            blocks: Vec::new(),
            target_chrom_sizes: HashMap::new(),
            source_chrom_sizes: HashMap::new(),
        }
    }
}

impl Default for ChainFile {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse a chain file from a reader
///
/// This function handles the core parsing logic, supporting any `BufRead` source.
pub fn parse_chain_reader<R: BufRead>(reader: R) -> Result<ChainFile, ChainParseError> {
    let mut result = ChainFile::new();
    let mut current_header: Option<ChainHeader> = None;
    let mut source_pos: u64 = 0;
    let mut target_pos: u64 = 0;
    let mut line_number: usize = 0;
    
    for line_result in reader.lines() {
        line_number += 1;
        let line = line_result?;
        let trimmed = line.trim();
        
        // Skip empty lines and comments
        if trimmed.is_empty() || trimmed.starts_with('#') {
            // Empty line marks end of chain block
            if current_header.is_some() {
                current_header = None;
            }
            continue;
        }
        
        if trimmed.starts_with("chain") {
            // Parse header line
            let header = ChainHeader::parse(trimmed, line_number)?;
            
            // Store chromosome sizes
            result.target_chrom_sizes.insert(header.target_name.clone(), header.target_size);
            result.source_chrom_sizes.insert(header.source_name.clone(), header.source_size);
            
            // Initialize positions
            source_pos = header.source_start;
            target_pos = header.target_start;
            current_header = Some(header);
        } else if let Some(ref header) = current_header {
            // Parse data line
            let data = DataLine::parse(trimmed, line_number)?;
            
            // Calculate target coordinates based on strand
            let (block_target_start, block_target_end) = if header.target_strand == Strand::Plus {
                // Positive strand: direct mapping
                (target_pos, target_pos + data.size)
            } else {
                // Negative strand: flip coordinates
                // target_start = target_size - (target_pos + size)
                // target_end = target_size - target_pos
                let flipped_start = header.target_size - (target_pos + data.size);
                let flipped_end = header.target_size - target_pos;
                (flipped_start, flipped_end)
            };
            
            // Calculate source coordinates based on strand
            let (block_source_start, block_source_end) = if header.source_strand == Strand::Plus {
                // Positive strand: direct mapping
                (source_pos, source_pos + data.size)
            } else {
                // Negative strand: flip coordinates
                let flipped_start = header.source_size - (source_pos + data.size);
                let flipped_end = header.source_size - source_pos;
                (flipped_start, flipped_end)
            };
            
            // Create alignment block
            let block = ChainBlock {
                source_chrom: header.source_name.clone(),
                source_start: block_source_start,
                source_end: block_source_end,
                target_chrom: header.target_name.clone(),
                target_start: block_target_start,
                target_end: block_target_end,
                target_strand: header.target_strand,
            };
            
            result.blocks.push(block);
            
            // Update positions for next block
            source_pos += data.size + data.source_gap;
            target_pos += data.size + data.target_gap;
        }
    }
    
    Ok(result)
}


/// Parse a chain file from a path
///
/// Automatically detects and handles compression:
/// - .gz extension or gzip magic bytes (1f 8b)
/// - .bz2 extension or bzip2 magic bytes (42 5a 68)
/// - Plain text otherwise
pub fn parse_chain_file(path: &Path) -> Result<ChainFile, ChainParseError> {
    use std::fs::File;
    use std::io::Read;
    
    let mut file = File::open(path)?;
    let extension = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    
    // Read first few bytes to detect compression format
    let mut magic = [0u8; 3];
    let bytes_read = file.read(&mut magic)?;
    
    // Reset file position
    drop(file);
    let file = File::open(path)?;
    
    // Detect format by extension or magic bytes
    let format = if extension == "gz" || (bytes_read >= 2 && magic[0] == 0x1f && magic[1] == 0x8b) {
        CompressionFormat::Gzip
    } else if extension == "bz2" || (bytes_read >= 3 && magic[0] == 0x42 && magic[1] == 0x5a && magic[2] == 0x68) {
        // BZ2 magic: "BZh" (0x42 0x5a 0x68)
        CompressionFormat::Bzip2
    } else {
        CompressionFormat::Plain
    };
    
    match format {
        CompressionFormat::Gzip => {
            let decoder = flate2::read::GzDecoder::new(file);
            let reader = BufReader::with_capacity(128 * 1024, decoder);
            parse_chain_reader(reader)
        }
        CompressionFormat::Bzip2 => {
            let decoder = bzip2::read::BzDecoder::new(file);
            let reader = BufReader::with_capacity(128 * 1024, decoder);
            parse_chain_reader(reader)
        }
        CompressionFormat::Plain => {
            let reader = BufReader::with_capacity(128 * 1024, file);
            parse_chain_reader(reader)
        }
    }
}

/// Compression format for chain files
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompressionFormat {
    /// Plain text (uncompressed)
    Plain,
    /// Gzip compressed (.gz)
    Gzip,
    /// Bzip2 compressed (.bz2)
    Bzip2,
}

/// Detect compression format from file path and/or content
pub fn detect_compression(path: &Path) -> Result<CompressionFormat, ChainParseError> {
    use std::fs::File;
    use std::io::Read;
    
    let extension = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    
    // First check by extension
    if extension == "gz" {
        return Ok(CompressionFormat::Gzip);
    }
    if extension == "bz2" {
        return Ok(CompressionFormat::Bzip2);
    }
    
    // Then check by magic bytes
    let mut file = File::open(path)?;
    let mut magic = [0u8; 3];
    let bytes_read = file.read(&mut magic)?;
    
    if bytes_read >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        return Ok(CompressionFormat::Gzip);
    }
    if bytes_read >= 3 && magic[0] == 0x42 && magic[1] == 0x5a && magic[2] == 0x68 {
        return Ok(CompressionFormat::Bzip2);
    }
    
    Ok(CompressionFormat::Plain)
}

/// Parse a chain file from bytes (for testing)
pub fn parse_chain_bytes(data: &[u8]) -> Result<ChainFile, ChainParseError> {
    let reader = BufReader::new(data);
    parse_chain_reader(reader)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_parse_chain_header() {
        // UCSC chain format: chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
        // Our mapping: UCSC target (t) = our source, UCSC query (q) = our target
        let line = "chain 1000 chr1 248956422 + 10000 20000 chr1 249250621 + 10500 20500 1";
        let header = ChainHeader::parse(line, 1).unwrap();
        
        assert_eq!(header.score, 1000);
        // UCSC target fields -> our source fields
        assert_eq!(header.source_name, "chr1");
        assert_eq!(header.source_size, 248956422);
        assert_eq!(header.source_strand, Strand::Plus);
        assert_eq!(header.source_start, 10000);
        assert_eq!(header.source_end, 20000);
        // UCSC query fields -> our target fields
        assert_eq!(header.target_name, "chr1");
        assert_eq!(header.target_size, 249250621);
        assert_eq!(header.target_strand, Strand::Plus);
        assert_eq!(header.target_start, 10500);
        assert_eq!(header.target_end, 20500);
        assert_eq!(header.chain_id, "1");
    }
    
    #[test]
    fn test_parse_chain_header_negative_strand() {
        // UCSC: tStrand='-', qStrand='+'
        // Our mapping: source_strand=Minus (from t), target_strand=Plus (from q)
        let line = "chain 500 chr2 242193529 - 5000 15000 chr2 243199373 + 5500 15500 2";
        let header = ChainHeader::parse(line, 1).unwrap();
        
        assert_eq!(header.source_strand, Strand::Minus);
        assert_eq!(header.target_strand, Strand::Plus);
    }
    
    #[test]
    fn test_parse_chain_header_both_negative_strand() {
        // Both strands can be negative in real chain files
        // UCSC: tStrand='-', qStrand='-'
        // Our mapping: source_strand=Minus (from t), target_strand=Minus (from q)
        let line = "chain 500 chr2 242193529 - 5000 15000 chr2 243199373 - 5500 15500 2";
        let header = ChainHeader::parse(line, 1).unwrap();
        
        assert_eq!(header.source_strand, Strand::Minus);
        assert_eq!(header.target_strand, Strand::Minus);
    }
    
    #[test]
    fn test_parse_data_line_three_fields() {
        // UCSC format: size dt dq
        // Our mapping: dt = source_gap, dq = target_gap
        let line = "100 50 30";
        let data = DataLine::parse(line, 1).unwrap();
        
        assert_eq!(data.size, 100);
        assert_eq!(data.source_gap, 50);  // dt
        assert_eq!(data.target_gap, 30);  // dq
    }
    
    #[test]
    fn test_parse_data_line_one_field() {
        let line = "200";
        let data = DataLine::parse(line, 1).unwrap();
        
        assert_eq!(data.size, 200);
        assert_eq!(data.target_gap, 0);
        assert_eq!(data.source_gap, 0);
    }

    
    #[test]
    fn test_parse_simple_chain() {
        let chain_data = b"\
chain 1000 chr1 1000 + 100 400 chr1 1000 + 100 400 1
100 50 50
100 50 50
100
";
        let result = parse_chain_bytes(chain_data).unwrap();
        
        assert_eq!(result.blocks.len(), 3);
        
        // First block
        assert_eq!(result.blocks[0].source_start, 100);
        assert_eq!(result.blocks[0].source_end, 200);
        assert_eq!(result.blocks[0].target_start, 100);
        assert_eq!(result.blocks[0].target_end, 200);
        
        // Second block (after gaps)
        assert_eq!(result.blocks[1].source_start, 250); // 100 + 100 + 50
        assert_eq!(result.blocks[1].source_end, 350);
        assert_eq!(result.blocks[1].target_start, 250); // 100 + 100 + 50
        assert_eq!(result.blocks[1].target_end, 350);
        
        // Third block
        assert_eq!(result.blocks[2].source_start, 400); // 250 + 100 + 50
        assert_eq!(result.blocks[2].source_end, 500);
    }
    
    #[test]
    fn test_parse_negative_strand_chain() {
        // UCSC: tStrand='-' means the source (user query) is on negative strand
        // Coordinates should be flipped for the SOURCE, not target
        // But wait - in UCSC, tStrand is always '+', qStrand can be '-'
        // Let's use qStrand='-' which means target is on negative strand
        let chain_data = b"\
chain 500 chr2 1000 + 100 300 chr2 1000 - 100 300 2
100 50 50
100
";
        let result = parse_chain_bytes(chain_data).unwrap();
        
        assert_eq!(result.blocks.len(), 2);
        
        // Source (from UCSC target): positive strand, no flip
        // Target (from UCSC query): negative strand, coordinates flipped
        // First block: target_size=1000, target_pos=100, size=100
        // flipped_start = 1000 - (100 + 100) = 800
        // flipped_end = 1000 - 100 = 900
        assert_eq!(result.blocks[0].source_start, 100);
        assert_eq!(result.blocks[0].source_end, 200);
        assert_eq!(result.blocks[0].target_start, 800);
        assert_eq!(result.blocks[0].target_end, 900);
        assert_eq!(result.blocks[0].target_strand, Strand::Minus);
        
        // Second block: target_pos=250 (100+100+50), size=100
        // flipped_start = 1000 - (250 + 100) = 650
        // flipped_end = 1000 - 250 = 750
        assert_eq!(result.blocks[1].source_start, 250);
        assert_eq!(result.blocks[1].source_end, 350);
        assert_eq!(result.blocks[1].target_start, 650);
        assert_eq!(result.blocks[1].target_end, 750);
    }
    
    #[test]
    fn test_parse_multiple_chains() {
        let chain_data = b"\
chain 1000 chr1 1000 + 0 100 chr1 1000 + 0 100 1
100

chain 500 chr2 2000 + 0 50 chr2 2000 + 0 50 2
50
";
        let result = parse_chain_bytes(chain_data).unwrap();
        
        assert_eq!(result.blocks.len(), 2);
        assert_eq!(result.blocks[0].source_chrom, "chr1");
        assert_eq!(result.blocks[1].source_chrom, "chr2");
        
        // Check chromosome sizes
        assert_eq!(result.target_chrom_sizes.get("chr1"), Some(&1000));
        assert_eq!(result.target_chrom_sizes.get("chr2"), Some(&2000));
    }
    
    #[test]
    fn test_parse_chain_with_comments() {
        let chain_data = b"\
# This is a comment
chain 1000 chr1 1000 + 0 100 chr1 1000 + 0 100 1
100
";
        let result = parse_chain_bytes(chain_data).unwrap();
        assert_eq!(result.blocks.len(), 1);
    }
    
    #[test]
    fn test_parse_chain_error_line_number() {
        let chain_data = b"\
chain 1000 chr1 1000 + 0 100 chr1 1000 + 0 100 1
invalid_data
";
        let result = parse_chain_bytes(chain_data);
        assert!(result.is_err());
        
        let err = result.unwrap_err();
        assert_eq!(err.line_number, Some(2));
    }
    
    #[test]
    fn test_compression_format_enum() {
        assert_eq!(CompressionFormat::Plain, CompressionFormat::Plain);
        assert_eq!(CompressionFormat::Gzip, CompressionFormat::Gzip);
        assert_eq!(CompressionFormat::Bzip2, CompressionFormat::Bzip2);
        assert_ne!(CompressionFormat::Plain, CompressionFormat::Gzip);
    }
    
    #[test]
    fn test_error_invalid_strand() {
        let line = "chain 1000 chr1 1000 X 0 100 chr1 1000 + 0 100 1";
        let result = ChainHeader::parse(line, 5);
        assert!(result.is_err());
        
        let err = result.unwrap_err();
        assert_eq!(err.kind, ChainParseErrorKind::InvalidStrand);
        assert_eq!(err.line_number, Some(5));
        assert!(err.message.contains("'X'"));
    }
    
    #[test]
    fn test_error_invalid_number() {
        let line = "chain abc chr1 1000 + 0 100 chr1 1000 + 0 100 1";
        let result = ChainHeader::parse(line, 3);
        assert!(result.is_err());
        
        let err = result.unwrap_err();
        assert_eq!(err.kind, ChainParseErrorKind::InvalidNumber);
        assert_eq!(err.line_number, Some(3));
        assert!(err.message.contains("abc"));
    }
    
    #[test]
    fn test_error_invalid_coordinates_start_gt_end() {
        // target_start (200) > target_end (100)
        let line = "chain 1000 chr1 1000 + 200 100 chr1 1000 + 0 100 1";
        let result = ChainHeader::parse(line, 1);
        assert!(result.is_err());
        
        let err = result.unwrap_err();
        assert_eq!(err.kind, ChainParseErrorKind::InvalidCoordinates);
        assert!(err.message.contains("start") && err.message.contains("end"));
    }
    
    #[test]
    fn test_error_invalid_coordinates_end_gt_size() {
        // target_end (2000) > target_size (1000)
        let line = "chain 1000 chr1 1000 + 0 2000 chr1 1000 + 0 100 1";
        let result = ChainHeader::parse(line, 1);
        assert!(result.is_err());
        
        let err = result.unwrap_err();
        assert_eq!(err.kind, ChainParseErrorKind::InvalidCoordinates);
        assert!(err.message.contains("end") && err.message.contains("size"));
    }
    
    #[test]
    fn test_error_too_few_fields() {
        let line = "chain 1000 chr1 1000 + 0 100";
        let result = ChainHeader::parse(line, 1);
        assert!(result.is_err());
        
        let err = result.unwrap_err();
        assert_eq!(err.kind, ChainParseErrorKind::InvalidHeader);
        assert!(err.message.contains("12"));
    }
    
    #[test]
    fn test_error_invalid_data_line_fields() {
        let line = "100 50";  // 2 fields, should be 1 or 3
        let result = DataLine::parse(line, 10);
        assert!(result.is_err());
        
        let err = result.unwrap_err();
        assert_eq!(err.kind, ChainParseErrorKind::InvalidDataLine);
        assert_eq!(err.line_number, Some(10));
    }
    
    #[test]
    fn test_error_kind_check() {
        let err = ChainParseError::invalid_strand('?', 1);
        assert!(err.is_kind(ChainParseErrorKind::InvalidStrand));
        assert!(!err.is_kind(ChainParseErrorKind::InvalidNumber));
    }
    
    #[test]
    fn test_error_display() {
        let err = ChainParseError::with_line("Test error message", 42);
        let display = format!("{}", err);
        assert!(display.contains("Line 42"));
        assert!(display.contains("Test error message"));
    }
}


#[cfg(test)]
mod integration_tests {
    use super::*;
    use std::path::PathBuf;
    
    /// Test parsing a real gzip chain file (if available)
    #[test]
    fn test_parse_real_chain_file_gz() {
        let chain_path = PathBuf::from("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz");
        
        if !chain_path.exists() {
            eprintln!("Skipping test: chain file not found at {:?}", chain_path);
            return;
        }
        
        // Verify format detection
        let format = detect_compression(&chain_path).unwrap();
        assert_eq!(format, CompressionFormat::Gzip, "Should detect gzip format");
        
        let result = parse_chain_file(&chain_path);
        assert!(result.is_ok(), "Failed to parse chain file: {:?}", result.err());
        
        let chain_file = result.unwrap();
        
        // Verify we got some blocks
        assert!(!chain_file.blocks.is_empty(), "No blocks parsed from chain file");
        
        // Verify chromosome sizes were captured
        assert!(!chain_file.target_chrom_sizes.is_empty(), "No target chrom sizes");
        assert!(!chain_file.source_chrom_sizes.is_empty(), "No source chrom sizes");
        
        // Print some stats
        eprintln!("Parsed {} alignment blocks", chain_file.blocks.len());
        eprintln!("Target chromosomes: {}", chain_file.target_chrom_sizes.len());
        eprintln!("Source chromosomes: {}", chain_file.source_chrom_sizes.len());
        
        // Verify block structure
        for block in chain_file.blocks.iter().take(5) {
            assert!(block.source_end > block.source_start, "Invalid source range");
            assert!(block.target_end > block.target_start, "Invalid target range");
            assert!(!block.source_chrom.is_empty(), "Empty source chrom");
            assert!(!block.target_chrom.is_empty(), "Empty target chrom");
        }
    }
    
    /// Test that gzip and plain text parsing produce identical results
    #[test]
    fn test_gz_plain_equivalence() {
        use std::io::Write;
        use flate2::write::GzEncoder;
        use flate2::Compression;
        
        let chain_data = b"\
chain 1000 chr1 1000 + 100 400 chr1 1000 + 100 400 1
100 50 50
100
";
        
        // Parse plain text
        let plain_result = parse_chain_bytes(chain_data).unwrap();
        
        // Create gzip compressed version
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(chain_data).unwrap();
        let gz_data = encoder.finish().unwrap();
        
        // Write to temp file and parse
        let temp_dir = std::env::temp_dir();
        let gz_path = temp_dir.join("test_chain.chain.gz");
        std::fs::write(&gz_path, &gz_data).unwrap();
        
        let gz_result = parse_chain_file(&gz_path).unwrap();
        
        // Clean up
        let _ = std::fs::remove_file(&gz_path);
        
        // Compare results
        assert_eq!(plain_result.blocks.len(), gz_result.blocks.len());
        for (plain_block, gz_block) in plain_result.blocks.iter().zip(gz_result.blocks.iter()) {
            assert_eq!(plain_block, gz_block);
        }
    }
    
    /// Test that bzip2 and plain text parsing produce identical results
    #[test]
    fn test_bz2_plain_equivalence() {
        use std::io::Write;
        use bzip2::write::BzEncoder;
        use bzip2::Compression;
        
        let chain_data = b"\
chain 1000 chr1 1000 + 100 400 chr1 1000 + 100 400 1
100 50 50
100
";
        
        // Parse plain text
        let plain_result = parse_chain_bytes(chain_data).unwrap();
        
        // Create bzip2 compressed version
        let mut encoder = BzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(chain_data).unwrap();
        let bz2_data = encoder.finish().unwrap();
        
        // Write to temp file and parse
        let temp_dir = std::env::temp_dir();
        let bz2_path = temp_dir.join("test_chain.chain.bz2");
        std::fs::write(&bz2_path, &bz2_data).unwrap();
        
        let bz2_result = parse_chain_file(&bz2_path).unwrap();
        
        // Clean up
        let _ = std::fs::remove_file(&bz2_path);
        
        // Compare results
        assert_eq!(plain_result.blocks.len(), bz2_result.blocks.len());
        for (plain_block, bz2_block) in plain_result.blocks.iter().zip(bz2_result.blocks.iter()) {
            assert_eq!(plain_block, bz2_block);
        }
    }
    
    /// Test format detection by magic bytes (without extension)
    #[test]
    fn test_format_detection_by_magic() {
        use std::io::Write;
        use flate2::write::GzEncoder;
        use flate2::Compression;
        
        let chain_data = b"chain 1000 chr1 1000 + 0 100 chr1 1000 + 0 100 1\n100\n";
        
        // Create gzip file without .gz extension
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(chain_data).unwrap();
        let gz_data = encoder.finish().unwrap();
        
        let temp_dir = std::env::temp_dir();
        let path_no_ext = temp_dir.join("test_chain_no_ext");
        std::fs::write(&path_no_ext, &gz_data).unwrap();
        
        // Should detect gzip by magic bytes
        let format = detect_compression(&path_no_ext).unwrap();
        assert_eq!(format, CompressionFormat::Gzip, "Should detect gzip by magic bytes");
        
        // Should parse correctly
        let result = parse_chain_file(&path_no_ext);
        assert!(result.is_ok(), "Should parse gzip file without extension");
        
        // Clean up
        let _ = std::fs::remove_file(&path_no_ext);
    }
}
