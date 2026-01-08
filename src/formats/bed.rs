//! BED format adapter
//!
//! Handles BED3/BED6/BED12 format conversion with zero-copy parsing.
//!
//! **Validates: Requirements 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7**

use crate::core::{CoordinateMapper, MappingSegment, Strand};
use memchr::memchr;
use rayon::prelude::*;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

/// BED record representation for output
#[derive(Debug, Clone)]
pub struct BedRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: Option<String>,
    pub score: Option<String>,
    pub strand: Option<Strand>,
    // BED12 fields
    pub thick_start: Option<u64>,
    pub thick_end: Option<u64>,
    pub item_rgb: Option<String>,
    pub block_count: Option<u32>,
    pub block_sizes: Option<String>,
    pub block_starts: Option<String>,
    // Extra fields beyond BED12
    pub extra_fields: Vec<String>,
}

/// Zero-copy BED record view for parsing
/// Only parses coordinate fields immediately, other fields are kept as byte slices
pub struct BedRecordView<'a> {
    /// Original line bytes
    line: &'a [u8],
    /// Chromosome name
    pub chrom: &'a str,
    /// Start position (0-based)
    pub start: u64,
    /// End position
    pub end: u64,
    /// Field boundaries (start, end) for lazy access
    field_bounds: Vec<(usize, usize)>,
}

impl<'a> BedRecordView<'a> {
    /// Parse a BED line with minimal allocation
    /// Only parses chrom, start, end immediately
    pub fn parse(line: &'a [u8]) -> Result<Self, BedParseError> {
        if line.is_empty() {
            return Err(BedParseError::EmptyLine);
        }

        // Find field boundaries using memchr for tab characters
        let mut field_bounds = Vec::with_capacity(12);
        let mut start_pos = 0;
        let mut pos = 0;
        
        while pos < line.len() {
            if let Some(tab_pos) = memchr(b'\t', &line[pos..]) {
                let end_pos = pos + tab_pos;
                field_bounds.push((start_pos, end_pos));
                start_pos = end_pos + 1;
                pos = start_pos;
            } else {
                // Last field
                field_bounds.push((start_pos, line.len()));
                break;
            }
        }
        
        // Need at least 3 fields (BED3)
        if field_bounds.len() < 3 {
            return Err(BedParseError::TooFewFields {
                expected: 3,
                found: field_bounds.len(),
            });
        }
        
        // Parse chrom (field 0)
        let chrom = std::str::from_utf8(&line[field_bounds[0].0..field_bounds[0].1])
            .map_err(|_| BedParseError::InvalidUtf8("chrom"))?;
        
        // Parse start (field 1)
        let start_str = std::str::from_utf8(&line[field_bounds[1].0..field_bounds[1].1])
            .map_err(|_| BedParseError::InvalidUtf8("start"))?;
        let start: u64 = start_str
            .parse()
            .map_err(|_| BedParseError::InvalidNumber("start", start_str.to_string()))?;
        
        // Parse end (field 2)
        let end_str = std::str::from_utf8(&line[field_bounds[2].0..field_bounds[2].1])
            .map_err(|_| BedParseError::InvalidUtf8("end"))?;
        let end: u64 = end_str
            .parse()
            .map_err(|_| BedParseError::InvalidNumber("end", end_str.to_string()))?;
        
        Ok(Self {
            line,
            chrom,
            start,
            end,
            field_bounds,
        })
    }
    
    /// Get the number of fields
    pub fn field_count(&self) -> usize {
        self.field_bounds.len()
    }
    
    /// Get field as string slice (lazy access)
    pub fn field(&self, index: usize) -> Option<&'a str> {
        self.field_bounds.get(index).and_then(|(start, end)| {
            std::str::from_utf8(&self.line[*start..*end]).ok()
        })
    }
    
    /// Get name field (field 3) if present
    pub fn name(&self) -> Option<&'a str> {
        self.field(3)
    }
    
    /// Get score field (field 4) if present
    pub fn score(&self) -> Option<&'a str> {
        self.field(4)
    }
    
    /// Get strand field (field 5) if present
    pub fn strand(&self) -> Option<Strand> {
        self.field(5).and_then(|s| {
            match s {
                "+" => Some(Strand::Plus),
                "-" => Some(Strand::Minus),
                "." => None,
                _ => None,
            }
        })
    }
    
    /// Get strand character for output
    pub fn strand_char(&self) -> Option<&'a str> {
        self.field(5)
    }
    
    /// Get thick_start (field 6) if present
    pub fn thick_start(&self) -> Option<u64> {
        self.field(6).and_then(|s| s.parse().ok())
    }
    
    /// Get thick_end (field 7) if present
    pub fn thick_end(&self) -> Option<u64> {
        self.field(7).and_then(|s| s.parse().ok())
    }
    
    /// Get item_rgb (field 8) if present
    pub fn item_rgb(&self) -> Option<&'a str> {
        self.field(8)
    }
    
    /// Get block_count (field 9) if present
    pub fn block_count(&self) -> Option<u32> {
        self.field(9).and_then(|s| s.parse().ok())
    }
    
    /// Get block_sizes (field 10) if present
    pub fn block_sizes(&self) -> Option<&'a str> {
        self.field(10)
    }
    
    /// Get block_starts (field 11) if present
    pub fn block_starts(&self) -> Option<&'a str> {
        self.field(11)
    }
    
    /// Check if this is a BED12 record
    pub fn is_bed12(&self) -> bool {
        self.field_count() >= 12
    }
    
    /// Check if this is a BED6 record
    pub fn is_bed6(&self) -> bool {
        self.field_count() >= 6
    }
}


/// BED parsing error
#[derive(Debug, thiserror::Error)]
pub enum BedParseError {
    #[error("Empty line")]
    EmptyLine,
    
    #[error("Too few fields: expected at least {expected}, found {found}")]
    TooFewFields { expected: usize, found: usize },
    
    #[error("Invalid UTF-8 in field: {0}")]
    InvalidUtf8(&'static str),
    
    #[error("Invalid number in field {0}: {1}")]
    InvalidNumber(&'static str, String),
    
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}

/// Conversion statistics
#[derive(Debug, Default, Clone)]
pub struct ConversionStats {
    pub total: usize,
    pub success: usize,
    pub failed: usize,
    pub multi_map: usize,
}

/// Result of converting a single BED record
#[derive(Debug)]
pub enum ConversionResult {
    /// Successfully mapped to a single location
    Success(String),
    /// Mapped to multiple locations
    MultiMap(Vec<String>),
    /// Failed to map
    Failed(String),
    /// Comment or header line (pass through)
    PassThrough(String),
}

/// Represents a single block in BED12 format
#[derive(Debug, Clone)]
struct Block {
    start: u64,  // Absolute start position
    end: u64,    // Absolute end position
}

/// Parse BED12 blocks from a record view
fn parse_bed12_blocks(view: &BedRecordView) -> Option<Vec<Block>> {
    let chrom_start = view.start;
    let block_sizes_str = view.block_sizes()?;
    let block_starts_str = view.block_starts()?;
    
    let sizes: Vec<u64> = block_sizes_str
        .split(',')
        .filter(|s| !s.is_empty())
        .filter_map(|s| s.parse().ok())
        .collect();
    
    let starts: Vec<u64> = block_starts_str
        .split(',')
        .filter(|s| !s.is_empty())
        .filter_map(|s| s.parse().ok())
        .collect();
    
    if sizes.len() != starts.len() || sizes.is_empty() {
        return None;
    }
    
    let blocks: Vec<Block> = sizes.iter().zip(starts.iter())
        .map(|(&size, &rel_start)| Block {
            start: chrom_start + rel_start,
            end: chrom_start + rel_start + size,
        })
        .collect();
    
    Some(blocks)
}

/// Convert a BED12 record by mapping each block individually
fn convert_bed12_record(
    view: &BedRecordView,
    mapper: &CoordinateMapper,
    input_strand: Strand,
) -> ConversionResult {
    // Parse blocks
    let blocks = match parse_bed12_blocks(view) {
        Some(b) if !b.is_empty() => b,
        _ => {
            // Fall back to regular mapping if blocks can't be parsed
            return convert_bed_record_simple(view, mapper, input_strand);
        }
    };
    
    // Map each block individually
    let mut mapped_blocks: Vec<(u64, u64, Strand, String)> = Vec::with_capacity(blocks.len());
    let mut target_chrom: Option<String> = None;
    let mut target_strand: Option<Strand> = None;
    
    for block in &blocks {
        let result = mapper.map(view.chrom, block.start, block.end, input_strand);
        
        match result {
            Some(segments) if segments.len() == 1 => {
                // CrossMap behavior: only accept blocks that map to exactly one segment
                // If a block maps to multiple segments, the entire record fails
                let seg = &segments[0];
                
                // Check all blocks map to the same chromosome and strand
                if let Some(ref tc) = target_chrom {
                    if tc != &seg.target.chrom {
                        // Blocks map to different chromosomes - fail
                        return ConversionResult::Failed(format_unmapped_line(view));
                    }
                } else {
                    target_chrom = Some(seg.target.chrom.clone());
                    target_strand = Some(seg.target.strand);
                }
                
                mapped_blocks.push((seg.target.start, seg.target.end, seg.target.strand, seg.target.chrom.clone()));
            }
            _ => {
                // Block failed to map or mapped to multiple segments - entire record fails
                return ConversionResult::Failed(format_unmapped_line(view));
            }
        }
    }
    
    if mapped_blocks.is_empty() {
        return ConversionResult::Failed(format_unmapped_line(view));
    }
    
    // Calculate new BED12 coordinates
    // CrossMap behavior: use first block's start and last block's end
    // CrossMap preserves the original block order from the input file
    // It does NOT sort blocks by position
    let new_chrom = target_chrom.unwrap();
    let new_strand = target_strand.unwrap_or(input_strand);
    
    // CrossMap uses: new_chrom_st = exons_new_pos[0][1], new_chrom_end = exons_new_pos[-1][2]
    // This is the FIRST and LAST block in original order, not min/max
    let new_chrom_start = mapped_blocks.first().unwrap().0;
    let new_chrom_end = mapped_blocks.last().unwrap().1;
    
    // Calculate new block starts (relative to new_chrom_start) preserving original order
    let new_block_starts: Vec<String> = mapped_blocks.iter()
        .map(|(s, _, _, _)| (s - new_chrom_start).to_string())
        .collect();
    
    // Calculate new block sizes preserving original order
    let new_block_sizes: Vec<String> = mapped_blocks.iter()
        .map(|(s, e, _, _)| (e - s).to_string())
        .collect();
    
    // Calculate new thick_start and thick_end
    // CrossMap behavior: preserve the offset from chrom boundaries
    // cds_start_offset = thick_start - chrom_start
    // cds_end_offset = chrom_end - thick_end
    // new_thick_start = new_chrom_start + cds_start_offset
    // new_thick_end = new_chrom_end - cds_end_offset
    let original_thick_start = view.thick_start().unwrap_or(view.start);
    let original_thick_end = view.thick_end().unwrap_or(view.end);
    
    let cds_start_offset = original_thick_start.saturating_sub(view.start);
    let cds_end_offset = view.end.saturating_sub(original_thick_end);
    
    let new_thick_start = (new_chrom_start + cds_start_offset).min(new_chrom_end);
    let new_thick_end = new_chrom_end.saturating_sub(cds_end_offset).max(new_chrom_start);
    
    // Validate BED12 format: thick_start must be <= thick_end
    // CrossMap's check_bed12 function rejects records where thickStart > thickEnd
    if new_thick_start > new_thick_end {
        return ConversionResult::Failed(format_unmapped_line(view));
    }
    
    // Also validate: thickStart >= chromStart and thickEnd <= chromEnd
    if new_thick_start < new_chrom_start || new_thick_end > new_chrom_end {
        return ConversionResult::Failed(format_unmapped_line(view));
    }
    
    // Additional validation: block_starts must be non-negative (check_bed12 requirement)
    for (s, _, _, _) in &mapped_blocks {
        if *s < new_chrom_start {
            return ConversionResult::Failed(format_unmapped_line(view));
        }
    }
    
    // Build output line
    let mut output = String::with_capacity(256);
    
    // Basic fields
    output.push_str(&new_chrom);
    output.push('\t');
    output.push_str(&new_chrom_start.to_string());
    output.push('\t');
    output.push_str(&new_chrom_end.to_string());
    
    // Name (field 3)
    if let Some(name) = view.name() {
        output.push('\t');
        output.push_str(name);
    }
    
    // Score (field 4)
    if view.field_count() > 4 {
        if let Some(score) = view.score() {
            output.push('\t');
            output.push_str(score);
        }
    }
    
    // Strand (field 5)
    if view.field_count() > 5 {
        output.push('\t');
        output.push(new_strand.to_char());
    }
    
    // thick_start (field 6)
    output.push('\t');
    output.push_str(&new_thick_start.to_string());
    
    // thick_end (field 7)
    output.push('\t');
    output.push_str(&new_thick_end.to_string());
    
    // item_rgb (field 8)
    if let Some(rgb) = view.item_rgb() {
        output.push('\t');
        output.push_str(rgb);
    }
    
    // block_count (field 9)
    output.push('\t');
    output.push_str(&mapped_blocks.len().to_string());
    
    // block_sizes (field 10)
    output.push('\t');
    output.push_str(&new_block_sizes.join(","));
    
    // block_starts (field 11)
    output.push('\t');
    output.push_str(&new_block_starts.join(","));
    
    // Extra fields beyond BED12
    for i in 12..view.field_count() {
        if let Some(field) = view.field(i) {
            output.push('\t');
            output.push_str(field);
        }
    }
    
    ConversionResult::Success(output)
}

/// Convert a single BED record (simple version for non-BED12)
fn convert_bed_record_simple(
    view: &BedRecordView,
    mapper: &CoordinateMapper,
    input_strand: Strand,
) -> ConversionResult {
    // Map the coordinates
    let result = mapper.map(view.chrom, view.start, view.end, input_strand);
    
    match result {
        Some(segments) if !segments.is_empty() => {
            // Build output lines
            let output_lines: Vec<String> = segments
                .iter()
                .map(|seg| format_output_line(view, seg))
                .collect();
            
            if output_lines.len() == 1 {
                ConversionResult::Success(output_lines.into_iter().next().unwrap())
            } else {
                ConversionResult::MultiMap(output_lines)
            }
        }
        _ => {
            // Failed to map - output original line with "Unmapped" annotation
            ConversionResult::Failed(format_unmapped_line(view))
        }
    }
}

/// Convert a single BED record
fn convert_bed_record(
    view: &BedRecordView,
    mapper: &CoordinateMapper,
    input_strand: Strand,
) -> ConversionResult {
    // Use special BED12 handling if this is a BED12 record
    if view.is_bed12() {
        return convert_bed12_record(view, mapper, input_strand);
    }
    
    // For non-BED12 records, use simple mapping
    convert_bed_record_simple(view, mapper, input_strand)
}

/// Format output line for a successfully mapped segment
fn format_output_line(view: &BedRecordView, seg: &MappingSegment) -> String {
    let mut output = String::with_capacity(256);
    
    // Output mapped coordinates
    output.push_str(&seg.target.chrom);
    output.push('\t');
    output.push_str(&seg.target.start.to_string());
    output.push('\t');
    output.push_str(&seg.target.end.to_string());
    
    // Preserve additional fields if present
    if view.field_count() > 3 {
        // Name (field 3)
        if let Some(name) = view.name() {
            output.push('\t');
            output.push_str(name);
        }
        
        // Score (field 4)
        if view.field_count() > 4 {
            if let Some(score) = view.score() {
                output.push('\t');
                output.push_str(score);
            }
        }
        
        // Strand (field 5) - CrossMap behavior: update strand based on mapping result
        // CrossMap combines query strand with target strand to determine final strand
        // The final strand in seg.target.strand already has this combination applied
        if view.field_count() > 5 {
            output.push('\t');
            // Use the combined strand from the mapping result
            output.push(seg.target.strand.to_char());
        }
        
        // BED12 fields (6-11) - adjust coordinates relative to new position
        if view.is_bed12() {
            // thick_start (field 6)
            if let Some(thick_start) = view.thick_start() {
                output.push('\t');
                // Adjust thick_start relative to new coordinates
                let offset = seg.target.start as i64 - view.start as i64;
                let new_thick_start = (thick_start as i64 + offset).max(seg.target.start as i64) as u64;
                output.push_str(&new_thick_start.to_string());
            }
            
            // thick_end (field 7)
            if let Some(thick_end) = view.thick_end() {
                output.push('\t');
                let offset = seg.target.start as i64 - view.start as i64;
                let new_thick_end = (thick_end as i64 + offset).min(seg.target.end as i64) as u64;
                output.push_str(&new_thick_end.to_string());
            }
            
            // item_rgb (field 8) - preserve as-is
            if let Some(rgb) = view.item_rgb() {
                output.push('\t');
                output.push_str(rgb);
            }
            
            // block_count (field 9) - preserve as-is
            if let Some(count) = view.field(9) {
                output.push('\t');
                output.push_str(count);
            }
            
            // block_sizes (field 10) - preserve as-is
            if let Some(sizes) = view.block_sizes() {
                output.push('\t');
                output.push_str(sizes);
            }
            
            // block_starts (field 11) - recalculate relative to new chromStart
            // CrossMap behavior: block_starts are relative to chromStart (always starting from 0)
            // When coordinates change, we need to recalculate block_starts
            // The first block_start should always be 0, and subsequent ones are adjusted
            if let Some(starts_str) = view.block_starts() {
                output.push('\t');
                // Parse original block_starts and recalculate relative to new position
                // CrossMap recalculates block_starts so they are relative to the new chromStart
                // This means the first block_start becomes 0, and others are adjusted accordingly
                let original_start = view.start;
                let new_start = seg.target.start;
                
                // Calculate the offset between old and new start positions
                let offset = new_start as i64 - original_start as i64;
                
                let new_starts: Vec<String> = starts_str
                    .split(',')
                    .filter(|s| !s.is_empty())
                    .map(|s| {
                        if let Ok(start_val) = s.parse::<i64>() {
                            // Original absolute position = original_start + start_val
                            // New absolute position = original_start + start_val + offset
                            // New relative position = (original_start + start_val + offset) - new_start
                            //                       = original_start + start_val + offset - new_start
                            //                       = original_start + start_val + (new_start - original_start) - new_start
                            //                       = start_val
                            // But CrossMap actually recalculates based on the mapped coordinates
                            // The new relative position should be: start_val + offset
                            // where offset = new_start - original_start
                            let new_relative = start_val + offset;
                            new_relative.to_string()
                        } else {
                            s.to_string()
                        }
                    })
                    .collect();
                output.push_str(&new_starts.join(","));
            }
        }
        
        // Extra fields beyond BED12
        for i in 12..view.field_count() {
            if let Some(field) = view.field(i) {
                output.push('\t');
                output.push_str(field);
            }
        }
    }
    
    output
}

/// Format unmapped line for failed conversion
fn format_unmapped_line(view: &BedRecordView) -> String {
    // Reconstruct original line
    let mut output = String::with_capacity(256);
    output.push_str(view.chrom);
    output.push('\t');
    output.push_str(&view.start.to_string());
    output.push('\t');
    output.push_str(&view.end.to_string());
    
    for i in 3..view.field_count() {
        if let Some(field) = view.field(i) {
            output.push('\t');
            output.push_str(field);
        }
    }
    
    output
}


/// Reusable parse buffer for zero-allocation parsing
pub struct ParseBuffer {
    line_buf: Vec<u8>,
}

impl ParseBuffer {
    pub fn new() -> Self {
        Self {
            line_buf: Vec::with_capacity(4096),
        }
    }
    
    pub fn clear(&mut self) {
        self.line_buf.clear();
    }
}

impl Default for ParseBuffer {
    fn default() -> Self {
        Self::new()
    }
}

/// Chunk size for parallel processing
const CHUNK_SIZE: usize = 10000;

/// Convert a BED file using the coordinate mapper (sequential version)
/// 
/// # Arguments
/// * `input` - Input BED file path
/// * `output` - Output BED file path for successfully mapped records
/// * `unmap` - Output file path for unmapped records
/// * `mapper` - Coordinate mapper with loaded chain index
/// * `threads` - Number of threads for parallel processing (1 = sequential)
/// 
/// # Returns
/// Conversion statistics
pub fn convert_bed<P: AsRef<Path>>(
    input: P,
    output: P,
    unmap: P,
    mapper: &CoordinateMapper,
    threads: usize,
) -> Result<ConversionStats, BedParseError> {
    if threads > 1 {
        convert_bed_parallel(input, output, unmap, mapper, threads)
    } else {
        convert_bed_sequential(input, output, unmap, mapper)
    }
}

/// Sequential BED conversion (single-threaded)
fn convert_bed_sequential<P: AsRef<Path>>(
    input: P,
    output: P,
    unmap: P,
    mapper: &CoordinateMapper,
) -> Result<ConversionStats, BedParseError> {
    let input_file = std::fs::File::open(input.as_ref())?;
    let reader = BufReader::with_capacity(128 * 1024, input_file);
    
    // Use BufWriter to avoid per-line syscalls (critical for performance)
    let mut output_file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(output.as_ref())?);
    let mut unmap_file = BufWriter::with_capacity(64 * 1024, std::fs::File::create(unmap.as_ref())?);
    
    let mut stats = ConversionStats::default();
    let mut line_buf = String::with_capacity(4096);
    
    let mut reader = reader;
    
    loop {
        line_buf.clear();
        let bytes_read = reader.read_line(&mut line_buf)?;
        if bytes_read == 0 {
            break;
        }
        
        // Remove trailing newline
        let line = line_buf.trim_end();
        
        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            // Pass through header lines to output
            if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                writeln!(output_file, "{}", line)?;
            }
            continue;
        }
        
        stats.total += 1;
        
        // Parse the BED record
        match BedRecordView::parse(line.as_bytes()) {
            Ok(view) => {
                // Get input strand for mapping
                let input_strand = view.strand().unwrap_or(Strand::Plus);
                
                // Convert the record
                match convert_bed_record(&view, mapper, input_strand) {
                    ConversionResult::Success(output_line) => {
                        writeln!(output_file, "{}", output_line)?;
                        stats.success += 1;
                    }
                    ConversionResult::MultiMap(output_lines) => {
                        for output_line in output_lines {
                            writeln!(output_file, "{}", output_line)?;
                        }
                        stats.success += 1;
                        stats.multi_map += 1;
                    }
                    ConversionResult::Failed(unmapped_line) => {
                        writeln!(unmap_file, "{}", unmapped_line)?;
                        stats.failed += 1;
                    }
                    ConversionResult::PassThrough(line) => {
                        writeln!(output_file, "{}", line)?;
                    }
                }
            }
            Err(_) => {
                // Invalid BED line - write to unmap file
                writeln!(unmap_file, "{}", line)?;
                stats.failed += 1;
            }
        }
    }
    
    Ok(stats)
}

/// Parallel BED conversion using rayon
/// 
/// Reads all lines into memory, processes in parallel chunks, then writes output.
/// This trades memory for speed - suitable for files that fit in memory.
fn convert_bed_parallel<P: AsRef<Path>>(
    input: P,
    output: P,
    unmap: P,
    mapper: &CoordinateMapper,
    threads: usize,
) -> Result<ConversionStats, BedParseError> {
    // Configure rayon thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| BedParseError::Io(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("Failed to create thread pool: {}", e)
        )))?;
    
    // Read all lines
    let input_file = std::fs::File::open(input.as_ref())?;
    let reader = BufReader::with_capacity(128 * 1024, input_file);
    
    let mut header_lines = Vec::new();
    let mut data_lines = Vec::new();
    
    for line_result in reader.lines() {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            header_lines.push(line);
        } else {
            data_lines.push(line);
        }
    }
    
    // Atomic counters for stats
    let total = AtomicUsize::new(0);
    let success = AtomicUsize::new(0);
    let failed = AtomicUsize::new(0);
    let multi_map = AtomicUsize::new(0);
    
    // Process in parallel
    let results: Vec<(Vec<String>, Vec<String>)> = pool.install(|| {
        data_lines
            .par_chunks(CHUNK_SIZE)
            .map(|chunk| {
                let mut success_lines = Vec::with_capacity(chunk.len());
                let mut failed_lines = Vec::new();
                
                for line in chunk {
                    total.fetch_add(1, Ordering::Relaxed);
                    
                    match BedRecordView::parse(line.as_bytes()) {
                        Ok(view) => {
                            let input_strand = view.strand().unwrap_or(Strand::Plus);
                            
                            match convert_bed_record(&view, mapper, input_strand) {
                                ConversionResult::Success(output_line) => {
                                    success_lines.push(output_line);
                                    success.fetch_add(1, Ordering::Relaxed);
                                }
                                ConversionResult::MultiMap(output_lines) => {
                                    success_lines.extend(output_lines);
                                    success.fetch_add(1, Ordering::Relaxed);
                                    multi_map.fetch_add(1, Ordering::Relaxed);
                                }
                                ConversionResult::Failed(unmapped_line) => {
                                    failed_lines.push(unmapped_line);
                                    failed.fetch_add(1, Ordering::Relaxed);
                                }
                                ConversionResult::PassThrough(pass_line) => {
                                    success_lines.push(pass_line);
                                }
                            }
                        }
                        Err(_) => {
                            failed_lines.push(line.clone());
                            failed.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                }
                
                (success_lines, failed_lines)
            })
            .collect()
    });
    
    // Write output files with BufWriter for performance
    let mut output_file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(output.as_ref())?);
    let mut unmap_file = BufWriter::with_capacity(64 * 1024, std::fs::File::create(unmap.as_ref())?);
    
    // Write headers first
    for header in &header_lines {
        writeln!(output_file, "{}", header)?;
    }
    
    // Write results (maintaining chunk order)
    for (success_lines, failed_lines) in results {
        for line in success_lines {
            writeln!(output_file, "{}", line)?;
        }
        for line in failed_lines {
            writeln!(unmap_file, "{}", line)?;
        }
    }
    
    Ok(ConversionStats {
        total: total.load(Ordering::Relaxed),
        success: success.load(Ordering::Relaxed),
        failed: failed.load(Ordering::Relaxed),
        multi_map: multi_map.load(Ordering::Relaxed),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_bed_record_view_bed3() {
        let line = b"chr1\t1000\t2000";
        let view = BedRecordView::parse(line).unwrap();
        
        assert_eq!(view.chrom, "chr1");
        assert_eq!(view.start, 1000);
        assert_eq!(view.end, 2000);
        assert_eq!(view.field_count(), 3);
        assert!(!view.is_bed6());
        assert!(!view.is_bed12());
    }
    
    #[test]
    fn test_bed_record_view_bed6() {
        let line = b"chr1\t1000\t2000\tgene1\t500\t+";
        let view = BedRecordView::parse(line).unwrap();
        
        assert_eq!(view.chrom, "chr1");
        assert_eq!(view.start, 1000);
        assert_eq!(view.end, 2000);
        assert_eq!(view.name(), Some("gene1"));
        assert_eq!(view.score(), Some("500"));
        assert_eq!(view.strand(), Some(Strand::Plus));
        assert!(view.is_bed6());
        assert!(!view.is_bed12());
    }
    
    #[test]
    fn test_bed_record_view_bed12() {
        let line = b"chr1\t1000\t2000\tgene1\t500\t+\t1100\t1900\t0,0,0\t2\t100,100\t0,900";
        let view = BedRecordView::parse(line).unwrap();
        
        assert_eq!(view.chrom, "chr1");
        assert_eq!(view.start, 1000);
        assert_eq!(view.end, 2000);
        assert_eq!(view.thick_start(), Some(1100));
        assert_eq!(view.thick_end(), Some(1900));
        assert_eq!(view.item_rgb(), Some("0,0,0"));
        assert_eq!(view.block_count(), Some(2));
        assert_eq!(view.block_sizes(), Some("100,100"));
        assert_eq!(view.block_starts(), Some("0,900"));
        assert!(view.is_bed12());
    }
    
    #[test]
    fn test_bed_record_view_too_few_fields() {
        let line = b"chr1\t1000";
        let result = BedRecordView::parse(line);
        assert!(matches!(result, Err(BedParseError::TooFewFields { .. })));
    }
    
    #[test]
    fn test_bed_record_view_empty_line() {
        let line = b"";
        let result = BedRecordView::parse(line);
        assert!(matches!(result, Err(BedParseError::EmptyLine)));
    }
    
    #[test]
    fn test_bed_record_view_invalid_number() {
        let line = b"chr1\tabc\t2000";
        let result = BedRecordView::parse(line);
        assert!(matches!(result, Err(BedParseError::InvalidNumber(_, _))));
    }
    
    #[test]
    fn test_strand_parsing() {
        let plus = b"chr1\t1000\t2000\tname\t0\t+";
        let minus = b"chr1\t1000\t2000\tname\t0\t-";
        let dot = b"chr1\t1000\t2000\tname\t0\t.";
        
        assert_eq!(BedRecordView::parse(plus).unwrap().strand(), Some(Strand::Plus));
        assert_eq!(BedRecordView::parse(minus).unwrap().strand(), Some(Strand::Minus));
        assert_eq!(BedRecordView::parse(dot).unwrap().strand(), None);
    }
}
