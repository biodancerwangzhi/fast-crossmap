//! GFF/GTF format adapter
//!
//! Handles GFF3 and GTF format conversion with zero-copy parsing.
//! GFF uses 1-based coordinates (unlike BED which is 0-based).
//!
//! **Validates: Requirements 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7**

use crate::core::{CoordinateMapper, Strand};
use memchr::memchr;
use rayon::prelude::*;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

/// GFF/GTF parse error
#[derive(Debug, Clone)]
pub enum GffParseError {
    EmptyLine,
    TooFewFields { expected: usize, found: usize },
    InvalidUtf8(&'static str),
    InvalidNumber(&'static str, String),
    InvalidStrand(String),
}

impl std::fmt::Display for GffParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GffParseError::EmptyLine => write!(f, "Empty line"),
            GffParseError::TooFewFields { expected, found } => {
                write!(f, "Too few fields: expected {}, found {}", expected, found)
            }
            GffParseError::InvalidUtf8(field) => write!(f, "Invalid UTF-8 in field: {}", field),
            GffParseError::InvalidNumber(field, value) => {
                write!(f, "Invalid number in field {}: {}", field, value)
            }
            GffParseError::InvalidStrand(s) => write!(f, "Invalid strand: {}", s),
        }
    }
}

impl std::error::Error for GffParseError {}

/// Zero-copy GFF/GTF record view for parsing
/// GFF format: seqname, source, feature, start, end, score, strand, frame, attributes
/// All coordinates are 1-based, closed interval [start, end]
pub struct GffRecordView<'a> {
    /// Original line bytes (kept for potential future use)
    #[allow(dead_code)]
    line: &'a [u8],
    /// Sequence name (chromosome)
    pub seqname: &'a str,
    /// Source field
    pub source: &'a str,
    /// Feature type
    pub feature: &'a str,
    /// Start position (1-based)
    pub start: u64,
    /// End position (1-based, inclusive)
    pub end: u64,
    /// Score field (as string, may be ".")
    pub score: &'a str,
    /// Strand
    pub strand: Option<Strand>,
    /// Strand character for output
    pub strand_char: &'a str,
    /// Frame field
    pub frame: &'a str,
    /// Attributes field
    pub attributes: &'a str,
}


impl<'a> GffRecordView<'a> {
    /// Parse a GFF/GTF line with minimal allocation
    /// GFF has exactly 9 tab-separated fields
    pub fn parse(line: &'a [u8]) -> Result<Self, GffParseError> {
        if line.is_empty() {
            return Err(GffParseError::EmptyLine);
        }

        // Find field boundaries using memchr for tab characters
        let mut field_bounds = Vec::with_capacity(9);
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
        
        // GFF requires exactly 9 fields
        if field_bounds.len() < 9 {
            return Err(GffParseError::TooFewFields {
                expected: 9,
                found: field_bounds.len(),
            });
        }
        
        // Helper to get field as str
        let get_field = |idx: usize, name: &'static str| -> Result<&'a str, GffParseError> {
            let (start, end) = field_bounds[idx];
            std::str::from_utf8(&line[start..end])
                .map_err(|_| GffParseError::InvalidUtf8(name))
        };
        
        // Parse all fields
        let seqname = get_field(0, "seqname")?;
        let source = get_field(1, "source")?;
        let feature = get_field(2, "feature")?;
        
        // Parse start (1-based)
        let start_str = get_field(3, "start")?;
        let start: u64 = start_str
            .parse()
            .map_err(|_| GffParseError::InvalidNumber("start", start_str.to_string()))?;
        
        // Parse end (1-based, inclusive)
        let end_str = get_field(4, "end")?;
        let end: u64 = end_str
            .parse()
            .map_err(|_| GffParseError::InvalidNumber("end", end_str.to_string()))?;
        
        let score = get_field(5, "score")?;
        let strand_char = get_field(6, "strand")?;
        let frame = get_field(7, "frame")?;
        let attributes = get_field(8, "attributes")?;
        
        // Parse strand
        let strand = match strand_char {
            "+" => Some(Strand::Plus),
            "-" => Some(Strand::Minus),
            "." => None,
            _ => return Err(GffParseError::InvalidStrand(strand_char.to_string())),
        };
        
        Ok(Self {
            line,
            seqname,
            source,
            feature,
            start,
            end,
            score,
            strand,
            strand_char,
            frame,
            attributes,
        })
    }
    
    /// Get the feature size (end - start + 1 for 1-based coordinates)
    pub fn size(&self) -> u64 {
        self.end - self.start + 1
    }
}


/// Conversion statistics
#[derive(Debug, Clone, Default)]
pub struct ConversionStats {
    pub total: usize,
    pub success: usize,
    pub failed: usize,
    pub comments: usize,
}

/// Convert a single GFF record
/// Returns None if conversion fails (unmapped, size changed, or multiple mappings)
fn convert_gff_record(
    view: &GffRecordView,
    mapper: &CoordinateMapper,
) -> Option<String> {
    // Get query strand (use Plus if unstranded)
    let query_strand = view.strand.unwrap_or(Strand::Plus);
    
    // Convert 1-based GFF coordinates to 0-based for mapping
    // GFF: [start, end] 1-based inclusive
    // Internal: [start, end) 0-based half-open
    let start_0based = view.start - 1;
    let end_0based = view.end; // end is exclusive in 0-based
    
    // Map coordinates
    let segments = mapper.map(view.seqname, start_0based, end_0based, query_strand)?;
    
    // GFF requires exact match: single segment, no size change
    if segments.is_empty() {
        return None;
    }
    
    // Multiple mappings = fail
    if segments.len() > 1 {
        return None;
    }
    
    let seg = &segments[0];
    
    // Check size preservation (exact match required)
    let original_size = view.size();
    let mapped_size = seg.target.end - seg.target.start;
    if mapped_size != original_size {
        return None;
    }
    
    // Convert back to 1-based coordinates for GFF output
    let new_start = seg.target.start + 1;
    let new_end = seg.target.end;
    
    // Determine output strand
    // CrossMap behavior: use the strand from the mapping result
    // fields[6] = a[1][3] in CrossMap's mapgff.py
    let output_strand = seg.target.strand.to_char();
    
    // Build output line
    Some(format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        seg.target.chrom,
        view.source,
        view.feature,
        new_start,
        new_end,
        view.score,
        output_strand,
        view.frame,
        view.attributes
    ))
}


/// Chunk size for parallel processing
const CHUNK_SIZE: usize = 10000;

/// Convert a GFF/GTF file
///
/// # Arguments
/// * `input` - Input GFF/GTF file path
/// * `output` - Output GFF/GTF file path
/// * `mapper` - Coordinate mapper
/// * `threads` - Number of threads (1 = sequential)
///
/// # Returns
/// Conversion statistics
pub fn convert_gff<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    threads: usize,
) -> Result<ConversionStats, std::io::Error> {
    let input_file = std::fs::File::open(input.as_ref())?;
    let reader = BufReader::with_capacity(128 * 1024, input_file);
    
    // Prepare output files with BufWriter for performance
    let output_path = output.as_ref();
    let unmap_path = output_path.with_extension("gff.unmap");
    
    let mut output_file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(output_path)?);
    let mut unmap_file = BufWriter::with_capacity(64 * 1024, std::fs::File::create(&unmap_path)?);
    
    // Atomic counters for parallel processing
    let total = AtomicUsize::new(0);
    let success = AtomicUsize::new(0);
    let failed = AtomicUsize::new(0);
    let comments = AtomicUsize::new(0);
    
    // Collect lines for processing
    let lines: Vec<String> = reader.lines().filter_map(|l| l.ok()).collect();
    
    if threads <= 1 {
        // Sequential processing
        for line in &lines {
            // Skip empty lines
            if line.is_empty() {
                continue;
            }
            
            // Pass through comment lines (starting with #)
            if line.starts_with('#') {
                writeln!(output_file, "{}", line)?;
                comments.fetch_add(1, Ordering::Relaxed);
                continue;
            }
            
            total.fetch_add(1, Ordering::Relaxed);
            
            // Parse and convert
            match GffRecordView::parse(line.as_bytes()) {
                Ok(view) => {
                    if let Some(converted) = convert_gff_record(&view, mapper) {
                        writeln!(output_file, "{}", converted)?;
                        success.fetch_add(1, Ordering::Relaxed);
                    } else {
                        writeln!(unmap_file, "{}", line)?;
                        failed.fetch_add(1, Ordering::Relaxed);
                    }
                }
                Err(_) => {
                    // Parse error - write to unmap
                    writeln!(unmap_file, "{}", line)?;
                    failed.fetch_add(1, Ordering::Relaxed);
                }
            }
        }
    } else {
        // Parallel processing
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                // First pass: write comments (must be sequential to preserve order)
                let mut data_lines: Vec<(usize, &String)> = Vec::new();
                
                for (idx, line) in lines.iter().enumerate() {
                    if line.is_empty() {
                        continue;
                    }
                    
                    if line.starts_with('#') {
                        writeln!(output_file, "{}", line).ok();
                        comments.fetch_add(1, Ordering::Relaxed);
                    } else {
                        data_lines.push((idx, line));
                    }
                }
                
                // Process data lines in parallel
                let results: Vec<(usize, Option<String>, &String)> = data_lines
                    .par_chunks(CHUNK_SIZE)
                    .flat_map(|chunk| {
                        chunk.iter().map(|(idx, line)| {
                            let result = GffRecordView::parse(line.as_bytes())
                                .ok()
                                .and_then(|view| convert_gff_record(&view, mapper));
                            (*idx, result, *line)
                        }).collect::<Vec<_>>()
                    })
                    .collect();
                
                // Write results (sequential to maintain order)
                for (_idx, result, original) in results {
                    total.fetch_add(1, Ordering::Relaxed);
                    match result {
                        Some(converted) => {
                            writeln!(output_file, "{}", converted).ok();
                            success.fetch_add(1, Ordering::Relaxed);
                        }
                        None => {
                            writeln!(unmap_file, "{}", original).ok();
                            failed.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                }
            });
    }
    
    Ok(ConversionStats {
        total: total.load(Ordering::Relaxed),
        success: success.load(Ordering::Relaxed),
        failed: failed.load(Ordering::Relaxed),
        comments: comments.load(Ordering::Relaxed),
    })
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gff_record_view_basic() {
        let line = b"chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tgene_id \"ENSG00000001\"";
        let view = GffRecordView::parse(line).unwrap();
        
        assert_eq!(view.seqname, "chr1");
        assert_eq!(view.source, "ensembl");
        assert_eq!(view.feature, "gene");
        assert_eq!(view.start, 1000);
        assert_eq!(view.end, 2000);
        assert_eq!(view.score, ".");
        assert_eq!(view.strand, Some(Strand::Plus));
        assert_eq!(view.strand_char, "+");
        assert_eq!(view.frame, ".");
        assert_eq!(view.attributes, "gene_id \"ENSG00000001\"");
        assert_eq!(view.size(), 1001);
    }

    #[test]
    fn test_gff_record_view_negative_strand() {
        let line = b"chr2\trefseq\texon\t5000\t5500\t100\t-\t0\tID=exon1";
        let view = GffRecordView::parse(line).unwrap();
        
        assert_eq!(view.seqname, "chr2");
        assert_eq!(view.strand, Some(Strand::Minus));
        assert_eq!(view.strand_char, "-");
        assert_eq!(view.score, "100");
        assert_eq!(view.frame, "0");
        assert_eq!(view.size(), 501);
    }

    #[test]
    fn test_gff_record_view_unstranded() {
        let line = b"chrX\t.\tregion\t100\t200\t.\t.\t.\t.";
        let view = GffRecordView::parse(line).unwrap();
        
        assert_eq!(view.strand, None);
        assert_eq!(view.strand_char, ".");
    }

    #[test]
    fn test_gff_record_view_too_few_fields() {
        let line = b"chr1\tensembl\tgene\t1000\t2000";
        let result = GffRecordView::parse(line);
        assert!(matches!(result, Err(GffParseError::TooFewFields { .. })));
    }

    #[test]
    fn test_gff_record_view_empty_line() {
        let line = b"";
        let result = GffRecordView::parse(line);
        assert!(matches!(result, Err(GffParseError::EmptyLine)));
    }

    #[test]
    fn test_gff_record_view_invalid_strand() {
        let line = b"chr1\t.\tgene\t1000\t2000\t.\tX\t.\t.";
        let result = GffRecordView::parse(line);
        assert!(matches!(result, Err(GffParseError::InvalidStrand(_))));
    }

    #[test]
    fn test_gff_record_view_gtf_format() {
        // GTF format with gene_id and transcript_id
        let line = b"chr1\thavana\ttranscript\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";";
        let view = GffRecordView::parse(line).unwrap();
        
        assert_eq!(view.seqname, "chr1");
        assert_eq!(view.source, "havana");
        assert_eq!(view.feature, "transcript");
        assert_eq!(view.start, 11869);
        assert_eq!(view.end, 14409);
        assert!(view.attributes.contains("gene_id"));
        assert!(view.attributes.contains("transcript_id"));
    }
}
