//! GVCF format adapter
//!
//! Handles GVCF (Genomic VCF) format conversion with support for non-variant regions.
//! GVCF extends VCF with END= INFO field for non-variant blocks.
//!
//! **Validates: Requirements 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7**

use crate::core::{dna, CoordinateMapper, Strand};
use memchr::memchr;
use std::cell::{Cell, RefCell};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

/// GVCF parsing error
#[derive(Debug, Clone)]
pub enum GvcfParseError {
    EmptyLine,
    TooFewFields { expected: usize, found: usize },
    InvalidUtf8(&'static str),
    InvalidNumber(&'static str, String),
}

impl std::fmt::Display for GvcfParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GvcfParseError::EmptyLine => write!(f, "Empty line"),
            GvcfParseError::TooFewFields { expected, found } => {
                write!(f, "Too few fields: expected {}, found {}", expected, found)
            }
            GvcfParseError::InvalidUtf8(field) => write!(f, "Invalid UTF-8 in field: {}", field),
            GvcfParseError::InvalidNumber(field, value) => {
                write!(f, "Invalid number in field {}: {}", field, value)
            }
        }
    }
}

impl std::error::Error for GvcfParseError {}

/// Zero-copy GVCF record view for parsing
/// Extends VCF parsing with END= support for non-variant regions
pub struct GvcfRecordView<'a> {
    /// Original line bytes
    #[allow(dead_code)]
    line: &'a [u8],
    /// Chromosome name
    pub chrom: &'a str,
    /// Position (1-based)
    pub pos: u64,
    /// Field boundaries (start, end) for lazy access
    field_bounds: Vec<(usize, usize)>,
    /// Cached INFO parsing
    info_parsed: Cell<bool>,
    info_cache: RefCell<Option<HashMap<String, String>>>,
}


impl<'a> GvcfRecordView<'a> {
    /// Parse a GVCF line with minimal allocation
    pub fn parse(line: &'a [u8]) -> Result<Self, GvcfParseError> {
        if line.is_empty() {
            return Err(GvcfParseError::EmptyLine);
        }

        // Find field boundaries using memchr for tab characters
        let mut field_bounds = Vec::with_capacity(10);
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
        
        // VCF requires at least 8 fields
        if field_bounds.len() < 8 {
            return Err(GvcfParseError::TooFewFields {
                expected: 8,
                found: field_bounds.len(),
            });
        }
        
        // Parse CHROM (field 0)
        let chrom = std::str::from_utf8(&line[field_bounds[0].0..field_bounds[0].1])
            .map_err(|_| GvcfParseError::InvalidUtf8("CHROM"))?;
        
        // Parse POS (field 1)
        let pos_str = std::str::from_utf8(&line[field_bounds[1].0..field_bounds[1].1])
            .map_err(|_| GvcfParseError::InvalidUtf8("POS"))?;
        let pos: u64 = pos_str
            .parse()
            .map_err(|_| GvcfParseError::InvalidNumber("POS", pos_str.to_string()))?;
        
        Ok(Self {
            line,
            chrom,
            pos,
            field_bounds,
            info_parsed: Cell::new(false),
            info_cache: RefCell::new(None),
        })
    }
    
    /// Get the number of fields
    pub fn field_count(&self) -> usize {
        self.field_bounds.len()
    }
    
    /// Get field as string slice (lazy access)
    fn field(&self, index: usize) -> Option<&'a str> {
        self.field_bounds.get(index).and_then(|(start, end)| {
            std::str::from_utf8(&self.line[*start..*end]).ok()
        })
    }
    
    /// Get ID field (field 2)
    pub fn id(&self) -> Option<&'a str> {
        self.field(2)
    }
    
    /// Get REF allele (field 3)
    pub fn ref_allele(&self) -> Option<&'a str> {
        self.field(3)
    }
    
    /// Get ALT alleles (field 4)
    pub fn alt_alleles(&self) -> Option<&'a str> {
        self.field(4)
    }
    
    /// Get QUAL field (field 5)
    pub fn qual(&self) -> Option<&'a str> {
        self.field(5)
    }
    
    /// Get FILTER field (field 6)
    pub fn filter(&self) -> Option<&'a str> {
        self.field(6)
    }
    
    /// Get INFO field (field 7)
    pub fn info(&self) -> Option<&'a str> {
        self.field(7)
    }
    
    /// Get FORMAT field (field 8) if present
    pub fn format(&self) -> Option<&'a str> {
        self.field(8)
    }
    
    /// Get sample fields (fields 9+)
    pub fn samples(&self) -> Vec<&'a str> {
        (9..self.field_bounds.len())
            .filter_map(|i| self.field(i))
            .collect()
    }
    
    /// Parse INFO field into key-value pairs (cached)
    pub fn parse_info(&self) -> HashMap<String, String> {
        if self.info_parsed.get() {
            return self.info_cache.borrow().clone().unwrap_or_default();
        }
        
        let mut info_map = HashMap::new();
        if let Some(info_str) = self.info() {
            if info_str != "." {
                for item in info_str.split(';') {
                    if let Some(eq_pos) = item.find('=') {
                        let key = &item[..eq_pos];
                        let value = &item[eq_pos + 1..];
                        info_map.insert(key.to_string(), value.to_string());
                    } else {
                        // Flag (no value)
                        info_map.insert(item.to_string(), String::new());
                    }
                }
            }
        }
        
        self.info_parsed.set(true);
        *self.info_cache.borrow_mut() = Some(info_map.clone());
        info_map
    }
    
    /// Get END position from INFO field (for non-variant blocks)
    /// Returns None if END is not present
    pub fn end_position(&self) -> Option<u64> {
        let info = self.parse_info();
        info.get("END").and_then(|v| v.parse().ok())
    }
    
    /// Check if this is a non-variant block (has END= in INFO)
    pub fn is_non_variant_block(&self) -> bool {
        self.end_position().is_some()
    }
    
    /// Check if ALT is <NON_REF> or <*> (GVCF non-variant marker)
    pub fn is_gvcf_non_ref(&self) -> bool {
        if let Some(alt) = self.alt_alleles() {
            alt == "<NON_REF>" || alt == "<*>" || alt == "."
        } else {
            false
        }
    }
}


/// Stub for FASTA reader (reference genome access)
/// In production, this would use rust-htslib or similar
pub mod fasta_stub {
    use std::path::Path;
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader};
    
    /// Simple FASTA reader for reference genome
    pub struct FastaReader {
        /// Chromosome sequences (loaded on demand)
        sequences: HashMap<String, Vec<u8>>,
        /// Path to FASTA file
        path: std::path::PathBuf,
    }
    
    impl FastaReader {
        /// Open a FASTA file
        pub fn open<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
            Ok(Self {
                sequences: HashMap::new(),
                path: path.as_ref().to_path_buf(),
            })
        }
        
        /// Load a chromosome sequence
        fn load_chrom(&mut self, chrom: &str) -> std::io::Result<()> {
            if self.sequences.contains_key(chrom) {
                return Ok(());
            }
            
            let file = std::fs::File::open(&self.path)?;
            let reader = BufReader::new(file);
            
            let mut current_chrom: Option<String> = None;
            let mut current_seq = Vec::new();
            let mut found = false;
            
            for line in reader.lines() {
                let line = line?;
                if line.starts_with('>') {
                    // Save previous chromosome if it matches
                    if let Some(ref c) = current_chrom {
                        if c == chrom {
                            found = true;
                            break;
                        }
                    }
                    // Parse new chromosome name
                    let name = line[1..].split_whitespace().next().unwrap_or("");
                    // Handle both "chr1" and "1" formats
                    let normalized = if name.starts_with("chr") {
                        name.to_string()
                    } else {
                        format!("chr{}", name)
                    };
                    
                    if normalized == chrom || name == chrom {
                        current_chrom = Some(chrom.to_string());
                        current_seq.clear();
                    } else {
                        current_chrom = None;
                    }
                } else if current_chrom.is_some() {
                    current_seq.extend(line.as_bytes().iter().filter(|b| !b.is_ascii_whitespace()));
                }
            }
            
            if found || current_chrom.is_some() {
                self.sequences.insert(chrom.to_string(), current_seq);
            }
            
            Ok(())
        }
        
        /// Fetch sequence at given position (0-based, half-open)
        pub fn fetch(&mut self, chrom: &str, start: u64, end: u64) -> Option<String> {
            // Try to load chromosome if not already loaded
            if !self.sequences.contains_key(chrom) {
                self.load_chrom(chrom).ok()?;
            }
            
            let seq = self.sequences.get(chrom)?;
            let start = start as usize;
            let end = end as usize;
            
            if start >= seq.len() || end > seq.len() {
                return None;
            }
            
            String::from_utf8(seq[start..end].to_vec()).ok()
        }
    }
}


/// Conversion statistics
#[derive(Debug, Clone, Default)]
pub struct ConversionStats {
    pub total: usize,
    pub success: usize,
    pub failed: usize,
    pub headers: usize,
}

/// Result of converting a single GVCF record
#[allow(dead_code)]
enum ConversionResult {
    /// Successfully mapped
    Success(String),
    /// Failed to map
    Failed(String, String),
    /// Header line (pass through, reserved for future use)
    Header(String),
}

/// Reconstruct a GVCF line from view
fn reconstruct_line(view: &GvcfRecordView) -> String {
    let mut parts = Vec::new();
    for i in 0..view.field_count() {
        if let Some(field) = view.field(i) {
            parts.push(field.to_string());
        }
    }
    parts.join("\t")
}

/// Update INFO field with new END value
fn update_info_end(info: &str, new_end: u64) -> String {
    let parts: Vec<&str> = info.split(';').collect();
    
    let mut result = Vec::new();
    let mut found = false;
    for part in &parts {
        if part.starts_with("END=") {
            result.push(format!("END={}", new_end));
            found = true;
        } else {
            result.push(part.to_string());
        }
    }
    
    if !found {
        result.push(format!("END={}", new_end));
    }
    
    result.join(";")
}

/// Convert a single GVCF record
fn convert_gvcf_record(
    view: &GvcfRecordView,
    mapper: &CoordinateMapper,
    ref_genome: Option<&mut fasta_stub::FastaReader>,
    no_comp_allele: bool,
) -> ConversionResult {
    // Check if this is a non-variant block (has END=)
    let is_block = view.is_non_variant_block();
    let end_pos = view.end_position();
    
    // Map the region
    let start = view.pos - 1; // Convert to 0-based
    let end = if let Some(e) = end_pos {
        e // END is already 1-based, use as-is for end (exclusive in 0-based)
    } else {
        // For variant records, map just the first position
        start + 1
    };
    
    let result = mapper.map(view.chrom, start, end, Strand::Plus);
    
    match result {
        Some(segments) if segments.len() == 1 => {
            let seg = &segments[0];
            let target_chrom = &seg.target.chrom;
            let target_start = seg.target.start;
            let target_end = seg.target.end;
            let target_strand = seg.target.strand;
            
            // Get original fields
            let ref_allele = view.ref_allele().unwrap_or("N");
            let alt_alleles_str = view.alt_alleles().unwrap_or(".");
            
            // Calculate new position (1-based)
            let new_pos = target_start + 1;
            
            // Get new REF from target reference if available
            let new_ref = if let Some(ref_reader) = ref_genome {
                match ref_reader.fetch(target_chrom, target_start, target_start + 1) {
                    Some(seq) => seq.to_uppercase(),
                    None => ref_allele.to_string(),
                }
            } else {
                ref_allele.to_string()
            };
            
            // Process ALT alleles
            let new_alt = if view.is_gvcf_non_ref() {
                // Keep <NON_REF> or <*> as-is
                alt_alleles_str.to_string()
            } else if dna::is_dna(alt_alleles_str) {
                // Process DNA alleles
                let mut alt_parts = Vec::new();
                for alt in alt_alleles_str.split(',') {
                    if dna::is_dna(alt) {
                        let updated = if target_strand == Strand::Minus {
                            dna::revcomp(alt)
                        } else {
                            alt.to_string()
                        };
                        
                        // Check REF == ALT filter
                        if !no_comp_allele && updated == new_ref {
                            // Skip this allele
                            continue;
                        }
                        alt_parts.push(updated);
                    } else {
                        alt_parts.push(alt.to_string());
                    }
                }
                
                if alt_parts.is_empty() {
                    // All alleles filtered out
                    return ConversionResult::Failed(
                        reconstruct_line(view),
                        "Fail(REF==ALT)".to_string(),
                    );
                }
                alt_parts.join(",")
            } else {
                alt_alleles_str.to_string()
            };
            
            // Update INFO field for non-variant blocks
            let new_info = if is_block {
                // Update END= to new target end position (1-based)
                let new_end = target_end; // target_end is already the correct 1-based end
                update_info_end(view.info().unwrap_or("."), new_end)
            } else {
                view.info().unwrap_or(".").to_string()
            };
            
            // Build output line
            let mut output = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                target_chrom,
                new_pos,
                view.id().unwrap_or("."),
                new_ref,
                new_alt,
                view.qual().unwrap_or("."),
                view.filter().unwrap_or("."),
                new_info
            );
            
            // Add FORMAT and samples if present
            if let Some(format) = view.format() {
                output.push('\t');
                output.push_str(format);
                for sample in view.samples() {
                    output.push('\t');
                    output.push_str(sample);
                }
            }
            
            ConversionResult::Success(output)
        }
        Some(segments) if segments.is_empty() => {
            ConversionResult::Failed(reconstruct_line(view), "Fail(Unmapped)".to_string())
        }
        Some(_) => {
            // Multiple mappings
            ConversionResult::Failed(reconstruct_line(view), "Fail(Multiple)".to_string())
        }
        None => {
            ConversionResult::Failed(reconstruct_line(view), "Fail(Unmapped)".to_string())
        }
    }
}


/// Update contig header with target assembly information
fn update_contig_header(line: &str, mapper: &CoordinateMapper) -> String {
    // Parse contig header: ##contig=<ID=chr1,length=248956422>
    if !line.starts_with("##contig=") {
        return line.to_string();
    }
    
    // Extract ID from header
    let id_start = line.find("ID=").map(|i| i + 3);
    let id_end = id_start.and_then(|s| {
        line[s..].find(',').or_else(|| line[s..].find('>')).map(|e| s + e)
    });
    
    if let (Some(start), Some(end)) = (id_start, id_end) {
        let chrom = &line[start..end];
        
        // Get target size for this chromosome
        if let Some(size) = mapper.index().target_chrom_size(chrom) {
            return format!("##contig=<ID={},length={}>", chrom, size);
        }
    }
    
    line.to_string()
}

/// Chunk size for parallel processing (reserved for future use)
#[allow(dead_code)]
const CHUNK_SIZE: usize = 10000;

/// Convert a GVCF file
///
/// # Arguments
/// * `input` - Input GVCF file path
/// * `output` - Output GVCF file path
/// * `mapper` - Coordinate mapper
/// * `ref_genome` - Optional path to target reference genome (FASTA)
/// * `no_comp_allele` - If true, don't filter REF==ALT
/// * `_threads` - Number of threads (reserved for future parallel processing)
///
/// # Returns
/// Conversion statistics
pub fn convert_gvcf<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    ref_genome: Option<P>,
    no_comp_allele: bool,
    _threads: usize,
) -> Result<ConversionStats, std::io::Error> {
    let input_file = std::fs::File::open(input.as_ref())?;
    let reader = BufReader::with_capacity(128 * 1024, input_file);
    
    // Prepare output files with BufWriter for performance
    let output_path = output.as_ref();
    let unmap_path = output_path.with_extension("gvcf.unmap");
    
    let mut output_file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(output_path)?);
    let mut unmap_file = BufWriter::with_capacity(64 * 1024, std::fs::File::create(&unmap_path)?);
    
    // Open reference genome if provided
    let mut ref_reader = ref_genome
        .map(|p| fasta_stub::FastaReader::open(p.as_ref()))
        .transpose()?;
    
    // Atomic counters
    let total = AtomicUsize::new(0);
    let success = AtomicUsize::new(0);
    let failed = AtomicUsize::new(0);
    let headers = AtomicUsize::new(0);
    
    // Collect lines
    let lines: Vec<String> = reader.lines().filter_map(|l| l.ok()).collect();
    
    // Process sequentially (GVCF often needs reference genome access which isn't thread-safe)
    for line in &lines {
        if line.is_empty() {
            continue;
        }
        
        // Handle header lines
        if line.starts_with('#') {
            if line.starts_with("##contig=") {
                // Update contig headers
                let updated = update_contig_header(line, mapper);
                writeln!(output_file, "{}", updated)?;
            } else {
                writeln!(output_file, "{}", line)?;
            }
            headers.fetch_add(1, Ordering::Relaxed);
            continue;
        }
        
        total.fetch_add(1, Ordering::Relaxed);
        
        // Parse and convert
        match GvcfRecordView::parse(line.as_bytes()) {
            Ok(view) => {
                let result = convert_gvcf_record(&view, mapper, ref_reader.as_mut(), no_comp_allele);
                match result {
                    ConversionResult::Success(converted) => {
                        writeln!(output_file, "{}", converted)?;
                        success.fetch_add(1, Ordering::Relaxed);
                    }
                    ConversionResult::Failed(original, _reason) => {
                        writeln!(unmap_file, "{}", original)?;
                        failed.fetch_add(1, Ordering::Relaxed);
                    }
                    ConversionResult::Header(h) => {
                        writeln!(output_file, "{}", h)?;
                        headers.fetch_add(1, Ordering::Relaxed);
                    }
                }
            }
            Err(_) => {
                writeln!(unmap_file, "{}", line)?;
                failed.fetch_add(1, Ordering::Relaxed);
            }
        }
    }
    
    Ok(ConversionStats {
        total: total.load(Ordering::Relaxed),
        success: success.load(Ordering::Relaxed),
        failed: failed.load(Ordering::Relaxed),
        headers: headers.load(Ordering::Relaxed),
    })
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gvcf_record_view_basic() {
        let line = b"chr1\t100\t.\tA\tG\t30\tPASS\tDP=100";
        let view = GvcfRecordView::parse(line).unwrap();
        
        assert_eq!(view.chrom, "chr1");
        assert_eq!(view.pos, 100);
        assert_eq!(view.ref_allele(), Some("A"));
        assert_eq!(view.alt_alleles(), Some("G"));
        assert!(!view.is_non_variant_block());
        assert!(!view.is_gvcf_non_ref());
    }

    #[test]
    fn test_gvcf_record_view_non_variant_block() {
        let line = b"chr1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=200;DP=50";
        let view = GvcfRecordView::parse(line).unwrap();
        
        assert_eq!(view.chrom, "chr1");
        assert_eq!(view.pos, 100);
        assert!(view.is_non_variant_block());
        assert_eq!(view.end_position(), Some(200));
        assert!(view.is_gvcf_non_ref());
    }

    #[test]
    fn test_gvcf_record_view_star_allele() {
        let line = b"chr1\t100\t.\tA\t<*>\t.\t.\tEND=150";
        let view = GvcfRecordView::parse(line).unwrap();
        
        assert!(view.is_gvcf_non_ref());
        assert!(view.is_non_variant_block());
        assert_eq!(view.end_position(), Some(150));
    }

    #[test]
    fn test_gvcf_record_view_with_samples() {
        let line = b"chr1\t100\t.\tA\tG\t30\tPASS\tDP=100\tGT:DP\t0/1:30";
        let view = GvcfRecordView::parse(line).unwrap();
        
        assert_eq!(view.format(), Some("GT:DP"));
        assert_eq!(view.samples(), vec!["0/1:30"]);
    }

    #[test]
    fn test_gvcf_info_parsing() {
        let line = b"chr1\t100\t.\tA\t<NON_REF>\t.\t.\tEND=200;DP=50;MQ=60";
        let view = GvcfRecordView::parse(line).unwrap();
        let info = view.parse_info();
        
        assert_eq!(info.get("END"), Some(&"200".to_string()));
        assert_eq!(info.get("DP"), Some(&"50".to_string()));
        assert_eq!(info.get("MQ"), Some(&"60".to_string()));
    }

    #[test]
    fn test_update_info_end() {
        let info = "END=100;DP=50";
        let updated = update_info_end(info, 200);
        assert!(updated.contains("END=200"));
        assert!(updated.contains("DP=50"));
        
        let info2 = "DP=50";
        let updated2 = update_info_end(info2, 300);
        assert!(updated2.contains("END=300"));
        assert!(updated2.contains("DP=50"));
    }

    #[test]
    fn test_gvcf_record_view_empty_line() {
        let line = b"";
        let result = GvcfRecordView::parse(line);
        assert!(matches!(result, Err(GvcfParseError::EmptyLine)));
    }

    #[test]
    fn test_gvcf_record_view_too_few_fields() {
        let line = b"chr1\t100\t.\tA";
        let result = GvcfRecordView::parse(line);
        assert!(matches!(result, Err(GvcfParseError::TooFewFields { .. })));
    }
}
