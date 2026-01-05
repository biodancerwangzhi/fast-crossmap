//! VCF format adapter
//!
//! Handles VCF format conversion with zero-copy parsing.
//!
//! **Validates: Requirements 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7**

use crate::core::{dna, CoordinateMapper, Strand};
use memchr::memchr;
use rayon::prelude::*;
use std::cell::{Cell, RefCell};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

/// VCF record representation for output
#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u64,
    pub id: String,
    pub ref_allele: String,
    pub alt_alleles: Vec<String>,
    pub qual: String,
    pub filter: String,
    pub info: String,
    pub format: Option<String>,
    pub samples: Vec<String>,
}

/// Zero-copy VCF record view for parsing
/// Only parses CHROM and POS immediately, other fields are kept as byte slices
pub struct VcfRecordView<'a> {
    /// Original line bytes
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

impl<'a> VcfRecordView<'a> {
    /// Parse a VCF line with minimal allocation
    /// Only parses CHROM and POS immediately
    pub fn parse(line: &'a [u8]) -> Result<Self, VcfParseError> {
        if line.is_empty() {
            return Err(VcfParseError::EmptyLine);
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
        
        // VCF requires at least 8 fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
        if field_bounds.len() < 8 {
            return Err(VcfParseError::TooFewFields {
                expected: 8,
                found: field_bounds.len(),
            });
        }
        
        // Parse CHROM (field 0)
        let chrom = std::str::from_utf8(&line[field_bounds[0].0..field_bounds[0].1])
            .map_err(|_| VcfParseError::InvalidUtf8("CHROM"))?;
        
        // Parse POS (field 1)
        let pos_str = std::str::from_utf8(&line[field_bounds[1].0..field_bounds[1].1])
            .map_err(|_| VcfParseError::InvalidUtf8("POS"))?;
        let pos: u64 = pos_str
            .parse()
            .map_err(|_| VcfParseError::InvalidNumber("POS", pos_str.to_string()))?;
        
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
    pub fn field(&self, index: usize) -> Option<&'a str> {
        self.field_bounds.get(index).and_then(|(start, end)| {
            std::str::from_utf8(&self.line[*start..*end]).ok()
        })
    }
    
    /// Get ID field (field 2)
    pub fn id(&self) -> Option<&'a str> {
        self.field(2)
    }
    
    /// Get REF field (field 3)
    pub fn ref_allele(&self) -> Option<&'a str> {
        self.field(3)
    }
    
    /// Get ALT field (field 4)
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
        (9..self.field_count())
            .filter_map(|i| self.field(i))
            .collect()
    }
    
    /// Parse INFO field lazily (only when needed)
    pub fn parse_info(&self) -> HashMap<String, String> {
        if !self.info_parsed.get() {
            let info_str = self.info().unwrap_or(".");
            let mut map = HashMap::new();
            
            if info_str != "." {
                for item in info_str.split(';') {
                    if let Some(eq_pos) = item.find('=') {
                        let key = item[..eq_pos].to_string();
                        let value = item[eq_pos + 1..].to_string();
                        map.insert(key, value);
                    } else {
                        // Flag without value
                        map.insert(item.to_string(), String::new());
                    }
                }
            }
            
            *self.info_cache.borrow_mut() = Some(map.clone());
            self.info_parsed.set(true);
            map
        } else {
            self.info_cache.borrow().clone().unwrap_or_default()
        }
    }
    
    /// Get variant type based on REF and ALT lengths
    pub fn variant_type(&self) -> VariantType {
        let ref_len = self.ref_allele().map(|s| s.len()).unwrap_or(0);
        let alt = self.alt_alleles().unwrap_or(".");
        
        // Get first ALT allele for type determination
        let first_alt = alt.split(',').next().unwrap_or(".");
        let alt_len = first_alt.len();
        
        if ref_len == alt_len {
            VariantType::Substitution
        } else if alt_len > ref_len {
            VariantType::Insertion
        } else {
            VariantType::Deletion
        }
    }
}

/// Variant type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    Substitution,
    Insertion,
    Deletion,
}

/// VCF parsing error
#[derive(Debug, thiserror::Error)]
pub enum VcfParseError {
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
}

/// Result of converting a single VCF record
#[derive(Debug)]
pub enum ConversionResult {
    /// Successfully mapped
    Success(String),
    /// Failed to map with reason
    Failed(String, String),
    /// Header line (pass through to output)
    Header(String),
    /// Header line (pass through to unmap)
    UnmapHeader(String),
    /// Contig header (needs special handling)
    ContigHeader(String),
}

/// Convert a single VCF record
fn convert_vcf_record(
    view: &VcfRecordView,
    mapper: &CoordinateMapper,
    ref_genome: Option<&pysam_stub::FastaReader>,
    no_comp_allele: bool,
) -> ConversionResult {
    // Map the first position of REF allele (VCF is 1-based)
    let start = view.pos - 1; // Convert to 0-based
    let end = start + 1; // Map only the first position
    
    let result = mapper.map(view.chrom, start, end, Strand::Plus);
    
    match result {
        Some(segments) if segments.len() == 1 => {
            let seg = &segments[0];
            let target_chrom = &seg.target.chrom;
            let target_start = seg.target.start;
            let target_strand = seg.target.strand; // Fixed: access strand from target MapResult
            
            // Get original fields
            let ref_allele = view.ref_allele().unwrap_or("N");
            let alt_alleles_str = view.alt_alleles().unwrap_or(".");
            let _ref_allele_size = ref_allele.len();
            
            // Determine variant type
            let _v_type = view.variant_type();
            
            // Calculate new REF position based on strand and variant type
            let (new_pos, new_ref) = if let Some(ref_reader) = ref_genome {
                // Get REF from target reference genome
                let ref_start = target_start;
                let ref_end = ref_start + 1;
                
                match ref_reader.fetch(target_chrom, ref_start, ref_end) {
                    Some(seq) => (target_start + 1, seq.to_uppercase()),
                    None => {
                        return ConversionResult::Failed(
                            reconstruct_line(view),
                            "Fail(KeyError)".to_string(),
                        );
                    }
                }
            } else {
                // No reference genome provided, keep original REF
                (target_start + 1, ref_allele.to_string())
            };
            
            if new_ref.is_empty() {
                return ConversionResult::Failed(
                    reconstruct_line(view),
                    "Fail(KeyError)".to_string(),
                );
            }
            
            // Process ALT alleles
            let mut alt_alleles_updated = Vec::new();
            for alt_allele in alt_alleles_str.split(',') {
                if dna::is_dna(alt_allele) {
                    let updated = if ref_allele.len() != alt_allele.len() {
                        // Indel: replace first nucleotide with new REF, handle rest
                        if target_strand == Strand::Minus {
                            // Reverse complement the rest (after first nucleotide)
                            let first_char = new_ref.chars().next().unwrap_or('N');
                            if alt_allele.len() > 1 {
                                format!("{}{}", first_char, dna::revcomp(&alt_allele[1..]))
                            } else {
                                first_char.to_string()
                            }
                        } else {
                            // Forward strand: replace first nucleotide only
                            let first_char = new_ref.chars().next().unwrap_or('N');
                            if alt_allele.len() > 1 {
                                format!("{}{}", first_char, &alt_allele[1..])
                            } else {
                                first_char.to_string()
                            }
                        }
                    } else {
                        // Substitution
                        if target_strand == Strand::Minus {
                            dna::revcomp(alt_allele)
                        } else {
                            alt_allele.to_string()
                        }
                    };
                    
                    // Only add if different from REF
                    if updated != new_ref {
                        alt_alleles_updated.push(updated);
                    }
                } else {
                    // Non-DNA allele (e.g., <DEL>, <INS>), keep as-is
                    alt_alleles_updated.push(alt_allele.to_string());
                }
            }
            
            // Check if all ALT alleles were filtered out
            if alt_alleles_updated.is_empty() {
                return ConversionResult::Failed(
                    reconstruct_line(view),
                    "Fail(REF==ALT)".to_string(),
                );
            }
            
            // Check REF == ALT (unless noCompAllele is set)
            if !no_comp_allele && alt_alleles_updated.len() == 1 && alt_alleles_updated[0] == new_ref {
                return ConversionResult::Failed(
                    reconstruct_line(view),
                    "Fail(REF==ALT)".to_string(),
                );
            }
            
            // Build output line
            let output = format_output_line(
                view,
                target_chrom,
                new_pos,
                &new_ref,
                &alt_alleles_updated,
            );
            
            ConversionResult::Success(output)
        }
        Some(segments) if segments.len() > 1 => {
            // Multiple mappings
            ConversionResult::Failed(
                reconstruct_line(view),
                "Fail(Multiple_hits)".to_string(),
            )
        }
        _ => {
            // No mapping found
            ConversionResult::Failed(
                reconstruct_line(view),
                "Fail(Unmap)".to_string(),
            )
        }
    }
}

/// Format output line for a successfully mapped VCF record
fn format_output_line(
    view: &VcfRecordView,
    chrom: &str,
    pos: u64,
    ref_allele: &str,
    alt_alleles: &[String],
) -> String {
    let mut output = String::with_capacity(512);
    
    // CHROM
    output.push_str(chrom);
    output.push('\t');
    
    // POS
    output.push_str(&pos.to_string());
    output.push('\t');
    
    // ID
    output.push_str(view.id().unwrap_or("."));
    output.push('\t');
    
    // REF
    output.push_str(ref_allele);
    output.push('\t');
    
    // ALT
    output.push_str(&alt_alleles.join(","));
    output.push('\t');
    
    // QUAL
    output.push_str(view.qual().unwrap_or("."));
    output.push('\t');
    
    // FILTER
    output.push_str(view.filter().unwrap_or("."));
    output.push('\t');
    
    // INFO - update END if present
    let info = view.info().unwrap_or(".");
    // TODO: Update END field if present
    output.push_str(info);
    
    // FORMAT and samples
    if let Some(format) = view.format() {
        output.push('\t');
        output.push_str(format);
        
        for sample in view.samples() {
            output.push('\t');
            output.push_str(sample);
        }
    }
    
    output
}

/// Reconstruct original line from view
fn reconstruct_line(view: &VcfRecordView) -> String {
    let mut output = String::with_capacity(512);
    
    for i in 0..view.field_count() {
        if i > 0 {
            output.push('\t');
        }
        if let Some(field) = view.field(i) {
            output.push_str(field);
        }
    }
    
    output
}


/// Stub module for FASTA reading (placeholder for pysam-like functionality)
pub mod pysam_stub {
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader};
    use std::path::Path;
    
    /// Simple FASTA reader for reference genome
    pub struct FastaReader {
        sequences: HashMap<String, Vec<u8>>,
    }
    
    impl FastaReader {
        /// Open a FASTA file
        pub fn open<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
            let file = std::fs::File::open(path)?;
            let reader = BufReader::new(file);
            let mut sequences = HashMap::new();
            let mut current_name = String::new();
            let mut current_seq = Vec::new();
            
            for line in reader.lines() {
                let line = line?;
                if line.starts_with('>') {
                    if !current_name.is_empty() {
                        sequences.insert(current_name.clone(), current_seq.clone());
                    }
                    current_name = line[1..].split_whitespace().next().unwrap_or("").to_string();
                    current_seq.clear();
                } else {
                    current_seq.extend(line.trim().bytes());
                }
            }
            
            if !current_name.is_empty() {
                sequences.insert(current_name, current_seq);
            }
            
            Ok(Self { sequences })
        }
        
        /// Fetch a region from the reference
        pub fn fetch(&self, chrom: &str, start: u64, end: u64) -> Option<String> {
            // Try with and without chr prefix
            let seq = self.sequences.get(chrom)
                .or_else(|| {
                    if chrom.starts_with("chr") {
                        self.sequences.get(&chrom[3..])
                    } else {
                        self.sequences.get(&format!("chr{}", chrom))
                    }
                })?;
            
            let start = start as usize;
            let end = (end as usize).min(seq.len());
            
            if start >= seq.len() {
                return None;
            }
            
            Some(String::from_utf8_lossy(&seq[start..end]).to_string())
        }
        
        /// Get chromosome names
        pub fn references(&self) -> Vec<&str> {
            self.sequences.keys().map(|s| s.as_str()).collect()
        }
        
        /// Get chromosome lengths
        pub fn lengths(&self) -> Vec<usize> {
            self.sequences.values().map(|s| s.len()).collect()
        }
    }
}

/// Chunk size for parallel processing
const CHUNK_SIZE: usize = 10000;

/// Convert a VCF file using the coordinate mapper
/// 
/// # Arguments
/// * `input` - Input VCF file path
/// * `output` - Output VCF file path for successfully mapped records
/// * `unmap` - Output file path for unmapped records (will be output.unmap)
/// * `mapper` - Coordinate mapper with loaded chain index
/// * `ref_genome` - Optional path to target reference genome FASTA
/// * `no_comp_allele` - If true, keep variants where REF==ALT
/// * `threads` - Number of threads for parallel processing (1 = sequential)
/// 
/// # Returns
/// Conversion statistics
pub fn convert_vcf<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    ref_genome: Option<P>,
    no_comp_allele: bool,
    threads: usize,
) -> Result<ConversionStats, VcfParseError> {
    if threads > 1 {
        convert_vcf_parallel(input, output, mapper, ref_genome, no_comp_allele, threads)
    } else {
        convert_vcf_sequential(input, output, mapper, ref_genome, no_comp_allele)
    }
}

/// Sequential VCF conversion (single-threaded)
fn convert_vcf_sequential<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    ref_genome: Option<P>,
    no_comp_allele: bool,
) -> Result<ConversionStats, VcfParseError> {
    let input_file = std::fs::File::open(input.as_ref())?;
    let reader = BufReader::with_capacity(128 * 1024, input_file);
    
    let output_path = output.as_ref();
    let unmap_path = output_path.with_extension("vcf.unmap");
    
    // Use BufWriter for performance
    let mut output_file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(output_path)?);
    let mut unmap_file = BufWriter::with_capacity(64 * 1024, std::fs::File::create(&unmap_path)?);
    
    // Load reference genome if provided
    let ref_reader = ref_genome
        .map(|p| pysam_stub::FastaReader::open(p.as_ref()))
        .transpose()?;
    
    let mut stats = ConversionStats::default();
    let mut line_buf = String::with_capacity(4096);
    let mut reader = reader;
    
    // Track if we've seen the #CHROM header
    let mut _seen_chrom_header = false;
    
    loop {
        line_buf.clear();
        let bytes_read = reader.read_line(&mut line_buf)?;
        if bytes_read == 0 {
            break;
        }
        
        let line = line_buf.trim_end();
        
        if line.is_empty() {
            continue;
        }
        
        // Handle header lines
        if line.starts_with('#') {
            if line.starts_with("##fileformat") 
                || line.starts_with("##INFO")
                || line.starts_with("##FILTER")
                || line.starts_with("##FORMAT")
                || line.starts_with("##ALT")
                || line.starts_with("##SAMPLE")
                || line.starts_with("##PEDIGREE")
            {
                // Write to both files
                writeln!(output_file, "{}", line)?;
                writeln!(unmap_file, "{}", line)?;
            } else if line.starts_with("##assembly") || line.starts_with("##contig") {
                // Write only to unmap file
                writeln!(unmap_file, "{}", line)?;
            } else if line.starts_with("#CHROM") {
                _seen_chrom_header = true;
                // Write contig headers for target assembly
                if let Some(ref reader) = ref_reader {
                    for (chrom, len) in reader.references().iter().zip(reader.lengths()) {
                        writeln!(output_file, "##contig=<ID={},length={}>", chrom, len)?;
                    }
                }
                // Write liftover metadata
                writeln!(output_file, "##liftOverProgram=FastCrossMap")?;
                // Write column header to both files
                writeln!(output_file, "{}", line)?;
                writeln!(unmap_file, "{}", line)?;
            } else {
                // Other header lines - write to output only
                writeln!(output_file, "{}", line)?;
            }
            continue;
        }
        
        stats.total += 1;
        
        // Parse the VCF record
        match VcfRecordView::parse(line.as_bytes()) {
            Ok(view) => {
                match convert_vcf_record(&view, mapper, ref_reader.as_ref(), no_comp_allele) {
                    ConversionResult::Success(output_line) => {
                        writeln!(output_file, "{}", output_line)?;
                        stats.success += 1;
                    }
                    ConversionResult::Failed(original, reason) => {
                        writeln!(unmap_file, "{}\t{}", original, reason)?;
                        stats.failed += 1;
                    }
                    ConversionResult::Header(h) => {
                        writeln!(output_file, "{}", h)?;
                    }
                    ConversionResult::UnmapHeader(h) => {
                        writeln!(unmap_file, "{}", h)?;
                    }
                    ConversionResult::ContigHeader(_) => {
                        // Already handled above
                    }
                }
            }
            Err(_) => {
                // Invalid VCF line - write to unmap file
                writeln!(unmap_file, "{}\tFail(ParseError)", line)?;
                stats.failed += 1;
            }
        }
    }
    
    Ok(stats)
}

/// Parallel VCF conversion using rayon
fn convert_vcf_parallel<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    ref_genome: Option<P>,
    no_comp_allele: bool,
    threads: usize,
) -> Result<ConversionStats, VcfParseError> {
    // Configure rayon thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| VcfParseError::Io(std::io::Error::new(
            std::io::ErrorKind::Other,
            format!("Failed to create thread pool: {}", e)
        )))?;
    
    // Read all lines
    let input_file = std::fs::File::open(input.as_ref())?;
    let reader = BufReader::with_capacity(128 * 1024, input_file);
    
    let mut header_lines_output = Vec::new();
    let mut header_lines_unmap = Vec::new();
    let mut data_lines = Vec::new();
    
    // Load reference genome if provided
    let ref_reader = ref_genome
        .map(|p| pysam_stub::FastaReader::open(p.as_ref()))
        .transpose()?;
    
    for line_result in reader.lines() {
        let line = line_result?;
        if line.is_empty() {
            continue;
        }
        
        if line.starts_with('#') {
            if line.starts_with("##fileformat") 
                || line.starts_with("##INFO")
                || line.starts_with("##FILTER")
                || line.starts_with("##FORMAT")
                || line.starts_with("##ALT")
                || line.starts_with("##SAMPLE")
                || line.starts_with("##PEDIGREE")
            {
                header_lines_output.push(line.clone());
                header_lines_unmap.push(line);
            } else if line.starts_with("##assembly") || line.starts_with("##contig") {
                header_lines_unmap.push(line);
            } else if line.starts_with("#CHROM") {
                // Add contig headers for target assembly
                if let Some(ref reader) = ref_reader {
                    for (chrom, len) in reader.references().iter().zip(reader.lengths()) {
                        header_lines_output.push(format!("##contig=<ID={},length={}>", chrom, len));
                    }
                }
                header_lines_output.push("##liftOverProgram=FastCrossMap".to_string());
                header_lines_output.push(line.clone());
                header_lines_unmap.push(line);
            } else {
                header_lines_output.push(line);
            }
        } else {
            data_lines.push(line);
        }
    }
    
    // Atomic counters for stats
    let total = AtomicUsize::new(0);
    let success = AtomicUsize::new(0);
    let failed = AtomicUsize::new(0);
    
    // Process in parallel
    let results: Vec<(Vec<String>, Vec<String>)> = pool.install(|| {
        data_lines
            .par_chunks(CHUNK_SIZE)
            .map(|chunk| {
                let mut success_lines = Vec::with_capacity(chunk.len());
                let mut failed_lines = Vec::new();
                
                for line in chunk {
                    total.fetch_add(1, Ordering::Relaxed);
                    
                    match VcfRecordView::parse(line.as_bytes()) {
                        Ok(view) => {
                            match convert_vcf_record(&view, mapper, ref_reader.as_ref(), no_comp_allele) {
                                ConversionResult::Success(output_line) => {
                                    success_lines.push(output_line);
                                    success.fetch_add(1, Ordering::Relaxed);
                                }
                                ConversionResult::Failed(original, reason) => {
                                    failed_lines.push(format!("{}\t{}", original, reason));
                                    failed.fetch_add(1, Ordering::Relaxed);
                                }
                                _ => {}
                            }
                        }
                        Err(_) => {
                            failed_lines.push(format!("{}\tFail(ParseError)", line));
                            failed.fetch_add(1, Ordering::Relaxed);
                        }
                    }
                }
                
                (success_lines, failed_lines)
            })
            .collect()
    });
    
    // Write output files with BufWriter for performance
    let output_path = output.as_ref();
    let unmap_path = output_path.with_extension("vcf.unmap");
    
    let mut output_file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(output_path)?);
    let mut unmap_file = BufWriter::with_capacity(64 * 1024, std::fs::File::create(&unmap_path)?);
    
    // Write headers
    for header in &header_lines_output {
        writeln!(output_file, "{}", header)?;
    }
    for header in &header_lines_unmap {
        writeln!(unmap_file, "{}", header)?;
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
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_vcf_record_view_basic() {
        let line = b"chr1\t12345\trs123\tA\tG\t30\tPASS\tDP=100";
        let view = VcfRecordView::parse(line).unwrap();
        
        assert_eq!(view.chrom, "chr1");
        assert_eq!(view.pos, 12345);
        assert_eq!(view.id(), Some("rs123"));
        assert_eq!(view.ref_allele(), Some("A"));
        assert_eq!(view.alt_alleles(), Some("G"));
        assert_eq!(view.qual(), Some("30"));
        assert_eq!(view.filter(), Some("PASS"));
        assert_eq!(view.info(), Some("DP=100"));
    }
    
    #[test]
    fn test_vcf_record_view_with_samples() {
        let line = b"chr1\t12345\t.\tA\tG\t.\t.\t.\tGT:DP\t0/1:30\t1/1:25";
        let view = VcfRecordView::parse(line).unwrap();
        
        assert_eq!(view.chrom, "chr1");
        assert_eq!(view.pos, 12345);
        assert_eq!(view.format(), Some("GT:DP"));
        assert_eq!(view.samples(), vec!["0/1:30", "1/1:25"]);
    }
    
    #[test]
    fn test_vcf_record_view_too_few_fields() {
        let line = b"chr1\t12345\trs123";
        let result = VcfRecordView::parse(line);
        assert!(matches!(result, Err(VcfParseError::TooFewFields { .. })));
    }
    
    #[test]
    fn test_vcf_record_view_empty_line() {
        let line = b"";
        let result = VcfRecordView::parse(line);
        assert!(matches!(result, Err(VcfParseError::EmptyLine)));
    }
    
    #[test]
    fn test_variant_type_detection() {
        // Substitution
        let line = b"chr1\t100\t.\tA\tG\t.\t.\t.";
        let view = VcfRecordView::parse(line).unwrap();
        assert_eq!(view.variant_type(), VariantType::Substitution);
        
        // Insertion
        let line = b"chr1\t100\t.\tA\tAG\t.\t.\t.";
        let view = VcfRecordView::parse(line).unwrap();
        assert_eq!(view.variant_type(), VariantType::Insertion);
        
        // Deletion
        let line = b"chr1\t100\t.\tAG\tA\t.\t.\t.";
        let view = VcfRecordView::parse(line).unwrap();
        assert_eq!(view.variant_type(), VariantType::Deletion);
    }
    
    #[test]
    fn test_info_parsing() {
        let line = b"chr1\t100\t.\tA\tG\t.\t.\tDP=100;AF=0.5;DB";
        let view = VcfRecordView::parse(line).unwrap();
        let info = view.parse_info();
        
        assert_eq!(info.get("DP"), Some(&"100".to_string()));
        assert_eq!(info.get("AF"), Some(&"0.5".to_string()));
        assert_eq!(info.get("DB"), Some(&"".to_string())); // Flag
    }
    
    #[test]
    fn test_multi_allelic() {
        let line = b"chr1\t100\t.\tA\tG,T,C\t.\t.\t.";
        let view = VcfRecordView::parse(line).unwrap();
        
        assert_eq!(view.alt_alleles(), Some("G,T,C"));
    }
}
