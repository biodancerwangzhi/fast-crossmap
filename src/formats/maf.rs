//! MAF (Mutation Annotation Format) adapter
//!
//! Handles MAF format conversion for mutation annotation data.
//! MAF is a tab-delimited format used by TCGA and other cancer genomics projects.
//!
//! **Validates: Requirements 8.1, 8.2, 8.3, 8.4, 8.5, 8.6**

use crate::core::{dna, CoordinateMapper, Strand};
use memchr::memchr;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

/// MAF parsing error
#[derive(Debug, Clone)]
pub enum MafParseError {
    EmptyLine,
    TooFewFields { expected: usize, found: usize },
    InvalidUtf8(&'static str),
    InvalidNumber(&'static str, String),
    MissingColumn(String),
}

impl std::fmt::Display for MafParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MafParseError::EmptyLine => write!(f, "Empty line"),
            MafParseError::TooFewFields { expected, found } => {
                write!(f, "Too few fields: expected {}, found {}", expected, found)
            }
            MafParseError::InvalidUtf8(field) => write!(f, "Invalid UTF-8 in field: {}", field),
            MafParseError::InvalidNumber(field, value) => {
                write!(f, "Invalid number in field {}: {}", field, value)
            }
            MafParseError::MissingColumn(col) => write!(f, "Missing required column: {}", col),
        }
    }
}

impl std::error::Error for MafParseError {}

/// Standard MAF column names
pub const MAF_COLUMNS: &[&str] = &[
    "Hugo_Symbol",
    "Entrez_Gene_Id",
    "Center",
    "NCBI_Build",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "dbSNP_RS",
    "dbSNP_Val_Status",
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
];

/// Column indices for key fields
#[derive(Debug, Clone)]
pub struct MafColumnIndices {
    pub chromosome: usize,
    pub start_position: usize,
    pub end_position: usize,
    pub strand: usize,
    pub reference_allele: usize,
    pub ncbi_build: usize,
    pub hugo_symbol: Option<usize>,
}


impl MafColumnIndices {
    /// Parse column indices from header line
    pub fn from_header(header: &str) -> Result<Self, MafParseError> {
        let columns: Vec<&str> = header.split('\t').collect();
        
        let find_col = |name: &str| -> Result<usize, MafParseError> {
            columns.iter()
                .position(|&c| c == name)
                .ok_or_else(|| MafParseError::MissingColumn(name.to_string()))
        };
        
        Ok(Self {
            chromosome: find_col("Chromosome")?,
            start_position: find_col("Start_Position")?,
            end_position: find_col("End_Position")?,
            strand: find_col("Strand")?,
            reference_allele: find_col("Reference_Allele")?,
            ncbi_build: find_col("NCBI_Build")?,
            hugo_symbol: columns.iter().position(|&c| c == "Hugo_Symbol"),
        })
    }
}

/// Zero-copy MAF record view for parsing
pub struct MafRecordView<'a> {
    /// Original line bytes
    #[allow(dead_code)]
    line: &'a [u8],
    /// Field values
    fields: Vec<&'a str>,
    /// Column indices
    indices: &'a MafColumnIndices,
}

impl<'a> MafRecordView<'a> {
    /// Parse a MAF line with column indices
    pub fn parse(line: &'a [u8], indices: &'a MafColumnIndices) -> Result<Self, MafParseError> {
        if line.is_empty() {
            return Err(MafParseError::EmptyLine);
        }

        // Find field boundaries using memchr for tab characters
        let mut fields = Vec::with_capacity(20);
        let mut start_pos = 0;
        let mut pos = 0;
        
        while pos < line.len() {
            if let Some(tab_pos) = memchr(b'\t', &line[pos..]) {
                let end_pos = pos + tab_pos;
                let field = std::str::from_utf8(&line[start_pos..end_pos])
                    .map_err(|_| MafParseError::InvalidUtf8("field"))?;
                fields.push(field);
                start_pos = end_pos + 1;
                pos = start_pos;
            } else {
                // Last field
                let field = std::str::from_utf8(&line[start_pos..])
                    .map_err(|_| MafParseError::InvalidUtf8("field"))?;
                fields.push(field);
                break;
            }
        }
        
        // Verify we have enough fields
        let max_idx = [
            indices.chromosome,
            indices.start_position,
            indices.end_position,
            indices.strand,
            indices.reference_allele,
            indices.ncbi_build,
        ].into_iter().max().unwrap_or(0);
        
        if fields.len() <= max_idx {
            return Err(MafParseError::TooFewFields {
                expected: max_idx + 1,
                found: fields.len(),
            });
        }
        
        Ok(Self {
            line,
            fields,
            indices,
        })
    }
    
    /// Get chromosome
    pub fn chromosome(&self) -> &'a str {
        self.fields[self.indices.chromosome]
    }
    
    /// Get start position (1-based)
    pub fn start_position(&self) -> Result<u64, MafParseError> {
        self.fields[self.indices.start_position]
            .parse()
            .map_err(|_| MafParseError::InvalidNumber(
                "Start_Position",
                self.fields[self.indices.start_position].to_string()
            ))
    }
    
    /// Get end position (1-based)
    pub fn end_position(&self) -> Result<u64, MafParseError> {
        self.fields[self.indices.end_position]
            .parse()
            .map_err(|_| MafParseError::InvalidNumber(
                "End_Position",
                self.fields[self.indices.end_position].to_string()
            ))
    }
    
    /// Get strand
    pub fn strand(&self) -> Option<Strand> {
        match self.fields[self.indices.strand] {
            "+" => Some(Strand::Plus),
            "-" => Some(Strand::Minus),
            _ => None,
        }
    }
    
    /// Get reference allele
    pub fn reference_allele(&self) -> &'a str {
        self.fields[self.indices.reference_allele]
    }
    
    /// Get NCBI build
    pub fn ncbi_build(&self) -> &'a str {
        self.fields[self.indices.ncbi_build]
    }
    
    /// Get Hugo symbol if present
    pub fn hugo_symbol(&self) -> Option<&'a str> {
        self.indices.hugo_symbol.map(|i| self.fields[i])
    }
    
    /// Get all fields
    pub fn fields(&self) -> &[&'a str] {
        &self.fields
    }
    
    /// Get field count
    pub fn field_count(&self) -> usize {
        self.fields.len()
    }
}


/// Stub for FASTA reader (reference genome access)
pub mod fasta_stub {
    use std::path::Path;
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader};
    
    /// Simple FASTA reader for reference genome
    pub struct FastaReader {
        sequences: HashMap<String, Vec<u8>>,
        path: std::path::PathBuf,
    }
    
    impl FastaReader {
        pub fn open<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
            Ok(Self {
                sequences: HashMap::new(),
                path: path.as_ref().to_path_buf(),
            })
        }
        
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
                    if let Some(ref c) = current_chrom {
                        if c == chrom {
                            found = true;
                            break;
                        }
                    }
                    let name = line[1..].split_whitespace().next().unwrap_or("");
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
        
        pub fn fetch(&mut self, chrom: &str, start: u64, end: u64) -> Option<String> {
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


/// Convert a single MAF record
fn convert_maf_record(
    view: &MafRecordView,
    mapper: &CoordinateMapper,
    ref_genome: Option<&mut fasta_stub::FastaReader>,
    target_build: &str,
) -> Option<String> {
    // Get coordinates (MAF uses 1-based coordinates)
    let start = view.start_position().ok()?;
    let end = view.end_position().ok()?;
    let chrom = view.chromosome();
    
    // Get query strand
    let query_strand = view.strand().unwrap_or(Strand::Plus);
    
    // Convert to 0-based for mapping
    let start_0based = start - 1;
    let end_0based = end; // end is exclusive in 0-based
    
    // Map coordinates
    let segments = mapper.map(chrom, start_0based, end_0based, query_strand)?;
    
    // Require single mapping
    if segments.len() != 1 {
        return None;
    }
    
    let seg = &segments[0];
    let target_chrom = &seg.target.chrom;
    let target_start = seg.target.start + 1; // Convert back to 1-based
    let target_end = seg.target.end;
    let target_strand = seg.target.strand;
    
    // Get new reference allele from target genome if available
    let new_ref = if let Some(ref_reader) = ref_genome {
        match ref_reader.fetch(target_chrom, seg.target.start, seg.target.end) {
            Some(seq) => {
                let seq_upper = seq.to_uppercase();
                // Reverse complement if target strand is negative
                if target_strand == Strand::Minus {
                    dna::revcomp(&seq_upper)
                } else {
                    seq_upper
                }
            }
            None => view.reference_allele().to_string(),
        }
    } else {
        // No reference genome, keep original or reverse complement
        let ref_allele = view.reference_allele();
        if target_strand == Strand::Minus && dna::is_dna(ref_allele) {
            dna::revcomp(ref_allele)
        } else {
            ref_allele.to_string()
        }
    };
    
    // Build output line with updated fields
    let mut output_fields: Vec<String> = view.fields().iter().map(|s| s.to_string()).collect();
    
    // Update Chromosome
    output_fields[view.indices.chromosome] = target_chrom.clone();
    
    // Update Start_Position
    output_fields[view.indices.start_position] = target_start.to_string();
    
    // Update End_Position
    output_fields[view.indices.end_position] = target_end.to_string();
    
    // Update Reference_Allele
    output_fields[view.indices.reference_allele] = new_ref;
    
    // Update NCBI_Build
    output_fields[view.indices.ncbi_build] = target_build.to_string();
    
    // Update Strand if target strand differs
    if target_strand == Strand::Minus {
        let current_strand = view.strand();
        let new_strand = match current_strand {
            Some(Strand::Plus) => "-",
            Some(Strand::Minus) => "+",
            None => ".",
        };
        output_fields[view.indices.strand] = new_strand.to_string();
    }
    
    Some(output_fields.join("\t"))
}

/// Convert a MAF file
///
/// # Arguments
/// * `input` - Input MAF file path
/// * `output` - Output MAF file path
/// * `mapper` - Coordinate mapper
/// * `ref_genome` - Optional path to target reference genome (FASTA)
/// * `target_build` - Target assembly name (e.g., "GRCh38")
///
/// # Returns
/// Conversion statistics
pub fn convert_maf<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    ref_genome: Option<P>,
    target_build: &str,
) -> Result<ConversionStats, std::io::Error> {
    let input_file = std::fs::File::open(input.as_ref())?;
    let reader = BufReader::with_capacity(128 * 1024, input_file);
    
    // Prepare output files with BufWriter for performance
    let output_path = output.as_ref();
    let unmap_path = output_path.with_extension("maf.unmap");
    
    let mut output_file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(output_path)?);
    let mut unmap_file = BufWriter::with_capacity(64 * 1024, std::fs::File::create(&unmap_path)?);
    
    // Open reference genome if provided
    let mut ref_reader = ref_genome
        .map(|p| fasta_stub::FastaReader::open(p.as_ref()))
        .transpose()?;
    
    // Counters
    let total = AtomicUsize::new(0);
    let success = AtomicUsize::new(0);
    let failed = AtomicUsize::new(0);
    let headers = AtomicUsize::new(0);
    
    let mut column_indices: Option<MafColumnIndices> = None;
    
    for line in reader.lines() {
        let line = line?;
        
        if line.is_empty() {
            continue;
        }
        
        // Handle header/comment lines
        if line.starts_with('#') {
            writeln!(output_file, "{}", line)?;
            headers.fetch_add(1, Ordering::Relaxed);
            continue;
        }
        
        // First non-comment line should be the column header
        if column_indices.is_none() {
            match MafColumnIndices::from_header(&line) {
                Ok(indices) => {
                    column_indices = Some(indices);
                    writeln!(output_file, "{}", line)?;
                    headers.fetch_add(1, Ordering::Relaxed);
                    continue;
                }
                Err(_) => {
                    // Not a valid header, treat as data
                }
            }
        }
        
        // Need column indices to process data
        let indices = match &column_indices {
            Some(i) => i,
            None => {
                writeln!(unmap_file, "{}", line)?;
                failed.fetch_add(1, Ordering::Relaxed);
                continue;
            }
        };
        
        total.fetch_add(1, Ordering::Relaxed);
        
        // Parse and convert
        match MafRecordView::parse(line.as_bytes(), indices) {
            Ok(view) => {
                if let Some(converted) = convert_maf_record(&view, mapper, ref_reader.as_mut(), target_build) {
                    writeln!(output_file, "{}", converted)?;
                    success.fetch_add(1, Ordering::Relaxed);
                } else {
                    writeln!(unmap_file, "{}", line)?;
                    failed.fetch_add(1, Ordering::Relaxed);
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

    fn create_test_indices() -> MafColumnIndices {
        MafColumnIndices {
            chromosome: 4,
            start_position: 5,
            end_position: 6,
            strand: 7,
            reference_allele: 10,
            ncbi_build: 3,
            hugo_symbol: Some(0),
        }
    }

    #[test]
    fn test_maf_column_indices_from_header() {
        let header = "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele";
        let indices = MafColumnIndices::from_header(header).unwrap();
        
        assert_eq!(indices.hugo_symbol, Some(0));
        assert_eq!(indices.ncbi_build, 3);
        assert_eq!(indices.chromosome, 4);
        assert_eq!(indices.start_position, 5);
        assert_eq!(indices.end_position, 6);
        assert_eq!(indices.strand, 7);
        assert_eq!(indices.reference_allele, 10);
    }

    #[test]
    fn test_maf_record_view_basic() {
        let indices = create_test_indices();
        let line = b"TP53\t7157\tBCM\tGRCh37\tchr17\t7577120\t7577120\t+\tMissense_Mutation\tSNP\tG\tG\tA\trs121912651\t.\tTCGA-A1-A0SK-01A";
        
        let view = MafRecordView::parse(line, &indices).unwrap();
        
        assert_eq!(view.hugo_symbol(), Some("TP53"));
        assert_eq!(view.chromosome(), "chr17");
        assert_eq!(view.start_position().unwrap(), 7577120);
        assert_eq!(view.end_position().unwrap(), 7577120);
        assert_eq!(view.strand(), Some(Strand::Plus));
        assert_eq!(view.reference_allele(), "G");
        assert_eq!(view.ncbi_build(), "GRCh37");
    }

    #[test]
    fn test_maf_record_view_negative_strand() {
        let indices = create_test_indices();
        let line = b"BRCA1\t672\tBCM\tGRCh37\tchr17\t41276044\t41276044\t-\tMissense_Mutation\tSNP\tC\tC\tT\t.\t.\tTCGA-A1-A0SK-01A";
        
        let view = MafRecordView::parse(line, &indices).unwrap();
        
        assert_eq!(view.strand(), Some(Strand::Minus));
        assert_eq!(view.reference_allele(), "C");
    }

    #[test]
    fn test_maf_column_indices_missing_column() {
        let header = "Hugo_Symbol\tEntrez_Gene_Id\tCenter";
        let result = MafColumnIndices::from_header(header);
        assert!(result.is_err());
    }

    #[test]
    fn test_maf_record_view_empty_line() {
        let indices = create_test_indices();
        let line = b"";
        let result = MafRecordView::parse(line, &indices);
        assert!(matches!(result, Err(MafParseError::EmptyLine)));
    }
}
