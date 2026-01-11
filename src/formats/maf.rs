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
    /// Loads all sequences into memory at once for fast access
    pub struct FastaReader {
        /// Chromosome sequences
        sequences: HashMap<String, Vec<u8>>,
    }
    
    impl FastaReader {
        /// Open a FASTA file and load all sequences into memory
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
        
        /// Fetch sequence at given position (0-based, half-open)
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
    ref_genome: Option<&fasta_stub::FastaReader>,
    target_build: &str,
) -> Option<String> {
    // Get coordinates (MAF uses 1-based coordinates)
    let start = view.start_position().ok()?;
    let end = view.end_position().ok()?;
    let chrom = view.chromosome();
    
    // Convert to 0-based for mapping
    // CrossMap uses: start = int(fields[5])-1, end = int(fields[6])
    let start_0based = start - 1;
    let end_0based = end; // end is exclusive in 0-based
    
    // Map coordinates - CrossMap always uses '+' strand for mapping
    let segments = mapper.map(chrom, start_0based, end_0based, Strand::Plus)?;
    
    // Require single mapping (len(a) == 2 in CrossMap means one mapping)
    if segments.len() != 1 {
        return None;
    }
    
    let seg = &segments[0];
    let target_chrom = &seg.target.chrom;
    let target_start = seg.target.start + 1; // Convert back to 1-based
    let target_end = seg.target.end;
    let target_strand = seg.target.strand;
    
    // Get new reference allele from target genome if available
    // CrossMap: fields[10] = refFasta.fetch(target_chr, target_start, target_end).upper()
    // Then: if a[1][3] == '-': fields[10] = revcomp_DNA(fields[10], True)
    let new_ref = if let Some(ref_reader) = ref_genome {
        match ref_reader.fetch(target_chrom, seg.target.start, seg.target.end) {
            Some(seq) => {
                let seq_upper = seq.to_uppercase();
                // Reverse complement if target strand is negative (matches CrossMap)
                if target_strand == Strand::Minus {
                    dna::revcomp(&seq_upper)
                } else {
                    seq_upper
                }
            }
            None => return None, // CrossMap fails if fetch fails
        }
    } else {
        // No reference genome - this shouldn't happen for MAF
        view.reference_allele().to_string()
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
    
    // Update NCBI_Build (CrossMap: fields[3] = ref_name)
    output_fields[view.indices.ncbi_build] = target_build.to_string();
    
    // NOTE: CrossMap does NOT update the Strand field, only the Reference_Allele
    // So we should NOT flip the strand here to match CrossMap behavior
    
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
    let ref_reader = ref_genome
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
                if let Some(converted) = convert_maf_record(&view, mapper, ref_reader.as_ref(), target_build) {
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
