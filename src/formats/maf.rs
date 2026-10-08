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


/// Conversion statistics
#[derive(Debug, Clone, Default)]
pub struct ConversionStats {
    pub total: usize,
    pub success: usize,
    pub failed: usize,
    pub headers: usize,
}


/// Intermediate result from coordinate mapping (no ref genome needed)
struct MappedRecord {
    fields: Vec<String>,
    indices: MafColumnIndices,
    target_chrom: String,
    target_start_0based: u64,
    target_end_0based: u64,
    target_start_1based: u64,
    target_end_1based: u64,
    target_strand: Strand,
}

/// A pending item in a processing chunk: either mapped or failed
enum ChunkItem {
    Mapped(MappedRecord),
    Failed(String),
}

const CHUNK_SIZE: usize = 10_000;

fn flush_chunk(
    chunk: &mut Vec<ChunkItem>,
    ref_reader: Option<&crate::core::fasta::FastaReader>,
    output_file: &mut BufWriter<std::fs::File>,
    unmap_file: &mut BufWriter<std::fs::File>,
    target_build: &str,
    success: &mut usize,
    failed: &mut usize,
) -> std::io::Result<()> {
    if chunk.is_empty() {
        return Ok(());
    }

    // Collect mapped indices and build batch requests
    let mapped_indices: Vec<usize> = chunk.iter().enumerate()
        .filter_map(|(i, item)| matches!(item, ChunkItem::Mapped(_)).then_some(i))
        .collect();

    let ref_seqs = if let Some(reader) = ref_reader {
        let requests: Vec<(String, u64, u64)> = mapped_indices.iter()
            .map(|&i| {
                if let ChunkItem::Mapped(rec) = &chunk[i] {
                    (rec.target_chrom.clone(), rec.target_start_0based, rec.target_end_0based)
                } else {
                    unreachable!()
                }
            })
            .collect();
        Some(reader.batch_fetch(&requests))
    } else {
        None
    };

    let mut fetch_idx = 0;
    for item in chunk.drain(..) {
        match item {
            ChunkItem::Mapped(rec) => {
                let new_ref = if let Some(ref seqs) = ref_seqs {
                    match &seqs[fetch_idx] {
                        Some(seq) => {
                            let seq_upper = seq.to_uppercase();
                            Some(if rec.target_strand == Strand::Minus {
                                dna::revcomp(&seq_upper)
                            } else {
                                seq_upper
                            })
                        }
                        None => None,
                    }
                } else {
                    Some(rec.fields[rec.indices.reference_allele].clone())
                };

                if let Some(new_ref) = new_ref {
                    let mut out_fields = rec.fields;
                    out_fields[rec.indices.chromosome] = rec.target_chrom;
                    out_fields[rec.indices.start_position] = rec.target_start_1based.to_string();
                    out_fields[rec.indices.end_position] = rec.target_end_1based.to_string();
                    out_fields[rec.indices.reference_allele] = new_ref;
                    out_fields[rec.indices.ncbi_build] = target_build.to_string();
                    writeln!(output_file, "{}", out_fields.join("\t"))?;
                    *success += 1;
                } else {
                    writeln!(unmap_file, "{}", rec.fields.join("\t"))?;
                    *failed += 1;
                }

                if ref_seqs.is_some() {
                    fetch_idx += 1;
                }
            }
            ChunkItem::Failed(line) => {
                writeln!(unmap_file, "{}", line)?;
            }
        }
    }

    Ok(())
}

/// Convert a MAF file using chunked batch-fetch for low memory.
pub fn convert_maf<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    ref_genome: Option<P>,
    target_build: &str,
) -> Result<ConversionStats, std::io::Error> {
    let input_file = std::fs::File::open(input.as_ref())?;
    let reader = BufReader::with_capacity(128 * 1024, input_file);

    let output_path = output.as_ref();
    let unmap_path = output_path.with_extension("maf.unmap");

    let mut output_file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(output_path)?);
    let mut unmap_file = BufWriter::with_capacity(64 * 1024, std::fs::File::create(&unmap_path)?);

    let ref_reader = ref_genome
        .map(|p| crate::core::fasta::FastaReader::open(p.as_ref()))
        .transpose()?;

    let mut total: usize = 0;
    let mut success: usize = 0;
    let mut failed: usize = 0;
    let mut headers: usize = 0;

    let mut column_indices: Option<MafColumnIndices> = None;
    let mut chunk: Vec<ChunkItem> = Vec::with_capacity(CHUNK_SIZE);

    for line in reader.lines() {
        let line = line?;

        if line.is_empty() {
            continue;
        }

        if line.starts_with('#') {
            writeln!(output_file, "{}", line)?;
            headers += 1;
            continue;
        }

        if column_indices.is_none() {
            match MafColumnIndices::from_header(&line) {
                Ok(indices) => {
                    column_indices = Some(indices);
                    writeln!(output_file, "{}", line)?;
                    headers += 1;
                    continue;
                }
                Err(_) => {}
            }
        }

        let indices = match &column_indices {
            Some(i) => i,
            None => {
                chunk.push(ChunkItem::Failed(line));
                failed += 1;
                if chunk.len() >= CHUNK_SIZE {
                    flush_chunk(&mut chunk, ref_reader.as_ref(), &mut output_file, &mut unmap_file, target_build, &mut success, &mut failed)?;
                }
                continue;
            }
        };

        total += 1;

        let parsed = MafRecordView::parse(line.as_bytes(), indices);
        let mapped = parsed.ok().and_then(|view| {
            let start = view.start_position().ok()?;
            let end = view.end_position().ok()?;
            let chrom = view.chromosome();
            let start_0based = start - 1;
            let end_0based = end;

            let segments = mapper.map(chrom, start_0based, end_0based, Strand::Plus)?;
            if segments.len() != 1 {
                return None;
            }

            let seg = &segments[0];
            Some(MappedRecord {
                fields: view.fields().iter().map(|s| s.to_string()).collect(),
                indices: indices.clone(),
                target_chrom: seg.target.chrom.clone(),
                target_start_0based: seg.target.start,
                target_end_0based: seg.target.end,
                target_start_1based: seg.target.start + 1,
                target_end_1based: seg.target.end,
                target_strand: seg.target.strand,
            })
        });

        match mapped {
            Some(rec) => chunk.push(ChunkItem::Mapped(rec)),
            None => {
                chunk.push(ChunkItem::Failed(line));
                failed += 1;
            }
        }

        if chunk.len() >= CHUNK_SIZE {
            flush_chunk(&mut chunk, ref_reader.as_ref(), &mut output_file, &mut unmap_file, target_build, &mut success, &mut failed)?;
        }
    }

    // Flush remaining
    flush_chunk(&mut chunk, ref_reader.as_ref(), &mut output_file, &mut unmap_file, target_build, &mut success, &mut failed)?;

    Ok(ConversionStats {
        total,
        success,
        failed,
        headers,
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
