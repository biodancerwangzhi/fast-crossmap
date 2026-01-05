//! Region format adapter
//!
//! Handles large genomic region conversion with partial mapping support.
//! Unlike BED conversion which requires exact mapping, region conversion
//! allows partial matches based on a minimum mapping ratio.
//!
//! **Validates: Requirements 11.1, 11.2, 11.3, 11.4, 11.5, 11.6**

use crate::core::{CoordinateMapper, Strand};
use std::collections::HashSet;
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::fs::File;
use std::path::Path;

/// Region conversion error
#[derive(Debug)]
pub enum RegionError {
    IoError(std::io::Error),
    InvalidFormat(String),
}

impl std::fmt::Display for RegionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RegionError::IoError(e) => write!(f, "IO error: {}", e),
            RegionError::InvalidFormat(msg) => write!(f, "Invalid format: {}", msg),
        }
    }
}

impl std::error::Error for RegionError {}

impl From<std::io::Error> for RegionError {
    fn from(e: std::io::Error) -> Self {
        RegionError::IoError(e)
    }
}

/// Conversion statistics
#[derive(Debug, Clone, Default)]
pub struct ConversionStats {
    pub total: usize,
    pub success: usize,
    pub failed: usize,
    pub cross_chrom: usize,
    pub low_ratio: usize,
    pub unmapped: usize,
}

/// Region mapping result
#[derive(Debug, Clone)]
pub struct RegionResult {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub map_ratio: f64,
}

/// Failure reason for region mapping
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FailureReason {
    Unmapped,
    CrossChrom,
    LowRatio,
    InvalidFormat,
}

impl FailureReason {
    pub fn as_str(&self) -> &'static str {
        match self {
            FailureReason::Unmapped => "Unmap",
            FailureReason::CrossChrom => "CrossChroms",
            FailureReason::LowRatio => "LowRatio",
            FailureReason::InvalidFormat => "InvalidFormat",
        }
    }
}

/// Map a single region with partial mapping support
///
/// Returns Ok(RegionResult) if mapping succeeds with ratio >= min_ratio
/// Returns Err(FailureReason) if mapping fails
pub fn map_region(
    mapper: &CoordinateMapper,
    chrom: &str,
    start: u64,
    end: u64,
    strand: Strand,
    min_ratio: f64,
) -> Result<RegionResult, FailureReason> {
    let total_query_length = end - start;
    if total_query_length == 0 {
        return Err(FailureReason::InvalidFormat);
    }
    
    // Get all mapping segments
    let segments = mapper.map(chrom, start, end, strand);
    
    if segments.is_none() || segments.as_ref().map(|s| s.is_empty()).unwrap_or(true) {
        return Err(FailureReason::Unmapped);
    }
    
    let segments = segments.unwrap();
    
    // Single segment = 100% match
    if segments.len() == 1 {
        let seg = &segments[0];
        return Ok(RegionResult {
            chrom: seg.target.chrom.clone(),
            start: seg.target.start,
            end: seg.target.end,
            strand: seg.target.strand,
            map_ratio: 1.0,
        });
    }
    
    // Multiple segments - calculate mapping ratio
    let mut mapped_bases: u64 = 0;
    let mut target_chroms = HashSet::new();
    let mut target_starts = Vec::new();
    let mut target_ends = Vec::new();
    let mut target_strand = Strand::Plus;
    
    for seg in &segments {
        // Count mapped bases from source
        mapped_bases += seg.source.end - seg.source.start;
        
        // Collect target info
        target_chroms.insert(seg.target.chrom.clone());
        target_starts.push(seg.target.start);
        target_ends.push(seg.target.end);
        target_strand = seg.target.strand;
    }
    
    let map_ratio = mapped_bases as f64 / total_query_length as f64;
    
    // Check if mapping crosses chromosomes
    if target_chroms.len() > 1 {
        return Err(FailureReason::CrossChrom);
    }
    
    // Check if mapping ratio meets threshold
    if map_ratio < min_ratio {
        return Err(FailureReason::LowRatio);
    }
    
    // Merge all target segments into one region
    let target_chrom = target_chroms.into_iter().next().unwrap();
    let merged_start = *target_starts.iter().min().unwrap();
    let merged_end = *target_ends.iter().max().unwrap();
    
    Ok(RegionResult {
        chrom: target_chrom,
        start: merged_start,
        end: merged_end,
        strand: target_strand,
        map_ratio,
    })
}

/// Parse a BED line and extract region info
pub fn parse_bed_line(line: &str) -> Result<(String, u64, u64, Strand, Vec<&str>), FailureReason> {
    let fields: Vec<&str> = line.split('\t').collect();
    
    if fields.len() < 3 {
        return Err(FailureReason::InvalidFormat);
    }
    
    let chrom = fields[0].to_string();
    let start: u64 = fields[1].parse().map_err(|_| FailureReason::InvalidFormat)?;
    let end: u64 = fields[2].parse().map_err(|_| FailureReason::InvalidFormat)?;
    
    if start > end {
        return Err(FailureReason::InvalidFormat);
    }
    
    // Try to find strand in fields
    let mut strand = Strand::Plus;
    for field in &fields {
        if *field == "+" {
            strand = Strand::Plus;
            break;
        } else if *field == "-" {
            strand = Strand::Minus;
            break;
        }
    }
    
    Ok((chrom, start, end, strand, fields))
}

/// Convert a region BED file
///
/// # Arguments
/// * `input` - Input BED file path
/// * `output` - Output BED file path (mapped regions)
/// * `mapper` - Coordinate mapper
/// * `min_ratio` - Minimum mapping ratio (default 0.85)
///
/// # Returns
/// Conversion statistics
pub fn convert_region<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    min_ratio: f64,
) -> Result<ConversionStats, RegionError> {
    let input_file = File::open(input.as_ref())?;
    let reader = BufReader::new(input_file);
    
    let output_file = File::create(output.as_ref())?;
    let mut writer = BufWriter::new(output_file);
    
    // Create unmap file
    let unmap_path = format!("{}.unmap", output.as_ref().display());
    let unmap_file = File::create(&unmap_path)?;
    let mut unmap_writer = BufWriter::new(unmap_file);
    
    let mut stats = ConversionStats::default();
    
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        
        // Skip comments and empty lines
        if trimmed.is_empty() 
            || trimmed.starts_with('#') 
            || trimmed.starts_with("track") 
            || trimmed.starts_with("browser") 
        {
            continue;
        }
        
        stats.total += 1;
        
        // Parse BED line
        let parsed = parse_bed_line(trimmed);
        if parsed.is_err() {
            writeln!(unmap_writer, "{}\tFail\tInvalidFormat", trimmed)?;
            stats.failed += 1;
            continue;
        }
        
        let (chrom, start, end, strand, fields) = parsed.unwrap();
        
        // Map the region
        match map_region(mapper, &chrom, start, end, strand, min_ratio) {
            Ok(result) => {
                // Build output line with updated coordinates
                let mut out_fields: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
                out_fields[0] = result.chrom;
                out_fields[1] = result.start.to_string();
                out_fields[2] = result.end.to_string();
                
                // Update strand if present
                for field in &mut out_fields {
                    if *field == "+" || *field == "-" {
                        *field = if result.strand == Strand::Plus { "+".to_string() } else { "-".to_string() };
                    }
                }
                
                writeln!(writer, "{}\tmap_ratio={:.4}", out_fields.join("\t"), result.map_ratio)?;
                stats.success += 1;
            }
            Err(reason) => {
                match reason {
                    FailureReason::Unmapped => {
                        writeln!(unmap_writer, "{}\tFail\t{}", trimmed, reason.as_str())?;
                        stats.unmapped += 1;
                    }
                    FailureReason::CrossChrom => {
                        writeln!(unmap_writer, "{}\tFail\t{}", trimmed, reason.as_str())?;
                        stats.cross_chrom += 1;
                    }
                    FailureReason::LowRatio => {
                        // For low ratio, we still want to show the ratio
                        // Need to recalculate to get the actual ratio
                        let total_len = end - start;
                        if let Some(segments) = mapper.map(&chrom, start, end, strand) {
                            let mapped: u64 = segments.iter()
                                .map(|s| s.source.end - s.source.start)
                                .sum();
                            let ratio = mapped as f64 / total_len as f64;
                            writeln!(unmap_writer, "{}\tFail\tmap_ratio={:.4}", trimmed, ratio)?;
                        } else {
                            writeln!(unmap_writer, "{}\tFail\t{}", trimmed, reason.as_str())?;
                        }
                        stats.low_ratio += 1;
                    }
                    FailureReason::InvalidFormat => {
                        writeln!(unmap_writer, "{}\tFail\t{}", trimmed, reason.as_str())?;
                        stats.failed += 1;
                    }
                }
                stats.failed += 1;
            }
        }
    }
    
    Ok(stats)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_bed_line_basic() {
        let (chrom, start, end, strand, fields) = parse_bed_line("chr1\t100\t200").unwrap();
        assert_eq!(chrom, "chr1");
        assert_eq!(start, 100);
        assert_eq!(end, 200);
        assert_eq!(strand, Strand::Plus);
        assert_eq!(fields.len(), 3);
    }

    #[test]
    fn test_parse_bed_line_with_strand() {
        let (chrom, start, end, strand, _) = parse_bed_line("chr1\t100\t200\tname\t0\t-").unwrap();
        assert_eq!(chrom, "chr1");
        assert_eq!(start, 100);
        assert_eq!(end, 200);
        assert_eq!(strand, Strand::Minus);
    }

    #[test]
    fn test_parse_bed_line_invalid() {
        assert!(parse_bed_line("chr1\t100").is_err());
        assert!(parse_bed_line("chr1\tabc\t200").is_err());
        assert!(parse_bed_line("chr1\t200\t100").is_err()); // start > end
    }

    #[test]
    fn test_failure_reason_as_str() {
        assert_eq!(FailureReason::Unmapped.as_str(), "Unmap");
        assert_eq!(FailureReason::CrossChrom.as_str(), "CrossChroms");
        assert_eq!(FailureReason::LowRatio.as_str(), "LowRatio");
        assert_eq!(FailureReason::InvalidFormat.as_str(), "InvalidFormat");
    }
}
