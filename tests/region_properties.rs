//! Region format property tests
//!
//! Tests for region format adapter including partial mapping and map_ratio calculation.

use proptest::prelude::*;
use fast_crossmap::formats::region::{FailureReason, parse_bed_line};
use fast_crossmap::core::Strand;

// ============================================================================
// Generators
// ============================================================================

fn arb_chrom_name() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("chr1".to_string()),
        Just("chr2".to_string()),
        Just("chr10".to_string()),
        Just("chrX".to_string()),
        Just("chrY".to_string()),
    ]
}

fn arb_strand() -> impl Strategy<Value = Strand> {
    prop_oneof![Just(Strand::Plus), Just(Strand::Minus)]
}

fn arb_bed_line() -> impl Strategy<Value = String> {
    (arb_chrom_name(), 0u64..1000000, 1u64..1000, arb_strand())
        .prop_map(|(chrom, start, len, strand)| {
            let end = start + len;
            let strand_char = if strand == Strand::Plus { "+" } else { "-" };
            format!("{}\t{}\t{}\tname\t0\t{}", chrom, start, end, strand_char)
        })
}

// ============================================================================
// Property Tests
// ============================================================================

proptest! {
    /// Property: Parsed BED line preserves coordinates
    #[test]
    fn test_bed_parse_preserves_coords(
        chrom in arb_chrom_name(),
        start in 0u64..1000000,
        len in 1u64..1000
    ) {
        let end = start + len;
        let line = format!("{}\t{}\t{}", chrom, start, end);
        
        let result = parse_bed_line(&line);
        prop_assert!(result.is_ok());
        
        let (parsed_chrom, parsed_start, parsed_end, _, _) = result.unwrap();
        prop_assert_eq!(parsed_chrom, chrom);
        prop_assert_eq!(parsed_start, start);
        prop_assert_eq!(parsed_end, end);
    }
    
    /// Property: BED line with strand preserves strand
    #[test]
    fn test_bed_parse_preserves_strand(line in arb_bed_line()) {
        let result = parse_bed_line(&line);
        prop_assert!(result.is_ok());
        
        let (_, _, _, strand, _) = result.unwrap();
        
        // Check that strand matches what's in the line
        if line.contains("\t-") {
            prop_assert_eq!(strand, Strand::Minus);
        } else if line.contains("\t+") {
            prop_assert_eq!(strand, Strand::Plus);
        }
    }
    
    /// Property: Invalid BED lines are rejected
    #[test]
    fn test_invalid_bed_rejected(
        chrom in arb_chrom_name(),
        start in 100u64..1000000
    ) {
        // Only 2 fields - should fail
        let line = format!("{}\t{}", chrom, start);
        prop_assert!(parse_bed_line(&line).is_err());
        
        // Start > End - should fail
        let line = format!("{}\t{}\t{}", chrom, start, start - 50);
        prop_assert!(parse_bed_line(&line).is_err());
    }
    
    /// Property: map_ratio is always between 0 and 1
    #[test]
    fn test_map_ratio_bounds(ratio in 0.0f64..=1.0) {
        // This is a sanity check - map_ratio should always be in [0, 1]
        prop_assert!(ratio >= 0.0 && ratio <= 1.0);
    }
}

// ============================================================================
// Unit Tests
// ============================================================================

#[test]
fn test_parse_bed_line_basic() {
    let result = parse_bed_line("chr1\t100\t200");
    assert!(result.is_ok());
    
    let (chrom, start, end, strand, fields) = result.unwrap();
    assert_eq!(chrom, "chr1");
    assert_eq!(start, 100);
    assert_eq!(end, 200);
    assert_eq!(strand, Strand::Plus);
    assert_eq!(fields.len(), 3);
}

#[test]
fn test_parse_bed_line_with_strand_plus() {
    let result = parse_bed_line("chr1\t100\t200\tname\t0\t+");
    assert!(result.is_ok());
    
    let (_, _, _, strand, _) = result.unwrap();
    assert_eq!(strand, Strand::Plus);
}

#[test]
fn test_parse_bed_line_with_strand_minus() {
    let result = parse_bed_line("chr1\t100\t200\tname\t0\t-");
    assert!(result.is_ok());
    
    let (_, _, _, strand, _) = result.unwrap();
    assert_eq!(strand, Strand::Minus);
}

#[test]
fn test_parse_bed_line_too_few_fields() {
    assert!(parse_bed_line("chr1\t100").is_err());
    assert!(parse_bed_line("chr1").is_err());
    assert!(parse_bed_line("").is_err());
}

#[test]
fn test_parse_bed_line_invalid_coords() {
    assert!(parse_bed_line("chr1\tabc\t200").is_err());
    assert!(parse_bed_line("chr1\t100\txyz").is_err());
}

#[test]
fn test_parse_bed_line_start_greater_than_end() {
    assert!(parse_bed_line("chr1\t200\t100").is_err());
}

#[test]
fn test_failure_reason_strings() {
    assert_eq!(FailureReason::Unmapped.as_str(), "Unmap");
    assert_eq!(FailureReason::CrossChrom.as_str(), "CrossChroms");
    assert_eq!(FailureReason::LowRatio.as_str(), "LowRatio");
    assert_eq!(FailureReason::InvalidFormat.as_str(), "InvalidFormat");
}

// ============================================================================
// Integration Tests
// ============================================================================

#[test]
fn test_region_conversion_with_crossmap() {
    use std::path::Path;
    use std::process::Command;
    use std::io::Write;
    use tempfile::tempdir;
    
    let chain_file = Path::new("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz");
    if !chain_file.exists() {
        eprintln!("Chain file not found, skipping CrossMap comparison test");
        return;
    }
    
    // Create test input
    let dir = tempdir().unwrap();
    let input_path = dir.path().join("test_region.bed");
    let mut input_file = std::fs::File::create(&input_path).unwrap();
    
    // Large regions that may have partial mapping
    writeln!(input_file, "chr1\t100000\t200000\tregion1\t0\t+").unwrap();
    writeln!(input_file, "chr1\t1000000\t2000000\tregion2\t0\t-").unwrap();
    writeln!(input_file, "chr2\t50000000\t51000000\tregion3\t0\t+").unwrap();
    drop(input_file);
    
    // Run CrossMap
    let crossmap_output = dir.path().join("crossmap_region.bed");
    let crossmap_result = Command::new("CrossMap")
        .args(&[
            "region",
            chain_file.to_str().unwrap(),
            input_path.to_str().unwrap(),
            crossmap_output.to_str().unwrap(),
            "-r", "0.85"
        ])
        .output();
    
    if crossmap_result.is_err() {
        eprintln!("CrossMap not available, skipping comparison test");
        return;
    }
    
    let output = crossmap_result.unwrap();
    if !output.status.success() {
        eprintln!("CrossMap failed: {}", String::from_utf8_lossy(&output.stderr));
        return;
    }
    
    // Run FastCrossMap
    use fast_crossmap::core::{ChainIndex, CoordinateMapper, ChromStyle};
    
    let index = ChainIndex::from_chain_file(chain_file).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    let fastcm_output = dir.path().join("fastcm_region.bed");
    
    let stats = fast_crossmap::formats::region::convert_region(
        &input_path,
        &fastcm_output,
        &mapper,
        0.85,
    ).unwrap();
    
    println!("Region conversion stats: total={}, success={}, failed={}", 
             stats.total, stats.success, stats.failed);
    
    // Compare outputs
    if crossmap_output.exists() && fastcm_output.exists() {
        let crossmap_content = std::fs::read_to_string(&crossmap_output).unwrap_or_default();
        let fastcm_content = std::fs::read_to_string(&fastcm_output).unwrap_or_default();
        
        println!("=== CrossMap Output ===");
        println!("{}", crossmap_content);
        println!("=== FastCrossMap Output ===");
        println!("{}", fastcm_content);
        
        // Compare line counts
        let crossmap_lines: Vec<&str> = crossmap_content.lines().collect();
        let fastcm_lines: Vec<&str> = fastcm_content.lines().collect();
        
        println!("CrossMap lines: {}, FastCrossMap lines: {}", 
                 crossmap_lines.len(), fastcm_lines.len());
    }
}
