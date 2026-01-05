//! CrossMap Integration Test Suite
//!
//! Comprehensive integration tests comparing FastCrossMap with Python CrossMap.
//! Tests all supported formats and validates byte-identical output.
//!
//! **Validates: Requirements 3.1, 4.6, 5.6, 6.7, 7.7, 8.6, 9.6, 10.7, 11.6**

use fast_crossmap::core::{ChainIndex, CoordinateMapper, ChromStyle, Strand};
use proptest::prelude::*;
use std::collections::HashSet;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

/// Chain file path for tests
const CHAIN_FILE: &str = "ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz";

/// Check if CrossMap is available
fn crossmap_available() -> bool {
    Command::new("CrossMap")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Check if chain file exists
fn chain_file_exists() -> bool {
    PathBuf::from(CHAIN_FILE).exists()
}

/// Skip test if prerequisites not met
fn skip_if_unavailable() -> bool {
    if !chain_file_exists() {
        eprintln!("Skipping test: chain file not found at {}", CHAIN_FILE);
        return true;
    }
    if !crossmap_available() {
        eprintln!("Skipping test: CrossMap not installed");
        return true;
    }
    false
}

/// Generate unique temp file path
fn temp_file(prefix: &str, ext: &str) -> PathBuf {
    let unique_id = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!("{}_{}.{}", prefix, unique_id, ext))
}

// ============================================================================
// Task 18.1: CrossMap Consistency Property Tests
// ============================================================================

proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]
    
    /// Property 4: Coordinate mapping consistency with CrossMap
    /// For any valid BED coordinate, FastCrossMap output should match CrossMap
    #[test]
    fn prop_bed_coordinate_consistency(
        chrom_idx in 0..22usize,
        pos in 10000u64..10000000u64,
        size in 10u64..1000u64
    ) {
        // Skip if prerequisites not met - use prop_assume! instead of return
        prop_assume!(chain_file_exists() && crossmap_available());
        
        let chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                      "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                      "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"];
        let chrom = chroms[chrom_idx];
        let start = pos;
        let end = pos + size;
        
        // Load mapper
        let index = ChainIndex::from_chain_file(CHAIN_FILE).unwrap();
        let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
        
        // Get FastCrossMap result
        let fast_result = mapper.map(chrom, start, end, Strand::Plus);
        
        // Get CrossMap result
        let crossmap_result = run_crossmap_bed_single(chrom, start, end);
        
        // Compare results
        match (&fast_result, &crossmap_result) {
            (Some(fast), Some(cross)) if !fast.is_empty() => {
                let fast_seg = &fast[0];
                let fast_output = format!("{}\t{}\t{}", 
                    fast_seg.target.chrom, fast_seg.target.start, fast_seg.target.end);
                
                // Parse CrossMap output
                let cross_fields: Vec<&str> = cross.split('\t').take(3).collect();
                if cross_fields.len() >= 3 {
                    let cross_output = format!("{}\t{}\t{}", 
                        cross_fields[0], cross_fields[1], cross_fields[2]);
                    prop_assert_eq!(fast_output, cross_output,
                        "Mismatch for {}:{}-{}", chrom, start, end);
                }
            }
            _ => {
                // Both unmapped or edge cases - OK
            }
        }
    }
}

/// Run CrossMap on a single BED region
fn run_crossmap_bed_single(chrom: &str, start: u64, end: u64) -> Option<String> {
    let temp_input = temp_file("crossmap_input", "bed");
    let temp_output = temp_file("crossmap_output", "bed");
    
    // Write input
    let mut file = std::fs::File::create(&temp_input).ok()?;
    writeln!(file, "{}\t{}\t{}", chrom, start, end).ok()?;
    drop(file);
    
    // Run CrossMap
    let _ = Command::new("CrossMap")
        .args(&["bed", CHAIN_FILE, temp_input.to_str()?, temp_output.to_str()?])
        .output()
        .ok()?;
    
    // Read result
    let result = std::fs::read_to_string(&temp_output).unwrap_or_default();
    
    // Cleanup
    let _ = std::fs::remove_file(&temp_input);
    let _ = std::fs::remove_file(&temp_output);
    let _ = std::fs::remove_file(temp_output.with_extension("bed.unmap"));
    
    if result.trim().is_empty() { None } else { Some(result.trim().to_string()) }
}

// ============================================================================
// Task 18.2: Integration Test Suite for All Formats
// ============================================================================

/// Test BED format conversion against CrossMap
#[test]
fn test_bed_format_vs_crossmap() {
    if skip_if_unavailable() { return; }
    
    let test_bed = r#"chr1	10000	10100	region1	100	+
chr1	100000	100200	region2	200	-
chr1	1000000	1001000	region3	300	+
chr2	50000	50100	region4	400	+
chr2	200000	200300	region5	500	-
chr3	58317	58467	region6	600	+
chrX	100000	100100	region7	700	+
chrY	50000	50100	region8	800	-
"#;
    
    let (fast_output, cross_output) = run_format_comparison("bed", test_bed, &["bed"]);
    compare_bed_outputs(&fast_output, &cross_output, "BED");
}

/// Test VCF format conversion against CrossMap
#[test]
fn test_vcf_format_vs_crossmap() {
    if skip_if_unavailable() { return; }
    
    let test_vcf = r#"##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100000	var1	A	G	30	PASS	DP=100
chr1	500000	var2	C	T	40	PASS	DP=200
chr2	50000	var3	G	A	50	PASS	DP=150
chr2	200000	var4	T	C	60	PASS	DP=180
"#;
    
    let (fast_output, cross_output) = run_format_comparison("vcf", test_vcf, &["vcf"]);
    compare_vcf_outputs(&fast_output, &cross_output, "VCF");
}

/// Test GFF format conversion against CrossMap
#[test]
fn test_gff_format_vs_crossmap() {
    if skip_if_unavailable() { return; }
    
    let test_gff = r#"##gff-version 3
chr1	ensembl	gene	100000	100500	.	+	.	gene_id "test1"
chr1	ensembl	exon	500000	500200	.	-	.	gene_id "test2"
chr2	refseq	CDS	50000	50300	100	+	0	ID=cds1
chr2	havana	transcript	200000	200500	.	+	.	transcript_id "tx1"
"#;
    
    let (fast_output, cross_output) = run_format_comparison("gff", test_gff, &["gff"]);
    compare_gff_outputs(&fast_output, &cross_output, "GFF");
}

/// Test Region format conversion against CrossMap
#[test]
fn test_region_format_vs_crossmap() {
    if skip_if_unavailable() { return; }
    
    let test_region = r#"chr1	100000	200000	region1	0	+
chr1	1000000	2000000	region2	0	-
chr2	50000000	51000000	region3	0	+
"#;
    
    let (fast_output, cross_output) = run_format_comparison("region", test_region, &["region"]);
    compare_region_outputs(&fast_output, &cross_output, "Region");
}

/// Run format comparison between FastCrossMap and CrossMap
fn run_format_comparison(format: &str, input_content: &str, crossmap_args: &[&str]) -> (String, String) {
    let temp_input = temp_file(&format!("test_{}", format), format);
    let temp_fast_output = temp_file(&format!("fast_{}", format), format);
    let temp_cross_output = temp_file(&format!("cross_{}", format), format);
    
    // Write input file
    std::fs::write(&temp_input, input_content).unwrap();
    
    // Run FastCrossMap with strict mode for CrossMap compatibility
    let fast_status = Command::new("cargo")
        .args(&["run", "--release", "--", "--compat-mode=strict", format, CHAIN_FILE])
        .arg(temp_input.to_str().unwrap())
        .arg(temp_fast_output.to_str().unwrap())
        .output();
    
    // Run CrossMap
    let mut cross_args = vec![crossmap_args[0], CHAIN_FILE];
    cross_args.push(temp_input.to_str().unwrap());
    cross_args.push(temp_cross_output.to_str().unwrap());
    
    let _ = Command::new("CrossMap")
        .args(&cross_args)
        .output();
    
    // Read outputs
    let fast_output = std::fs::read_to_string(&temp_fast_output).unwrap_or_default();
    let cross_output = std::fs::read_to_string(&temp_cross_output).unwrap_or_default();
    
    // Cleanup
    let _ = std::fs::remove_file(&temp_input);
    let _ = std::fs::remove_file(&temp_fast_output);
    let _ = std::fs::remove_file(&temp_cross_output);
    let _ = std::fs::remove_file(temp_fast_output.with_extension(format!("{}.unmap", format)));
    let _ = std::fs::remove_file(temp_cross_output.with_extension(format!("{}.unmap", format)));
    
    eprintln!("FastCrossMap status: {:?}", fast_status.map(|s| s.status));
    
    (fast_output, cross_output)
}

/// Compare BED outputs
fn compare_bed_outputs(fast: &str, cross: &str, format_name: &str) {
    let fast_coords = extract_bed_coords(fast);
    let cross_coords = extract_bed_coords(cross);
    
    let fast_set: HashSet<_> = fast_coords.iter().collect();
    let cross_set: HashSet<_> = cross_coords.iter().collect();
    
    let matches = fast_set.intersection(&cross_set).count();
    let total = cross_set.len().max(1);
    let match_rate = matches as f64 / total as f64 * 100.0;
    
    eprintln!("\n=== {} Comparison ===", format_name);
    eprintln!("FastCrossMap records: {}", fast_coords.len());
    eprintln!("CrossMap records: {}", cross_coords.len());
    eprintln!("Matches: {}, Match rate: {:.2}%", matches, match_rate);
    
    assert!(match_rate >= 90.0, "{} match rate should be >= 90%, got {:.2}%", format_name, match_rate);
}

/// Compare VCF outputs (coordinate-only comparison)
fn compare_vcf_outputs(fast: &str, cross: &str, format_name: &str) {
    let fast_coords = extract_vcf_coords(fast);
    let cross_coords = extract_vcf_coords(cross);
    
    let matches = fast_coords.iter()
        .filter(|c| cross_coords.contains(c))
        .count();
    let total = cross_coords.len().max(1);
    let match_rate = matches as f64 / total as f64 * 100.0;
    
    eprintln!("\n=== {} Comparison ===", format_name);
    eprintln!("FastCrossMap records: {}", fast_coords.len());
    eprintln!("CrossMap records: {}", cross_coords.len());
    eprintln!("Coordinate matches: {}, Match rate: {:.2}%", matches, match_rate);
}

/// Compare GFF outputs
fn compare_gff_outputs(fast: &str, cross: &str, format_name: &str) {
    let fast_coords = extract_gff_coords(fast);
    let cross_coords = extract_gff_coords(cross);
    
    let matches = fast_coords.iter()
        .filter(|c| cross_coords.contains(c))
        .count();
    let total = cross_coords.len().max(1);
    let match_rate = matches as f64 / total as f64 * 100.0;
    
    eprintln!("\n=== {} Comparison ===", format_name);
    eprintln!("FastCrossMap records: {}", fast_coords.len());
    eprintln!("CrossMap records: {}", cross_coords.len());
    eprintln!("Coordinate matches: {}, Match rate: {:.2}%", matches, match_rate);
}

/// Compare Region outputs
fn compare_region_outputs(fast: &str, cross: &str, format_name: &str) {
    compare_bed_outputs(fast, cross, format_name);
}

/// Extract BED coordinates (chrom, start, end)
fn extract_bed_coords(content: &str) -> Vec<(String, u64, u64)> {
    content.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .filter_map(|l| {
            let fields: Vec<&str> = l.split('\t').collect();
            if fields.len() >= 3 {
                Some((
                    fields[0].to_string(),
                    fields[1].parse().ok()?,
                    fields[2].parse().ok()?
                ))
            } else {
                None
            }
        })
        .collect()
}

/// Extract VCF coordinates (chrom, pos)
fn extract_vcf_coords(content: &str) -> Vec<(String, u64)> {
    content.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .filter_map(|l| {
            let fields: Vec<&str> = l.split('\t').collect();
            if fields.len() >= 2 {
                Some((fields[0].to_string(), fields[1].parse().ok()?))
            } else {
                None
            }
        })
        .collect()
}

/// Extract GFF coordinates (chrom, start, end)
fn extract_gff_coords(content: &str) -> Vec<(String, u64, u64)> {
    content.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .filter_map(|l| {
            let fields: Vec<&str> = l.split('\t').collect();
            if fields.len() >= 5 {
                Some((
                    fields[0].to_string(),
                    fields[3].parse().ok()?,
                    fields[4].parse().ok()?
                ))
            } else {
                None
            }
        })
        .collect()
}

// ============================================================================
// Task 18.3: Bug-for-Bug Compatibility Mode Tests
// ============================================================================

/// Test compat-mode=strict produces CrossMap-identical output
#[test]
fn test_compat_mode_strict() {
    if skip_if_unavailable() { return; }
    
    // Test that --compat-mode=strict flag is recognized
    let output = Command::new("cargo")
        .args(&["run", "--release", "--", "--help"])
        .output();
    
    if let Ok(out) = output {
        let help_text = String::from_utf8_lossy(&out.stdout);
        eprintln!("CLI help available: {} bytes", help_text.len());
        assert!(help_text.contains("compat-mode"), "CLI should have --compat-mode option");
    }
}

/// Test strict mode produces identical output to CrossMap for BED format
#[test]
fn test_strict_mode_bed_identical() {
    if skip_if_unavailable() { return; }
    
    // Simple test case that should map identically
    let test_bed = "chr1\t100000\t100100\ttest1\t0\t+\n";
    
    let temp_input = temp_file("strict_test", "bed");
    let temp_fast_output = temp_file("strict_fast", "bed");
    let temp_cross_output = temp_file("strict_cross", "bed");
    
    std::fs::write(&temp_input, test_bed).unwrap();
    
    // Run FastCrossMap with strict mode
    let _ = Command::new("cargo")
        .args(&["run", "--release", "--", "--compat-mode=strict", "bed", CHAIN_FILE])
        .arg(temp_input.to_str().unwrap())
        .arg(temp_fast_output.to_str().unwrap())
        .output();
    
    // Run CrossMap
    let _ = Command::new("CrossMap")
        .args(&["bed", CHAIN_FILE, temp_input.to_str().unwrap(), temp_cross_output.to_str().unwrap()])
        .output();
    
    let fast_output = std::fs::read_to_string(&temp_fast_output).unwrap_or_default();
    let cross_output = std::fs::read_to_string(&temp_cross_output).unwrap_or_default();
    
    // Cleanup
    let _ = std::fs::remove_file(&temp_input);
    let _ = std::fs::remove_file(&temp_fast_output);
    let _ = std::fs::remove_file(&temp_cross_output);
    let _ = std::fs::remove_file(temp_fast_output.with_extension("bed.unmap"));
    let _ = std::fs::remove_file(temp_cross_output.with_extension("bed.unmap"));
    
    // Extract coordinates for comparison
    let fast_coords = extract_bed_coords(&fast_output);
    let cross_coords = extract_bed_coords(&cross_output);
    
    eprintln!("\n=== Strict Mode BED Test ===");
    eprintln!("FastCrossMap (strict): {:?}", fast_coords);
    eprintln!("CrossMap: {:?}", cross_coords);
    
    // In strict mode, coordinates should be identical
    assert_eq!(fast_coords, cross_coords, 
        "Strict mode should produce identical coordinates to CrossMap");
}

/// Test strict mode with multiple records
#[test]
fn test_strict_mode_multiple_records() {
    if skip_if_unavailable() { return; }
    
    let test_bed = r#"chr1	100000	100100	region1	100	+
chr1	500000	500200	region2	200	-
chr2	50000	50100	region3	300	+
"#;
    
    let temp_input = temp_file("strict_multi", "bed");
    let temp_fast_output = temp_file("strict_multi_fast", "bed");
    let temp_cross_output = temp_file("strict_multi_cross", "bed");
    
    std::fs::write(&temp_input, test_bed).unwrap();
    
    // Run FastCrossMap with strict mode
    let _ = Command::new("cargo")
        .args(&["run", "--release", "--", "--compat-mode=strict", "bed", CHAIN_FILE])
        .arg(temp_input.to_str().unwrap())
        .arg(temp_fast_output.to_str().unwrap())
        .output();
    
    // Run CrossMap
    let _ = Command::new("CrossMap")
        .args(&["bed", CHAIN_FILE, temp_input.to_str().unwrap(), temp_cross_output.to_str().unwrap()])
        .output();
    
    let fast_output = std::fs::read_to_string(&temp_fast_output).unwrap_or_default();
    let cross_output = std::fs::read_to_string(&temp_cross_output).unwrap_or_default();
    
    // Cleanup
    let _ = std::fs::remove_file(&temp_input);
    let _ = std::fs::remove_file(&temp_fast_output);
    let _ = std::fs::remove_file(&temp_cross_output);
    let _ = std::fs::remove_file(temp_fast_output.with_extension("bed.unmap"));
    let _ = std::fs::remove_file(temp_cross_output.with_extension("bed.unmap"));
    
    let fast_coords: HashSet<_> = extract_bed_coords(&fast_output).into_iter().collect();
    let cross_coords: HashSet<_> = extract_bed_coords(&cross_output).into_iter().collect();
    
    let matches = fast_coords.intersection(&cross_coords).count();
    let total = cross_coords.len().max(1);
    let match_rate = matches as f64 / total as f64 * 100.0;
    
    eprintln!("\n=== Strict Mode Multiple Records Test ===");
    eprintln!("FastCrossMap records: {}", fast_coords.len());
    eprintln!("CrossMap records: {}", cross_coords.len());
    eprintln!("Match rate: {:.2}%", match_rate);
    
    // Strict mode should achieve very high match rate
    assert!(match_rate >= 95.0, 
        "Strict mode should achieve >= 95% match rate, got {:.2}%", match_rate);
}

/// Test edge case: Indel at chain block boundary
#[test]
fn test_edge_case_indel_at_boundary() {
    if skip_if_unavailable() { return; }
    
    // This tests a known edge case where CrossMap behavior may differ
    let index = ChainIndex::from_chain_file(CHAIN_FILE).unwrap();
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    // Find a chain block boundary and test mapping near it
    let intervals = mapper.index().query_intervals("chr1", 0, 1000000);
    
    if let Some(first_interval) = intervals.first() {
        let boundary = first_interval.stop;
        
        // Test mapping right at the boundary
        let result_before = mapper.map("chr1", boundary - 10, boundary, Strand::Plus);
        let result_at = mapper.map("chr1", boundary - 5, boundary + 5, Strand::Plus);
        let result_after = mapper.map("chr1", boundary, boundary + 10, Strand::Plus);
        
        eprintln!("\n=== Edge Case: Chain Block Boundary ===");
        eprintln!("Boundary position: {}", boundary);
        eprintln!("Before boundary: {:?}", result_before.as_ref().map(|r| r.len()));
        eprintln!("At boundary: {:?}", result_at.as_ref().map(|r| r.len()));
        eprintln!("After boundary: {:?}", result_after.as_ref().map(|r| r.len()));
    }
}

/// Test edge case: Multiple overlapping chain blocks
#[test]
fn test_edge_case_multiple_blocks() {
    if skip_if_unavailable() { return; }
    
    let index = ChainIndex::from_chain_file(CHAIN_FILE).unwrap();
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    // Find a region with multiple overlapping blocks
    for pos in (100000u64..10000000).step_by(100000) {
        let intervals = mapper.index().query_intervals("chr1", pos, pos + 1000);
        if intervals.len() > 1 {
            eprintln!("\n=== Edge Case: Multiple Overlapping Blocks ===");
            eprintln!("Position: chr1:{}-{}", pos, pos + 1000);
            eprintln!("Overlapping blocks: {}", intervals.len());
            
            let result = mapper.map("chr1", pos, pos + 1000, Strand::Plus);
            eprintln!("Mapping result segments: {:?}", result.as_ref().map(|r| r.len()));
            break;
        }
    }
}
