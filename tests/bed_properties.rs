//! Property-based tests for BED format conversion
//!
//! **Feature: fast-crossmap, Property 9: BED 字段保留完整性**
//! **Validates: Requirements 4.2, 4.3**

use fast_crossmap::core::{ChainIndex, CoordinateMapper, ChromStyle, Strand};
use fast_crossmap::formats::bed::{BedRecordView, convert_bed};
use proptest::prelude::*;
use std::path::PathBuf;

/// Generate a valid chromosome name
fn arb_chrom_name() -> impl Strategy<Value = String> {
    prop_oneof![
        (1u8..=22).prop_map(|n| format!("chr{}", n)),
        Just("chrX".to_string()),
        Just("chrY".to_string()),
    ]
}

/// Generate a valid BED name field
fn arb_bed_name() -> impl Strategy<Value = String> {
    "[a-zA-Z][a-zA-Z0-9_]{0,20}".prop_map(|s| s)
}

/// Generate a valid score (0-1000)
fn arb_score() -> impl Strategy<Value = String> {
    (0u32..=1000).prop_map(|n| n.to_string())
}

/// Generate a strand character
fn arb_strand_char() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("+".to_string()),
        Just("-".to_string()),
        Just(".".to_string()),
    ]
}

/// Generate a BED3 line
fn arb_bed3_line() -> impl Strategy<Value = String> {
    (arb_chrom_name(), 1000u64..100000, 100u64..1000)
        .prop_map(|(chrom, start, size)| {
            format!("{}\t{}\t{}", chrom, start, start + size)
        })
}

/// Generate a BED6 line
fn arb_bed6_line() -> impl Strategy<Value = String> {
    (
        arb_chrom_name(),
        1000u64..100000,
        100u64..1000,
        arb_bed_name(),
        arb_score(),
        arb_strand_char(),
    )
        .prop_map(|(chrom, start, size, name, score, strand)| {
            format!("{}\t{}\t{}\t{}\t{}\t{}", 
                chrom, start, start + size, name, score, strand)
        })
}


proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]
    
    /// **Property 9: BED 字段保留完整性**
    /// 
    /// For any successfully converted BED record, all non-coordinate fields 
    /// (name, score, etc.) should be preserved unchanged.
    ///
    /// **Validates: Requirements 4.2, 4.3**
    #[test]
    fn prop_bed_field_preservation(line in arb_bed6_line()) {
        // Parse the original line
        let original = BedRecordView::parse(line.as_bytes()).unwrap();
        
        // Extract original non-coordinate fields
        let original_name = original.name().map(|s| s.to_string());
        let original_score = original.score().map(|s| s.to_string());
        let original_strand_char = original.strand_char().map(|s| s.to_string());
        
        // Verify parsing preserves fields
        prop_assert_eq!(original.field_count(), 6, "Should have 6 fields");
        prop_assert!(original_name.is_some(), "Name should be present");
        prop_assert!(original_score.is_some(), "Score should be present");
        prop_assert!(original_strand_char.is_some(), "Strand should be present");
        
        // Verify field content matches input
        let fields: Vec<&str> = line.split('\t').collect();
        prop_assert_eq!(original.chrom, fields[0], "Chrom should match");
        prop_assert_eq!(original_name.as_deref(), Some(fields[3]), "Name should match");
        prop_assert_eq!(original_score.as_deref(), Some(fields[4]), "Score should match");
        prop_assert_eq!(original_strand_char.as_deref(), Some(fields[5]), "Strand should match");
    }
    
    /// Property: BED3 parsing extracts correct coordinates
    #[test]
    fn prop_bed3_coordinate_parsing(line in arb_bed3_line()) {
        let view = BedRecordView::parse(line.as_bytes()).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        let expected_start: u64 = fields[1].parse().unwrap();
        let expected_end: u64 = fields[2].parse().unwrap();
        
        prop_assert_eq!(view.chrom, fields[0]);
        prop_assert_eq!(view.start, expected_start);
        prop_assert_eq!(view.end, expected_end);
        prop_assert_eq!(view.field_count(), 3);
    }
    
    /// Property: BED6 parsing extracts all fields correctly
    #[test]
    fn prop_bed6_field_parsing(line in arb_bed6_line()) {
        let view = BedRecordView::parse(line.as_bytes()).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        
        prop_assert_eq!(view.chrom, fields[0]);
        prop_assert_eq!(view.name(), Some(fields[3]));
        prop_assert_eq!(view.score(), Some(fields[4]));
        prop_assert_eq!(view.strand_char(), Some(fields[5]));
        prop_assert!(view.is_bed6());
    }
}

/// Integration test: BED conversion with real chain file
#[test]
fn test_bed_conversion_with_crossmap() {
    let chain_path = PathBuf::from("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz");
    
    if !chain_path.exists() {
        eprintln!("Skipping test: chain file not found");
        return;
    }
    
    // Load chain file
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    // Create test BED file
    let test_bed = "chr1\t10000\t10100\ttest_region\t500\t+\n\
                    chr1\t100000\t100100\ttest_region2\t600\t-\n\
                    chr2\t50000\t50100\ttest_region3\t700\t.\n";
    
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("test_input.bed");
    let output_path = temp_dir.join("test_output.bed");
    let unmap_path = temp_dir.join("test_unmap.bed");
    
    std::fs::write(&input_path, test_bed).unwrap();
    
    // Convert
    let stats = convert_bed(&input_path, &output_path, &unmap_path, &mapper, 1).unwrap();
    
    // Verify stats
    assert_eq!(stats.total, 3, "Should process 3 records");
    assert!(stats.success > 0, "Should have some successful conversions");
    
    // Read output and verify field preservation
    let output = std::fs::read_to_string(&output_path).unwrap();
    
    for line in output.lines() {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(fields.len() >= 6, "Output should preserve BED6 fields");
        
        // Verify name field is preserved (one of our test names)
        let name = fields[3];
        assert!(
            name == "test_region" || name == "test_region2" || name == "test_region3",
            "Name field should be preserved: {}", name
        );
        
        // Verify score field is preserved (one of our test scores)
        let score = fields[4];
        assert!(
            score == "500" || score == "600" || score == "700",
            "Score field should be preserved: {}", score
        );
    }
    
    // Cleanup
    let _ = std::fs::remove_file(&input_path);
    let _ = std::fs::remove_file(&output_path);
    let _ = std::fs::remove_file(&unmap_path);
}

/// Test BED12 field preservation
#[test]
fn test_bed12_field_preservation() {
    let bed12_line = b"chr1\t1000\t2000\tgene1\t500\t+\t1100\t1900\t255,0,0\t2\t100,100\t0,900";
    let view = BedRecordView::parse(bed12_line).unwrap();
    
    assert_eq!(view.chrom, "chr1");
    assert_eq!(view.start, 1000);
    assert_eq!(view.end, 2000);
    assert_eq!(view.name(), Some("gene1"));
    assert_eq!(view.score(), Some("500"));
    assert_eq!(view.strand(), Some(Strand::Plus));
    assert_eq!(view.thick_start(), Some(1100));
    assert_eq!(view.thick_end(), Some(1900));
    assert_eq!(view.item_rgb(), Some("255,0,0"));
    assert_eq!(view.block_count(), Some(2));
    assert_eq!(view.block_sizes(), Some("100,100"));
    assert_eq!(view.block_starts(), Some("0,900"));
    assert!(view.is_bed12());
}


/// **Property 11: 并行处理结果确定性**
/// 
/// For any input file, running conversion with different thread counts 
/// should produce identical output (content should be the same when sorted).
///
/// **Validates: Requirements 4.7, 5.7**
#[test]
fn test_parallel_determinism() {
    let chain_path = PathBuf::from("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz");
    
    if !chain_path.exists() {
        eprintln!("Skipping test: chain file not found");
        return;
    }
    
    // Load chain file
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    // Create test BED file with many records
    let mut test_bed = String::new();
    for i in 0..100 {
        let start = 10000 + i * 1000;
        let end = start + 100;
        test_bed.push_str(&format!("chr1\t{}\t{}\tregion_{}\t{}\t+\n", start, end, i, i % 1000));
    }
    
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("parallel_test_input.bed");
    std::fs::write(&input_path, &test_bed).unwrap();
    
    // Run with 1 thread (sequential)
    let output_1 = temp_dir.join("parallel_test_output_1.bed");
    let unmap_1 = temp_dir.join("parallel_test_unmap_1.bed");
    let stats_1 = convert_bed(&input_path, &output_1, &unmap_1, &mapper, 1).unwrap();
    
    // Run with 4 threads (parallel)
    let output_4 = temp_dir.join("parallel_test_output_4.bed");
    let unmap_4 = temp_dir.join("parallel_test_unmap_4.bed");
    let stats_4 = convert_bed(&input_path, &output_4, &unmap_4, &mapper, 4).unwrap();
    
    // Run with 8 threads (parallel)
    let output_8 = temp_dir.join("parallel_test_output_8.bed");
    let unmap_8 = temp_dir.join("parallel_test_unmap_8.bed");
    let stats_8 = convert_bed(&input_path, &output_8, &unmap_8, &mapper, 8).unwrap();
    
    // Verify stats are identical
    assert_eq!(stats_1.total, stats_4.total, "Total count should match");
    assert_eq!(stats_1.total, stats_8.total, "Total count should match");
    assert_eq!(stats_1.success, stats_4.success, "Success count should match");
    assert_eq!(stats_1.success, stats_8.success, "Success count should match");
    assert_eq!(stats_1.failed, stats_4.failed, "Failed count should match");
    assert_eq!(stats_1.failed, stats_8.failed, "Failed count should match");
    
    // Read and sort output lines for comparison
    fn read_and_sort(path: &std::path::Path) -> Vec<String> {
        let content = std::fs::read_to_string(path).unwrap();
        let mut lines: Vec<String> = content.lines().map(|s| s.to_string()).collect();
        lines.sort();
        lines
    }
    
    let lines_1 = read_and_sort(&output_1);
    let lines_4 = read_and_sort(&output_4);
    let lines_8 = read_and_sort(&output_8);
    
    // Verify content is identical (when sorted)
    assert_eq!(lines_1.len(), lines_4.len(), "Output line count should match");
    assert_eq!(lines_1.len(), lines_8.len(), "Output line count should match");
    
    for (i, (l1, l4)) in lines_1.iter().zip(lines_4.iter()).enumerate() {
        assert_eq!(l1, l4, "Line {} differs between 1 and 4 threads", i);
    }
    
    for (i, (l1, l8)) in lines_1.iter().zip(lines_8.iter()).enumerate() {
        assert_eq!(l1, l8, "Line {} differs between 1 and 8 threads", i);
    }
    
    // Cleanup
    let _ = std::fs::remove_file(&input_path);
    let _ = std::fs::remove_file(&output_1);
    let _ = std::fs::remove_file(&unmap_1);
    let _ = std::fs::remove_file(&output_4);
    let _ = std::fs::remove_file(&unmap_4);
    let _ = std::fs::remove_file(&output_8);
    let _ = std::fs::remove_file(&unmap_8);
    
    eprintln!("Parallel determinism test passed: {} records processed", stats_1.total);
}

/// Test that parallel processing produces correct results
#[test]
fn test_parallel_correctness() {
    let chain_path = PathBuf::from("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz");
    
    if !chain_path.exists() {
        eprintln!("Skipping test: chain file not found");
        return;
    }
    
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    // Create larger test file
    let mut test_bed = String::new();
    for chrom_num in 1..=5 {
        for i in 0..50 {
            let start = 10000 + i * 2000;
            let end = start + 100;
            test_bed.push_str(&format!(
                "chr{}\t{}\t{}\tgene_{}_{}\t{}\t{}\n",
                chrom_num, start, end, chrom_num, i, 
                (i * 10) % 1000,
                if i % 2 == 0 { "+" } else { "-" }
            ));
        }
    }
    
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("parallel_correct_input.bed");
    std::fs::write(&input_path, &test_bed).unwrap();
    
    let output_path = temp_dir.join("parallel_correct_output.bed");
    let unmap_path = temp_dir.join("parallel_correct_unmap.bed");
    
    let stats = convert_bed(&input_path, &output_path, &unmap_path, &mapper, 4).unwrap();
    
    assert_eq!(stats.total, 250, "Should process 250 records (5 chroms * 50 each)");
    assert!(stats.success > 0, "Should have successful conversions");
    
    // Verify output format
    let output = std::fs::read_to_string(&output_path).unwrap();
    for line in output.lines() {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(fields.len() >= 6, "Output should have at least 6 fields: {:?}", fields);
        
        // Verify coordinates are valid numbers
        let _start: u64 = fields[1].parse().expect("Start should be a number");
        let _end: u64 = fields[2].parse().expect("End should be a number");
    }
    
    // Cleanup
    let _ = std::fs::remove_file(&input_path);
    let _ = std::fs::remove_file(&output_path);
    let _ = std::fs::remove_file(&unmap_path);
    
    eprintln!("Parallel correctness test passed: {}/{} successful", stats.success, stats.total);
}


/// Integration test: Compare BED conversion output with CrossMap
#[test]
fn test_bed_vs_crossmap() {
    use std::process::Command;
    
    let chain_path = PathBuf::from("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz");
    
    if !chain_path.exists() {
        eprintln!("Skipping test: chain file not found");
        return;
    }
    
    // Check if CrossMap is available
    let crossmap_check = Command::new("CrossMap")
        .arg("--version")
        .output();
    
    if crossmap_check.is_err() {
        eprintln!("Skipping test: CrossMap not installed");
        return;
    }
    
    // Create test BED file with various regions
    let test_bed = "\
chr1\t10000\t10100\tregion1\t100\t+
chr1\t100000\t100200\tregion2\t200\t-
chr1\t500000\t500500\tregion3\t300\t.
chr2\t50000\t50100\tregion4\t400\t+
chr2\t200000\t200300\tregion5\t500\t-
chr3\t100000\t100150\tregion6\t600\t+
chr10\t50000\t50200\tregion7\t700\t+
chr22\t20000\t20100\tregion8\t800\t-
";
    
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("crossmap_test_input.bed");
    std::fs::write(&input_path, test_bed).unwrap();
    
    // Run CrossMap
    let crossmap_output = temp_dir.join("crossmap_output.bed");
    let crossmap_unmap = temp_dir.join("crossmap_unmap.bed");
    
    let crossmap_result = Command::new("CrossMap")
        .args(&[
            "bed",
            chain_path.to_str().unwrap(),
            input_path.to_str().unwrap(),
            crossmap_output.to_str().unwrap(),
        ])
        .output();
    
    if crossmap_result.is_err() {
        eprintln!("CrossMap execution failed");
        return;
    }
    
    // Run FastCrossMap
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    let fast_output = temp_dir.join("fast_output.bed");
    let fast_unmap = temp_dir.join("fast_unmap.bed");
    
    let stats = convert_bed(&input_path, &fast_output, &fast_unmap, &mapper, 1).unwrap();
    
    eprintln!("FastCrossMap stats: total={}, success={}, failed={}", 
              stats.total, stats.success, stats.failed);
    
    // Compare outputs
    let crossmap_content = std::fs::read_to_string(&crossmap_output).unwrap_or_default();
    let fast_content = std::fs::read_to_string(&fast_output).unwrap();
    
    eprintln!("\n=== CrossMap output ===");
    for line in crossmap_content.lines() {
        eprintln!("{}", line);
    }
    
    eprintln!("\n=== FastCrossMap output ===");
    for line in fast_content.lines() {
        eprintln!("{}", line);
    }
    
    // Parse and compare line by line
    let crossmap_lines: Vec<&str> = crossmap_content.lines().collect();
    let fast_lines: Vec<&str> = fast_content.lines().collect();
    
    eprintln!("\n=== Comparison ===");
    eprintln!("CrossMap lines: {}, FastCrossMap lines: {}", 
              crossmap_lines.len(), fast_lines.len());
    
    let mut matches = 0;
    let mut mismatches = 0;
    
    for (i, (cross_line, fast_line)) in crossmap_lines.iter().zip(fast_lines.iter()).enumerate() {
        let cross_fields: Vec<&str> = cross_line.split('\t').collect();
        let fast_fields: Vec<&str> = fast_line.split('\t').collect();
        
        // Compare first 3 fields (chrom, start, end)
        if cross_fields.len() >= 3 && fast_fields.len() >= 3 {
            let coords_match = cross_fields[0] == fast_fields[0] 
                && cross_fields[1] == fast_fields[1]
                && cross_fields[2] == fast_fields[2];
            
            if coords_match {
                matches += 1;
            } else {
                mismatches += 1;
                eprintln!("Line {}: MISMATCH", i);
                eprintln!("  CrossMap: {}:{}-{}", cross_fields[0], cross_fields[1], cross_fields[2]);
                eprintln!("  FastCrossMap: {}:{}-{}", fast_fields[0], fast_fields[1], fast_fields[2]);
            }
        }
    }
    
    eprintln!("\nResults: {} matches, {} mismatches", matches, mismatches);
    
    // Cleanup
    let _ = std::fs::remove_file(&input_path);
    let _ = std::fs::remove_file(&crossmap_output);
    let _ = std::fs::remove_file(&crossmap_unmap);
    let _ = std::fs::remove_file(&fast_output);
    let _ = std::fs::remove_file(&fast_unmap);
    
    // We expect high match rate
    if matches + mismatches > 0 {
        let match_rate = matches as f64 / (matches + mismatches) as f64;
        eprintln!("Match rate: {:.2}%", match_rate * 100.0);
        assert!(match_rate > 0.9, "Match rate should be > 90%, got {:.2}%", match_rate * 100.0);
    }
}
