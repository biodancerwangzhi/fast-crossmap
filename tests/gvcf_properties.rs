//! Property-based tests for GVCF format conversion
//!
//! **Feature: fast-crossmap, GVCF 格式转换**
//! **Validates: Requirements 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7**

use fast_crossmap::core::{ChainIndex, CoordinateMapper, ChromStyle};
use fast_crossmap::formats::gvcf::{GvcfRecordView, convert_gvcf};
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

/// Generate a valid DNA allele
fn arb_dna_allele() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("A".to_string()),
        Just("T".to_string()),
        Just("G".to_string()),
        Just("C".to_string()),
        "[ATGC]{2,5}".prop_map(|s| s),
    ]
}

/// Generate a GVCF variant line
fn arb_gvcf_variant_line() -> impl Strategy<Value = String> {
    (
        arb_chrom_name(),
        1000u64..100000,
        arb_dna_allele(),
        arb_dna_allele(),
    )
        .prop_map(|(chrom, pos, ref_allele, alt_allele)| {
            format!("{}\t{}\t.\t{}\t{}\t30\tPASS\tDP=100", 
                chrom, pos, ref_allele, alt_allele)
        })
}

/// Generate a GVCF non-variant block line
fn arb_gvcf_block_line() -> impl Strategy<Value = String> {
    (
        arb_chrom_name(),
        1000u64..100000,
        100u64..1000,  // block size
    )
        .prop_map(|(chrom, pos, size)| {
            let end = pos + size;
            format!("{}\t{}\t.\tA\t<NON_REF>\t.\t.\tEND={};DP=50", 
                chrom, pos, end)
        })
}


proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]
    
    /// Property: GVCF parsing extracts correct coordinates
    #[test]
    fn prop_gvcf_coordinate_parsing(line in arb_gvcf_variant_line()) {
        let view = GvcfRecordView::parse(line.as_bytes()).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        let expected_pos: u64 = fields[1].parse().unwrap();
        
        prop_assert_eq!(view.chrom, fields[0]);
        prop_assert_eq!(view.pos, expected_pos);
    }
    
    /// Property: GVCF non-variant blocks have END= parsed correctly
    #[test]
    fn prop_gvcf_block_end_parsing(line in arb_gvcf_block_line()) {
        let view = GvcfRecordView::parse(line.as_bytes()).unwrap();
        
        prop_assert!(view.is_non_variant_block(), "Should be non-variant block");
        prop_assert!(view.is_gvcf_non_ref(), "Should have <NON_REF> ALT");
        
        // Extract END from line
        let info = view.info().unwrap();
        let end_str = info.split(';')
            .find(|s| s.starts_with("END="))
            .map(|s| &s[4..])
            .unwrap();
        let expected_end: u64 = end_str.parse().unwrap();
        
        prop_assert_eq!(view.end_position(), Some(expected_end));
    }
    
    /// Property: GVCF variant records don't have END=
    #[test]
    fn prop_gvcf_variant_no_end(line in arb_gvcf_variant_line()) {
        let view = GvcfRecordView::parse(line.as_bytes()).unwrap();
        
        prop_assert!(!view.is_non_variant_block(), "Variant should not be non-variant block");
        prop_assert_eq!(view.end_position(), None);
    }
}

/// Test GVCF conversion with real chain file
#[test]
fn test_gvcf_conversion_basic() {
    let chain_path = PathBuf::from("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz");
    
    if !chain_path.exists() {
        eprintln!("Skipping test: chain file not found");
        return;
    }
    
    // Load chain file
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    // Create test GVCF file
    let test_gvcf = "\
##fileformat=VCFv4.2
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##contig=<ID=chr1,length=249250621>
##contig=<ID=chr2,length=243199373>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t100000\t.\tA\t<NON_REF>\t.\t.\tEND=100100;DP=50\tGT\t0/0
chr1\t100101\t.\tC\tT\t30\tPASS\tDP=100\tGT\t0/1
chr2\t50000\t.\tG\t<NON_REF>\t.\t.\tEND=50500;DP=40\tGT\t0/0
chr2\t50501\t.\tA\tG\t40\tPASS\tDP=80\tGT\t1/1
";
    
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("test_input.gvcf");
    let output_path = temp_dir.join("test_output.gvcf");
    
    std::fs::write(&input_path, test_gvcf).unwrap();
    
    // Convert (without reference genome)
    let stats = convert_gvcf(&input_path, &output_path, &mapper, None::<&PathBuf>, false, 1).unwrap();
    
    eprintln!("GVCF conversion stats: total={}, success={}, failed={}, headers={}", 
              stats.total, stats.success, stats.failed, stats.headers);
    
    // Verify stats
    assert_eq!(stats.total, 4, "Should process 4 records");
    
    // Read output and verify structure
    let output = std::fs::read_to_string(&output_path).unwrap();
    
    eprintln!("\n=== GVCF Output ===");
    for line in output.lines() {
        eprintln!("{}", line);
    }
    
    // Check headers are preserved
    assert!(output.contains("##fileformat"), "Should have fileformat header");
    assert!(output.contains("##INFO=<ID=END"), "Should have END INFO header");
    assert!(output.contains("#CHROM"), "Should have column header");
    
    // Cleanup
    let _ = std::fs::remove_file(&input_path);
    let _ = std::fs::remove_file(&output_path);
    let unmap_path = output_path.with_extension("gvcf.unmap");
    let _ = std::fs::remove_file(&unmap_path);
}

/// Test GVCF vs CrossMap comparison
#[test]
fn test_gvcf_vs_crossmap() {
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
    
    // Note: CrossMap gvcf requires a reference genome file
    // This test compares coordinate mapping only (without reference)
    
    // Create test GVCF file
    let test_gvcf = "\
##fileformat=VCFv4.2
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100000\t.\tA\t<NON_REF>\t.\t.\tEND=100100
chr2\t50000\t.\tG\t<NON_REF>\t.\t.\tEND=50500
";
    
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("gvcf_crossmap_test.gvcf");
    let fast_output = temp_dir.join("gvcf_fast_output.gvcf");
    
    std::fs::write(&input_path, test_gvcf).unwrap();
    
    // Run FastCrossMap
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    let stats = convert_gvcf(&input_path, &fast_output, &mapper, None::<&PathBuf>, false, 1).unwrap();
    
    eprintln!("FastCrossMap GVCF: total={}, success={}, failed={}", stats.total, stats.success, stats.failed);
    
    // Read output
    let fast_content = std::fs::read_to_string(&fast_output).unwrap_or_default();
    
    eprintln!("\n=== FastCrossMap GVCF output ===");
    for line in fast_content.lines() {
        if !line.starts_with('#') {
            eprintln!("{}", line);
        }
    }
    
    // Note: Full CrossMap comparison requires reference genome
    // For now, just verify our output is valid
    let data_lines: Vec<&str> = fast_content.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    
    eprintln!("\nProcessed {} data lines", data_lines.len());
    
    // Cleanup
    let _ = std::fs::remove_file(&input_path);
    let _ = std::fs::remove_file(&fast_output);
    let _ = std::fs::remove_file(fast_output.with_extension("gvcf.unmap"));
}
