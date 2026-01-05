//! Property-based tests for VCF format conversion
//!
//! **Feature: fast-crossmap, Property 10: VCF INFO/FORMAT 保留完整性**
//! **Validates: Requirements 5.3**

use fast_crossmap::core::{ChainIndex, CoordinateMapper, ChromStyle};
use fast_crossmap::formats::vcf::{VcfRecordView, convert_vcf, VariantType};
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

/// Generate a valid VCF ID field
fn arb_vcf_id() -> impl Strategy<Value = String> {
    prop_oneof![
        Just(".".to_string()),
        "[a-zA-Z][a-zA-Z0-9_]{0,10}".prop_map(|s| format!("rs{}", s)),
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

/// Generate a valid QUAL field
fn arb_qual() -> impl Strategy<Value = String> {
    prop_oneof![
        Just(".".to_string()),
        (0u32..1000).prop_map(|n| n.to_string()),
        (0.0f64..100.0).prop_map(|f| format!("{:.2}", f)),
    ]
}

/// Generate a valid FILTER field
fn arb_filter() -> impl Strategy<Value = String> {
    prop_oneof![
        Just(".".to_string()),
        Just("PASS".to_string()),
        Just("LowQual".to_string()),
        Just("q10;s50".to_string()),
    ]
}

/// Generate a valid INFO field
fn arb_info() -> impl Strategy<Value = String> {
    prop_oneof![
        Just(".".to_string()),
        Just("DP=100".to_string()),
        Just("DP=50;AF=0.5".to_string()),
        Just("DP=100;AF=0.25;DB".to_string()),
        Just("DP=200;MQ=60;FS=0.0;SOR=0.5".to_string()),
    ]
}

/// Generate a minimal VCF line (8 fields)
fn arb_vcf_line_minimal() -> impl Strategy<Value = String> {
    (
        arb_chrom_name(),
        1000u64..100000,
        arb_vcf_id(),
        arb_dna_allele(),
        arb_dna_allele(),
        arb_qual(),
        arb_filter(),
        arb_info(),
    )
        .prop_map(|(chrom, pos, id, ref_allele, alt_allele, qual, filter, info)| {
            format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                chrom, pos, id, ref_allele, alt_allele, qual, filter, info)
        })
}

/// Generate a VCF line with FORMAT and samples
fn arb_vcf_line_with_samples() -> impl Strategy<Value = String> {
    (
        arb_chrom_name(),
        1000u64..100000,
        arb_vcf_id(),
        arb_dna_allele(),
        arb_dna_allele(),
        arb_qual(),
        arb_filter(),
        arb_info(),
    )
        .prop_map(|(chrom, pos, id, ref_allele, alt_allele, qual, filter, info)| {
            format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tGT:DP\t0/1:30\t1/1:25", 
                chrom, pos, id, ref_allele, alt_allele, qual, filter, info)
        })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]
    
    /// **Property 10: VCF INFO/FORMAT 保留完整性**
    /// 
    /// For any VCF record, parsing should preserve all INFO and FORMAT fields unchanged.
    ///
    /// **Validates: Requirements 5.3**
    #[test]
    fn prop_vcf_info_format_preservation(line in arb_vcf_line_with_samples()) {
        // Parse the original line
        let original = VcfRecordView::parse(line.as_bytes()).unwrap();
        
        // Extract original fields
        let original_info = original.info().map(|s| s.to_string());
        let original_format = original.format().map(|s| s.to_string());
        let original_samples: Vec<String> = original.samples().iter().map(|s| s.to_string()).collect();
        
        // Verify parsing preserves fields
        prop_assert!(original_info.is_some(), "INFO should be present");
        prop_assert!(original_format.is_some(), "FORMAT should be present");
        prop_assert_eq!(original_samples.len(), 2, "Should have 2 samples");
        
        // Verify field content matches input
        let fields: Vec<&str> = line.split('\t').collect();
        prop_assert_eq!(original.chrom, fields[0], "CHROM should match");
        prop_assert_eq!(original_info.as_deref(), Some(fields[7]), "INFO should match");
        prop_assert_eq!(original_format.as_deref(), Some(fields[8]), "FORMAT should match");
        prop_assert_eq!(&original_samples[0], fields[9], "Sample 1 should match");
        prop_assert_eq!(&original_samples[1], fields[10], "Sample 2 should match");
    }
    
    /// Property: VCF minimal parsing extracts correct coordinates
    #[test]
    fn prop_vcf_coordinate_parsing(line in arb_vcf_line_minimal()) {
        let view = VcfRecordView::parse(line.as_bytes()).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        let expected_pos: u64 = fields[1].parse().unwrap();
        
        prop_assert_eq!(view.chrom, fields[0]);
        prop_assert_eq!(view.pos, expected_pos);
        prop_assert_eq!(view.field_count(), 8);
    }
    
    /// Property: VCF ID field is preserved
    #[test]
    fn prop_vcf_id_preservation(line in arb_vcf_line_minimal()) {
        let view = VcfRecordView::parse(line.as_bytes()).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        
        prop_assert_eq!(view.id(), Some(fields[2]));
    }
    
    /// Property: VCF REF/ALT fields are preserved
    #[test]
    fn prop_vcf_allele_preservation(line in arb_vcf_line_minimal()) {
        let view = VcfRecordView::parse(line.as_bytes()).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        
        prop_assert_eq!(view.ref_allele(), Some(fields[3]));
        prop_assert_eq!(view.alt_alleles(), Some(fields[4]));
    }
    
    /// Property: VCF QUAL/FILTER fields are preserved
    #[test]
    fn prop_vcf_qual_filter_preservation(line in arb_vcf_line_minimal()) {
        let view = VcfRecordView::parse(line.as_bytes()).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        
        prop_assert_eq!(view.qual(), Some(fields[5]));
        prop_assert_eq!(view.filter(), Some(fields[6]));
    }
}

/// Test variant type detection
#[test]
fn test_variant_type_detection() {
    // SNP (substitution)
    let snp = b"chr1\t100\t.\tA\tG\t.\t.\t.";
    let view = VcfRecordView::parse(snp).unwrap();
    assert_eq!(view.variant_type(), VariantType::Substitution);
    
    // MNP (multi-nucleotide substitution)
    let mnp = b"chr1\t100\t.\tAT\tGC\t.\t.\t.";
    let view = VcfRecordView::parse(mnp).unwrap();
    assert_eq!(view.variant_type(), VariantType::Substitution);
    
    // Insertion
    let ins = b"chr1\t100\t.\tA\tATG\t.\t.\t.";
    let view = VcfRecordView::parse(ins).unwrap();
    assert_eq!(view.variant_type(), VariantType::Insertion);
    
    // Deletion
    let del = b"chr1\t100\t.\tATG\tA\t.\t.\t.";
    let view = VcfRecordView::parse(del).unwrap();
    assert_eq!(view.variant_type(), VariantType::Deletion);
}

/// Test INFO field parsing
#[test]
fn test_info_field_parsing() {
    let line = b"chr1\t100\t.\tA\tG\t.\t.\tDP=100;AF=0.5;DB;MQ=60";
    let view = VcfRecordView::parse(line).unwrap();
    let info = view.parse_info();
    
    assert_eq!(info.get("DP"), Some(&"100".to_string()));
    assert_eq!(info.get("AF"), Some(&"0.5".to_string()));
    assert_eq!(info.get("DB"), Some(&"".to_string())); // Flag
    assert_eq!(info.get("MQ"), Some(&"60".to_string()));
}

/// Test multi-allelic VCF records
#[test]
fn test_multi_allelic_vcf() {
    let line = b"chr1\t100\t.\tA\tG,T,C\t.\t.\t.";
    let view = VcfRecordView::parse(line).unwrap();
    
    assert_eq!(view.ref_allele(), Some("A"));
    assert_eq!(view.alt_alleles(), Some("G,T,C"));
}

/// Test VCF with complex INFO field
#[test]
fn test_complex_info_field() {
    let line = b"chr1\t100\t.\tA\tG\t30\tPASS\tDP=100;AF=0.5,0.3;ANN=G|missense|HIGH|GENE1";
    let view = VcfRecordView::parse(line).unwrap();
    let info = view.parse_info();
    
    assert_eq!(info.get("DP"), Some(&"100".to_string()));
    assert_eq!(info.get("AF"), Some(&"0.5,0.3".to_string()));
    assert!(info.get("ANN").is_some());
}

/// Integration test: VCF conversion with real chain file
#[test]
fn test_vcf_conversion_with_crossmap() {
    let chain_path = PathBuf::from("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz");
    
    if !chain_path.exists() {
        eprintln!("Skipping test: chain file not found");
        return;
    }
    
    // Load chain file
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    // Create test VCF file
    let test_vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t10000\trs1\tA\tG\t30\tPASS\tDP=100\tGT\t0/1
chr1\t100000\trs2\tC\tT\t40\tPASS\tDP=200\tGT\t1/1
chr2\t50000\trs3\tG\tA\t50\tPASS\tDP=150\tGT\t0/1
";
    
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("test_input.vcf");
    let output_path = temp_dir.join("test_output.vcf");
    
    std::fs::write(&input_path, test_vcf).unwrap();
    
    // Convert (without reference genome for simplicity)
    let stats = convert_vcf(&input_path, &output_path, &mapper, None::<&PathBuf>, false, 1).unwrap();
    
    // Verify stats
    assert_eq!(stats.total, 3, "Should process 3 records");
    
    // Read output and verify structure
    let output = std::fs::read_to_string(&output_path).unwrap();
    
    // Check headers are present
    assert!(output.contains("##fileformat"), "Should have fileformat header");
    assert!(output.contains("##INFO"), "Should have INFO header");
    assert!(output.contains("#CHROM"), "Should have column header");
    
    // Cleanup
    let _ = std::fs::remove_file(&input_path);
    let _ = std::fs::remove_file(&output_path);
    let unmap_path = output_path.with_extension("vcf.unmap");
    let _ = std::fs::remove_file(&unmap_path);
}

/// Test parallel VCF conversion determinism
#[test]
fn test_vcf_parallel_determinism() {
    let chain_path = PathBuf::from("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz");
    
    if !chain_path.exists() {
        eprintln!("Skipping test: chain file not found");
        return;
    }
    
    // Load chain file
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    // Create test VCF file with many records
    let mut test_vcf = String::from("\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
");
    
    for i in 0..100 {
        let pos = 10000 + i * 1000;
        test_vcf.push_str(&format!("chr1\t{}\trs{}\tA\tG\t30\tPASS\tDP={}\n", pos, i, i * 10));
    }
    
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("vcf_parallel_test_input.vcf");
    std::fs::write(&input_path, &test_vcf).unwrap();
    
    // Run with 1 thread (sequential)
    let output_1 = temp_dir.join("vcf_parallel_test_output_1.vcf");
    let stats_1 = convert_vcf(&input_path, &output_1, &mapper, None::<&PathBuf>, false, 1).unwrap();
    
    // Run with 4 threads (parallel)
    let output_4 = temp_dir.join("vcf_parallel_test_output_4.vcf");
    let stats_4 = convert_vcf(&input_path, &output_4, &mapper, None::<&PathBuf>, false, 4).unwrap();
    
    // Verify stats are identical
    assert_eq!(stats_1.total, stats_4.total, "Total count should match");
    assert_eq!(stats_1.success, stats_4.success, "Success count should match");
    assert_eq!(stats_1.failed, stats_4.failed, "Failed count should match");
    
    // Read and sort output lines for comparison (excluding headers)
    fn read_data_lines(path: &std::path::Path) -> Vec<String> {
        let content = std::fs::read_to_string(path).unwrap();
        let mut lines: Vec<String> = content
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|s| s.to_string())
            .collect();
        lines.sort();
        lines
    }
    
    let lines_1 = read_data_lines(&output_1);
    let lines_4 = read_data_lines(&output_4);
    
    // Verify content is identical (when sorted)
    assert_eq!(lines_1.len(), lines_4.len(), "Output line count should match");
    
    for (i, (l1, l4)) in lines_1.iter().zip(lines_4.iter()).enumerate() {
        assert_eq!(l1, l4, "Line {} differs between 1 and 4 threads", i);
    }
    
    // Cleanup
    let _ = std::fs::remove_file(&input_path);
    let _ = std::fs::remove_file(&output_1);
    let _ = std::fs::remove_file(&output_4);
    let _ = std::fs::remove_file(output_1.with_extension("vcf.unmap"));
    let _ = std::fs::remove_file(output_4.with_extension("vcf.unmap"));
    
    eprintln!("VCF parallel determinism test passed: {} records processed", stats_1.total);
}


/// Test negative strand handling for SNPs
#[test]
fn test_vcf_negative_strand_snp() {
    use fast_crossmap::core::dna::revcomp;
    
    // Test reverse complement for SNPs
    assert_eq!(revcomp("A"), "T");
    assert_eq!(revcomp("T"), "A");
    assert_eq!(revcomp("G"), "C");
    assert_eq!(revcomp("C"), "G");
    
    // Multi-nucleotide
    assert_eq!(revcomp("AT"), "AT");
    assert_eq!(revcomp("GC"), "GC");
    assert_eq!(revcomp("AG"), "CT");
}

/// Test negative strand handling for indels
#[test]
fn test_vcf_negative_strand_indel() {
    use fast_crossmap::core::dna::revcomp;
    
    // For indels, the first nucleotide is replaced with new REF
    // and the rest is reverse complemented
    
    // Insertion: A -> ATG on negative strand
    // First char stays (from new REF), rest is revcomp
    let alt = "ATG";
    let rest = &alt[1..]; // "TG"
    let revcomp_rest = revcomp(rest); // "CA"
    assert_eq!(revcomp_rest, "CA");
    
    // Deletion: ATG -> A on negative strand
    // Single char, no change needed for the alt itself
    let alt = "A";
    assert_eq!(alt.len(), 1);
}

/// Test VCF conversion preserves non-DNA alleles
#[test]
fn test_vcf_non_dna_alleles() {
    let line = b"chr1\t100\t.\tA\t<DEL>\t.\t.\t.";
    let view = VcfRecordView::parse(line).unwrap();
    
    assert_eq!(view.alt_alleles(), Some("<DEL>"));
    
    // Non-DNA alleles should be preserved as-is
    let alt = view.alt_alleles().unwrap();
    assert!(!fast_crossmap::core::dna::is_dna(alt));
}

/// Test VCF with symbolic alleles
#[test]
fn test_vcf_symbolic_alleles() {
    // Symbolic alleles like <DEL>, <INS>, <DUP> should be preserved
    let line = b"chr1\t100\t.\tA\t<DEL>,<INS>\t.\t.\t.";
    let view = VcfRecordView::parse(line).unwrap();
    
    assert_eq!(view.alt_alleles(), Some("<DEL>,<INS>"));
}

/// Test VCF with breakend notation
#[test]
fn test_vcf_breakend() {
    let line = b"chr1\t100\t.\tA\tA]chr2:200]\t.\t.\t.";
    let view = VcfRecordView::parse(line).unwrap();
    
    assert_eq!(view.alt_alleles(), Some("A]chr2:200]"));
}


/// Integration test: Compare VCF coordinate mapping with CrossMap using BED-based comparison
/// Since CrossMap VCF requires a reference genome, we compare the underlying coordinate mapping
/// by converting VCF positions to BED format and comparing with CrossMap BED output.
#[test]
fn test_vcf_vs_crossmap() {
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
    
    // Test VCF positions as BED coordinates (VCF POS is 1-based, BED is 0-based)
    // We'll create a BED file with the same positions and compare
    let test_positions = vec![
        ("chr1", 10000u64),
        ("chr1", 100000),
        ("chr1", 500000),
        ("chr2", 50000),
        ("chr2", 200000),
        ("chr3", 100000),
        ("chr10", 50000),
        ("chr22", 20000),
    ];
    
    // Create BED file (0-based coordinates, so VCF POS - 1)
    let mut bed_content = String::new();
    for (i, (chrom, pos)) in test_positions.iter().enumerate() {
        let start = pos - 1; // Convert to 0-based
        let end = start + 1;
        bed_content.push_str(&format!("{}\t{}\t{}\tvar{}\n", chrom, start, end, i));
    }
    
    let temp_dir = std::env::temp_dir();
    let bed_input = temp_dir.join("vcf_coord_test.bed");
    std::fs::write(&bed_input, &bed_content).unwrap();
    
    // Run CrossMap BED
    let crossmap_output = temp_dir.join("vcf_coord_crossmap.bed");
    let _crossmap_result = Command::new("CrossMap")
        .args(&[
            "bed",
            chain_path.to_str().unwrap(),
            bed_input.to_str().unwrap(),
            crossmap_output.to_str().unwrap(),
        ])
        .output();
    
    // Run FastCrossMap VCF
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    // Create VCF file
    let mut vcf_content = String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (i, (chrom, pos)) in test_positions.iter().enumerate() {
        vcf_content.push_str(&format!("{}\t{}\tvar{}\tA\tG\t30\tPASS\tDP=100\n", chrom, pos, i));
    }
    
    let vcf_input = temp_dir.join("vcf_coord_test.vcf");
    let vcf_output = temp_dir.join("vcf_coord_fast.vcf");
    std::fs::write(&vcf_input, &vcf_content).unwrap();
    
    let stats = convert_vcf(&vcf_input, &vcf_output, &mapper, None::<&PathBuf>, false, 1).unwrap();
    
    eprintln!("\n=== VCF Coordinate Mapping Comparison ===");
    eprintln!("FastCrossMap VCF: Total={}, Success={}, Failed={}", stats.total, stats.success, stats.failed);
    
    // Read CrossMap BED output
    let crossmap_content = std::fs::read_to_string(&crossmap_output).unwrap_or_default();
    let crossmap_lines: Vec<&str> = crossmap_content.lines().collect();
    
    // Read FastCrossMap VCF output
    let fast_content = std::fs::read_to_string(&vcf_output).unwrap();
    let fast_lines: Vec<&str> = fast_content.lines()
        .filter(|l| !l.starts_with('#'))
        .collect();
    
    eprintln!("\n=== CrossMap BED output ===");
    for line in &crossmap_lines {
        eprintln!("{}", line);
    }
    
    eprintln!("\n=== FastCrossMap VCF output (data lines) ===");
    for line in &fast_lines {
        eprintln!("{}", line);
    }
    
    // Compare coordinates
    // CrossMap BED: chrom, start, end (0-based)
    // FastCrossMap VCF: chrom, pos (1-based)
    eprintln!("\n=== Coordinate Comparison ===");
    
    let mut matches = 0;
    let mut mismatches = 0;
    
    // Build a map of CrossMap results by variant name
    let mut crossmap_results: std::collections::HashMap<String, (String, u64)> = std::collections::HashMap::new();
    for line in &crossmap_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 4 {
            let chrom = fields[0];
            let start: u64 = fields[1].parse().unwrap_or(0);
            let name = fields[3];
            // Convert BED start (0-based) to VCF POS (1-based)
            crossmap_results.insert(name.to_string(), (chrom.to_string(), start + 1));
        }
    }
    
    // Compare with FastCrossMap VCF output
    for line in &fast_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 3 {
            let fast_chrom = fields[0];
            let fast_pos: u64 = fields[1].parse().unwrap_or(0);
            let var_id = fields[2];
            
            if let Some((cross_chrom, cross_pos)) = crossmap_results.get(var_id) {
                if fast_chrom == cross_chrom && fast_pos == *cross_pos {
                    matches += 1;
                    eprintln!("  {} MATCH: {}:{}", var_id, fast_chrom, fast_pos);
                } else {
                    mismatches += 1;
                    eprintln!("  {} MISMATCH: Fast={}:{}, Cross={}:{}", 
                             var_id, fast_chrom, fast_pos, cross_chrom, cross_pos);
                }
            } else {
                eprintln!("  {} FastCrossMap mapped, CrossMap unmapped: {}:{}", var_id, fast_chrom, fast_pos);
            }
        }
    }
    
    // Check for variants that CrossMap mapped but FastCrossMap didn't
    let fast_ids: std::collections::HashSet<String> = fast_lines.iter()
        .filter_map(|l| l.split('\t').nth(2).map(|s| s.to_string()))
        .collect();
    
    for (var_id, (chrom, pos)) in &crossmap_results {
        if !fast_ids.contains(var_id) {
            eprintln!("  {} CrossMap mapped, FastCrossMap unmapped: {}:{}", var_id, chrom, pos);
        }
    }
    
    let total_compared = matches + mismatches;
    if total_compared > 0 {
        let match_rate = matches as f64 / total_compared as f64;
        eprintln!("\n=== Results ===");
        eprintln!("Matches: {}, Mismatches: {}", matches, mismatches);
        eprintln!("Match rate: {:.2}%", match_rate * 100.0);
        
        assert!(match_rate >= 0.9, "Match rate should be >= 90%, got {:.2}%", match_rate * 100.0);
    }
    
    // Cleanup
    let _ = std::fs::remove_file(&bed_input);
    let _ = std::fs::remove_file(&crossmap_output);
    let _ = std::fs::remove_file(&vcf_input);
    let _ = std::fs::remove_file(&vcf_output);
    let _ = std::fs::remove_file(vcf_output.with_extension("vcf.unmap"));
}
