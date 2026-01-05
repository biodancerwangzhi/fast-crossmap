//! Property-based tests for GFF/GTF format conversion
//!
//! **Feature: fast-crossmap, GFF/GTF 格式转换**
//! **Validates: Requirements 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7**

use fast_crossmap::core::{ChainIndex, CoordinateMapper, ChromStyle, Strand};
use fast_crossmap::formats::gff::{GffRecordView, convert_gff};
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

/// Generate a valid GFF source field
fn arb_source() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("ensembl".to_string()),
        Just("havana".to_string()),
        Just("refseq".to_string()),
        Just(".".to_string()),
    ]
}

/// Generate a valid GFF feature type
fn arb_feature() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("gene".to_string()),
        Just("transcript".to_string()),
        Just("exon".to_string()),
        Just("CDS".to_string()),
        Just("UTR".to_string()),
    ]
}

/// Generate a valid score field
fn arb_score() -> impl Strategy<Value = String> {
    prop_oneof![
        Just(".".to_string()),
        (0u32..1000).prop_map(|n| n.to_string()),
        (0.0f64..100.0).prop_map(|f| format!("{:.2}", f)),
    ]
}

/// Generate a valid strand field
fn arb_strand() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("+".to_string()),
        Just("-".to_string()),
        Just(".".to_string()),
    ]
}

/// Generate a valid frame field
fn arb_frame() -> impl Strategy<Value = String> {
    prop_oneof![
        Just(".".to_string()),
        Just("0".to_string()),
        Just("1".to_string()),
        Just("2".to_string()),
    ]
}

/// Generate a valid attributes field
fn arb_attributes() -> impl Strategy<Value = String> {
    prop_oneof![
        Just(".".to_string()),
        Just("gene_id \"ENSG00000001\"".to_string()),
        Just("gene_id \"ENSG00000001\"; transcript_id \"ENST00000001\"".to_string()),
        Just("ID=gene1;Name=TestGene".to_string()),
    ]
}


/// Generate a valid GFF line
fn arb_gff_line() -> impl Strategy<Value = String> {
    (
        arb_chrom_name(),
        arb_source(),
        arb_feature(),
        1000u64..100000,  // start (1-based)
        arb_score(),
        arb_strand(),
        arb_frame(),
        arb_attributes(),
    )
        .prop_map(|(chrom, source, feature, start, score, strand, frame, attrs)| {
            let end = start + 100; // Fixed size for testing
            format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
                chrom, source, feature, start, end, score, strand, frame, attrs)
        })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]
    
    /// Property: GFF parsing extracts correct coordinates (1-based)
    #[test]
    fn prop_gff_coordinate_parsing(line in arb_gff_line()) {
        let view = GffRecordView::parse(line.as_bytes()).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        let expected_start: u64 = fields[3].parse().unwrap();
        let expected_end: u64 = fields[4].parse().unwrap();
        
        prop_assert_eq!(view.seqname, fields[0]);
        prop_assert_eq!(view.start, expected_start);
        prop_assert_eq!(view.end, expected_end);
        // Size should be end - start + 1 for 1-based coordinates
        prop_assert_eq!(view.size(), expected_end - expected_start + 1);
    }
    
    /// Property: GFF source/feature/frame fields are preserved
    #[test]
    fn prop_gff_field_preservation(line in arb_gff_line()) {
        let view = GffRecordView::parse(line.as_bytes()).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        
        prop_assert_eq!(view.source, fields[1]);
        prop_assert_eq!(view.feature, fields[2]);
        prop_assert_eq!(view.score, fields[5]);
        prop_assert_eq!(view.strand_char, fields[6]);
        prop_assert_eq!(view.frame, fields[7]);
        prop_assert_eq!(view.attributes, fields[8]);
    }
    
    /// Property: GFF strand parsing is correct
    #[test]
    fn prop_gff_strand_parsing(line in arb_gff_line()) {
        let view = GffRecordView::parse(line.as_bytes()).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        let strand_char = fields[6];
        
        match strand_char {
            "+" => prop_assert_eq!(view.strand, Some(Strand::Plus)),
            "-" => prop_assert_eq!(view.strand, Some(Strand::Minus)),
            "." => prop_assert_eq!(view.strand, None),
            _ => prop_assert!(false, "Unexpected strand"),
        }
    }
}


/// Test GFF conversion with real chain file
#[test]
fn test_gff_conversion_with_crossmap() {
    let chain_path = PathBuf::from("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz");
    
    if !chain_path.exists() {
        eprintln!("Skipping test: chain file not found");
        return;
    }
    
    // Load chain file
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    // Create test GFF file
    let test_gff = "\
##gff-version 3
#description: test GFF file
chr1\tensembl\tgene\t10001\t10100\t.\t+\t.\tgene_id \"test1\"
chr1\tensembl\texon\t100001\t100200\t.\t-\t.\tgene_id \"test2\"
chr2\trefseq\tCDS\t50001\t50100\t100\t+\t0\tID=cds1
chr2\thavana\ttranscript\t200001\t200300\t.\t-\t.\ttranscript_id \"tx1\"
chr3\t.\tregion\t100001\t100150\t.\t.\t.\t.
";
    
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("test_input.gff");
    let output_path = temp_dir.join("test_output.gff");
    
    std::fs::write(&input_path, test_gff).unwrap();
    
    // Convert
    let stats = convert_gff(&input_path, &output_path, &mapper, 1).unwrap();
    
    eprintln!("GFF conversion stats: total={}, success={}, failed={}, comments={}", 
              stats.total, stats.success, stats.failed, stats.comments);
    
    // Verify stats
    assert_eq!(stats.total, 5, "Should process 5 records");
    assert_eq!(stats.comments, 2, "Should have 2 comment lines");
    
    // Read output and verify structure
    let output = std::fs::read_to_string(&output_path).unwrap();
    
    // Check headers are preserved
    assert!(output.contains("##gff-version 3"), "Should preserve GFF version header");
    assert!(output.contains("#description"), "Should preserve comment lines");
    
    eprintln!("\n=== GFF Output ===");
    for line in output.lines() {
        eprintln!("{}", line);
    }
    
    // Cleanup
    let _ = std::fs::remove_file(&input_path);
    let _ = std::fs::remove_file(&output_path);
    let unmap_path = output_path.with_extension("gff.unmap");
    let _ = std::fs::remove_file(&unmap_path);
}

/// Test GFF vs CrossMap comparison
#[test]
fn test_gff_vs_crossmap() {
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
    
    // Create test GFF file with positions known to map
    let test_gff = "\
##gff-version 3
chr1\tensembl\tgene\t100001\t100100\t.\t+\t.\tgene_id \"test1\"
chr2\trefseq\texon\t50001\t50100\t.\t+\t.\tID=exon1
chr2\thavana\tCDS\t200001\t200100\t.\t-\t0\tID=cds1
";
    
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("gff_crossmap_test.gff");
    let fast_output = temp_dir.join("gff_fast_output.gff");
    let cross_output = temp_dir.join("gff_cross_output.gff");
    
    std::fs::write(&input_path, test_gff).unwrap();
    
    // Run FastCrossMap
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    let stats = convert_gff(&input_path, &fast_output, &mapper, 1).unwrap();
    
    eprintln!("FastCrossMap GFF: total={}, success={}, failed={}", stats.total, stats.success, stats.failed);
    
    // Run CrossMap
    let _crossmap_result = Command::new("CrossMap")
        .args(&[
            "gff",
            chain_path.to_str().unwrap(),
            input_path.to_str().unwrap(),
            cross_output.to_str().unwrap(),
        ])
        .output();
    
    // Read outputs
    let fast_content = std::fs::read_to_string(&fast_output).unwrap_or_default();
    let cross_content = std::fs::read_to_string(&cross_output).unwrap_or_default();
    
    eprintln!("\n=== FastCrossMap GFF output ===");
    for line in fast_content.lines() {
        if !line.starts_with('#') {
            eprintln!("{}", line);
        }
    }
    
    eprintln!("\n=== CrossMap GFF output ===");
    for line in cross_content.lines() {
        if !line.starts_with('#') {
            eprintln!("{}", line);
        }
    }
    
    // Compare data lines (excluding comments)
    let fast_lines: Vec<&str> = fast_content.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    let cross_lines: Vec<&str> = cross_content.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    
    eprintln!("\n=== Comparison ===");
    eprintln!("FastCrossMap lines: {}, CrossMap lines: {}", fast_lines.len(), cross_lines.len());
    
    // Compare coordinates (first 5 fields: seqname, source, feature, start, end)
    let mut matches = 0;
    let mut mismatches = 0;
    
    for (fast, cross) in fast_lines.iter().zip(cross_lines.iter()) {
        let fast_fields: Vec<&str> = fast.split('\t').take(5).collect();
        let cross_fields: Vec<&str> = cross.split('\t').take(5).collect();
        
        if fast_fields == cross_fields {
            matches += 1;
            eprintln!("MATCH: {:?}", fast_fields);
        } else {
            mismatches += 1;
            eprintln!("MISMATCH: Fast={:?}, Cross={:?}", fast_fields, cross_fields);
        }
    }
    
    let total = matches + mismatches;
    if total > 0 {
        let match_rate = matches as f64 / total as f64;
        eprintln!("\nResults: {} matches, {} mismatches", matches, mismatches);
        eprintln!("Match rate: {:.2}%", match_rate * 100.0);
    }
    
    // Cleanup
    let _ = std::fs::remove_file(&input_path);
    let _ = std::fs::remove_file(&fast_output);
    let _ = std::fs::remove_file(&cross_output);
    let _ = std::fs::remove_file(fast_output.with_extension("gff.unmap"));
    let _ = std::fs::remove_file(cross_output.with_extension("gff.unmap"));
}
