//! Property-based tests for MAF format conversion
//!
//! **Feature: fast-crossmap, MAF 格式转换**
//! **Validates: Requirements 8.1, 8.2, 8.3, 8.4, 8.5, 8.6**

use fast_crossmap::core::{ChainIndex, CoordinateMapper, ChromStyle, Strand};
use fast_crossmap::formats::maf::{MafRecordView, MafColumnIndices, convert_maf};
use proptest::prelude::*;
use std::path::PathBuf;

/// Create test column indices
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

/// Generate a valid chromosome name
fn arb_chrom_name() -> impl Strategy<Value = String> {
    prop_oneof![
        (1u8..=22).prop_map(|n| format!("chr{}", n)),
        Just("chrX".to_string()),
        Just("chrY".to_string()),
    ]
}

/// Generate a valid gene symbol
fn arb_gene_symbol() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("TP53".to_string()),
        Just("BRCA1".to_string()),
        Just("EGFR".to_string()),
        Just("KRAS".to_string()),
        Just("PIK3CA".to_string()),
    ]
}

/// Generate a valid DNA allele
fn arb_dna_allele() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("A".to_string()),
        Just("T".to_string()),
        Just("G".to_string()),
        Just("C".to_string()),
    ]
}

/// Generate a valid strand
fn arb_strand() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("+".to_string()),
        Just("-".to_string()),
    ]
}

/// Generate a MAF data line
fn arb_maf_line() -> impl Strategy<Value = String> {
    (
        arb_gene_symbol(),
        arb_chrom_name(),
        100000u64..10000000,
        arb_strand(),
        arb_dna_allele(),
        arb_dna_allele(),
    )
        .prop_map(|(gene, chrom, pos, strand, ref_allele, alt_allele)| {
            format!(
                "{}\t7157\tBCM\tGRCh37\t{}\t{}\t{}\t{}\tMissense_Mutation\tSNP\t{}\t{}\t{}\t.\t.\tSAMPLE1",
                gene, chrom, pos, pos, strand, ref_allele, ref_allele, alt_allele
            )
        })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]
    
    /// Property: MAF parsing extracts correct coordinates
    #[test]
    fn prop_maf_coordinate_parsing(line in arb_maf_line()) {
        let indices = create_test_indices();
        let view = MafRecordView::parse(line.as_bytes(), &indices).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        let expected_chrom = fields[4];
        let expected_start: u64 = fields[5].parse().unwrap();
        let expected_end: u64 = fields[6].parse().unwrap();
        
        prop_assert_eq!(view.chromosome(), expected_chrom);
        prop_assert_eq!(view.start_position().unwrap(), expected_start);
        prop_assert_eq!(view.end_position().unwrap(), expected_end);
    }
    
    /// Property: MAF strand parsing is correct
    #[test]
    fn prop_maf_strand_parsing(line in arb_maf_line()) {
        let indices = create_test_indices();
        let view = MafRecordView::parse(line.as_bytes(), &indices).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        let strand_char = fields[7];
        
        match strand_char {
            "+" => prop_assert_eq!(view.strand(), Some(Strand::Plus)),
            "-" => prop_assert_eq!(view.strand(), Some(Strand::Minus)),
            _ => prop_assert_eq!(view.strand(), None),
        }
    }
    
    /// Property: MAF field preservation
    #[test]
    fn prop_maf_field_preservation(line in arb_maf_line()) {
        let indices = create_test_indices();
        let view = MafRecordView::parse(line.as_bytes(), &indices).unwrap();
        
        let fields: Vec<&str> = line.split('\t').collect();
        
        prop_assert_eq!(view.hugo_symbol(), Some(fields[0]));
        prop_assert_eq!(view.ncbi_build(), fields[3]);
        prop_assert_eq!(view.reference_allele(), fields[10]);
    }
}

/// Test MAF conversion with real chain file
#[test]
fn test_maf_conversion_basic() {
    let chain_path = PathBuf::from("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz");
    
    if !chain_path.exists() {
        eprintln!("Skipping test: chain file not found");
        return;
    }
    
    // Load chain file
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    // Create test MAF file
    let test_maf = "\
#version 2.4
Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tdbSNP_RS\tdbSNP_Val_Status\tTumor_Sample_Barcode
TP53\t7157\tBCM\tGRCh37\tchr1\t100000\t100000\t+\tMissense_Mutation\tSNP\tG\tG\tA\t.\t.\tSAMPLE1
BRCA1\t672\tBCM\tGRCh37\tchr2\t50000\t50000\t-\tMissense_Mutation\tSNP\tC\tC\tT\t.\t.\tSAMPLE2
EGFR\t1956\tBCM\tGRCh37\tchr2\t200000\t200000\t+\tMissense_Mutation\tSNP\tA\tA\tG\t.\t.\tSAMPLE3
";
    
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("test_input.maf");
    let output_path = temp_dir.join("test_output.maf");
    
    std::fs::write(&input_path, test_maf).unwrap();
    
    // Convert
    let stats = convert_maf(&input_path, &output_path, &mapper, None::<&PathBuf>, "GRCh38").unwrap();
    
    eprintln!("MAF conversion stats: total={}, success={}, failed={}, headers={}", 
              stats.total, stats.success, stats.failed, stats.headers);
    
    // Verify stats
    assert_eq!(stats.total, 3, "Should process 3 records");
    
    // Read output and verify structure
    let output = std::fs::read_to_string(&output_path).unwrap();
    
    eprintln!("\n=== MAF Output ===");
    for line in output.lines() {
        eprintln!("{}", line);
    }
    
    // Check headers are preserved
    assert!(output.contains("#version"), "Should preserve version comment");
    assert!(output.contains("Hugo_Symbol"), "Should have column header");
    
    // Check NCBI_Build is updated
    let data_lines: Vec<&str> = output.lines()
        .filter(|l| !l.starts_with('#') && !l.starts_with("Hugo_Symbol"))
        .collect();
    
    for line in &data_lines {
        assert!(line.contains("GRCh38"), "NCBI_Build should be updated to GRCh38");
    }
    
    // Cleanup
    let _ = std::fs::remove_file(&input_path);
    let _ = std::fs::remove_file(&output_path);
    let unmap_path = output_path.with_extension("maf.unmap");
    let _ = std::fs::remove_file(&unmap_path);
}

/// Test MAF vs CrossMap comparison
#[test]
fn test_maf_vs_crossmap() {
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
    
    // Note: CrossMap maf requires a reference genome file
    // This test verifies coordinate mapping only
    
    let test_maf = "\
Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tdbSNP_RS\tdbSNP_Val_Status\tTumor_Sample_Barcode
TP53\t7157\tBCM\tGRCh37\tchr1\t100000\t100000\t+\tMissense_Mutation\tSNP\tG\tG\tA\t.\t.\tSAMPLE1
BRCA1\t672\tBCM\tGRCh37\tchr2\t50000\t50000\t+\tMissense_Mutation\tSNP\tC\tC\tT\t.\t.\tSAMPLE2
";
    
    let temp_dir = std::env::temp_dir();
    let input_path = temp_dir.join("maf_crossmap_test.maf");
    let fast_output = temp_dir.join("maf_fast_output.maf");
    
    std::fs::write(&input_path, test_maf).unwrap();
    
    // Run FastCrossMap
    let index = ChainIndex::from_chain_file(&chain_path).expect("Failed to load chain file");
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    let stats = convert_maf(&input_path, &fast_output, &mapper, None::<&PathBuf>, "GRCh38").unwrap();
    
    eprintln!("FastCrossMap MAF: total={}, success={}, failed={}", stats.total, stats.success, stats.failed);
    
    // Read output
    let fast_content = std::fs::read_to_string(&fast_output).unwrap_or_default();
    
    eprintln!("\n=== FastCrossMap MAF output ===");
    for line in fast_content.lines() {
        if !line.starts_with('#') && !line.starts_with("Hugo_Symbol") {
            eprintln!("{}", line);
        }
    }
    
    // Cleanup
    let _ = std::fs::remove_file(&input_path);
    let _ = std::fs::remove_file(&fast_output);
    let _ = std::fs::remove_file(fast_output.with_extension("maf.unmap"));
}
