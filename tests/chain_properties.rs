//! Property-based tests for Chain file parsing
//!
//! **Feature: fast-crossmap, Property 1: Chain 解析往返一致性**
//! **Validates: Requirements 1.1, 1.4, 1.5**

use fast_crossmap::core::parse_chain_bytes;
use fast_crossmap::Strand;
use proptest::prelude::*;

/// Generate a valid chromosome name
fn arb_chrom_name() -> impl Strategy<Value = String> {
    prop_oneof![
        // Standard chromosomes with chr prefix
        (1u8..=22).prop_map(|n| format!("chr{}", n)),
        Just("chrX".to_string()),
        Just("chrY".to_string()),
        Just("chrM".to_string()),
        // Standard chromosomes without chr prefix
        (1u8..=22).prop_map(|n| format!("{}", n)),
        Just("X".to_string()),
        Just("Y".to_string()),
        Just("MT".to_string()),
    ]
}

/// Generate a valid strand
fn arb_strand() -> impl Strategy<Value = Strand> {
    prop_oneof![Just(Strand::Plus), Just(Strand::Minus)]
}

/// Generate a valid chain header with consistent coordinates
fn arb_chain_header() -> impl Strategy<Value = ChainHeaderData> {
    (
        1u64..1000000,           // score
        arb_chrom_name(),        // target_name
        100000u64..500000000,    // target_size
        arb_strand(),            // target_strand
        arb_chrom_name(),        // source_name
        100000u64..500000000,    // source_size
        arb_strand(),            // source_strand
        1u64..1000,              // chain_id
    )
        .prop_flat_map(|(score, t_name, t_size, t_strand, s_name, s_size, s_strand, chain_id)| {
            // Generate valid start/end within size bounds
            let t_max_start = t_size.saturating_sub(10000).max(1);
            let s_max_start = s_size.saturating_sub(10000).max(1);
            
            (
                Just(score),
                Just(t_name),
                Just(t_size),
                Just(t_strand),
                0u64..t_max_start,  // target_start
                Just(s_name),
                Just(s_size),
                Just(s_strand),
                0u64..s_max_start,  // source_start
                Just(chain_id),
            )
        })
        .prop_flat_map(|(score, t_name, t_size, t_strand, t_start, s_name, s_size, s_strand, s_start, chain_id)| {
            // Generate end positions that are greater than start
            let t_max_end = t_size.min(t_start + 50000);
            let s_max_end = s_size.min(s_start + 50000);
            
            (
                Just(score),
                Just(t_name),
                Just(t_size),
                Just(t_strand),
                Just(t_start),
                (t_start + 100)..=t_max_end,  // target_end
                Just(s_name),
                Just(s_size),
                Just(s_strand),
                Just(s_start),
                (s_start + 100)..=s_max_end,  // source_end
                Just(chain_id),
            )
        })
        .prop_map(|(score, t_name, t_size, t_strand, t_start, t_end, s_name, s_size, s_strand, s_start, s_end, chain_id)| {
            ChainHeaderData {
                score,
                target_name: t_name,
                target_size: t_size,
                target_strand: t_strand,
                target_start: t_start,
                target_end: t_end,
                source_name: s_name,
                source_size: s_size,
                source_strand: s_strand,
                source_start: s_start,
                source_end: s_end,
                chain_id,
            }
        })
}


/// Data for generating chain headers
#[derive(Debug, Clone)]
struct ChainHeaderData {
    score: u64,
    target_name: String,
    target_size: u64,
    target_strand: Strand,
    target_start: u64,
    target_end: u64,
    source_name: String,
    source_size: u64,
    source_strand: Strand,
    source_start: u64,
    source_end: u64,
    chain_id: u64,
}

impl ChainHeaderData {
    fn to_header_line(&self) -> String {
        // UCSC chain format: chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
        // Our mapping: UCSC target (t) = our source, UCSC query (q) = our target
        // So we write: source fields first (as UCSC target), then target fields (as UCSC query)
        format!(
            "chain {} {} {} {} {} {} {} {} {} {} {} {}",
            self.score,
            self.source_name,      // UCSC tName = our source
            self.source_size,      // UCSC tSize = our source
            self.source_strand.to_char(), // UCSC tStrand = our source
            self.source_start,     // UCSC tStart = our source
            self.source_end,       // UCSC tEnd = our source
            self.target_name,      // UCSC qName = our target
            self.target_size,      // UCSC qSize = our target
            self.target_strand.to_char(), // UCSC qStrand = our target
            self.target_start,     // UCSC qStart = our target
            self.target_end,       // UCSC qEnd = our target
            self.chain_id,
        )
    }
}

/// Generate data lines for a chain block
/// Returns (data_lines, expected_blocks)
fn generate_data_lines(header: &ChainHeaderData) -> (Vec<String>, Vec<ExpectedBlock>) {
    let mut lines = Vec::new();
    let mut blocks = Vec::new();
    
    // Calculate available space
    let source_available = header.source_end.saturating_sub(header.source_start);
    let target_available = header.target_end.saturating_sub(header.target_start);
    let available = source_available.min(target_available);
    
    if available < 10 {
        // Not enough space for any blocks
        return (lines, blocks);
    }
    
    // Generate 1-3 alignment blocks
    let block_size = (available / 3).max(10).min(1000);
    let gap_size = 10u64;
    
    let mut source_pos = header.source_start;
    let mut target_pos = header.target_start;
    
    // First block with gaps
    if source_pos + block_size + gap_size < header.source_end 
        && target_pos + block_size + gap_size < header.target_end 
    {
        lines.push(format!("{} {} {}", block_size, gap_size, gap_size));
        blocks.push(ExpectedBlock {
            source_start: source_pos,
            source_end: source_pos + block_size,
            target_pos,
            size: block_size,
            header: header.clone(),
        });
        source_pos += block_size + gap_size;
        target_pos += block_size + gap_size;
    }
    
    // Second block (last, no gaps)
    let remaining_source = header.source_end.saturating_sub(source_pos);
    let remaining_target = header.target_end.saturating_sub(target_pos);
    let final_size = remaining_source.min(remaining_target).min(block_size);
    
    if final_size >= 10 {
        lines.push(format!("{}", final_size));
        blocks.push(ExpectedBlock {
            source_start: source_pos,
            source_end: source_pos + final_size,
            target_pos,
            size: final_size,
            header: header.clone(),
        });
    }
    
    (lines, blocks)
}

#[derive(Debug, Clone)]
struct ExpectedBlock {
    source_start: u64,
    source_end: u64,
    target_pos: u64,
    size: u64,
    header: ChainHeaderData,
}

impl ExpectedBlock {
    /// Calculate expected target coordinates after strand flipping
    fn expected_target_coords(&self) -> (u64, u64) {
        if self.header.target_strand == Strand::Plus {
            (self.target_pos, self.target_pos + self.size)
        } else {
            // Negative strand: flip coordinates
            let flipped_start = self.header.target_size - (self.target_pos + self.size);
            let flipped_end = self.header.target_size - self.target_pos;
            (flipped_start, flipped_end)
        }
    }
    
    /// Calculate expected source coordinates after strand flipping
    fn expected_source_coords(&self) -> (u64, u64) {
        if self.header.source_strand == Strand::Plus {
            (self.source_start, self.source_end)
        } else {
            // Negative strand: flip coordinates
            let flipped_start = self.header.source_size - self.source_end;
            let flipped_end = self.header.source_size - self.source_start;
            (flipped_start, flipped_end)
        }
    }
}


proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]
    
    /// **Property 1: Chain 解析往返一致性**
    /// 
    /// For any valid chain file, parsing it should produce ChainBlocks with
    /// correct source and target coordinates, properly handling strand flipping.
    ///
    /// **Validates: Requirements 1.1, 1.4, 1.5**
    #[test]
    fn prop_chain_parse_roundtrip(header in arb_chain_header()) {
        let (data_lines, expected_blocks) = generate_data_lines(&header);
        
        // Skip if no valid blocks could be generated
        if expected_blocks.is_empty() {
            return Ok(());
        }
        
        // Build chain file content
        let mut content = header.to_header_line();
        content.push('\n');
        for line in &data_lines {
            content.push_str(line);
            content.push('\n');
        }
        
        // Parse the chain file
        let result = parse_chain_bytes(content.as_bytes());
        prop_assert!(result.is_ok(), "Failed to parse chain: {:?}\nContent:\n{}", result.err(), content);
        
        let chain_file = result.unwrap();
        
        // Verify number of blocks
        prop_assert_eq!(
            chain_file.blocks.len(), 
            expected_blocks.len(),
            "Block count mismatch. Content:\n{}", content
        );
        
        // Verify each block
        for (i, (actual, expected)) in chain_file.blocks.iter().zip(expected_blocks.iter()).enumerate() {
            let (exp_source_start, exp_source_end) = expected.expected_source_coords();
            let (exp_target_start, exp_target_end) = expected.expected_target_coords();
            
            prop_assert_eq!(
                &actual.source_chrom, &expected.header.source_name,
                "Block {}: source_chrom mismatch", i
            );
            prop_assert_eq!(
                actual.source_start, exp_source_start,
                "Block {}: source_start mismatch. Expected {}, got {}. Header: {:?}", 
                i, exp_source_start, actual.source_start, expected.header
            );
            prop_assert_eq!(
                actual.source_end, exp_source_end,
                "Block {}: source_end mismatch. Expected {}, got {}. Header: {:?}", 
                i, exp_source_end, actual.source_end, expected.header
            );
            
            prop_assert_eq!(
                &actual.target_chrom, &expected.header.target_name,
                "Block {}: target_chrom mismatch", i
            );
            prop_assert_eq!(
                actual.target_start, exp_target_start,
                "Block {}: target_start mismatch. Expected {}, got {}. Header: {:?}", 
                i, exp_target_start, actual.target_start, expected.header
            );
            prop_assert_eq!(
                actual.target_end, exp_target_end,
                "Block {}: target_end mismatch. Expected {}, got {}. Header: {:?}", 
                i, exp_target_end, actual.target_end, expected.header
            );
            
            prop_assert_eq!(
                actual.target_strand, expected.header.target_strand,
                "Block {}: target_strand mismatch", i
            );
        }
        
        // Verify chromosome sizes were captured
        prop_assert!(
            chain_file.target_chrom_sizes.contains_key(&header.target_name),
            "Target chrom size not recorded"
        );
        prop_assert!(
            chain_file.source_chrom_sizes.contains_key(&header.source_name),
            "Source chrom size not recorded"
        );
        prop_assert_eq!(
            chain_file.target_chrom_sizes.get(&header.target_name),
            Some(&header.target_size),
            "Target chrom size mismatch"
        );
        prop_assert_eq!(
            chain_file.source_chrom_sizes.get(&header.source_name),
            Some(&header.source_size),
            "Source chrom size mismatch"
        );
    }
    
    /// Property: Block coordinates are always valid (start < end)
    #[test]
    fn prop_chain_block_coords_valid(header in arb_chain_header()) {
        let (data_lines, _) = generate_data_lines(&header);
        
        // Build chain file content
        let mut content = header.to_header_line();
        content.push('\n');
        for line in &data_lines {
            content.push_str(line);
            content.push('\n');
        }
        
        let result = parse_chain_bytes(content.as_bytes());
        if let Ok(chain_file) = result {
            for block in &chain_file.blocks {
                prop_assert!(
                    block.source_start < block.source_end,
                    "Invalid source range: {} >= {}", block.source_start, block.source_end
                );
                prop_assert!(
                    block.target_start < block.target_end,
                    "Invalid target range: {} >= {}", block.target_start, block.target_end
                );
            }
        }
    }
    
    /// Property: Multiple chains can be parsed correctly
    #[test]
    fn prop_multiple_chains_parse(
        header1 in arb_chain_header(),
        header2 in arb_chain_header()
    ) {
        let (data_lines1, blocks1) = generate_data_lines(&header1);
        let (data_lines2, blocks2) = generate_data_lines(&header2);
        
        // Build chain file with two chains
        let mut content = String::new();
        
        // First chain
        content.push_str(&header1.to_header_line());
        content.push('\n');
        for line in &data_lines1 {
            content.push_str(line);
            content.push('\n');
        }
        content.push('\n'); // Empty line separates chains
        
        // Second chain
        content.push_str(&header2.to_header_line());
        content.push('\n');
        for line in &data_lines2 {
            content.push_str(line);
            content.push('\n');
        }
        
        let result = parse_chain_bytes(content.as_bytes());
        prop_assert!(result.is_ok(), "Failed to parse multiple chains: {:?}", result.err());
        
        let chain_file = result.unwrap();
        let expected_total = blocks1.len() + blocks2.len();
        
        prop_assert_eq!(
            chain_file.blocks.len(),
            expected_total,
            "Expected {} blocks, got {}", expected_total, chain_file.blocks.len()
        );
    }
}


/// **Property 12: 压缩文件透明处理**
/// 
/// For any chain file, parsing the plain text version and the gzip/bzip2
/// compressed version should produce identical ChainFile structures.
///
/// **Validates: Requirements 1.2, 1.3**
mod compression_properties {
    use super::*;
    use std::io::Write;
    use flate2::write::GzEncoder;
    use flate2::Compression as GzCompression;
    use bzip2::write::BzEncoder;
    use bzip2::Compression as Bz2Compression;
    
    proptest! {
        #![proptest_config(ProptestConfig::with_cases(50))]
        
        /// Property 12a: Gzip compression transparency
        /// Parsing plain text and gzip-compressed chain files should produce identical results
        #[test]
        fn prop_gzip_transparency(header in arb_chain_header()) {
            let (data_lines, expected_blocks) = generate_data_lines(&header);
            
            if expected_blocks.is_empty() {
                return Ok(());
            }
            
            // Build chain file content
            let mut content = header.to_header_line();
            content.push('\n');
            for line in &data_lines {
                content.push_str(line);
                content.push('\n');
            }
            let content_bytes = content.as_bytes();
            
            // Parse plain text
            let plain_result = parse_chain_bytes(content_bytes);
            prop_assert!(plain_result.is_ok(), "Failed to parse plain text");
            let plain_chain = plain_result.unwrap();
            
            // Create gzip compressed version
            let mut encoder = GzEncoder::new(Vec::new(), GzCompression::default());
            encoder.write_all(content_bytes).unwrap();
            let gz_data = encoder.finish().unwrap();
            
            // Write to temp file and parse
            let temp_dir = std::env::temp_dir();
            let gz_path = temp_dir.join(format!("prop_test_{}.chain.gz", header.chain_id));
            std::fs::write(&gz_path, &gz_data).unwrap();
            
            let gz_result = fast_crossmap::parse_chain_file(&gz_path);
            let _ = std::fs::remove_file(&gz_path);
            
            prop_assert!(gz_result.is_ok(), "Failed to parse gzip file");
            let gz_chain = gz_result.unwrap();
            
            // Compare results
            prop_assert_eq!(
                plain_chain.blocks.len(), 
                gz_chain.blocks.len(),
                "Block count mismatch between plain and gzip"
            );
            
            for (i, (plain_block, gz_block)) in plain_chain.blocks.iter().zip(gz_chain.blocks.iter()).enumerate() {
                prop_assert_eq!(
                    plain_block, gz_block,
                    "Block {} differs between plain and gzip", i
                );
            }
            
            prop_assert_eq!(
                plain_chain.target_chrom_sizes,
                gz_chain.target_chrom_sizes,
                "Target chrom sizes differ"
            );
            prop_assert_eq!(
                plain_chain.source_chrom_sizes,
                gz_chain.source_chrom_sizes,
                "Source chrom sizes differ"
            );
        }
        
        /// Property 12b: Bzip2 compression transparency
        /// Parsing plain text and bzip2-compressed chain files should produce identical results
        #[test]
        fn prop_bzip2_transparency(header in arb_chain_header()) {
            let (data_lines, expected_blocks) = generate_data_lines(&header);
            
            if expected_blocks.is_empty() {
                return Ok(());
            }
            
            // Build chain file content
            let mut content = header.to_header_line();
            content.push('\n');
            for line in &data_lines {
                content.push_str(line);
                content.push('\n');
            }
            let content_bytes = content.as_bytes();
            
            // Parse plain text
            let plain_result = parse_chain_bytes(content_bytes);
            prop_assert!(plain_result.is_ok(), "Failed to parse plain text");
            let plain_chain = plain_result.unwrap();
            
            // Create bzip2 compressed version
            let mut encoder = BzEncoder::new(Vec::new(), Bz2Compression::default());
            encoder.write_all(content_bytes).unwrap();
            let bz2_data = encoder.finish().unwrap();
            
            // Write to temp file and parse
            let temp_dir = std::env::temp_dir();
            let bz2_path = temp_dir.join(format!("prop_test_{}.chain.bz2", header.chain_id));
            std::fs::write(&bz2_path, &bz2_data).unwrap();
            
            let bz2_result = fast_crossmap::parse_chain_file(&bz2_path);
            let _ = std::fs::remove_file(&bz2_path);
            
            prop_assert!(bz2_result.is_ok(), "Failed to parse bzip2 file");
            let bz2_chain = bz2_result.unwrap();
            
            // Compare results
            prop_assert_eq!(
                plain_chain.blocks.len(), 
                bz2_chain.blocks.len(),
                "Block count mismatch between plain and bzip2"
            );
            
            for (i, (plain_block, bz2_block)) in plain_chain.blocks.iter().zip(bz2_chain.blocks.iter()).enumerate() {
                prop_assert_eq!(
                    plain_block, bz2_block,
                    "Block {} differs between plain and bzip2", i
                );
            }
        }
        
        /// Property 12c: All compression formats produce identical results
        /// Plain, gzip, and bzip2 should all produce the same ChainFile
        #[test]
        fn prop_all_formats_equivalent(header in arb_chain_header()) {
            let (data_lines, expected_blocks) = generate_data_lines(&header);
            
            if expected_blocks.is_empty() {
                return Ok(());
            }
            
            // Build chain file content
            let mut content = header.to_header_line();
            content.push('\n');
            for line in &data_lines {
                content.push_str(line);
                content.push('\n');
            }
            let content_bytes = content.as_bytes();
            
            // Parse plain text
            let plain_chain = parse_chain_bytes(content_bytes).unwrap();
            
            // Create and parse gzip
            let mut gz_encoder = GzEncoder::new(Vec::new(), GzCompression::default());
            gz_encoder.write_all(content_bytes).unwrap();
            let gz_data = gz_encoder.finish().unwrap();
            
            let temp_dir = std::env::temp_dir();
            let gz_path = temp_dir.join(format!("prop_all_{}.chain.gz", header.chain_id));
            std::fs::write(&gz_path, &gz_data).unwrap();
            let gz_chain = fast_crossmap::parse_chain_file(&gz_path).unwrap();
            let _ = std::fs::remove_file(&gz_path);
            
            // Create and parse bzip2
            let mut bz2_encoder = BzEncoder::new(Vec::new(), Bz2Compression::default());
            bz2_encoder.write_all(content_bytes).unwrap();
            let bz2_data = bz2_encoder.finish().unwrap();
            
            let bz2_path = temp_dir.join(format!("prop_all_{}.chain.bz2", header.chain_id));
            std::fs::write(&bz2_path, &bz2_data).unwrap();
            let bz2_chain = fast_crossmap::parse_chain_file(&bz2_path).unwrap();
            let _ = std::fs::remove_file(&bz2_path);
            
            // All three should be identical
            prop_assert_eq!(plain_chain.blocks.len(), gz_chain.blocks.len());
            prop_assert_eq!(plain_chain.blocks.len(), bz2_chain.blocks.len());
            
            for i in 0..plain_chain.blocks.len() {
                prop_assert_eq!(&plain_chain.blocks[i], &gz_chain.blocks[i]);
                prop_assert_eq!(&plain_chain.blocks[i], &bz2_chain.blocks[i]);
            }
        }
    }
}
