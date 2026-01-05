//! Property-based tests for CoordinateMapper
//!
//! **Feature: fast-crossmap, Property 5, 6, 7: 坐标映射正确性**
//! **Validates: Requirements 3.4, 3.5, 3.6**

use fast_crossmap::core::{
    ChainIndex, ChainFile, ChainBlock, CoordinateMapper, ChromStyle, Strand,
    intersect_intervals,
};
use proptest::prelude::*;
use std::collections::HashMap;

/// Create a ChainFile with a single block for testing
fn create_single_block_chain(
    source_chrom: &str,
    source_start: u64,
    source_end: u64,
    target_chrom: &str,
    target_start: u64,
    target_end: u64,
    target_strand: Strand,
) -> ChainFile {
    let chrom_size = 500_000_000u64;
    
    let mut target_chrom_sizes = HashMap::new();
    let mut source_chrom_sizes = HashMap::new();
    target_chrom_sizes.insert(target_chrom.to_string(), chrom_size);
    source_chrom_sizes.insert(source_chrom.to_string(), chrom_size);
    
    ChainFile {
        blocks: vec![ChainBlock {
            source_chrom: source_chrom.to_string(),
            source_start,
            source_end,
            target_chrom: target_chrom.to_string(),
            target_start,
            target_end,
            target_strand,
        }],
        target_chrom_sizes,
        source_chrom_sizes,
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]
    
    /// **Property 5: 正链偏移计算正确性**
    /// 
    /// For any mapping to a positive strand target, the target start should equal
    /// `t_start + left_offset` where `left_offset = query_start - source_block_start`.
    ///
    /// **Validates: Requirements 3.4**
    #[test]
    fn prop_positive_strand_offset(
        s_start in 0u64..100000,
        block_size in 100u64..10000,
        t_start in 0u64..100000,
        query_offset in 0u64..100,
        query_len in 1u64..100,
    ) {
        let s_end = s_start + block_size;
        let t_end = t_start + block_size;
        
        // Ensure query is within block
        let q_start = s_start + query_offset.min(block_size - 1);
        let q_end = (q_start + query_len).min(s_end);
        
        if q_start >= q_end {
            return Ok(());
        }
        
        // Create chain with positive strand target
        let chain_file = create_single_block_chain(
            "chr1", s_start, s_end,
            "chr1", t_start, t_end,
            Strand::Plus,
        );
        let index = ChainIndex::from_chain_data(chain_file);
        let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
        
        // Map the query
        let results = mapper.map("chr1", q_start, q_end, Strand::Plus);
        prop_assert!(results.is_some(), "Mapping should succeed");
        let results = results.unwrap();
        prop_assert_eq!(results.len(), 1, "Should have exactly one result");
        
        let segment = &results[0];
        
        // Calculate expected values
        let left_offset = q_start - s_start;
        let expected_target_start = t_start + left_offset;
        let size = q_end - q_start;
        let expected_target_end = expected_target_start + size;
        
        prop_assert_eq!(
            segment.target.start, expected_target_start,
            "Target start should be t_start + left_offset: {} + {} = {}",
            t_start, left_offset, expected_target_start
        );
        prop_assert_eq!(
            segment.target.end, expected_target_end,
            "Target end should be target_start + size"
        );
    }

    
    /// **Property 6: 负链偏移计算正确性**
    /// 
    /// For any mapping to a negative strand target, the target start should equal
    /// `t_end - left_offset - size` where size is the mapped region size.
    ///
    /// **Validates: Requirements 3.5**
    #[test]
    fn prop_negative_strand_offset(
        s_start in 0u64..100000,
        block_size in 100u64..10000,
        t_start in 0u64..100000,
        query_offset in 0u64..100,
        query_len in 1u64..100,
    ) {
        let s_end = s_start + block_size;
        let t_end = t_start + block_size;
        
        // Ensure query is within block
        let q_start = s_start + query_offset.min(block_size - 1);
        let q_end = (q_start + query_len).min(s_end);
        
        if q_start >= q_end {
            return Ok(());
        }
        
        // Create chain with negative strand target
        let chain_file = create_single_block_chain(
            "chr1", s_start, s_end,
            "chr1", t_start, t_end,
            Strand::Minus,
        );
        let index = ChainIndex::from_chain_data(chain_file);
        let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
        
        // Map the query
        let results = mapper.map("chr1", q_start, q_end, Strand::Plus);
        prop_assert!(results.is_some(), "Mapping should succeed");
        let results = results.unwrap();
        prop_assert_eq!(results.len(), 1, "Should have exactly one result");
        
        let segment = &results[0];
        
        // Calculate expected values for negative strand
        let left_offset = q_start - s_start;
        let size = q_end - q_start;
        let expected_target_start = t_end - left_offset - size;
        let expected_target_end = expected_target_start + size;
        
        prop_assert_eq!(
            segment.target.start, expected_target_start,
            "Target start should be t_end - left_offset - size: {} - {} - {} = {}",
            t_end, left_offset, size, expected_target_start
        );
        prop_assert_eq!(
            segment.target.end, expected_target_end,
            "Target end should be target_start + size"
        );
        
        // Strand should be Minus (Plus query + Minus target = Minus)
        prop_assert_eq!(segment.target.strand, Strand::Minus);
    }
    
    /// **Property 7: 链方向互补正确性**
    /// 
    /// For any mapping, the output strand should follow the XOR rule:
    /// - Plus + Plus = Plus
    /// - Plus + Minus = Minus
    /// - Minus + Plus = Minus
    /// - Minus + Minus = Plus
    ///
    /// **Validates: Requirements 3.6**
    #[test]
    fn prop_strand_combination(
        s_start in 0u64..100000,
        block_size in 100u64..10000,
        t_start in 0u64..100000,
        query_strand in prop_oneof![Just(Strand::Plus), Just(Strand::Minus)],
        target_strand in prop_oneof![Just(Strand::Plus), Just(Strand::Minus)],
    ) {
        let s_end = s_start + block_size;
        let t_end = t_start + block_size;
        
        // Query in the middle of the block
        let q_start = s_start + block_size / 4;
        let q_end = s_start + block_size / 2;
        
        // Create chain
        let chain_file = create_single_block_chain(
            "chr1", s_start, s_end,
            "chr1", t_start, t_end,
            target_strand,
        );
        let index = ChainIndex::from_chain_data(chain_file);
        let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
        
        // Map the query
        let results = mapper.map("chr1", q_start, q_end, query_strand);
        prop_assert!(results.is_some(), "Mapping should succeed");
        let results = results.unwrap();
        prop_assert_eq!(results.len(), 1, "Should have exactly one result");
        
        let segment = &results[0];
        
        // Expected strand follows XOR rule
        let expected_strand = query_strand.combine(target_strand);
        
        prop_assert_eq!(
            segment.target.strand, expected_strand,
            "Strand combination: {:?} + {:?} should be {:?}",
            query_strand, target_strand, expected_strand
        );
    }
    
    /// Property: Intersection calculation is correct
    #[test]
    fn prop_intersect_intervals_symmetric(
        start1 in 0u64..10000,
        len1 in 1u64..1000,
        start2 in 0u64..10000,
        len2 in 1u64..1000,
    ) {
        let end1 = start1 + len1;
        let end2 = start2 + len2;
        
        let result1 = intersect_intervals(start1, end1, start2, end2);
        let result2 = intersect_intervals(start2, end2, start1, end1);
        
        // Intersection should be symmetric
        prop_assert_eq!(result1, result2, "Intersection should be symmetric");
        
        // If there's an intersection, verify it's correct
        if let Some((i_start, i_end)) = result1 {
            prop_assert!(i_start >= start1.max(start2));
            prop_assert!(i_end <= end1.min(end2));
            prop_assert!(i_start < i_end, "Intersection should be non-empty");
        }
    }
    
    /// Property: Mapped region size equals intersection size
    #[test]
    fn prop_mapped_size_equals_intersection(
        s_start in 0u64..100000,
        block_size in 100u64..10000,
        t_start in 0u64..100000,
        query_offset in 0u64..100,
        query_len in 1u64..100,
    ) {
        let s_end = s_start + block_size;
        let t_end = t_start + block_size;
        
        // Ensure query is within block
        let q_start = s_start + query_offset.min(block_size - 1);
        let q_end = (q_start + query_len).min(s_end);
        
        if q_start >= q_end {
            return Ok(());
        }
        
        let chain_file = create_single_block_chain(
            "chr1", s_start, s_end,
            "chr1", t_start, t_end,
            Strand::Plus,
        );
        let index = ChainIndex::from_chain_data(chain_file);
        let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
        
        let results = mapper.map("chr1", q_start, q_end, Strand::Plus);
        prop_assert!(results.is_some());
        let results = results.unwrap();
        
        if !results.is_empty() {
            let segment = &results[0];
            let source_size = segment.source.end - segment.source.start;
            let target_size = segment.target.end - segment.target.start;
            
            prop_assert_eq!(
                source_size, target_size,
                "Source and target sizes should be equal"
            );
        }
    }
    
    /// Property: Source region in result is the intersection
    #[test]
    fn prop_source_region_is_intersection(
        s_start in 0u64..100000,
        block_size in 100u64..10000,
        t_start in 0u64..100000,
        q_start_offset in 0u64..200,
        q_len in 1u64..200,
    ) {
        let s_end = s_start + block_size;
        let t_end = t_start + block_size;
        
        // Query may partially overlap
        let q_start = if q_start_offset < 100 {
            s_start.saturating_sub(100 - q_start_offset)
        } else {
            s_start + (q_start_offset - 100)
        };
        let q_end = q_start + q_len;
        
        let chain_file = create_single_block_chain(
            "chr1", s_start, s_end,
            "chr1", t_start, t_end,
            Strand::Plus,
        );
        let index = ChainIndex::from_chain_data(chain_file);
        let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
        
        let results = mapper.map("chr1", q_start, q_end, Strand::Plus);
        prop_assert!(results.is_some());
        let results = results.unwrap();
        
        // Calculate expected intersection
        let expected_intersection = intersect_intervals(q_start, q_end, s_start, s_end);
        
        if let Some((exp_start, exp_end)) = expected_intersection {
            prop_assert_eq!(results.len(), 1);
            let segment = &results[0];
            prop_assert_eq!(segment.source.start, exp_start);
            prop_assert_eq!(segment.source.end, exp_end);
        } else {
            prop_assert!(results.is_empty(), "No intersection means no results");
        }
    }
}
