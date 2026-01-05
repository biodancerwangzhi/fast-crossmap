//! Property-based tests for ChainIndex interval queries
//!
//! **Feature: fast-crossmap, Property 2: 区间查询完整性**
//! **Validates: Requirements 2.1, 2.2**

use fast_crossmap::core::{ChainIndex, ChainFile, ChainBlock};
use fast_crossmap::Strand;
use proptest::prelude::*;
use std::collections::HashMap;

/// Generate a random chromosome name
fn arb_chrom() -> impl Strategy<Value = String> {
    prop_oneof![
        (1u8..=22).prop_map(|n| format!("chr{}", n)),
        Just("chrX".to_string()),
        Just("chrY".to_string()),
    ]
}

/// Generate a list of non-overlapping intervals for a chromosome
fn arb_intervals(_chrom: String, chrom_size: u64) -> impl Strategy<Value = Vec<(u64, u64)>> {
    // Generate 1-10 intervals
    (1usize..=10).prop_flat_map(move |count| {
        let size = chrom_size;
        prop::collection::vec(
            (0u64..size.saturating_sub(100), 10u64..100),
            count,
        ).prop_map(move |pairs| {
            // Sort and make non-overlapping
            let mut intervals: Vec<(u64, u64)> = pairs
                .into_iter()
                .map(|(start, len)| (start, (start + len).min(size)))
                .collect();
            intervals.sort_by_key(|&(s, _)| s);
            
            // Remove overlaps
            let mut result = Vec::new();
            let mut last_end = 0u64;
            for (start, end) in intervals {
                let actual_start = start.max(last_end);
                if actual_start < end {
                    result.push((actual_start, end));
                    last_end = end + 10; // Gap between intervals
                }
            }
            result
        })
    })
}

/// Create a ChainFile from intervals
fn create_chain_file(intervals_by_chrom: HashMap<String, Vec<(u64, u64)>>) -> ChainFile {
    let mut blocks = Vec::new();
    let mut target_chrom_sizes = HashMap::new();
    let mut source_chrom_sizes = HashMap::new();
    
    for (chrom, intervals) in &intervals_by_chrom {
        let chrom_size = 500_000_000u64; // Large enough
        target_chrom_sizes.insert(chrom.clone(), chrom_size);
        source_chrom_sizes.insert(chrom.clone(), chrom_size);
        
        for &(start, end) in intervals {
            blocks.push(ChainBlock {
                source_chrom: chrom.clone(),
                source_start: start,
                source_end: end,
                target_chrom: chrom.clone(),
                target_start: start,
                target_end: end,
                target_strand: Strand::Plus,
            });
        }
    }
    
    ChainFile {
        blocks,
        target_chrom_sizes,
        source_chrom_sizes,
    }
}


proptest! {
    #![proptest_config(ProptestConfig::with_cases(100))]
    
    /// **Property 2: 区间查询完整性**
    /// 
    /// For any coordinate range [start, end] on a chromosome, querying the index
    /// should return all and only the chain blocks that overlap with that range.
    ///
    /// **Validates: Requirements 2.1, 2.2**
    #[test]
    fn prop_query_returns_all_overlapping(
        _chrom in arb_chrom(),
        intervals in arb_intervals("chr1".to_string(), 10000),
        query_start in 0u64..10000,
        query_len in 10u64..1000,
    ) {
        let query_end = query_start + query_len;
        
        // Create index with intervals on chr1
        let mut intervals_map = HashMap::new();
        intervals_map.insert("chr1".to_string(), intervals.clone());
        let chain_file = create_chain_file(intervals_map);
        let index = ChainIndex::from_chain_data(chain_file);
        
        // Query the index
        let results = index.query("chr1", query_start, query_end);
        
        // Manually compute expected overlaps
        let expected: Vec<&(u64, u64)> = intervals.iter()
            .filter(|&&(s, e)| s < query_end && e > query_start)
            .collect();
        
        // Verify count matches
        prop_assert_eq!(
            results.len(), 
            expected.len(),
            "Query [{}, {}) returned {} results, expected {} overlapping intervals",
            query_start, query_end, results.len(), expected.len()
        );
        
        // Verify each result actually overlaps
        for result in &results {
            let overlaps = result.target_start < query_end && result.target_end > query_start;
            prop_assert!(
                overlaps,
                "Result [{}, {}) does not overlap query [{}, {})",
                result.target_start, result.target_end, query_start, query_end
            );
        }
    }
    
    /// Property: Query returns only intervals that overlap
    #[test]
    fn prop_query_returns_only_overlapping(
        intervals in arb_intervals("chr1".to_string(), 10000),
        query_start in 0u64..10000,
        query_len in 10u64..1000,
    ) {
        let query_end = query_start + query_len;
        
        let mut intervals_map = HashMap::new();
        intervals_map.insert("chr1".to_string(), intervals.clone());
        let chain_file = create_chain_file(intervals_map);
        let index = ChainIndex::from_chain_data(chain_file);
        
        let results = index.query("chr1", query_start, query_end);
        
        // Every result must overlap the query
        for result in results {
            let overlaps = result.target_start < query_end && result.target_end > query_start;
            prop_assert!(
                overlaps,
                "Returned interval [{}, {}) does not overlap query [{}, {})",
                result.target_start, result.target_end, query_start, query_end
            );
        }
    }
    
    /// Property: Query beyond all intervals returns empty results
    #[test]
    fn prop_query_beyond_intervals_returns_empty(
        intervals in arb_intervals("chr1".to_string(), 5000),
    ) {
        let mut intervals_map = HashMap::new();
        intervals_map.insert("chr1".to_string(), intervals);
        let chain_file = create_chain_file(intervals_map);
        let index = ChainIndex::from_chain_data(chain_file);
        
        // Query far beyond any possible interval (intervals are in 0..5000 range)
        let results = index.query("chr1", 100000, 100100);
        prop_assert!(
            results.is_empty(),
            "Query beyond all intervals should return no results, got {}",
            results.len()
        );
    }
    
    /// Property: Query on non-existent chromosome returns empty
    #[test]
    fn prop_nonexistent_chrom_returns_empty(
        intervals in arb_intervals("chr1".to_string(), 10000),
        query_start in 0u64..10000,
        query_len in 10u64..1000,
    ) {
        let query_end = query_start + query_len;
        
        let mut intervals_map = HashMap::new();
        intervals_map.insert("chr1".to_string(), intervals);
        let chain_file = create_chain_file(intervals_map);
        let index = ChainIndex::from_chain_data(chain_file);
        
        // Query non-existent chromosome
        let results = index.query("chrNONE", query_start, query_end);
        prop_assert!(
            results.is_empty(),
            "Query on non-existent chromosome should return empty"
        );
    }
    
    /// Property: Point query (single position) works correctly
    #[test]
    fn prop_point_query(
        intervals in arb_intervals("chr1".to_string(), 10000),
        pos in 0u64..10000,
    ) {
        let mut intervals_map = HashMap::new();
        intervals_map.insert("chr1".to_string(), intervals.clone());
        let chain_file = create_chain_file(intervals_map);
        let index = ChainIndex::from_chain_data(chain_file);
        
        // Point query: [pos, pos+1)
        let results = index.query("chr1", pos, pos + 1);
        
        // Count intervals containing this position
        let expected_count = intervals.iter()
            .filter(|&&(s, e)| s <= pos && pos < e)
            .count();
        
        prop_assert_eq!(
            results.len(),
            expected_count,
            "Point query at {} returned {} results, expected {}",
            pos, results.len(), expected_count
        );
    }
}


/// **Property 3: 染色体名称变体等价性**
/// 
/// For any chromosome name query, the index should return the same results
/// whether queried as "chr1", "1", or "CHR1" (case-insensitive, chr-prefix agnostic).
///
/// **Validates: Requirements 2.3**
mod chrom_name_properties {
    use super::*;
    
    /// Generate chromosome name variants
    fn chrom_variants(base: &str) -> Vec<String> {
        vec![
            format!("chr{}", base),
            format!("Chr{}", base),
            format!("CHR{}", base),
            base.to_string(),
            base.to_uppercase(),
            base.to_lowercase(),
        ]
    }
    
    proptest! {
        #![proptest_config(ProptestConfig::with_cases(50))]
        
        /// Property 3: Different chromosome name variants return same results
        #[test]
        fn prop_chrom_variants_equivalent(
            chrom_num in 1u8..=22,
            query_start in 0u64..10000,
            query_len in 10u64..1000,
        ) {
            let query_end = query_start + query_len;
            let base_chrom = format!("{}", chrom_num);
            
            // Create index with chr-prefixed chromosome
            let intervals = vec![(100u64, 200u64), (500u64, 600u64), (1000u64, 2000u64)];
            let mut intervals_map = HashMap::new();
            intervals_map.insert(format!("chr{}", chrom_num), intervals);
            let chain_file = create_chain_file(intervals_map);
            let index = ChainIndex::from_chain_data(chain_file);
            
            // Query with different variants
            let variants = chrom_variants(&base_chrom);
            let mut results: Vec<Vec<_>> = Vec::new();
            
            for variant in &variants {
                let r = index.query(variant, query_start, query_end);
                results.push(r.into_iter().map(|v| (v.target_start, v.target_end)).collect());
            }
            
            // All variants should return the same results
            let first = &results[0];
            for (i, result) in results.iter().enumerate().skip(1) {
                prop_assert_eq!(
                    first, result,
                    "Variant '{}' returned different results than '{}'",
                    variants[i], variants[0]
                );
            }
        }
        
        /// Property: X/Y/M chromosomes work with variants
        #[test]
        fn prop_special_chrom_variants(
            chrom in prop_oneof![Just("X"), Just("Y"), Just("M")],
            query_start in 0u64..10000,
            query_len in 10u64..1000,
        ) {
            let query_end = query_start + query_len;
            
            // Create index with chr-prefixed chromosome
            let intervals = vec![(100u64, 200u64), (500u64, 600u64)];
            let mut intervals_map = HashMap::new();
            intervals_map.insert(format!("chr{}", chrom), intervals);
            let chain_file = create_chain_file(intervals_map);
            let index = ChainIndex::from_chain_data(chain_file);
            
            // Query with different variants
            let variants = vec![
                format!("chr{}", chrom),
                format!("Chr{}", chrom),
                chrom.to_string(),
                chrom.to_lowercase(),
            ];
            
            let mut results: Vec<Vec<_>> = Vec::new();
            for variant in &variants {
                let r = index.query(variant, query_start, query_end);
                results.push(r.into_iter().map(|v| (v.target_start, v.target_end)).collect());
            }
            
            // All variants should return the same results
            let first = &results[0];
            for (i, result) in results.iter().enumerate().skip(1) {
                prop_assert_eq!(
                    first, result,
                    "Variant '{}' returned different results than '{}'",
                    variants[i], variants[0]
                );
            }
        }
        
        /// Property: has_chrom works with all variants
        #[test]
        fn prop_has_chrom_variants(chrom_num in 1u8..=22) {
            let base_chrom = format!("{}", chrom_num);
            
            // Create index with chr-prefixed chromosome
            let intervals = vec![(100u64, 200u64)];
            let mut intervals_map = HashMap::new();
            intervals_map.insert(format!("chr{}", chrom_num), intervals);
            let chain_file = create_chain_file(intervals_map);
            let index = ChainIndex::from_chain_data(chain_file);
            
            // All variants should be found
            let variants = chrom_variants(&base_chrom);
            for variant in &variants {
                prop_assert!(
                    index.has_chrom(variant),
                    "has_chrom('{}') should return true",
                    variant
                );
            }
            
            // Non-existent chromosome should not be found
            prop_assert!(!index.has_chrom("chrNONE"));
        }
    }
}
