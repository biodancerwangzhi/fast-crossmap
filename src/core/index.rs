//! Interval index for efficient coordinate queries
//!
//! Uses rust-lapper for O(log n + k) interval queries.

use crate::core::chain::{parse_chain_file, ChainFile, ChainParseError};
use crate::core::Strand;
use rust_lapper::{Interval, Lapper};
use std::collections::HashMap;
use std::path::Path;

/// Value stored in each interval - target mapping information
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntervalValue {
    /// Target chromosome name
    pub target_chrom: String,
    /// Target start position (0-based, already flipped for negative strand)
    pub target_start: u64,
    /// Target end position (exclusive)
    pub target_end: u64,
    /// Target strand direction
    pub target_strand: Strand,
    /// Source chromosome name (for reference)
    pub source_chrom: String,
}

/// Type alias for chain intervals
pub type ChainInterval = Interval<u64, IntervalValue>;

/// Interval index organized by source chromosome
/// 
/// Provides O(log n + k) interval queries where n is the number of
/// intervals and k is the number of overlapping results.
pub struct ChainIndex {
    /// Source chromosome -> interval tree (using Lapper)
    maps: HashMap<String, Lapper<u64, IntervalValue>>,
    /// Target chromosome sizes
    pub target_sizes: HashMap<String, u64>,
    /// Source chromosome sizes
    pub source_sizes: HashMap<String, u64>,
    /// Normalized chromosome name mapping (lowercase -> original)
    chrom_aliases: HashMap<String, String>,
}


impl ChainIndex {
    /// Build index from a chain file
    /// 
    /// Automatically handles gzip and bzip2 compression.
    /// 
    /// # Example
    /// ```ignore
    /// let index = ChainIndex::from_chain_file("hg19ToHg38.chain.gz")?;
    /// ```
    pub fn from_chain_file<P: AsRef<Path>>(path: P) -> Result<Self, ChainParseError> {
        let chain_file = parse_chain_file(path.as_ref())?;
        Ok(Self::from_chain_data(chain_file))
    }
    
    /// Build index from parsed chain data
    pub fn from_chain_data(chain_file: ChainFile) -> Self {
        // Group blocks by source chromosome
        let mut blocks_by_chrom: HashMap<String, Vec<ChainInterval>> = HashMap::new();
        
        for block in chain_file.blocks {
            let interval = Interval {
                start: block.source_start,
                stop: block.source_end,
                val: IntervalValue {
                    target_chrom: block.target_chrom,
                    target_start: block.target_start,
                    target_end: block.target_end,
                    target_strand: block.target_strand,
                    source_chrom: block.source_chrom.clone(),
                },
            };
            
            blocks_by_chrom
                .entry(block.source_chrom)
                .or_default()
                .push(interval);
        }
        
        // Build Lapper for each chromosome
        let mut maps = HashMap::new();
        let mut chrom_aliases = HashMap::new();
        
        for (chrom, intervals) in blocks_by_chrom {
            // Store chromosome aliases for flexible lookup
            let normalized = normalize_chrom_key(&chrom);
            chrom_aliases.insert(normalized, chrom.clone());
            
            // Build the interval tree
            maps.insert(chrom, Lapper::new(intervals));
        }
        
        Self {
            maps,
            target_sizes: chain_file.target_chrom_sizes,
            source_sizes: chain_file.source_chrom_sizes,
            chrom_aliases,
        }
    }
    
    /// Query intervals overlapping the given range
    /// 
    /// Automatically handles chromosome name variants (chr1, 1, CHR1).
    /// Returns references to IntervalValue for each overlapping block.
    pub fn query(&self, chrom: &str, start: u64, end: u64) -> Vec<&IntervalValue> {
        let lapper = self.find_lapper(chrom);
        
        match lapper {
            Some(l) => l.find(start, end).map(|iv| &iv.val).collect(),
            None => vec![],
        }
    }
    
    /// Query intervals and return full Interval structs
    pub fn query_intervals(&self, chrom: &str, start: u64, end: u64) -> Vec<&ChainInterval> {
        let lapper = self.find_lapper(chrom);
        
        match lapper {
            Some(l) => l.find(start, end).collect(),
            None => vec![],
        }
    }
    
    /// Find the Lapper for a chromosome, trying different naming styles
    fn find_lapper(&self, chrom: &str) -> Option<&Lapper<u64, IntervalValue>> {
        // Try exact match first
        if let Some(l) = self.maps.get(chrom) {
            return Some(l);
        }
        
        // Try normalized lookup
        let normalized = normalize_chrom_key(chrom);
        if let Some(original) = self.chrom_aliases.get(&normalized) {
            return self.maps.get(original);
        }
        
        // Try common variants
        let variants = [
            chrom.to_string(),
            chrom.replace("chr", ""),
            chrom.replace("Chr", ""),
            chrom.replace("CHR", ""),
            format!("chr{}", chrom),
            format!("Chr{}", chrom),
        ];
        
        for variant in &variants {
            if let Some(l) = self.maps.get(variant) {
                return Some(l);
            }
        }
        
        None
    }
    
    /// Get the canonical chromosome name used in the index
    pub fn get_canonical_chrom(&self, chrom: &str) -> Option<&str> {
        if self.maps.contains_key(chrom) {
            return self.maps.keys().find(|k| *k == chrom).map(|s| s.as_str());
        }
        
        let normalized = normalize_chrom_key(chrom);
        self.chrom_aliases.get(&normalized).map(|s| s.as_str())
    }
    
    /// Check if a chromosome exists in the index
    pub fn has_chrom(&self, chrom: &str) -> bool {
        self.find_lapper(chrom).is_some()
    }
    
    /// Get all source chromosome names
    pub fn source_chroms(&self) -> impl Iterator<Item = &str> {
        self.maps.keys().map(|s| s.as_str())
    }
    
    /// Get the number of intervals for a chromosome
    pub fn interval_count(&self, chrom: &str) -> usize {
        self.find_lapper(chrom).map(|l| l.len()).unwrap_or(0)
    }
    
    /// Get total number of intervals across all chromosomes
    pub fn total_intervals(&self) -> usize {
        self.maps.values().map(|l| l.len()).sum()
    }
    
    /// Get target chromosome size
    pub fn target_chrom_size(&self, chrom: &str) -> Option<u64> {
        self.target_sizes.get(chrom).copied()
            .or_else(|| self.target_sizes.get(&chrom.replace("chr", "")).copied())
            .or_else(|| self.target_sizes.get(&format!("chr{}", chrom)).copied())
    }
    
    /// Get source chromosome size
    pub fn source_chrom_size(&self, chrom: &str) -> Option<u64> {
        self.source_sizes.get(chrom).copied()
            .or_else(|| self.source_sizes.get(&chrom.replace("chr", "")).copied())
            .or_else(|| self.source_sizes.get(&format!("chr{}", chrom)).copied())
    }
}

/// Normalize chromosome name for flexible matching
/// 
/// Converts to lowercase and removes common prefixes.
fn normalize_chrom_key(chrom: &str) -> String {
    let lower = chrom.to_lowercase();
    if lower.starts_with("chr") {
        lower[3..].to_string()
    } else {
        lower
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::chain::parse_chain_bytes;
    
    fn create_test_index() -> ChainIndex {
        let chain_data = b"\
chain 1000 chr1 1000 + 100 500 chr1 1000 + 100 500 1
100 50 50
100 50 50
100

chain 500 chr2 2000 + 0 200 chr2 2000 + 0 200 2
100 50 50
50
";
        let chain_file = parse_chain_bytes(chain_data).unwrap();
        ChainIndex::from_chain_data(chain_file)
    }
    
    #[test]
    fn test_index_creation() {
        let index = create_test_index();
        
        assert!(index.has_chrom("chr1"));
        assert!(index.has_chrom("chr2"));
        assert!(!index.has_chrom("chr3"));
        
        assert_eq!(index.total_intervals(), 5); // 3 from chr1 + 2 from chr2
    }
    
    #[test]
    fn test_query_basic() {
        let index = create_test_index();
        
        // Query chr1 at position 150 (should hit first block: 100-200)
        let results = index.query("chr1", 150, 160);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].target_start, 100);
        assert_eq!(results[0].target_end, 200);
    }
    
    #[test]
    fn test_query_no_overlap() {
        let index = create_test_index();
        
        // Query chr1 at position 50 (before any blocks)
        let results = index.query("chr1", 50, 60);
        assert!(results.is_empty());
    }
    
    #[test]
    fn test_query_multiple_overlaps() {
        let index = create_test_index();
        
        // Query chr1 spanning multiple blocks
        let results = index.query("chr1", 100, 500);
        assert_eq!(results.len(), 3); // All 3 blocks from chr1
    }
    
    #[test]
    fn test_chrom_name_variants() {
        let index = create_test_index();
        
        // Should find chr1 with different naming styles
        assert!(index.has_chrom("chr1"));
        assert!(index.has_chrom("1"));
        assert!(index.has_chrom("CHR1"));
        assert!(index.has_chrom("Chr1"));
        
        // Query should work with variants
        let results1 = index.query("chr1", 150, 160);
        let results2 = index.query("1", 150, 160);
        assert_eq!(results1.len(), results2.len());
    }
    
    #[test]
    fn test_chrom_sizes() {
        let index = create_test_index();
        
        assert_eq!(index.target_chrom_size("chr1"), Some(1000));
        assert_eq!(index.target_chrom_size("chr2"), Some(2000));
        assert_eq!(index.target_chrom_size("chr3"), None);
        
        assert_eq!(index.source_chrom_size("chr1"), Some(1000));
        assert_eq!(index.source_chrom_size("chr2"), Some(2000));
    }
    
    #[test]
    fn test_canonical_chrom() {
        let index = create_test_index();
        
        assert_eq!(index.get_canonical_chrom("chr1"), Some("chr1"));
        assert_eq!(index.get_canonical_chrom("1"), Some("chr1"));
        assert_eq!(index.get_canonical_chrom("chr3"), None);
    }
    
    #[test]
    fn test_source_chroms() {
        let index = create_test_index();
        
        let chroms: Vec<&str> = index.source_chroms().collect();
        assert!(chroms.contains(&"chr1"));
        assert!(chroms.contains(&"chr2"));
        assert_eq!(chroms.len(), 2);
    }
    
    #[test]
    fn test_interval_count() {
        let index = create_test_index();
        
        assert_eq!(index.interval_count("chr1"), 3);
        assert_eq!(index.interval_count("chr2"), 2);
        assert_eq!(index.interval_count("chr3"), 0);
    }
}

#[cfg(test)]
mod integration_tests {
    use super::*;
    use std::path::PathBuf;
    
    #[test]
    fn test_load_real_chain_file() {
        let chain_path = PathBuf::from("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz");
        
        if !chain_path.exists() {
            eprintln!("Skipping test: chain file not found");
            return;
        }
        
        let start = std::time::Instant::now();
        let index = ChainIndex::from_chain_file(&chain_path);
        let elapsed = start.elapsed();
        
        assert!(index.is_ok(), "Failed to load chain file: {:?}", index.err());
        let index = index.unwrap();
        
        eprintln!("Loaded {} intervals in {:?}", index.total_intervals(), elapsed);
        eprintln!("Source chromosomes: {}", index.maps.len());
        
        // Should load in reasonable time (< 5 seconds)
        assert!(elapsed.as_secs() < 10, "Loading took too long: {:?}", elapsed);
        
        // Should have chr1
        assert!(index.has_chrom("chr1") || index.has_chrom("1"));
        
        // Test a query
        let results = index.query("chr1", 1000000, 1000100);
        eprintln!("Query chr1:1000000-1000100 returned {} results", results.len());
    }
}
