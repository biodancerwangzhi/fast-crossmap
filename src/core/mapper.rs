//! Coordinate mapping algorithm
//!
//! Maps coordinates from source to target genome assembly.
//! 
//! The mapping algorithm follows CrossMap's logic:
//! 1. Query the interval index for overlapping chain blocks
//! 2. For each overlapping block, compute the intersection
//! 3. Calculate target coordinates using offset formulas
//! 4. Handle strand direction combinations

use crate::core::index::IntervalValue;
use crate::core::ChainIndex;

/// Compatibility mode for CrossMap behavior
/// 
/// Controls how edge cases and ambiguous situations are handled during coordinate mapping.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CompatMode {
    /// Default mode: use FastCrossMap's improved logic
    /// - May produce slightly different results in edge cases
    /// - Optimized for performance
    #[default]
    Improved,
    /// Strict mode: exactly match CrossMap behavior
    /// - Bug-for-bug compatibility with Python CrossMap
    /// - Handles edge cases identically to CrossMap
    /// - Use for validation and comparison testing
    Strict,
}

impl CompatMode {
    /// Parse from string (for CLI argument)
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "improved" | "default" => Some(CompatMode::Improved),
            "strict" | "crossmap" => Some(CompatMode::Strict),
            _ => None,
        }
    }
    
    /// Check if strict mode is enabled
    pub fn is_strict(&self) -> bool {
        matches!(self, CompatMode::Strict)
    }
}

/// Strand orientation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Hash)]
pub enum Strand {
    #[default]
    Plus,
    Minus,
}

impl Strand {
    /// Get the complement strand
    /// 
    /// # Examples
    /// ```
    /// use fast_crossmap::core::Strand;
    /// assert_eq!(Strand::Plus.complement(), Strand::Minus);
    /// assert_eq!(Strand::Minus.complement(), Strand::Plus);
    /// ```
    pub fn complement(&self) -> Self {
        match self {
            Strand::Plus => Strand::Minus,
            Strand::Minus => Strand::Plus,
        }
    }

    /// Parse strand from char
    /// 
    /// # Examples
    /// ```
    /// use fast_crossmap::core::Strand;
    /// assert_eq!(Strand::from_char('+'), Some(Strand::Plus));
    /// assert_eq!(Strand::from_char('-'), Some(Strand::Minus));
    /// assert_eq!(Strand::from_char('.'), None);
    /// ```
    pub fn from_char(c: char) -> Option<Self> {
        match c {
            '+' => Some(Strand::Plus),
            '-' => Some(Strand::Minus),
            _ => None,
        }
    }

    /// Parse strand from byte
    pub fn from_byte(b: u8) -> Option<Self> {
        match b {
            b'+' => Some(Strand::Plus),
            b'-' => Some(Strand::Minus),
            _ => None,
        }
    }

    /// Convert to char
    pub fn to_char(&self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
        }
    }

    /// Convert to byte
    pub fn to_byte(&self) -> u8 {
        match self {
            Strand::Plus => b'+',
            Strand::Minus => b'-',
        }
    }

    /// Combine two strands (for query strand + target strand)
    /// 
    /// When mapping coordinates, the final strand is determined by:
    /// - Plus + Plus = Plus
    /// - Plus + Minus = Minus
    /// - Minus + Plus = Minus
    /// - Minus + Minus = Plus
    /// 
    /// # Examples
    /// ```
    /// use fast_crossmap::core::Strand;
    /// assert_eq!(Strand::Plus.combine(Strand::Plus), Strand::Plus);
    /// assert_eq!(Strand::Plus.combine(Strand::Minus), Strand::Minus);
    /// assert_eq!(Strand::Minus.combine(Strand::Plus), Strand::Minus);
    /// assert_eq!(Strand::Minus.combine(Strand::Minus), Strand::Plus);
    /// ```
    pub fn combine(&self, other: Strand) -> Strand {
        match (self, other) {
            (Strand::Plus, Strand::Plus) => Strand::Plus,
            (Strand::Plus, Strand::Minus) => Strand::Minus,
            (Strand::Minus, Strand::Plus) => Strand::Minus,
            (Strand::Minus, Strand::Minus) => Strand::Plus,
        }
    }
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

/// Chromosome ID style for output formatting
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ChromStyle {
    /// Keep same style as input query
    #[default]
    AsIs,
    /// Short style without "chr" prefix: "1", "2", "X", "Y", "M"
    Short,
    /// Long style with "chr" prefix: "chr1", "chr2", "chrX", "chrY", "chrM"
    Long,
}

impl ChromStyle {
    /// Parse from string (for CLI argument)
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "asis" | "as_is" | "as-is" => Some(ChromStyle::AsIs),
            "short" | "s" => Some(ChromStyle::Short),
            "long" | "l" => Some(ChromStyle::Long),
            _ => None,
        }
    }
}

/// Update chromosome ID according to the specified style
/// 
/// # Arguments
/// * `chrom` - Original chromosome name
/// * `style` - Target chromosome style
/// 
/// # Returns
/// Chromosome name formatted according to the style
/// 
/// # Examples
/// ```
/// use fast_crossmap::core::{ChromStyle, update_chrom_id};
/// 
/// // Short style removes "chr" prefix
/// assert_eq!(update_chrom_id("chr1", ChromStyle::Short), "1");
/// assert_eq!(update_chrom_id("chrX", ChromStyle::Short), "X");
/// assert_eq!(update_chrom_id("1", ChromStyle::Short), "1");
/// 
/// // Long style adds "chr" prefix
/// assert_eq!(update_chrom_id("1", ChromStyle::Long), "chr1");
/// assert_eq!(update_chrom_id("X", ChromStyle::Long), "chrX");
/// assert_eq!(update_chrom_id("chr1", ChromStyle::Long), "chr1");
/// 
/// // AsIs keeps original
/// assert_eq!(update_chrom_id("chr1", ChromStyle::AsIs), "chr1");
/// assert_eq!(update_chrom_id("1", ChromStyle::AsIs), "1");
/// ```
pub fn update_chrom_id(chrom: &str, style: ChromStyle) -> String {
    match style {
        ChromStyle::AsIs => chrom.to_string(),
        ChromStyle::Short => {
            // Remove "chr" prefix if present (case-insensitive)
            if chrom.len() > 3 && chrom[..3].eq_ignore_ascii_case("chr") {
                chrom[3..].to_string()
            } else {
                chrom.to_string()
            }
        }
        ChromStyle::Long => {
            // Add "chr" prefix if not present
            if chrom.len() > 3 && chrom[..3].eq_ignore_ascii_case("chr") {
                // Already has chr prefix, normalize to lowercase "chr"
                format!("chr{}", &chrom[3..])
            } else {
                format!("chr{}", chrom)
            }
        }
    }
}

/// Normalize chromosome name for lookup (handles chr1/1/CHR1 variants)
/// 
/// Returns a canonical form for comparison purposes.
/// 
/// # Examples
/// ```
/// use fast_crossmap::core::normalize_chrom;
/// 
/// // All these should normalize to the same value
/// assert_eq!(normalize_chrom("chr1"), normalize_chrom("1"));
/// assert_eq!(normalize_chrom("CHR1"), normalize_chrom("chr1"));
/// assert_eq!(normalize_chrom("Chr1"), normalize_chrom("1"));
/// 
/// // Special chromosomes
/// assert_eq!(normalize_chrom("chrX"), normalize_chrom("X"));
/// assert_eq!(normalize_chrom("chrM"), normalize_chrom("MT"));
/// ```
pub fn normalize_chrom(chrom: &str) -> String {
    // Remove "chr" prefix if present (case-insensitive)
    let without_prefix = if chrom.len() > 3 && chrom[..3].eq_ignore_ascii_case("chr") {
        &chrom[3..]
    } else {
        chrom
    };
    
    // Normalize to uppercase
    let upper = without_prefix.to_uppercase();
    
    // Handle MT/M equivalence
    if upper == "M" {
        "MT".to_string()
    } else if upper == "MT" {
        "MT".to_string()
    } else {
        upper
    }
}

/// Check if two chromosome names are equivalent
/// 
/// # Examples
/// ```
/// use fast_crossmap::core::chroms_equivalent;
/// 
/// assert!(chroms_equivalent("chr1", "1"));
/// assert!(chroms_equivalent("CHR1", "chr1"));
/// assert!(chroms_equivalent("chrM", "MT"));
/// assert!(!chroms_equivalent("chr1", "chr2"));
/// ```
pub fn chroms_equivalent(chrom1: &str, chrom2: &str) -> bool {
    normalize_chrom(chrom1) == normalize_chrom(chrom2)
}

/// Result of coordinate mapping
#[derive(Debug, Clone, PartialEq)]
pub struct MapResult {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
}

/// A single mapping segment (source region + target region)
#[derive(Debug, Clone, PartialEq)]
pub struct MappingSegment {
    /// Source region that was mapped
    pub source: MapResult,
    /// Target region after mapping
    pub target: MapResult,
}

/// Compute the intersection of two intervals on the same chromosome
/// 
/// Returns None if intervals don't overlap or are on different chromosomes.
/// 
/// # Arguments
/// * `start1`, `end1` - First interval [start1, end1)
/// * `start2`, `end2` - Second interval [start2, end2)
/// 
/// # Returns
/// The intersection interval (start, end) or None if no overlap
#[inline]
pub fn intersect_intervals(start1: u64, end1: u64, start2: u64, end2: u64) -> Option<(u64, u64)> {
    if start1 >= end2 || end1 <= start2 {
        return None;
    }
    Some((start1.max(start2), end1.min(end2)))
}

/// Coordinate mapper using chain index
pub struct CoordinateMapper {
    index: ChainIndex,
    chrom_style: ChromStyle,
    compat_mode: CompatMode,
}

impl CoordinateMapper {
    pub fn new(index: ChainIndex, chrom_style: ChromStyle) -> Self {
        Self { 
            index, 
            chrom_style,
            compat_mode: CompatMode::default(),
        }
    }
    
    /// Create a new mapper with specified compatibility mode
    pub fn with_compat_mode(index: ChainIndex, chrom_style: ChromStyle, compat_mode: CompatMode) -> Self {
        Self { 
            index, 
            chrom_style,
            compat_mode,
        }
    }
    
    /// Set the compatibility mode
    pub fn set_compat_mode(&mut self, mode: CompatMode) {
        self.compat_mode = mode;
    }
    
    /// Get the compatibility mode
    pub fn compat_mode(&self) -> CompatMode {
        self.compat_mode
    }

    /// Get the chromosome style
    pub fn chrom_style(&self) -> ChromStyle {
        self.chrom_style
    }
    
    /// Get a reference to the underlying index
    pub fn index(&self) -> &ChainIndex {
        &self.index
    }
    
    /// Get target chromosome sizes
    pub fn target_sizes(&self) -> &std::collections::HashMap<String, u64> {
        &self.index.target_sizes
    }

    /// Map coordinates from source to target assembly
    /// 
    /// Returns None if the chromosome is not found in the index.
    /// Returns an empty Vec if no overlapping chain blocks are found.
    /// Returns a Vec of MappingSegment for each overlapping block.
    /// 
    /// # Arguments
    /// * `chrom` - Source chromosome name
    /// * `start` - Start position (0-based, inclusive)
    /// * `end` - End position (0-based, exclusive)
    /// * `strand` - Query strand direction
    /// 
    /// # Algorithm
    /// For each overlapping chain block:
    /// 1. Compute intersection of query with source block
    /// 2. Calculate left_offset = intersection_start - source_block_start
    /// 3. Calculate size = intersection_end - intersection_start
    /// 4. For positive target strand: target_start = t_start + left_offset
    /// 5. For negative target strand: target_start = t_end - left_offset - size
    /// 6. Combine query strand with target strand for final strand
    pub fn map(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        strand: Strand,
    ) -> Option<Vec<MappingSegment>> {
        // Check if chromosome exists
        if !self.index.has_chrom(chrom) {
            return None;
        }
        
        // Query overlapping intervals
        let intervals = self.index.query_intervals(chrom, start, end);
        
        if intervals.is_empty() {
            return Some(vec![]);
        }
        
        let mut results = Vec::with_capacity(intervals.len());
        
        for interval in intervals {
            // Source block coordinates from the interval
            let s_start = interval.start;
            let s_end = interval.stop;
            let target_info = &interval.val;
            
            // Compute intersection of query with source block
            let (real_start, real_end) = match intersect_intervals(start, end, s_start, s_end) {
                Some(intersection) => intersection,
                None => continue, // Should not happen since we queried overlapping
            };
            
            // Calculate mapping parameters
            let left_offset = real_start - s_start;
            let size = real_end - real_start;
            
            // Calculate target coordinates based on target strand
            let (target_start, target_end) = self.calculate_target_coords(
                target_info,
                left_offset,
                size,
            );
            
            // Combine strands: query_strand XOR target_strand
            let final_strand = strand.combine(target_info.target_strand);
            
            // Format chromosome according to style
            // For ChromStyle::AsIs, we want to preserve the user's input style
            // So if user queries "chr1", output should also use "chr1" style
            let target_chrom = match self.chrom_style {
                ChromStyle::AsIs => {
                    // Preserve user's chromosome naming style
                    // If user used "chr" prefix, add it to target; otherwise remove it
                    let user_has_chr = chrom.len() > 3 && chrom[..3].eq_ignore_ascii_case("chr");
                    let target_has_chr = target_info.target_chrom.len() > 3 
                        && target_info.target_chrom[..3].eq_ignore_ascii_case("chr");
                    
                    if user_has_chr && !target_has_chr {
                        // User used "chr1", chain has "1" -> output "chr1"
                        format!("chr{}", target_info.target_chrom)
                    } else if !user_has_chr && target_has_chr {
                        // User used "1", chain has "chr1" -> output "1"
                        target_info.target_chrom[3..].to_string()
                    } else {
                        // Same style, use as-is
                        target_info.target_chrom.clone()
                    }
                }
                _ => update_chrom_id(&target_info.target_chrom, self.chrom_style),
            };
            let source_chrom = match self.chrom_style {
                ChromStyle::AsIs => chrom.to_string(),
                _ => update_chrom_id(chrom, self.chrom_style),
            };
            
            results.push(MappingSegment {
                source: MapResult {
                    chrom: source_chrom,
                    start: real_start,
                    end: real_end,
                    strand,
                },
                target: MapResult {
                    chrom: target_chrom,
                    start: target_start,
                    end: target_end,
                    strand: final_strand,
                },
            });
        }
        
        Some(results)
    }
    
    /// Calculate target coordinates based on strand direction
    /// 
    /// For positive strand: target_start = t_start + left_offset
    /// For negative strand: coordinates are already flipped in the index,
    ///   so we calculate from t_end backwards: target_start = t_end - left_offset - size
    #[inline]
    fn calculate_target_coords(
        &self,
        target_info: &IntervalValue,
        left_offset: u64,
        size: u64,
    ) -> (u64, u64) {
        match target_info.target_strand {
            Strand::Plus => {
                // Positive strand: simple offset from start
                let target_start = target_info.target_start + left_offset;
                let target_end = target_start + size;
                (target_start, target_end)
            }
            Strand::Minus => {
                // Negative strand: coordinates in index are already flipped
                // We need to map from the "end" of the block going backwards
                // The index stores flipped coordinates, so:
                // - target_start in index = original target_size - (original_target_pos + block_size)
                // - target_end in index = original target_size - original_target_pos
                // 
                // For a query at left_offset from source_start:
                // The corresponding position in target is at (target_end - left_offset - size)
                let target_start = target_info.target_end - left_offset - size;
                let target_end = target_start + size;
                (target_start, target_end)
            }
        }
    }
    
    /// Map a single position (useful for VCF)
    /// 
    /// Returns the first mapping result for a single base position.
    pub fn map_single(
        &self,
        chrom: &str,
        pos: u64,
        strand: Strand,
    ) -> Option<MappingSegment> {
        let results = self.map(chrom, pos, pos + 1, strand)?;
        results.into_iter().next()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::chain::parse_chain_bytes;
    use crate::core::ChainIndex;

    fn create_test_index() -> ChainIndex {
        // Create a simple chain file for testing
        // Source: chr1:100-500, Target: chr1:100-500 (positive strand)
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
    fn test_strand_complement() {
        assert_eq!(Strand::Plus.complement(), Strand::Minus);
        assert_eq!(Strand::Minus.complement(), Strand::Plus);
    }

    #[test]
    fn test_strand_complement_involution() {
        // complement(complement(x)) == x
        assert_eq!(Strand::Plus.complement().complement(), Strand::Plus);
        assert_eq!(Strand::Minus.complement().complement(), Strand::Minus);
    }

    #[test]
    fn test_strand_from_char() {
        assert_eq!(Strand::from_char('+'), Some(Strand::Plus));
        assert_eq!(Strand::from_char('-'), Some(Strand::Minus));
        assert_eq!(Strand::from_char('.'), None);
        assert_eq!(Strand::from_char('x'), None);
    }

    #[test]
    fn test_strand_to_char() {
        assert_eq!(Strand::Plus.to_char(), '+');
        assert_eq!(Strand::Minus.to_char(), '-');
    }

    #[test]
    fn test_strand_combine() {
        // Same strand = Plus
        assert_eq!(Strand::Plus.combine(Strand::Plus), Strand::Plus);
        assert_eq!(Strand::Minus.combine(Strand::Minus), Strand::Plus);
        
        // Different strand = Minus
        assert_eq!(Strand::Plus.combine(Strand::Minus), Strand::Minus);
        assert_eq!(Strand::Minus.combine(Strand::Plus), Strand::Minus);
    }

    #[test]
    fn test_strand_display() {
        assert_eq!(format!("{}", Strand::Plus), "+");
        assert_eq!(format!("{}", Strand::Minus), "-");
    }

    #[test]
    fn test_chrom_style_from_str() {
        assert_eq!(ChromStyle::from_str("asis"), Some(ChromStyle::AsIs));
        assert_eq!(ChromStyle::from_str("short"), Some(ChromStyle::Short));
        assert_eq!(ChromStyle::from_str("long"), Some(ChromStyle::Long));
        assert_eq!(ChromStyle::from_str("LONG"), Some(ChromStyle::Long));
        assert_eq!(ChromStyle::from_str("invalid"), None);
    }

    #[test]
    fn test_update_chrom_id_short() {
        assert_eq!(update_chrom_id("chr1", ChromStyle::Short), "1");
        assert_eq!(update_chrom_id("chrX", ChromStyle::Short), "X");
        assert_eq!(update_chrom_id("chrM", ChromStyle::Short), "M");
        assert_eq!(update_chrom_id("CHR1", ChromStyle::Short), "1");
        assert_eq!(update_chrom_id("1", ChromStyle::Short), "1");
        assert_eq!(update_chrom_id("X", ChromStyle::Short), "X");
    }

    #[test]
    fn test_update_chrom_id_long() {
        assert_eq!(update_chrom_id("1", ChromStyle::Long), "chr1");
        assert_eq!(update_chrom_id("X", ChromStyle::Long), "chrX");
        assert_eq!(update_chrom_id("M", ChromStyle::Long), "chrM");
        assert_eq!(update_chrom_id("chr1", ChromStyle::Long), "chr1");
        assert_eq!(update_chrom_id("CHR1", ChromStyle::Long), "chr1");
    }

    #[test]
    fn test_update_chrom_id_asis() {
        assert_eq!(update_chrom_id("chr1", ChromStyle::AsIs), "chr1");
        assert_eq!(update_chrom_id("1", ChromStyle::AsIs), "1");
        assert_eq!(update_chrom_id("CHR1", ChromStyle::AsIs), "CHR1");
    }

    #[test]
    fn test_normalize_chrom() {
        assert_eq!(normalize_chrom("chr1"), "1");
        assert_eq!(normalize_chrom("CHR1"), "1");
        assert_eq!(normalize_chrom("1"), "1");
        assert_eq!(normalize_chrom("chrX"), "X");
        assert_eq!(normalize_chrom("x"), "X");
        
        // MT/M equivalence
        assert_eq!(normalize_chrom("chrM"), "MT");
        assert_eq!(normalize_chrom("M"), "MT");
        assert_eq!(normalize_chrom("MT"), "MT");
        assert_eq!(normalize_chrom("chrMT"), "MT");
    }
    
    #[test]
    fn test_chroms_equivalent() {
        assert!(chroms_equivalent("chr1", "1"));
        assert!(chroms_equivalent("CHR1", "chr1"));
        assert!(chroms_equivalent("chrM", "MT"));
        assert!(chroms_equivalent("M", "chrMT"));
        assert!(!chroms_equivalent("chr1", "chr2"));
        assert!(!chroms_equivalent("chrX", "chrY"));
    }
    
    #[test]
    fn test_intersect_intervals() {
        // Overlapping intervals
        assert_eq!(intersect_intervals(0, 100, 50, 150), Some((50, 100)));
        assert_eq!(intersect_intervals(50, 150, 0, 100), Some((50, 100)));
        
        // One contains the other
        assert_eq!(intersect_intervals(0, 100, 25, 75), Some((25, 75)));
        assert_eq!(intersect_intervals(25, 75, 0, 100), Some((25, 75)));
        
        // Exact match
        assert_eq!(intersect_intervals(0, 100, 0, 100), Some((0, 100)));
        
        // No overlap
        assert_eq!(intersect_intervals(0, 50, 50, 100), None);
        assert_eq!(intersect_intervals(0, 50, 100, 150), None);
        
        // Adjacent (no overlap)
        assert_eq!(intersect_intervals(0, 50, 50, 100), None);
    }
    
    #[test]
    fn test_map_basic_positive_strand() {
        let index = create_test_index();
        let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
        
        // Query within first block (100-200)
        let results = mapper.map("chr1", 120, 180, Strand::Plus);
        assert!(results.is_some());
        let results = results.unwrap();
        assert_eq!(results.len(), 1);
        
        let segment = &results[0];
        assert_eq!(segment.source.start, 120);
        assert_eq!(segment.source.end, 180);
        assert_eq!(segment.target.start, 120); // Same as source for identity mapping
        assert_eq!(segment.target.end, 180);
        assert_eq!(segment.target.strand, Strand::Plus);
    }
    
    #[test]
    fn test_map_nonexistent_chrom() {
        let index = create_test_index();
        let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
        
        let results = mapper.map("chrNONE", 0, 100, Strand::Plus);
        assert!(results.is_none());
    }
    
    #[test]
    fn test_map_no_overlap() {
        let index = create_test_index();
        let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
        
        // Query before any blocks
        let results = mapper.map("chr1", 0, 50, Strand::Plus);
        assert!(results.is_some());
        assert!(results.unwrap().is_empty());
    }
    
    #[test]
    fn test_map_partial_overlap() {
        let index = create_test_index();
        let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
        
        // Query overlapping start of first block (100-200)
        let results = mapper.map("chr1", 50, 150, Strand::Plus);
        assert!(results.is_some());
        let results = results.unwrap();
        assert_eq!(results.len(), 1);
        
        let segment = &results[0];
        // Intersection should be [100, 150)
        assert_eq!(segment.source.start, 100);
        assert_eq!(segment.source.end, 150);
        assert_eq!(segment.target.start, 100);
        assert_eq!(segment.target.end, 150);
    }
    
    #[test]
    fn test_map_single() {
        let index = create_test_index();
        let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
        
        // Map single position
        let result = mapper.map_single("chr1", 150, Strand::Plus);
        assert!(result.is_some());
        let segment = result.unwrap();
        assert_eq!(segment.source.start, 150);
        assert_eq!(segment.source.end, 151);
        assert_eq!(segment.target.start, 150);
        assert_eq!(segment.target.end, 151);
    }
    
    #[test]
    fn test_map_strand_combination() {
        let index = create_test_index();
        let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
        
        // Query with minus strand on positive target
        let results = mapper.map("chr1", 120, 180, Strand::Minus);
        assert!(results.is_some());
        let results = results.unwrap();
        assert_eq!(results.len(), 1);
        
        // Plus target + Minus query = Minus
        assert_eq!(results[0].target.strand, Strand::Minus);
    }
    
    #[test]
    fn test_map_chrom_style() {
        let index = create_test_index();
        
        // Test Short style
        let mapper = CoordinateMapper::new(index, ChromStyle::Short);
        let results = mapper.map("chr1", 120, 180, Strand::Plus);
        assert!(results.is_some());
        let results = results.unwrap();
        assert_eq!(results[0].target.chrom, "1");
        assert_eq!(results[0].source.chrom, "1");
    }
}
