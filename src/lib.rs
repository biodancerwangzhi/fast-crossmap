//! # FastCrossMap
//!
//! High-performance genome coordinate liftover tool, written in pure Rust.
//!
//! FastCrossMap converts genomic coordinates between reference assemblies
//! (e.g. hg19 → hg38) using UCSC chain files. It supports 8 formats:
//! BED, VCF, GFF, GVCF, MAF, WIG, BigWig, and BAM/SAM.
//!
//! ## Quick start
//!
//! ```no_run
//! use fast_crossmap::{ChainIndex, CoordinateMapper, ChromStyle, Strand};
//!
//! // Load chain file (supports .chain, .chain.gz, .chain.bz2)
//! let index = ChainIndex::from_chain_file("hg19ToHg38.over.chain.gz").unwrap();
//!
//! // Create mapper
//! let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
//!
//! // Map a single position
//! if let Some(seg) = mapper.map_single("chr1", 1000, Strand::Plus) {
//!     println!("{}:{}", seg.target.chrom, seg.target.start);
//! }
//!
//! // Map a region
//! if let Some(segments) = mapper.map("chr1", 1000, 2000, Strand::Plus) {
//!     for seg in &segments {
//!         println!("{}:{}-{}", seg.target.chrom, seg.target.start, seg.target.end);
//!     }
//! }
//! ```
//!
//! ## File-level conversion
//!
//! Each format module provides a `convert_*` function that reads an input file,
//! lifts over coordinates, and writes the result:
//!
//! ```no_run
//! use fast_crossmap::{ChainIndex, CoordinateMapper, ChromStyle, bed};
//!
//! let index = ChainIndex::from_chain_file("hg19ToHg38.over.chain.gz").unwrap();
//! let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
//! let stats = bed::convert_bed("input.bed", "output.bed", "unmapped.bed", &mapper, 4).unwrap();
//! println!("Converted {} records, {} failed", stats.success, stats.failed);
//! ```

pub mod core;
pub mod formats;

// Re-export core types
pub use core::{
    ChainBlock, ChainFile, ChainFileError, ChainHeader, ChainIndex, ChainParseError,
    ChromStyle, CompatMode, ConversionError, CoordinateMapper, FastCrossMapError,
    MapResult, MappingError, Strand, parse_chain_file, parse_chain_bytes,
};

// Re-export all format modules
pub use formats::{bed, vcf, gff, gvcf, maf, wig, region};
#[cfg(feature = "bam")]
pub use formats::bam;

// Re-export bigwig converter (submodule of wig)
pub use formats::wig::bigwig::convert_bigwig;
