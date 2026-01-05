//! FastCrossMap - High-performance genome coordinate liftover
//!
//! A Rust reimplementation of CrossMap focusing on BED and VCF conversion.
//!
//! # Features
//!
//! - 10x+ faster than Python CrossMap
//! - Zero-copy parsing for minimal memory allocation
//! - Parallel processing with rayon
//! - Support for compressed chain files (gzip, bzip2)
//!
//! # Example
//!
//! ```ignore
//! use fast_crossmap::{ChainIndex, CoordinateMapper, ChromStyle, Strand};
//!
//! // Load chain file
//! let index = ChainIndex::from_chain_file("hg19ToHg38.chain")?;
//!
//! // Create mapper
//! let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
//!
//! // Map coordinates
//! let result = mapper.map("chr1", 1000, 2000, Strand::Plus);
//! ```

pub mod core;
pub mod formats;

// Re-export commonly used types
pub use core::{
    ChainBlock, ChainFile, ChainFileError, ChainHeader, ChainIndex, ChainParseError, 
    ChromStyle, ConversionError, CoordinateMapper, FastCrossMapError, MapResult, 
    MappingError, Strand, parse_chain_file, parse_chain_bytes,
};
pub use formats::{bed, vcf};
