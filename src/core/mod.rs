//! Core coordinate mapping functionality
//!
//! This module contains the chain file parser, interval index,
//! and coordinate mapping algorithms.

mod chain;
pub mod dna;
mod error;
mod index;
pub mod io;
mod mapper;

pub use chain::{
    parse_chain_file, parse_chain_bytes, parse_chain_reader, 
    ChainBlock, ChainFile, ChainHeader, CompressionFormat,
    ChainParseError as ChainFileError, ChainParseErrorKind,
    detect_compression,
};
pub use error::{
    ChainParseError, ChainResult, ConversionError, ConversionResult,
    FastCrossMapError, MappingError, MappingResult, Result,
};
pub use index::{ChainIndex, ChainInterval, IntervalValue};
pub use io::{
    ByteLineIterator, IoStrategy, LineIterator, SmartReader,
    DEFAULT_BUFFER_SIZE, LARGE_BUFFER_SIZE, MMAP_THRESHOLD,
};
pub use mapper::{ChromStyle, CompatMode, CoordinateMapper, MapResult, MappingSegment, Strand, normalize_chrom, update_chrom_id, chroms_equivalent, intersect_intervals};
