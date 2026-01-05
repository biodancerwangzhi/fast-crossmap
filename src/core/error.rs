//! Error types for FastCrossMap
//!
//! Defines all error types used throughout the library.

use std::path::PathBuf;
use thiserror::Error;

/// Main error type for FastCrossMap operations
#[derive(Debug, Error)]
pub enum FastCrossMapError {
    /// Chain file parsing errors
    #[error("Chain parse error: {0}")]
    ChainParse(#[from] ChainParseError),

    /// Coordinate mapping errors
    #[error("Mapping error: {0}")]
    Mapping(#[from] MappingError),

    /// Format conversion errors
    #[error("Conversion error: {0}")]
    Conversion(#[from] ConversionError),

    /// I/O errors
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// Errors that can occur during chain file parsing
#[derive(Debug, Error)]
pub enum ChainParseError {
    /// Invalid chain header format
    #[error("Invalid chain header at line {line}: {message}")]
    InvalidHeader { line: usize, message: String },

    /// Invalid data line format
    #[error("Invalid data line at line {line}: {message}")]
    InvalidDataLine { line: usize, message: String },

    /// Source strand must be '+'
    #[error("Source strand must be '+', got '{strand}' at line {line}")]
    InvalidSourceStrand { line: usize, strand: String },

    /// Target strand must be '+' or '-'
    #[error("Target strand must be '+' or '-', got '{strand}' at line {line}")]
    InvalidTargetStrand { line: usize, strand: String },

    /// Failed to parse integer
    #[error("Failed to parse integer '{value}' at line {line}: {message}")]
    ParseInt {
        line: usize,
        value: String,
        message: String,
    },

    /// File not found
    #[error("Chain file not found: {0}")]
    FileNotFound(PathBuf),

    /// Unsupported compression format
    #[error("Unsupported compression format: {0}")]
    UnsupportedCompression(String),

    /// I/O error during parsing
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// Errors that can occur during coordinate mapping
#[derive(Debug, Error)]
pub enum MappingError {
    /// Chromosome not found in index
    #[error("Chromosome not found in index: {0}")]
    ChromosomeNotFound(String),

    /// Invalid coordinate range
    #[error("Invalid coordinate range: start ({start}) > end ({end})")]
    InvalidRange { start: u64, end: u64 },

    /// No mapping found
    #[error("No mapping found for {chrom}:{start}-{end}")]
    NoMapping {
        chrom: String,
        start: u64,
        end: u64,
    },
}

/// Errors that can occur during format conversion
#[derive(Debug, Error)]
pub enum ConversionError {
    /// Invalid BED format
    #[error("Invalid BED format at line {line}: {message}")]
    InvalidBed { line: usize, message: String },

    /// Invalid VCF format
    #[error("Invalid VCF format at line {line}: {message}")]
    InvalidVcf { line: usize, message: String },

    /// Invalid GFF format
    #[error("Invalid GFF format at line {line}: {message}")]
    InvalidGff { line: usize, message: String },

    /// Invalid BAM/SAM format
    #[error("Invalid BAM/SAM format: {0}")]
    InvalidBam(String),

    /// Reference genome error
    #[error("Reference genome error: {0}")]
    ReferenceGenome(String),

    /// Output write error
    #[error("Failed to write output: {0}")]
    WriteError(String),

    /// I/O error
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// Result type alias for FastCrossMap operations
pub type Result<T> = std::result::Result<T, FastCrossMapError>;

/// Result type alias for chain parsing operations
pub type ChainResult<T> = std::result::Result<T, ChainParseError>;

/// Result type alias for mapping operations
pub type MappingResult<T> = std::result::Result<T, MappingError>;

/// Result type alias for conversion operations
pub type ConversionResult<T> = std::result::Result<T, ConversionError>;
