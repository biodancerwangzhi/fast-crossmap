//! File format adapters
//!
//! Adapters for different genomic file formats (BED, VCF, GVCF, GFF/GTF, MAF, Wiggle/BigWig, BAM/SAM/CRAM, Region).

#[cfg(feature = "bam")]
pub mod bam;
pub mod bed;
pub mod gff;
pub mod gvcf;
pub mod maf;
pub mod region;
pub mod vcf;
pub mod wig;

#[cfg(feature = "bam")]
pub use bam::{BamError, AlignmentTag, CigarOp, CigarReconstructor, ConversionStats as BamConversionStats, convert_bam};
pub use bed::{BedRecordView, BedParseError, convert_bed, ConversionStats as BedConversionStats};
pub use gff::{GffRecordView, GffParseError, convert_gff, ConversionStats as GffConversionStats};
pub use gvcf::{GvcfRecordView, GvcfParseError, convert_gvcf, ConversionStats as GvcfConversionStats};
pub use maf::{MafRecordView, MafParseError, MafColumnIndices, convert_maf, ConversionStats as MafConversionStats};
pub use region::{RegionError, RegionResult, FailureReason, map_region, convert_region, parse_bed_line, ConversionStats as RegionConversionStats};
pub use vcf::{VcfRecordView, VcfParseError, convert_vcf, ConversionStats as VcfConversionStats};
pub use wig::{WigReader, WigDeclaration, WigFormat, WigDataPoint, BedGraphRecord, WigParseError, convert_wig, ConversionStats as WigConversionStats};
