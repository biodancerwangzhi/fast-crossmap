//! BAM/SAM/CRAM format adapter
//!
//! Handles BAM/SAM/CRAM format conversion for alignment data.
//! Uses rust-htslib for reading and writing BAM files.
//!
//! **Validates: Requirements 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7**

use crate::core::{CoordinateMapper, Strand};
use rust_htslib::bam::{self, Read, Record, Header, HeaderView};
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Cigar, CigarString};
use std::collections::HashMap;
use std::path::Path;

/// BAM conversion error
#[derive(Debug)]
pub enum BamError {
    HtslibError(rust_htslib::errors::Error),
    IoError(std::io::Error),
    InvalidCigar(String),
    MappingFailed(String),
}

impl std::fmt::Display for BamError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BamError::HtslibError(e) => write!(f, "HTSlib error: {}", e),
            BamError::IoError(e) => write!(f, "IO error: {}", e),
            BamError::InvalidCigar(msg) => write!(f, "Invalid CIGAR: {}", msg),
            BamError::MappingFailed(msg) => write!(f, "Mapping failed: {}", msg),
        }
    }
}

impl std::error::Error for BamError {}

impl From<rust_htslib::errors::Error> for BamError {
    fn from(e: rust_htslib::errors::Error) -> Self {
        BamError::HtslibError(e)
    }
}

impl From<std::io::Error> for BamError {
    fn from(e: std::io::Error) -> Self {
        BamError::IoError(e)
    }
}

/// Alignment tags for tracking mapping status
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentTag {
    QF, NN, NU, NM, UN, UU, UM, MN, MU, MM, SN, SM, SU,
}

impl AlignmentTag {
    pub fn as_str(&self) -> &'static str {
        match self {
            AlignmentTag::QF => "QF", AlignmentTag::NN => "NN",
            AlignmentTag::NU => "NU", AlignmentTag::NM => "NM",
            AlignmentTag::UN => "UN", AlignmentTag::UU => "UU",
            AlignmentTag::UM => "UM", AlignmentTag::MN => "MN",
            AlignmentTag::MU => "MU", AlignmentTag::MM => "MM",
            AlignmentTag::SN => "SN", AlignmentTag::SM => "SM",
            AlignmentTag::SU => "SU",
        }
    }
}

/// Conversion statistics
#[derive(Debug, Clone, Default)]
pub struct ConversionStats {
    pub total: usize,
    pub mapped: usize,
    pub unmapped: usize,
    pub failed: usize,
    pub paired: usize,
    pub single: usize,
}

/// CIGAR operation types
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    Match(u32), Insertion(u32), Deletion(u32), Skip(u32),
    SoftClip(u32), HardClip(u32), Padding(u32), Equal(u32), Diff(u32),
}

impl CigarOp {
    pub fn len(&self) -> u32 {
        match self {
            CigarOp::Match(n) | CigarOp::Insertion(n) | CigarOp::Deletion(n) |
            CigarOp::Skip(n) | CigarOp::SoftClip(n) | CigarOp::HardClip(n) |
            CigarOp::Padding(n) | CigarOp::Equal(n) | CigarOp::Diff(n) => *n,
        }
    }
    
    pub fn consumes_reference(&self) -> bool {
        matches!(self, CigarOp::Match(_) | CigarOp::Deletion(_) | 
                 CigarOp::Skip(_) | CigarOp::Equal(_) | CigarOp::Diff(_))
    }
    
    pub fn consumes_query(&self) -> bool {
        matches!(self, CigarOp::Match(_) | CigarOp::Insertion(_) | 
                 CigarOp::SoftClip(_) | CigarOp::Equal(_) | CigarOp::Diff(_))
    }
    
    pub fn from_htslib(cigar: &Cigar) -> Self {
        match cigar {
            Cigar::Match(n) => CigarOp::Match(*n),
            Cigar::Ins(n) => CigarOp::Insertion(*n),
            Cigar::Del(n) => CigarOp::Deletion(*n),
            Cigar::RefSkip(n) => CigarOp::Skip(*n),
            Cigar::SoftClip(n) => CigarOp::SoftClip(*n),
            Cigar::HardClip(n) => CigarOp::HardClip(*n),
            Cigar::Pad(n) => CigarOp::Padding(*n),
            Cigar::Equal(n) => CigarOp::Equal(*n),
            Cigar::Diff(n) => CigarOp::Diff(*n),
        }
    }
    
    pub fn to_htslib(&self) -> Cigar {
        match self {
            CigarOp::Match(n) => Cigar::Match(*n),
            CigarOp::Insertion(n) => Cigar::Ins(*n),
            CigarOp::Deletion(n) => Cigar::Del(*n),
            CigarOp::Skip(n) => Cigar::RefSkip(*n),
            CigarOp::SoftClip(n) => Cigar::SoftClip(*n),
            CigarOp::HardClip(n) => Cigar::HardClip(*n),
            CigarOp::Padding(n) => Cigar::Pad(*n),
            CigarOp::Equal(n) => Cigar::Equal(*n),
            CigarOp::Diff(n) => Cigar::Diff(*n),
        }
    }
}

/// CIGAR string reconstructor
pub struct CigarReconstructor;

impl CigarReconstructor {
    pub fn reverse(cigar: &[CigarOp]) -> Vec<CigarOp> {
        cigar.iter().rev().cloned().collect()
    }
    
    pub fn merge_adjacent(cigar: &[CigarOp]) -> Vec<CigarOp> {
        if cigar.is_empty() { return vec![]; }
        let mut result = Vec::with_capacity(cigar.len());
        let mut current = cigar[0];
        for op in cigar.iter().skip(1) {
            let can_merge = match (&current, op) {
                (CigarOp::Match(a), CigarOp::Match(b)) => Some(CigarOp::Match(a + b)),
                (CigarOp::Insertion(a), CigarOp::Insertion(b)) => Some(CigarOp::Insertion(a + b)),
                (CigarOp::Deletion(a), CigarOp::Deletion(b)) => Some(CigarOp::Deletion(a + b)),
                (CigarOp::Skip(a), CigarOp::Skip(b)) => Some(CigarOp::Skip(a + b)),
                (CigarOp::SoftClip(a), CigarOp::SoftClip(b)) => Some(CigarOp::SoftClip(a + b)),
                (CigarOp::HardClip(a), CigarOp::HardClip(b)) => Some(CigarOp::HardClip(a + b)),
                (CigarOp::Equal(a), CigarOp::Equal(b)) => Some(CigarOp::Equal(a + b)),
                (CigarOp::Diff(a), CigarOp::Diff(b)) => Some(CigarOp::Diff(a + b)),
                _ => None,
            };
            if let Some(merged) = can_merge { current = merged; }
            else { result.push(current); current = *op; }
        }
        result.push(current);
        result
    }
    
    pub fn query_length(cigar: &[CigarOp]) -> u32 {
        cigar.iter().filter(|op| op.consumes_query()).map(|op| op.len()).sum()
    }
    
    pub fn reference_length(cigar: &[CigarOp]) -> u32 {
        cigar.iter().filter(|op| op.consumes_reference()).map(|op| op.len()).sum()
    }
    
    pub fn validate(cigar: &[CigarOp], seq_len: u32) -> bool {
        Self::query_length(cigar) == seq_len
    }
    
    pub fn handle_break(cigar: &[CigarOp], break_pos: u32) -> (Vec<CigarOp>, Vec<CigarOp>) {
        let mut left = Vec::new();
        let mut right = Vec::new();
        let mut ref_pos = 0u32;
        let mut query_pos = 0u32;
        let mut in_left = true;
        
        for op in cigar {
            if in_left {
                let ref_consumed = if op.consumes_reference() { op.len() } else { 0 };
                if ref_pos + ref_consumed > break_pos {
                    let left_len = break_pos - ref_pos;
                    let right_len = ref_consumed - left_len;
                    if left_len > 0 {
                        left.push(match op {
                            CigarOp::Match(_) => CigarOp::Match(left_len),
                            CigarOp::Deletion(_) => CigarOp::Deletion(left_len),
                            CigarOp::Skip(_) => CigarOp::Skip(left_len),
                            CigarOp::Equal(_) => CigarOp::Equal(left_len),
                            CigarOp::Diff(_) => CigarOp::Diff(left_len),
                            _ => *op,
                        });
                    }
                    let query_consumed = if op.consumes_query() { op.len() } else { 0 };
                    if query_consumed > left_len {
                        left.push(CigarOp::SoftClip(query_consumed - left_len));
                    }
                    if right_len > 0 {
                        right.push(CigarOp::SoftClip(query_pos + left_len));
                        right.push(match op {
                            CigarOp::Match(_) => CigarOp::Match(right_len),
                            CigarOp::Deletion(_) => CigarOp::Deletion(right_len),
                            CigarOp::Skip(_) => CigarOp::Skip(right_len),
                            CigarOp::Equal(_) => CigarOp::Equal(right_len),
                            CigarOp::Diff(_) => CigarOp::Diff(right_len),
                            _ => *op,
                        });
                    }
                    in_left = false;
                } else { left.push(*op); }
                ref_pos += ref_consumed;
                if op.consumes_query() { query_pos += op.len(); }
            } else { right.push(*op); }
        }
        (Self::merge_adjacent(&left), Self::merge_adjacent(&right))
    }
}


fn ops_to_cigar_string(ops: &[CigarOp]) -> CigarString {
    CigarString(ops.iter().map(|op| op.to_htslib()).collect())
}

fn parse_cigar(record: &Record) -> Vec<CigarOp> {
    record.cigar().iter().map(|c| CigarOp::from_htslib(&c)).collect()
}

/// Build new BAM header with target chromosome sizes
fn build_target_header(
    original_header: &HeaderView,
    target_sizes: &HashMap<String, u64>,
) -> Header {
    let mut new_header = Header::new();
    
    // Add HD line
    let mut hd_record = HeaderRecord::new(b"HD");
    hd_record.push_tag(b"VN", "1.6");
    hd_record.push_tag(b"SO", "unsorted");
    new_header.push_record(&hd_record);
    
    // Add SQ lines for target chromosomes
    for (chrom, size) in target_sizes {
        let mut sq_record = HeaderRecord::new(b"SQ");
        sq_record.push_tag(b"SN", chrom);
        sq_record.push_tag(b"LN", &size.to_string());
        new_header.push_record(&sq_record);
    }
    
    // Copy PG lines from original header text
    let header_text = String::from_utf8_lossy(original_header.as_bytes());
    for line in header_text.lines() {
        if line.starts_with("@PG") || line.starts_with("@RG") || line.starts_with("@CO") {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() > 1 {
                let tag = &parts[0][1..];
                let mut record = HeaderRecord::new(tag.as_bytes());
                for part in &parts[1..] {
                    if let Some(idx) = part.find(':') {
                        let key = &part[..idx];
                        let value = &part[idx+1..];
                        record.push_tag(key.as_bytes(), value);
                    }
                }
                new_header.push_record(&record);
            }
        }
    }
    
    new_header
}

fn get_chrom_name(header: &HeaderView, tid: i32) -> Option<String> {
    if tid < 0 { return None; }
    let name = header.tid2name(tid as u32);
    Some(String::from_utf8_lossy(name).to_string())
}

fn get_tid(header: &HeaderView, chrom: &str) -> Option<i32> {
    header.tid(chrom.as_bytes()).map(|t| t as i32)
}

fn revcomp_seq(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' | b'a' => b'T', b'T' | b't' => b'A',
        b'C' | b'c' => b'G', b'G' | b'g' => b'C',
        _ => b'N',
    }).collect()
}

fn reverse_qual(qual: &[u8]) -> Vec<u8> {
    qual.iter().rev().cloned().collect()
}

fn convert_record(
    record: &Record,
    input_header: &HeaderView,
    output_header: &HeaderView,
    mapper: &CoordinateMapper,
) -> Option<(Record, AlignmentTag)> {
    if record.is_unmapped() { return None; }
    
    let tid = record.tid();
    let chrom = get_chrom_name(input_header, tid)?;
    let start = record.pos() as u64;
    let cigar_ops = parse_cigar(record);
    let ref_len = CigarReconstructor::reference_length(&cigar_ops) as u64;
    let end = start + ref_len;
    
    let query_strand = if record.is_reverse() { Strand::Minus } else { Strand::Plus };
    let segments = mapper.map(&chrom, start, end, query_strand)?;
    if segments.len() != 1 { return None; }
    
    let seg = &segments[0];
    let target_chrom = &seg.target.chrom;
    let target_start = seg.target.start;
    let target_strand = seg.target.strand;
    let target_tid = get_tid(output_header, target_chrom)?;
    
    let need_revcomp = target_strand == Strand::Minus && query_strand == Strand::Plus
        || target_strand == Strand::Plus && query_strand == Strand::Minus;
    
    let (new_seq, new_qual, new_cigar) = if need_revcomp {
        let seq = record.seq().as_bytes();
        let qual: Vec<u8> = record.qual().to_vec();
        (revcomp_seq(&seq), reverse_qual(&qual), CigarReconstructor::reverse(&cigar_ops))
    } else {
        (record.seq().as_bytes(), record.qual().to_vec(), cigar_ops)
    };
    
    let cigar_string = ops_to_cigar_string(&new_cigar);
    let mut new_record = Record::new();
    new_record.set(record.qname(), Some(&cigar_string), &new_seq, &new_qual);
    new_record.set_tid(target_tid);
    new_record.set_pos(target_start as i64);
    
    let mut flags = record.flags();
    if need_revcomp { flags ^= 0x10; }
    new_record.set_flags(flags);
    new_record.set_mapq(record.mapq());
    
    let tag = if record.is_paired() {
        if record.is_mate_unmapped() { AlignmentTag::MU } else { AlignmentTag::MM }
    } else { AlignmentTag::SM };
    
    // Copy auxiliary tags
    for aux in record.aux_iter() {
        if let Ok(aux) = aux {
            let tag_name = aux.0;
            if tag_name == b"QF" || tag_name == b"OC" || tag_name == b"OP" { continue; }
            if let Ok(value) = record.aux(tag_name) {
                let _ = new_record.push_aux(tag_name, value);
            }
        }
    }
    
    Some((new_record, tag))
}

/// Convert a BAM/SAM/CRAM file
pub fn convert_bam<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    threads: usize,
) -> Result<ConversionStats, BamError> {
    let mut reader = bam::Reader::from_path(input.as_ref())?;
    reader.set_threads(threads)?;
    let input_header = reader.header().clone();
    
    let target_sizes = mapper.target_sizes();
    let output_header = build_target_header(&input_header, target_sizes);
    
    let mut writer = bam::Writer::from_path(output.as_ref(), &output_header, bam::Format::Bam)?;
    writer.set_threads(threads)?;
    let output_header_view = writer.header().clone();
    
    let mut stats = ConversionStats::default();
    let mut record = Record::new();
    
    while reader.read(&mut record).is_some() {
        stats.total += 1;
        if record.is_paired() { stats.paired += 1; } else { stats.single += 1; }
        if record.is_unmapped() { stats.unmapped += 1; continue; }
        
        match convert_record(&record, &input_header, &output_header_view, mapper) {
            Some((new_record, _tag)) => { writer.write(&new_record)?; stats.mapped += 1; }
            None => { stats.failed += 1; }
        }
    }
    
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_op_from_htslib() {
        assert_eq!(CigarOp::from_htslib(&Cigar::Match(10)), CigarOp::Match(10));
        assert_eq!(CigarOp::from_htslib(&Cigar::Ins(5)), CigarOp::Insertion(5));
    }

    #[test]
    fn test_cigar_consumes_reference() {
        assert!(CigarOp::Match(10).consumes_reference());
        assert!(CigarOp::Deletion(5).consumes_reference());
        assert!(!CigarOp::Insertion(5).consumes_reference());
        assert!(!CigarOp::SoftClip(10).consumes_reference());
    }

    #[test]
    fn test_cigar_consumes_query() {
        assert!(CigarOp::Match(10).consumes_query());
        assert!(CigarOp::Insertion(5).consumes_query());
        assert!(!CigarOp::Deletion(5).consumes_query());
    }

    #[test]
    fn test_cigar_reverse() {
        let cigar = vec![CigarOp::Match(10), CigarOp::Insertion(2), CigarOp::Match(5)];
        let reversed = CigarReconstructor::reverse(&cigar);
        assert_eq!(reversed, vec![CigarOp::Match(5), CigarOp::Insertion(2), CigarOp::Match(10)]);
    }

    #[test]
    fn test_cigar_merge_adjacent() {
        let cigar = vec![CigarOp::Match(10), CigarOp::Match(5), CigarOp::Insertion(2)];
        let merged = CigarReconstructor::merge_adjacent(&cigar);
        assert_eq!(merged, vec![CigarOp::Match(15), CigarOp::Insertion(2)]);
    }

    #[test]
    fn test_cigar_query_length() {
        let cigar = vec![CigarOp::Match(10), CigarOp::Insertion(2), CigarOp::Deletion(3), CigarOp::Match(5)];
        assert_eq!(CigarReconstructor::query_length(&cigar), 17); // 10 + 2 + 5
    }

    #[test]
    fn test_cigar_reference_length() {
        let cigar = vec![CigarOp::Match(10), CigarOp::Insertion(2), CigarOp::Deletion(3), CigarOp::Match(5)];
        assert_eq!(CigarReconstructor::reference_length(&cigar), 18); // 10 + 3 + 5
    }

    #[test]
    fn test_revcomp_seq() {
        assert_eq!(revcomp_seq(b"ACGT"), b"ACGT");
        assert_eq!(revcomp_seq(b"AACGT"), b"ACGTT");
    }

    #[test]
    fn test_reverse_qual() {
        assert_eq!(reverse_qual(&[10, 20, 30, 40]), vec![40, 30, 20, 10]);
    }

    #[test]
    fn test_alignment_tag_as_str() {
        assert_eq!(AlignmentTag::QF.as_str(), "QF");
        assert_eq!(AlignmentTag::MM.as_str(), "MM");
    }
}
