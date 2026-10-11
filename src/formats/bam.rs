//! BAM/SAM format adapter
//!
//! Handles BAM/SAM format conversion for alignment data.
//! Uses noodles (pure Rust) for reading and writing BAM files.

use crate::core::{CoordinateMapper, MappingSegment, Strand};
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        record::Flags,
        record::MappingQuality,
        record::cigar::{op::Kind as CigarKind, Op as SamOp},
        record_buf::{
            self,
            data::field::Value as BufValue,
        },
        RecordBuf,
    },
    header::record::value::{
        map::{self, ReferenceSequence},
        Map,
    },
};
use noodles_sam::alignment::record::data::field::Tag as SamTag;
use std::collections::HashMap;
use std::io;
use std::num::NonZero;
use std::path::Path;

/// BAM conversion error
#[derive(Debug)]
pub enum BamError {
    IoError(std::io::Error),
    InvalidCigar(String),
    MappingFailed(String),
}

impl std::fmt::Display for BamError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BamError::IoError(e) => write!(f, "IO error: {}", e),
            BamError::InvalidCigar(msg) => write!(f, "Invalid CIGAR: {}", msg),
            BamError::MappingFailed(msg) => write!(f, "Mapping failed: {}", msg),
        }
    }
}

impl std::error::Error for BamError {}

impl From<std::io::Error> for BamError {
    fn from(e: std::io::Error) -> Self {
        BamError::IoError(e)
    }
}

/// PLACEHOLDER_CONTINUE

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

    /// The paired-end categories, in the order CrossMap prints them.
    pub const PAIRED: [AlignmentTag; 9] = [
        AlignmentTag::UU, AlignmentTag::UN, AlignmentTag::UM,
        AlignmentTag::MM, AlignmentTag::MU, AlignmentTag::MN,
        AlignmentTag::NN, AlignmentTag::NU, AlignmentTag::NM,
    ];

    /// The single-end categories, in the order CrossMap prints them.
    pub const SINGLE: [AlignmentTag; 3] = [
        AlignmentTag::SU, AlignmentTag::SM, AlignmentTag::SN,
    ];

    /// A dense index into `ConversionStats::counts`.
    pub fn index(&self) -> usize {
        match self {
            AlignmentTag::QF => 0, AlignmentTag::NN => 1,
            AlignmentTag::NU => 2, AlignmentTag::NM => 3,
            AlignmentTag::UN => 4, AlignmentTag::UU => 5,
            AlignmentTag::UM => 6, AlignmentTag::MN => 7,
            AlignmentTag::MU => 8, AlignmentTag::MM => 9,
            AlignmentTag::SN => 10, AlignmentTag::SM => 11,
            AlignmentTag::SU => 12,
        }
    }

    /// A human-readable label, matching CrossMap's summary wording.
    pub fn description(&self) -> &'static str {
        match self {
            AlignmentTag::QF => "QC failed",
            AlignmentTag::NN => "R1 unmap, R2 unmap",
            AlignmentTag::NU => "R1 unmap, R2 unique",
            AlignmentTag::NM => "R1 unmap, R2 multiple",
            AlignmentTag::UN => "R1 unique, R2 unmapp",
            AlignmentTag::UU => "R1 unique, R2 unique",
            AlignmentTag::UM => "R1 unique, R2 multiple",
            AlignmentTag::MN => "R1 multiple, R2 unmapped",
            AlignmentTag::MU => "R1 multiple, R2 unique",
            AlignmentTag::MM => "R1 multiple, R2 multiple",
            AlignmentTag::SN => "Unmapped",
            AlignmentTag::SM => "Multiple mapped",
            AlignmentTag::SU => "Uniquley mapped",
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
    /// How many records fell into each `AlignmentTag` category. CrossMap
    /// reports these itself, so they double as a compatibility check.
    pub counts: [usize; 13],
}

impl ConversionStats {
    pub fn count_of(&self, tag: AlignmentTag) -> usize {
        self.counts[tag.index()]
    }

    pub fn count(&mut self, tag: AlignmentTag) {
        self.counts[tag.index()] += 1;
    }
}

/// PLACEHOLDER_CIGAR

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

    pub fn from_noodles(op: SamOp) -> Self {
        let n = op.len() as u32;
        match op.kind() {
            CigarKind::Match => CigarOp::Match(n),
            CigarKind::Insertion => CigarOp::Insertion(n),
            CigarKind::Deletion => CigarOp::Deletion(n),
            CigarKind::Skip => CigarOp::Skip(n),
            CigarKind::SoftClip => CigarOp::SoftClip(n),
            CigarKind::HardClip => CigarOp::HardClip(n),
            CigarKind::Pad => CigarOp::Padding(n),
            CigarKind::SequenceMatch => CigarOp::Equal(n),
            CigarKind::SequenceMismatch => CigarOp::Diff(n),
        }
    }

    pub fn to_noodles(&self) -> SamOp {
        match self {
            CigarOp::Match(n) => SamOp::new(CigarKind::Match, *n as usize),
            CigarOp::Insertion(n) => SamOp::new(CigarKind::Insertion, *n as usize),
            CigarOp::Deletion(n) => SamOp::new(CigarKind::Deletion, *n as usize),
            CigarOp::Skip(n) => SamOp::new(CigarKind::Skip, *n as usize),
            CigarOp::SoftClip(n) => SamOp::new(CigarKind::SoftClip, *n as usize),
            CigarOp::HardClip(n) => SamOp::new(CigarKind::HardClip, *n as usize),
            CigarOp::Padding(n) => SamOp::new(CigarKind::Pad, *n as usize),
            CigarOp::Equal(n) => SamOp::new(CigarKind::SequenceMatch, *n as usize),
            CigarOp::Diff(n) => SamOp::new(CigarKind::SequenceMismatch, *n as usize),
        }
    }
}

/// PLACEHOLDER_RECONSTRUCTOR

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

/// PLACEHOLDER_HELPERS

fn ops_to_noodles_cigar(ops: &[CigarOp]) -> record_buf::Cigar {
    ops.iter().map(|op| op.to_noodles()).collect()
}

fn parse_cigar_from_bam(record: &bam::Record) -> Vec<CigarOp> {
    record.cigar().iter()
        .filter_map(|r| r.ok())
        .map(CigarOp::from_noodles)
        .collect()
}

fn get_chrom_name(header: &sam::Header, tid: usize) -> Option<String> {
    header.reference_sequences().get_index(tid)
        .map(|(name, _)| name.to_string())
}

fn get_tid(header: &sam::Header, chrom: &str) -> Option<usize> {
    if let Some(idx) = header.reference_sequences().get_index_of(chrom.as_bytes()) {
        return Some(idx);
    }
    let with_chr = format!("chr{}", chrom);
    if let Some(idx) = header.reference_sequences().get_index_of(with_chr.as_bytes()) {
        return Some(idx);
    }
    let without_chr = chrom.strip_prefix("chr").unwrap_or(chrom);
    header.reference_sequences().get_index_of(without_chr.as_bytes())
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

fn seq_bytes_from_record(record: &bam::Record) -> Vec<u8> {
    let seq = record.sequence();
    (0..seq.len()).filter_map(|i| seq.get(i)).collect()
}

fn qual_bytes_from_record(record: &bam::Record) -> Vec<u8> {
    let qs = record.quality_scores();
    let raw = qs.as_bytes();
    if raw.is_empty() || (raw.len() == 1 && raw[0] == 0xff) {
        vec![0xff; record.sequence().len()]
    } else {
        raw.to_vec()
    }
}

/// PLACEHOLDER_HEADER

fn build_target_header(
    original_header: &sam::Header,
    target_sizes: &HashMap<String, u64>,
) -> sam::Header {
    // CrossMap replaces @HD outright with a fixed `VN:1.0 SO:coordinate`
    // (cmmodule/sam_header.py), keeping nothing from the input's own @HD.
    // noodles' `Map::<map::Header>` only models VN, so SO is carried as an
    // "other" field -- which is exactly where its writer looks for it.
    let mut hd = Map::<map::Header>::new(map::header::Version::new(1, 0));
    hd.other_fields_mut()
        .insert(map::header::tag::SORT_ORDER, "coordinate".into());

    let mut builder = sam::Header::builder().set_header(hd);

    let mut sorted_chroms: Vec<_> = target_sizes.iter().collect();
    sorted_chroms.sort_by(|a, b| a.0.cmp(b.0));

    for (chrom, size) in sorted_chroms {
        let len = NonZero::try_from(*size as usize).unwrap_or(NonZero::new(1).unwrap());
        builder = builder.add_reference_sequence(
            chrom.as_str(),
            Map::<ReferenceSequence>::new(len),
        );
    }

    for (id, pg) in original_header.programs().roots() {
        builder = builder.add_program(id.to_vec(), pg.clone());
    }
    for (id, rg) in original_header.read_groups() {
        builder = builder.add_read_group(id.to_vec(), rg.clone());
    }
    for comment in original_header.comments() {
        builder = builder.add_comment(comment.to_vec());
    }

    builder.build()
}

/// PLACEHOLDER_CONVERT

fn copy_aux_data(src: &bam::Record) -> record_buf::Data {
    let qf = SamTag::new(b'Q', b'F');
    let oc = SamTag::new(b'O', b'C');
    let op = SamTag::new(b'O', b'P');
    let mut data = record_buf::Data::default();
    for result in src.data().iter() {
        if let Ok((tag, value)) = result {
            if tag == qf || tag == oc || tag == op { continue; }
            if let Some(buf_val) = value_to_buf(&value) {
                data.insert(tag, buf_val);
            }
        }
    }
    data
}

fn value_to_buf(v: &sam::alignment::record::data::field::Value<'_>) -> Option<BufValue> {
    use sam::alignment::record::data::field::Value;
    match v {
        Value::Character(c) => Some(BufValue::Character(*c)),
        Value::Int8(n) => Some(BufValue::Int8(*n)),
        Value::UInt8(n) => Some(BufValue::UInt8(*n)),
        Value::Int16(n) => Some(BufValue::Int16(*n)),
        Value::UInt16(n) => Some(BufValue::UInt16(*n)),
        Value::Int32(n) => Some(BufValue::Int32(*n)),
        Value::UInt32(n) => Some(BufValue::UInt32(*n)),
        Value::Float(n) => Some(BufValue::Float(*n)),
        Value::String(s) => Some(BufValue::String((*s).into())),
        Value::Hex(s) => Some(BufValue::Hex((*s).into())),
        Value::Array(_) => None,
    }
}

/// PLACEHOLDER_CONVERT_RECORD

/// Insert-size statistics used to re-derive the "proper pair" bit (`0x02`).
///
/// The defaults are CrossMap's (`cmmodule/mapbam.py`): a pair counts as proper
/// when its two ends land on opposite strands *and* the new template length
/// falls within `mean +/- fold * stdev`. With these defaults that window is
/// `[110, 290]`.
#[derive(Debug, Clone, Copy)]
pub struct InsertSizeParams {
    pub mean: f64,
    pub stdev: f64,
    pub fold: f64,
}

impl Default for InsertSizeParams {
    fn default() -> Self {
        Self { mean: 200.0, stdev: 30.0, fold: 3.0 }
    }
}

/// An input alignment with every format-specific detail already decoded, so
/// the BAM and SAM readers can share a single conversion routine.
struct DecodedRecord {
    flags: Flags,
    name: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
    mapq: Option<MappingQuality>,
    /// Source reference name and 0-based half-open interval. Absent for a
    /// record that was never placed, which is what `is_unmapped` means.
    chrom: Option<String>,
    start: u64,
    end: u64,
    /// Source CIGAR, in input order.
    cigar: Vec<CigarOp>,
    data: record_buf::Data,
    /// Mate reference name and 0-based alignment start, when the record
    /// carries them.
    mate_chrom: Option<String>,
    mate_start: Option<u64>,
    mate_reverse: bool,
}

impl DecodedRecord {
    fn query_strand(&self) -> Strand {
        if self.flags.is_reverse_complemented() { Strand::Minus } else { Strand::Plus }
    }

    fn mate_strand(&self) -> Strand {
        if self.mate_reverse { Strand::Minus } else { Strand::Plus }
    }

    /// The record's own liftover, or `None` when it cannot be placed -- either
    /// because it was already unmapped or because the chain does not cover it.
    fn self_mapping(&self, mapper: &CoordinateMapper) -> Option<Vec<MappingSegment>> {
        if self.flags.is_unmapped() { return None; }
        let chrom = self.chrom.as_deref()?;
        map_interval(mapper, chrom, self.start, self.end, self.query_strand())
    }

    /// The mate's liftover.
    ///
    /// CrossMap never looks at the mate record itself: it lifts over a single
    /// base at the mate's own alignment start (`mapbam.py`). The mate's CIGAR
    /// and length are therefore irrelevant.
    fn mate_mapping(&self, mapper: &CoordinateMapper) -> Option<Vec<MappingSegment>> {
        if self.flags.is_mate_unmapped() { return None; }
        let chrom = self.mate_chrom.as_deref()?;
        let start = self.mate_start?;
        map_interval(mapper, chrom, start, start + 1, self.mate_strand())
    }
}

/// Map one source interval, normalising "not mapped" to `None`.
///
/// This mirrors `cmmodule.utils.map_coordinates`, whose callers treat a
/// missing mapping, a single mapping and several mappings as three distinct
/// outcomes.
fn map_interval(
    mapper: &CoordinateMapper,
    chrom: &str,
    start: u64,
    end: u64,
    strand: Strand,
) -> Option<Vec<MappingSegment>> {
    match mapper.map(chrom, start, end, strand) {
        Some(v) if !v.is_empty() => Some(v),
        _ => None,
    }
}

/// Resolve a mapped interval to the output header's reference id, keeping the
/// segments themselves for the strand and position.
fn place(
    output_header: &sam::Header,
    segments: Option<Vec<MappingSegment>>,
) -> Option<(usize, Vec<MappingSegment>)> {
    let segments = segments?;
    let tid = get_tid(output_header, &segments[0].target.chrom)?;
    Some((tid, segments))
}

/// Everything the conversion decides about one output record, before it is
/// turned into a `RecordBuf`. Positions are 1-based.
struct OutSpec {
    flags: Flags,
    /// `None` writes `*` in the RNAME column.
    tid: Option<usize>,
    /// `None` writes `0` in the POS column. Unmapped output ignores this and
    /// writes 1 -- see `build_output`.
    pos: Option<usize>,
    /// `None` writes `255` ("unavailable") in the MAPQ column.
    mapq: Option<MappingQuality>,
    /// The CIGAR to write, in source order; empty writes `*`. It is reversed
    /// for the output when `reverse_complement` is set.
    cigar: Vec<CigarOp>,
    reverse_complement: bool,
    /// `None` writes `*` in the RNEXT column.
    mate_tid: Option<usize>,
    /// The PNEXT column. CrossMap leaves its zero-based field at zero and
    /// pysam renders that one-based, so its output always reads `1`; writing
    /// position 1 outright reproduces that on both the SAM and BAM paths.
    mate_pos: usize,
    template_length: i32,
}

/// Reproduce how CrossMap re-serialises the auxiliary fields.
///
/// CrossMap passes the tags through pysam (`set_tags`, then an explicit re-set
/// of the read group tag). That round trip does two things the SAM
/// specification does not ask for:
///
/// 1. every character (`A`) field comes back as a string (`Z`);
/// 2. `RG` is written last, so it always ends up at the end of the field list.
///
/// Everything else survives untouched: the order of the other fields, integer
/// widths (an `i` stays an `i` even at 2e9), floats, arrays, and one-character
/// string fields, which are *not* promoted to `A`.
fn crossmap_tags(data: &record_buf::Data) -> record_buf::Data {
    let mut out: Vec<(SamTag, BufValue)> = Vec::with_capacity(data.len());
    let mut read_group = None;

    for (tag, value) in data.iter() {
        let value = match value {
            BufValue::Character(c) => BufValue::String(vec![*c].into()),
            v => v.clone(),
        };

        if tag == SamTag::READ_GROUP {
            read_group = Some((tag, value));
        } else {
            out.push((tag, value));
        }
    }

    // `Data` is order preserving, so appending after the loop moves `RG` last.
    if let Some(entry) = read_group {
        out.push(entry);
    }

    out.into_iter().collect()
}

fn build_output(spec: OutSpec, rec: &DecodedRecord) -> RecordBuf {
    let (seq, qual) = if spec.reverse_complement {
        (revcomp_seq(&rec.seq), reverse_qual(&rec.qual))
    } else {
        (rec.seq.clone(), rec.qual.clone())
    };

    let cigar = if spec.reverse_complement {
        CigarReconstructor::reverse(&spec.cigar)
    } else {
        spec.cigar.clone()
    };

    let mut builder = RecordBuf::builder()
        .set_name(rec.name.as_slice())
        .set_flags(spec.flags)
        .set_cigar(ops_to_noodles_cigar(&cigar))
        .set_mate_alignment_start(Position::MIN)
        .set_template_length(spec.template_length)
        .set_sequence(record_buf::Sequence::from(seq.as_slice()))
        .set_quality_scores(record_buf::QualityScores::from(qual))
        .set_data(crossmap_tags(&rec.data));

    if let Some(tid) = spec.tid {
        builder = builder.set_reference_sequence_id(tid);
    }
    match spec.pos.and_then(Position::new) {
        Some(pos) => builder = builder.set_alignment_start(pos),
        // CrossMap leaves the field at zero for an unmapped read and pysam
        // renders it one-based, so its output reads `POS=1` under `RNAME=*`.
        // Writing position 1 reproduces that on both the SAM and BAM paths
        // (the BAM encoder subtracts one again).
        None if spec.tid.is_none() => {
            builder = builder.set_alignment_start(Position::MIN);
        }
        None => {}
    }
    if let Some(mapq) = spec.mapq {
        builder = builder.set_mapping_quality(mapq);
    }
    if let Some(mate_tid) = spec.mate_tid {
        builder = builder.set_mate_reference_sequence_id(mate_tid);
    }
    if let Some(mate_pos) = Position::new(spec.mate_pos) {
        builder = builder.set_mate_alignment_start(mate_pos);
    }

    builder.build()
}

fn maybe(flag: Flags, on: bool) -> Flags {
    if on { flag } else { Flags::default() }
}

/// Recognise a BAM alignment into the format-independent form.
fn decode_bam_record(record: &bam::Record, input_header: &sam::Header) -> DecodedRecord {
    let flags = record.flags();
    let chrom = record.reference_sequence_id()
        .and_then(|r| r.ok())
        .and_then(|tid| get_chrom_name(input_header, tid));
    let start = record.alignment_start()
        .and_then(|r| r.ok())
        .map(|p| p.get() as u64 - 1);
    let cigar = parse_cigar_from_bam(record);

    let mate_chrom = record.mate_reference_sequence_id()
        .and_then(|r| r.ok())
        .and_then(|tid| get_chrom_name(input_header, tid));
    let mate_start = record.mate_alignment_start()
        .and_then(|r| r.ok())
        .map(|p| p.get() as u64 - 1);

    let (start0, end) = match start {
        Some(s) => (s, s + CigarReconstructor::reference_length(&cigar) as u64),
        None => (0, 0),
    };

    DecodedRecord {
        flags,
        name: record.name().map(|n| n.to_vec()).unwrap_or_else(|| b"*".to_vec()),
        seq: seq_bytes_from_record(record),
        qual: qual_bytes_from_record(record),
        mapq: record.mapping_quality(),
        chrom,
        start: start0,
        end,
        cigar,
        data: copy_aux_data(record),
        mate_chrom,
        mate_start,
        mate_reverse: flags.is_mate_reverse_complemented(),
    }
}

/// Recognise a SAM alignment buffer into the format-independent form.
fn decode_sam_record(record_buf: &RecordBuf, input_header: &sam::Header) -> DecodedRecord {
    let flags = record_buf.flags();
    let chrom = record_buf.reference_sequence_id()
        .and_then(|tid| get_chrom_name(input_header, tid));
    let start = record_buf.alignment_start().map(|p| p.get() as u64 - 1);
    let cigar: Vec<CigarOp> = record_buf.cigar().as_ref().iter()
        .map(|op| CigarOp::from_noodles(*op))
        .collect();

    let mate_chrom = record_buf.mate_reference_sequence_id()
        .and_then(|tid| get_chrom_name(input_header, tid));
    let mate_start = record_buf.mate_alignment_start().map(|p| p.get() as u64 - 1);

    let (start0, end) = match start {
        Some(s) => (s, s + CigarReconstructor::reference_length(&cigar) as u64),
        None => (0, 0),
    };

    DecodedRecord {
        flags,
        name: record_buf.name().map(|n| n.to_vec()).unwrap_or_else(|| b"*".to_vec()),
        seq: record_buf.sequence().as_ref().to_vec(),
        qual: record_buf.quality_scores().as_ref().to_vec(),
        mapq: record_buf.mapping_quality(),
        chrom,
        start: start0,
        end,
        cigar,
        data: record_buf.data().clone(),
        mate_chrom,
        mate_start,
        mate_reverse: flags.is_mate_reverse_complemented(),
    }
}

/// Convert one decoded alignment, reproducing CrossMap's field rules.
///
/// Returns the record to emit, how to count it, and which CrossMap category it
/// fell into. A record is never dropped: a read that cannot be placed is
/// written out with no reference, as CrossMap does.
fn convert_decoded(
    rec: &DecodedRecord,
    output_header: &sam::Header,
    mapper: &CoordinateMapper,
    is: &InsertSizeParams,
) -> (BamClass, RecordBuf, AlignmentTag) {
    let (spec, tag) = if rec.flags.is_segmented() {
        convert_paired(rec, output_header, mapper, is)
    } else {
        convert_single(rec, output_header, mapper)
    };

    let class = if spec.tid.is_some() {
        BamClass::Mapped
    } else if rec.flags.is_unmapped() {
        BamClass::Unmapped
    } else {
        BamClass::Failed
    };

    (class, build_output(spec, rec), tag)
}

/// CrossMap's single-end branch (`mapbam.py`).
fn convert_single(
    rec: &DecodedRecord,
    output_header: &sam::Header,
    mapper: &CoordinateMapper,
) -> (OutSpec, AlignmentTag) {
    let mut spec = OutSpec {
        // The FLAG is rebuilt from zero, so anything the input carried beyond
        // the two bits re-derived below -- QC fail, duplicate, supplementary
        // -- is dropped, exactly as CrossMap drops it.
        flags: Flags::default(),
        tid: None,
        pos: None,
        mapq: None,
        cigar: rec.cigar.clone(),
        reverse_complement: false,
        mate_tid: None,
        mate_pos: 1,
        template_length: 0,
    };

    match place(output_header, rec.self_mapping(mapper)) {
        None => {
            spec.flags |= Flags::UNMAPPED;
            // A read that was already unmapped keeps its CIGAR; one that
            // merely failed to lift over loses it.
            if !rec.flags.is_unmapped() {
                spec.cigar = Vec::new();
            }
            (spec, AlignmentTag::SN)
        }
        Some((tid, segments)) => {
            let target = segments[0].target.clone();
            let multiple = segments.len() > 1;

            if target.strand == Strand::Minus {
                spec.flags |= Flags::REVERSE_COMPLEMENTED;
            }
            if multiple {
                spec.flags |= Flags::SECONDARY;
            }

            spec.tid = Some(tid);
            spec.pos = Some(target.start as usize + 1);
            spec.mapq = rec.mapq;
            spec.reverse_complement = target.strand != rec.query_strand();

            let tag = if multiple { AlignmentTag::SM } else { AlignmentTag::SU };
            (spec, tag)
        }
    }
}

/// Which of CrossMap's nine paired categories a pair falls into, from the shape
/// of each end's mapping: `None` if the end did not lift over, `Some(false)` if
/// it landed in exactly one place, `Some(true)` if it landed in several. This
/// is the `if`/`elif` ladder in `mapbam.py` collapsed into a match.
fn paired_category(read1: Option<bool>, read2: Option<bool>) -> AlignmentTag {
    match (read1, read2) {
        (None, None) => AlignmentTag::NN,
        (None, Some(false)) => AlignmentTag::NU,
        (None, Some(true)) => AlignmentTag::NM,
        (Some(false), None) => AlignmentTag::UN,
        (Some(true), None) => AlignmentTag::MN,
        (Some(false), Some(false)) => AlignmentTag::UU,
        (Some(false), Some(true)) => AlignmentTag::UM,
        (Some(true), Some(false)) => AlignmentTag::MU,
        (Some(true), Some(true)) => AlignmentTag::MM,
    }
}

/// CrossMap's paired-end branch (`mapbam.py`): nine categories, each with its
/// own rule for the FLAG, RNAME/POS, MAPQ, CIGAR, RNEXT/PNEXT and TLEN.
fn convert_paired(
    rec: &DecodedRecord,
    output_header: &sam::Header,
    mapper: &CoordinateMapper,
    is: &InsertSizeParams,
) -> (OutSpec, AlignmentTag) {
    // Only SEGMENTED and the first/last-segment bit survive from the input;
    // pysam's `is_read1`/`is_read2` are edge triggered, so a record with both
    // bits set is treated as read 1.
    let mut base = Flags::SEGMENTED;
    if rec.flags.is_first_segment() && !rec.flags.is_last_segment() {
        base |= Flags::FIRST_SEGMENT;
    } else if rec.flags.is_last_segment() && !rec.flags.is_first_segment() {
        base |= Flags::LAST_SEGMENT;
    }

    // QC-failed pairs are written out before anything is mapped.
    if rec.flags.is_qc_fail() {
        let spec = OutSpec {
            flags: base | Flags::QC_FAIL,
            tid: None,
            pos: None,
            mapq: None,
            cigar: rec.cigar.clone(),
            reverse_complement: false,
            mate_tid: None,
            mate_pos: 1,
            template_length: 0,
        };
        return (spec, AlignmentTag::QF);
    }

    let read1 = place(output_header, rec.self_mapping(mapper));
    let read2 = place(output_header, rec.mate_mapping(mapper));

    let seg1 = read1.as_ref().map(|(_, segments)| &segments[0]);
    let tid1 = read1.as_ref().map(|(tid, _)| *tid);
    let pos1 = seg1.map(|s| s.target.start as usize + 1);
    let flip1 = seg1.map_or(false, |s| s.target.strand != rec.query_strand());
    let minus1 = seg1.map_or(false, |s| s.target.strand == Strand::Minus);
    let multiple1 = read1.as_ref().map_or(false, |(_, segments)| segments.len() > 1);

    let seg2 = read2.as_ref().map(|(_, segments)| &segments[0]);
    let tid2 = read2.as_ref().map(|(tid, _)| *tid);
    let pos2 = seg2.map(|s| s.target.start as usize + 1);
    let minus2 = seg2.map_or(false, |s| s.target.strand == Strand::Minus);
    let multiple2 = read2.as_ref().map_or(false, |(_, segments)| segments.len() > 1);

    let tag = paired_category(tid1.map(|_| multiple1), tid2.map(|_| multiple2));

    let (flags, tid, pos, mapq, reverse_complement, mate_tid, mate_pos, template_length) =
        match tag {
            // Neither end placed: the read keeps its bases and nothing else.
            AlignmentTag::NN => (
                base | Flags::UNMAPPED,
                None, None, None, false, None, 1, 0,
            ),

            // Read 1 unplaced, read 2 placed: CrossMap parks the read on its
            // mate's new position and points RNEXT at the very same place.
            AlignmentTag::NU | AlignmentTag::NM => {
                let flags = base
                    | maybe(Flags::MATE_REVERSE_COMPLEMENTED, minus2)
                    | maybe(Flags::SECONDARY, tag == AlignmentTag::NM);
                // A read that was already unmapped keeps its own MAPQ; one
                // that merely failed to lift over has none.
                let mapq = if rec.flags.is_unmapped() { rec.mapq } else { None };
                (flags, tid2, pos2, mapq, false, tid2, pos2.unwrap_or(1), 0)
            }

            // Read 1 placed, read 2 unplaced.
            AlignmentTag::UN => (
                base | maybe(Flags::REVERSE_COMPLEMENTED, minus1) | Flags::MATE_UNMAPPED,
                tid1, pos1, rec.mapq, flip1, tid1, pos1.unwrap_or(1), 0,
            ),
            AlignmentTag::MN => (
                base | Flags::SECONDARY
                    | maybe(Flags::REVERSE_COMPLEMENTED, minus1)
                    | Flags::MATE_UNMAPPED,
                tid1, pos1, None, flip1, tid1, pos1.unwrap_or(1), 0,
            ),

            // Both ends placed. The only difference between the four is
            // whether read 1 and/or read 2 mapped to several places.
            AlignmentTag::UU | AlignmentTag::UM | AlignmentTag::MU | AlignmentTag::MM => {
                let mut flags = base
                    | maybe(Flags::REVERSE_COMPLEMENTED, minus1)
                    | maybe(Flags::MATE_REVERSE_COMPLEMENTED, minus2);
                if multiple1 {
                    flags |= Flags::SECONDARY;
                }
                if tag == AlignmentTag::UM {
                    flags |= Flags::SECONDARY;
                }

                let (mapq, template_length) = if tag == AlignmentTag::UU {
                    let start = pos1.unwrap_or(1) as i64;
                    let next = pos2.unwrap_or(1) as i64;
                    let tlen = (start - next).abs()
                        + CigarReconstructor::reference_length(&rec.cigar) as i64;
                    let opposite = seg1.map(|s| s.target.strand)
                        != seg2.map(|s| s.target.strand);
                    if opposite
                        && tlen as f64 <= is.mean + is.fold * is.stdev
                        && tlen as f64 >= is.mean - is.fold * is.stdev
                    {
                        flags |= Flags::PROPERLY_SEGMENTED;
                    }
                    (rec.mapq, tlen as i32)
                } else {
                    // A read that mapped to several places has no usable
                    // mapping quality of its own.
                    (if multiple1 { None } else { rec.mapq }, 0)
                };

                (flags, tid1, pos1, mapq, flip1, tid2, pos2.unwrap_or(1), template_length)
            }

            AlignmentTag::QF
            | AlignmentTag::SN
            | AlignmentTag::SM
            | AlignmentTag::SU => unreachable!("only the paired categories reach here"),
        };

    let spec = OutSpec {
        flags,
        tid,
        pos,
        mapq,
        cigar: rec.cigar.clone(),
        reverse_complement,
        mate_tid,
        mate_pos,
        template_length,
    };
    (spec, tag)
}


/// PLACEHOLDER_MAIN

/// Convert a BAM/SAM file
pub fn convert_bam<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    threads: usize,
) -> Result<ConversionStats, BamError> {
    let input_path = input.as_ref();
    let output_path = output.as_ref();

    let is_sam_input = input_path.extension().and_then(|e| e.to_str()) == Some("sam");
    let is_sam_output = output_path.extension().and_then(|e| e.to_str()) == Some("sam");

    let mut stats = ConversionStats::default();
    let target_sizes = mapper.target_sizes();

    // Two input paths (SAM text, BAM BGZF) × two output paths (plain SAM,
    // parallel BGZF). Each combination keeps its reader and writer concrete.
    if is_sam_input {
        let mut reader = std::fs::File::open(input_path)
            .map(std::io::BufReader::new)
            .map(sam::io::Reader::new)?;
        let input_header = reader.read_header()?;
        let output_header = build_target_header(&input_header, target_sizes);

        if is_sam_output {
            let mut writer = std::fs::File::create(output_path)
                .map(std::io::BufWriter::new)
                .map(sam::io::Writer::new)?;
            writer.write_header(&output_header)?;
            process_sam_records(&mut reader, &input_header, &output_header, mapper, &mut writer, threads, &mut stats)?;
        } else {
            let mut writer = new_bam_writer(output_path, threads)?;
            writer.write_header(&output_header)?;
            process_sam_records(&mut reader, &input_header, &output_header, mapper, &mut writer, threads, &mut stats)?;
            finish_bam_writer(writer)?;
        }
    } else {
        // With more than one worker the BGZF input is decompressed across the
        // workers too: `MultithreadedReader` yields decompressed bytes, so
        // `Reader::from` wraps it without adding a second decoder.
        let file = std::fs::File::open(input_path)?;
        let inner: Box<dyn io::Read + Send> = if threads > 1 {
            Box::new(noodles_bgzf::io::MultithreadedReader::with_worker_count(bgzf_workers(threads), file))
        } else {
            Box::new(noodles_bgzf::io::Reader::new(file))
        };
        let mut reader = bam::io::Reader::from(inner);
        let input_header = reader.read_header()?;
        let output_header = build_target_header(&input_header, target_sizes);

        if is_sam_output {
            let mut writer = std::fs::File::create(output_path)
                .map(std::io::BufWriter::new)
                .map(sam::io::Writer::new)?;
            writer.write_header(&output_header)?;
            process_bam_records(&mut reader, &input_header, &output_header, mapper, &mut writer, threads, &mut stats)?;
        } else {
            let mut writer = new_bam_writer(output_path, threads)?;
            writer.write_header(&output_header)?;
            process_bam_records(&mut reader, &input_header, &output_header, mapper, &mut writer, threads, &mut stats)?;
            finish_bam_writer(writer)?;
        }
    }

    Ok(stats)
}

/// BGZF worker count for a given thread budget.
///
/// Decompression and compression are memory-bandwidth bound and stop scaling
/// well before the conversion stage does, so their pools are capped. Handing
/// each of the three pools the full `threads` count oversubscribes the machine
/// and makes high `-t` values *slower* than moderate ones.
///
/// `FCM_BGZF_CAP` (reader) and `FCM_BGZF_WCAP` (writer) override the cap, for
/// measuring where the BGZF pools stop paying off.
fn bgzf_workers(threads: usize) -> NonZero<usize> {
    bgzf_workers_capped(threads, "FCM_BGZF_CAP")
}

fn bgzf_write_workers(threads: usize) -> NonZero<usize> {
    bgzf_workers_capped(threads, "FCM_BGZF_WCAP")
}

fn bgzf_workers_capped(threads: usize, env: &str) -> NonZero<usize> {
    let cap = std::env::var(env)
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(4)
        .max(1);
    NonZero::new(threads.clamp(1, cap)).expect("clamped to at least 1")
}

/// The BAM writer type: BGZF frames compressed on a pool, reassembled in order.
type BamWriter = bam::io::Writer<noodles_bgzf::io::MultithreadedWriter<std::fs::File>>;

/// Open a BAM output file with a BGZF compression pool.
///
/// Block boundaries depend only on the byte stream, not the worker count, so
/// the output is byte-identical for any `threads`.
fn new_bam_writer(output_path: &Path, threads: usize) -> Result<BamWriter, BamError> {
    let base = std::fs::File::create(output_path)?;
    Ok(bam::io::Writer::from(
        noodles_bgzf::io::MultithreadedWriter::with_worker_count(bgzf_write_workers(threads), base),
    ))
}

/// Flush the remaining BGZF blocks and append the EOF marker.
fn finish_bam_writer(writer: BamWriter) -> Result<(), BamError> {
    writer.into_inner().finish()?;
    Ok(())
}

/// PLACEHOLDER_PROCESS

/// How a converted alignment record is counted in the statistics.
#[derive(Default, Clone, Copy, PartialEq, Eq)]
enum BamClass {
    /// The record was placed at a reference position.
    #[default]
    Mapped,
    /// The record was already flagged unmapped on input and stayed that way.
    Unmapped,
    /// The record was mapped on input but could not be lifted over.
    Failed,
}

/// One converted alignment record.
#[derive(Default)]
struct BamOut {
    /// The record to emit. Every input record produces exactly one: a read
    /// that cannot be placed is written out with no reference, as CrossMap
    /// writes it.
    record: RecordBuf,
    /// How to count the record.
    class: BamClass,
    /// Whether the input record was part of a segment pair.
    segmented: bool,
    /// Which CrossMap category the record fell into, for the statistics.
    tag: Option<AlignmentTag>,
}

/// Number of alignment records converted as one unit.
///
/// `FCM_BATCH` overrides it, for measuring how batch size interacts with the
/// worker count. A batch that is too small spends its time in rayon's
/// split/collect hand-off rather than in conversion.
fn bam_batch_size() -> usize {
    crate::core::pipeline::batch_size_or(1024)
}

/// Accumulate one converted record into the statistics and write it out.
fn write_bam_out<W: sam::alignment::io::Write>(
    out: &BamOut,
    output_header: &sam::Header,
    writer: &mut W,
    stats: &mut ConversionStats,
) -> io::Result<()> {
    stats.total += 1;
    if out.segmented {
        stats.paired += 1;
    } else {
        stats.single += 1;
    }
    match out.class {
        BamClass::Mapped => stats.mapped += 1,
        BamClass::Unmapped => stats.unmapped += 1,
        BamClass::Failed => stats.failed += 1,
    }
    if let Some(tag) = out.tag {
        stats.count(tag);
    }
    writer.write_alignment_record(output_header, &out.record)?;
    Ok(())
}

fn process_bam_records<R: io::Read + Send, W: sam::alignment::io::Write>(
    reader: &mut bam::io::Reader<R>,
    input_header: &sam::Header,
    output_header: &sam::Header,
    mapper: &CoordinateMapper,
    writer: &mut W,
    threads: usize,
    stats: &mut ConversionStats,
) -> Result<(), BamError> {
    if threads <= 1 {
        return process_bam_records_sequential(
            reader, input_header, output_header, mapper, writer, stats,
        );
    }

    use crate::core::pipeline::{run_ordered_pipeline, ItemBatch};

    let is = InsertSizeParams::default();

    run_ordered_pipeline::<ItemBatch<bam::Record>, BamOut, _, _, _>(
        threads,
        |record: &bam::Record, out: &mut BamOut| {
            let decoded = decode_bam_record(record, input_header);
            let (class, converted, tag) =
                convert_decoded(&decoded, output_header, mapper, &is);
            out.segmented = decoded.flags.is_segmented();
            out.class = class;
            out.tag = Some(tag);
            out.record = converted;
            class == BamClass::Mapped
        },
        |batch: &mut ItemBatch<bam::Record>| {
            let batch_size = bam_batch_size();
            while batch.items.len() < batch_size {
                let mut record = bam::Record::default();
                if reader.read_record(&mut record)? == 0 {
                    break;
                }
                batch.items.push(record);
            }
            Ok(batch.items.len())
        },
        |_batch: &ItemBatch<bam::Record>, results: &[_]| {
            for result in results {
                write_bam_out(&result.output, output_header, writer, stats)?;
            }
            Ok(())
        },
    )?;

    Ok(())
}

fn process_bam_records_sequential<R: io::Read>(
    reader: &mut bam::io::Reader<R>,
    input_header: &sam::Header,
    output_header: &sam::Header,
    mapper: &CoordinateMapper,
    writer: &mut impl sam::alignment::io::Write,
    stats: &mut ConversionStats,
) -> Result<(), BamError> {
    let is = InsertSizeParams::default();

    for result in reader.records() {
        let record = result?;
        let decoded = decode_bam_record(&record, input_header);
        let (class, converted, tag) = convert_decoded(&decoded, output_header, mapper, &is);

        write_bam_out(
            &BamOut {
                record: converted,
                class,
                segmented: decoded.flags.is_segmented(),
                tag: Some(tag),
            },
            output_header,
            writer,
            stats,
        )?;
    }
    Ok(())
}

fn process_sam_records<R: io::BufRead + Send, W: sam::alignment::io::Write>(
    reader: &mut sam::io::Reader<R>,
    input_header: &sam::Header,
    output_header: &sam::Header,
    mapper: &CoordinateMapper,
    writer: &mut W,
    threads: usize,
    stats: &mut ConversionStats,
) -> Result<(), BamError> {
    if threads <= 1 {
        return process_sam_records_sequential(
            reader, input_header, output_header, mapper, writer, stats,
        );
    }

    use crate::core::pipeline::{run_ordered_pipeline, ItemBatch};

    let is = InsertSizeParams::default();

    run_ordered_pipeline::<ItemBatch<RecordBuf>, BamOut, _, _, _>(
        threads,
        |record_buf: &RecordBuf, out: &mut BamOut| {
            let decoded = decode_sam_record(record_buf, input_header);
            let (class, converted, tag) =
                convert_decoded(&decoded, output_header, mapper, &is);
            out.segmented = decoded.flags.is_segmented();
            out.class = class;
            out.tag = Some(tag);
            out.record = converted;
            class == BamClass::Mapped
        },
        |batch: &mut ItemBatch<RecordBuf>| {
            let batch_size = bam_batch_size();
            while batch.items.len() < batch_size {
                let mut record_buf = RecordBuf::default();
                if reader.read_record_buf(input_header, &mut record_buf)? == 0 {
                    break;
                }
                batch.items.push(record_buf);
            }
            Ok(batch.items.len())
        },
        |_batch: &ItemBatch<RecordBuf>, results: &[_]| {
            for result in results {
                write_bam_out(&result.output, output_header, writer, stats)?;
            }
            Ok(())
        },
    )?;

    Ok(())
}

fn process_sam_records_sequential<R: io::BufRead>(

    reader: &mut sam::io::Reader<R>,
    input_header: &sam::Header,
    output_header: &sam::Header,
    mapper: &CoordinateMapper,
    writer: &mut impl sam::alignment::io::Write,
    stats: &mut ConversionStats,
) -> Result<(), BamError> {
    let is = InsertSizeParams::default();

    for result in reader.record_bufs(input_header) {
        let record_buf = result?;
        let decoded = decode_sam_record(&record_buf, input_header);
        let (class, converted, tag) = convert_decoded(&decoded, output_header, mapper, &is);

        write_bam_out(
            &BamOut {
                record: converted,
                class,
                segmented: decoded.flags.is_segmented(),
                tag: Some(tag),
            },
            output_header,
            writer,
            stats,
        )?;
    }
    Ok(())
}

/// PLACEHOLDER_TESTS

#[cfg(test)]
mod tests {
    use super::*;

    /// Every shape of `(read1, read2)` mapping must land in its own category:
    /// unmapped (`None`), uniquely mapped (`Some(false)`), or mapped several
    /// places (`Some(true)`).
    #[test]
    fn test_paired_category_dispatch() {
        use AlignmentTag::*;
        assert_eq!(paired_category(None, None), NN);
        assert_eq!(paired_category(None, Some(false)), NU);
        assert_eq!(paired_category(None, Some(true)), NM);
        assert_eq!(paired_category(Some(false), None), UN);
        assert_eq!(paired_category(Some(true), None), MN);
        assert_eq!(paired_category(Some(false), Some(false)), UU);
        assert_eq!(paired_category(Some(false), Some(true)), UM);
        assert_eq!(paired_category(Some(true), Some(false)), MU);
        assert_eq!(paired_category(Some(true), Some(true)), MM);
    }

    #[test]
    fn test_crossmap_tags_promotes_char_and_moves_read_group_last() {
        let tag = |a: u8, b: u8| SamTag::new(a, b);
        let data: record_buf::Data = [
            (tag(b'X', b'T'), BufValue::Character(b'M')),
            (SamTag::READ_GROUP, BufValue::String(vec![b'G'].into())),
            (tag(b'N', b'M'), BufValue::Int32(8)),
            (tag(b'Z', b'C'), BufValue::String(vec![b'X'].into())),
        ]
        .into_iter()
        .collect();

        let got: Vec<_> = crossmap_tags(&data)
            .iter()
            .map(|(t, v)| (t, v.clone()))
            .collect();

        // Same order as the input, minus RG, which is appended last.
        assert_eq!(got[0], (tag(b'X', b'T'), BufValue::String(vec![b'M'].into())));
        assert_eq!(got[1], (tag(b'N', b'M'), BufValue::Int32(8)));
        // A one-character string field is left alone.
        assert_eq!(got[2], (tag(b'Z', b'C'), BufValue::String(vec![b'X'].into())));
        assert_eq!(got[3], (SamTag::READ_GROUP, BufValue::String(vec![b'G'].into())));
        assert_eq!(got.len(), 4);
    }

    #[test]
    fn test_crossmap_tags_without_read_group() {
        let tag = |a: u8, b: u8| SamTag::new(a, b);
        let data: record_buf::Data = [(tag(b'N', b'M'), BufValue::Int32(1))]
            .into_iter()
            .collect();

        let got: Vec<_> = crossmap_tags(&data).iter().map(|(t, _)| t).collect();
        assert_eq!(got, vec![tag(b'N', b'M')]);
    }

    #[test]
    fn test_cigar_op_from_noodles() {
        assert_eq!(CigarOp::from_noodles(SamOp::new(CigarKind::Match, 10)), CigarOp::Match(10));
        assert_eq!(CigarOp::from_noodles(SamOp::new(CigarKind::Insertion, 5)), CigarOp::Insertion(5));
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
        assert_eq!(CigarReconstructor::query_length(&cigar), 17);
    }

    #[test]
    fn test_cigar_reference_length() {
        let cigar = vec![CigarOp::Match(10), CigarOp::Insertion(2), CigarOp::Deletion(3), CigarOp::Match(5)];
        assert_eq!(CigarReconstructor::reference_length(&cigar), 18);
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