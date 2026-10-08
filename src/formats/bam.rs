//! BAM/SAM format adapter
//!
//! Handles BAM/SAM format conversion for alignment data.
//! Uses noodles (pure Rust) for reading and writing BAM files.

use crate::core::{CoordinateMapper, Strand};
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam::{
    self as sam,
    alignment::{
        record::Flags,
        record::cigar::{op::Kind as CigarKind, Op as SamOp},
        record_buf::{
            self,
            data::field::Value as BufValue,
        },
        RecordBuf,
    },
    header::record::value::{map::ReferenceSequence, Map},
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
    let mut builder = sam::Header::builder()
        .set_header(Default::default());

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

fn convert_record(
    record: &bam::Record,
    input_header: &sam::Header,
    output_header: &sam::Header,
    mapper: &CoordinateMapper,
) -> Option<(RecordBuf, AlignmentTag)> {
    let flags = record.flags();
    if flags.is_unmapped() { return None; }

    let tid = record.reference_sequence_id()?.ok()?;
    let chrom = get_chrom_name(input_header, tid)?;
    let start = record.alignment_start()?.ok()?.get() as u64 - 1;
    let cigar_ops = parse_cigar_from_bam(record);
    let ref_len = CigarReconstructor::reference_length(&cigar_ops) as u64;
    let end = start + ref_len;

    let query_strand = if flags.is_reverse_complemented() { Strand::Minus } else { Strand::Plus };
    let segments = mapper.map(&chrom, start, end, query_strand)?;
    if segments.is_empty() { return None; }

    let is_multiple = segments.len() > 1;
    let seg = &segments[0];
    let target_chrom = &seg.target.chrom;
    let target_start = seg.target.start;
    let target_strand = seg.target.strand;
    let target_tid = get_tid(output_header, target_chrom)?;

    let need_revcomp = target_strand == Strand::Minus && query_strand == Strand::Plus
        || target_strand == Strand::Plus && query_strand == Strand::Minus;

    let seq_bytes = seq_bytes_from_record(record);
    let qual_bytes = qual_bytes_from_record(record);

    let (new_seq, new_qual, new_cigar) = if need_revcomp {
        (revcomp_seq(&seq_bytes), reverse_qual(&qual_bytes), CigarReconstructor::reverse(&cigar_ops))
    } else {
        (seq_bytes, qual_bytes, cigar_ops)
    };

    let noodles_cigar = ops_to_noodles_cigar(&new_cigar);
    let pos = Position::try_from((target_start + 1) as usize).ok()?;

    let mut new_flags = flags;
    if need_revcomp { new_flags ^= Flags::REVERSE_COMPLEMENTED; }
    if is_multiple { new_flags |= Flags::SECONDARY; }

    let aux_data = copy_aux_data(record);

    let mut builder = RecordBuf::builder()
        .set_name(record.name().unwrap_or(b"*".into()))
        .set_flags(new_flags)
        .set_reference_sequence_id(target_tid)
        .set_alignment_start(pos)
        .set_cigar(noodles_cigar)
        .set_template_length(0)
        .set_sequence(record_buf::Sequence::from(new_seq.as_slice()))
        .set_quality_scores(record_buf::QualityScores::from(new_qual))
        .set_data(aux_data);
    if let Some(mq) = record.mapping_quality() {
        builder = builder.set_mapping_quality(mq);
    }
    let new_record = builder
        .build();

    let tag = if flags.is_segmented() {
        if flags.is_mate_unmapped() { AlignmentTag::MU } else { AlignmentTag::MM }
    } else if is_multiple { AlignmentTag::SM } else { AlignmentTag::SU };

    Some((new_record, tag))
}

/// PLACEHOLDER_MAIN

/// Build the placeholder record CrossMap writes for a read it could not place.
///
/// The read keeps its name, bases and qualities, but loses its reference,
/// position, CIGAR and template length. Leaving the mapping quality unset is
/// how 255 ("unavailable") is written.
///
/// `flags` is the caller's choice: CrossMap resets the flags to a bare
/// `FLAG_UNMAPPED` for a read that failed to lift over, but ORs it into the
/// incoming flags for a read that was already unmapped.
fn unmapped_placeholder(name: &[u8], seq: &[u8], qual: &[u8], flags: Flags) -> RecordBuf {
    RecordBuf::builder()
        .set_name(name)
        .set_flags(flags)
        .set_sequence(record_buf::Sequence::from(seq))
        .set_quality_scores(record_buf::QualityScores::from(qual.to_vec()))
        .build()
}

/// The placeholder for a record that was already flagged unmapped on input,
/// built straight from the input `bam::Record`.
fn already_unmapped_bam_placeholder(record: &bam::Record) -> RecordBuf {
    unmapped_placeholder(
        record.name().map(|n| &**n).unwrap_or(&b"*"[..]),
        &seq_bytes_from_record(record),
        &qual_bytes_from_record(record),
        record.flags() | Flags::UNMAPPED,
    )
}

/// The placeholder for a record that could not be converted, built straight
/// from the input `bam::Record`.
fn failed_bam_placeholder(record: &bam::Record) -> RecordBuf {
    unmapped_placeholder(
        record.name().map(|n| &**n).unwrap_or(&b"*"[..]),
        &seq_bytes_from_record(record),
        &qual_bytes_from_record(record),
        Flags::UNMAPPED,
    )
}

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
    /// The record mapped.
    #[default]
    Mapped,
    /// The record was already flagged unmapped on input.
    Unmapped,
    /// The record could not be mapped and was replaced (BAM) or dropped (SAM).
    Failed,
}

/// One converted alignment record.
#[derive(Default)]
struct BamOut {
    /// The record to emit, or `None` to drop it. The SAM path drops records
    /// whose coordinates could not be resolved; the BAM path always emits one.
    record: Option<RecordBuf>,
    /// How to count the record.
    class: BamClass,
    /// Whether the input record was part of a segment pair.
    segmented: bool,
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
    if let Some(record) = &out.record {
        writer.write_alignment_record(output_header, record)?;
    }
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

    run_ordered_pipeline::<ItemBatch<bam::Record>, BamOut, _, _, _>(
        threads,
        |record: &bam::Record, out: &mut BamOut| {
            let flags = record.flags();
            out.segmented = flags.is_segmented();
            if flags.is_unmapped() {
                out.class = BamClass::Unmapped;
                out.record = Some(already_unmapped_bam_placeholder(record));
                return false;
            }
            match convert_record(record, input_header, output_header, mapper) {
                Some((new_record, _tag)) => {
                    out.class = BamClass::Mapped;
                    out.record = Some(new_record);
                    true
                }
                None => {
                    out.class = BamClass::Failed;
                    out.record = Some(failed_bam_placeholder(record));
                    false
                }
            }
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
    for result in reader.records() {
        let record = result?;
        stats.total += 1;
        let flags = record.flags();
        if flags.is_segmented() { stats.paired += 1; } else { stats.single += 1; }

        if flags.is_unmapped() {
            stats.unmapped += 1;
            let new_record = already_unmapped_bam_placeholder(&record);
            writer.write_alignment_record(output_header, &new_record)?;
            continue;
        }

        match convert_record(&record, input_header, output_header, mapper) {
            Some((new_record, _tag)) => {
                writer.write_alignment_record(output_header, &new_record)?;
                stats.mapped += 1;
            }
            None => {
                stats.failed += 1;
                let new_record = failed_bam_placeholder(&record);
                writer.write_alignment_record(output_header, &new_record)?;
            }
        }
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

    run_ordered_pipeline::<ItemBatch<RecordBuf>, BamOut, _, _, _>(
        threads,
        |record_buf: &RecordBuf, out: &mut BamOut| {
            let flags = record_buf.flags();
            out.segmented = flags.is_segmented();

            // A record already flagged unmapped is passed through untouched.
            if flags.is_unmapped() {
                out.class = BamClass::Unmapped;
                out.record = Some(unmapped_placeholder(
                    record_buf.name().map(|n| &**n).unwrap_or(&b"*"[..]),
                    record_buf.sequence().as_ref(),
                    record_buf.quality_scores().as_ref(),
                    flags | Flags::UNMAPPED,
                ));
                return false;
            }

            match convert_sam_record(record_buf, input_header, output_header, mapper) {
                Some(record) => {
                    out.class = BamClass::Mapped;
                    out.record = Some(record);
                    true
                }
                None => {
                    // CrossMap keeps a read that failed to lift over, as an
                    // unmapped record with the original bases.
                    out.class = BamClass::Failed;
                    out.record = Some(unmapped_placeholder(
                        record_buf.name().map(|n| &**n).unwrap_or(&b"*"[..]),
                        record_buf.sequence().as_ref(),
                        record_buf.quality_scores().as_ref(),
                        Flags::UNMAPPED,
                    ));
                    false
                }
            }
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

/// Convert one SAM record buffer, returning `None` if it cannot be placed.
///
/// This mirrors the body of the sequential loop: every early `continue` in it
/// (missing reference id/name/start, empty mapping, unknown target) becomes a
/// `None` here.
fn convert_sam_record(
    record_buf: &RecordBuf,
    input_header: &sam::Header,
    output_header: &sam::Header,
    mapper: &CoordinateMapper,
) -> Option<RecordBuf> {
    let flags = record_buf.flags();
    let tid = record_buf.reference_sequence_id()?;
    let chrom = get_chrom_name(input_header, tid)?;
    let start_pos = record_buf.alignment_start()?.get() as u64 - 1;

    let cigar_ops: Vec<CigarOp> = record_buf.cigar().as_ref().iter().map(|op| CigarOp::from_noodles(*op)).collect();
    let ref_len = CigarReconstructor::reference_length(&cigar_ops) as u64;
    let end = start_pos + ref_len;

    let query_strand = if flags.is_reverse_complemented() { Strand::Minus } else { Strand::Plus };

    let segments = mapper.map(&chrom, start_pos, end, query_strand)?;
    if segments.is_empty() { return None; }

    let is_multiple = segments.len() > 1;
    let seg = &segments[0];
    let target_tid = get_tid(output_header, &seg.target.chrom)?;

    let target_strand = seg.target.strand;
    let need_revcomp = target_strand == Strand::Minus && query_strand == Strand::Plus
        || target_strand == Strand::Plus && query_strand == Strand::Minus;

    let seq_vec: Vec<u8> = record_buf.sequence().as_ref().to_vec();
    let qual_vec: Vec<u8> = record_buf.quality_scores().as_ref().to_vec();

    let (new_seq, new_qual, new_cigar) = if need_revcomp {
        (revcomp_seq(&seq_vec), reverse_qual(&qual_vec), CigarReconstructor::reverse(&cigar_ops))
    } else {
        (seq_vec, qual_vec, cigar_ops)
    };

    let noodles_cigar = ops_to_noodles_cigar(&new_cigar);
    let pos = Position::try_from((seg.target.start + 1) as usize).ok()?;

    let mut new_flags = flags;
    if need_revcomp { new_flags ^= Flags::REVERSE_COMPLEMENTED; }
    if is_multiple { new_flags |= Flags::SECONDARY; }

    let mut builder = RecordBuf::builder()
        .set_name(record_buf.name().unwrap_or(b"*".into()))
        .set_flags(new_flags)
        .set_reference_sequence_id(target_tid)
        .set_alignment_start(pos)
        .set_cigar(noodles_cigar)
        .set_template_length(0)
        .set_sequence(record_buf::Sequence::from(new_seq.as_slice()))
        .set_quality_scores(record_buf::QualityScores::from(new_qual))
        .set_data(record_buf.data().clone());
    if let Some(mq) = record_buf.mapping_quality() {
        builder = builder.set_mapping_quality(mq);
    }
    Some(builder.build())
}

fn process_sam_records_sequential<R: io::BufRead>(

    reader: &mut sam::io::Reader<R>,
    input_header: &sam::Header,
    output_header: &sam::Header,
    mapper: &CoordinateMapper,
    writer: &mut impl sam::alignment::io::Write,
    stats: &mut ConversionStats,
) -> Result<(), BamError> {
    for result in reader.record_bufs(input_header) {
        let record_buf = result?;
        stats.total += 1;
        let flags = record_buf.flags();
        if flags.is_segmented() { stats.paired += 1; } else { stats.single += 1; }

        // A read keeps its name, bases and qualities whether it was already
        // unmapped or could not be lifted over; only the flags differ.
        let name = record_buf.name().map(|n| &**n).unwrap_or(&b"*"[..]);
        let seq = record_buf.sequence().as_ref();
        let qual = record_buf.quality_scores().as_ref();

        if flags.is_unmapped() {
            stats.unmapped += 1;
            let new_record = unmapped_placeholder(name, seq, qual, flags | Flags::UNMAPPED);
            writer.write_alignment_record(output_header, &new_record)?;
            continue;
        }

        match convert_sam_record(&record_buf, input_header, output_header, mapper) {
            Some(new_record) => {
                writer.write_alignment_record(output_header, &new_record)?;
                stats.mapped += 1;
            }
            None => {
                // CrossMap keeps a read that failed to lift over, as an
                // unmapped record carrying the original bases.
                stats.failed += 1;
                let new_record = unmapped_placeholder(name, seq, qual, Flags::UNMAPPED);
                writer.write_alignment_record(output_header, &new_record)?;
            }
        }
    }
    Ok(())
}

/// PLACEHOLDER_TESTS

#[cfg(test)]
mod tests {
    use super::*;

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