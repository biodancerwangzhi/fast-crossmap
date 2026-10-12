//! Wiggle/BigWig format adapter
//!
//! Handles Wiggle (variableStep, fixedStep) and BigWig format conversion.
//! Outputs bedGraph (.bgr) and optionally BigWig (.bw) files.
//!
//! **Validates: Requirements 9.1, 9.2, 9.3, 9.4, 9.5, 9.6**

use crate::core::{CoordinateMapper, Strand};
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, Write};
use std::path::{Path, PathBuf};

/// Wiggle parsing error
#[derive(Debug, Clone)]
pub enum WigParseError {
    EmptyLine,
    InvalidFormat(String),
    InvalidNumber(String),
    MissingChrom,
    MissingSpan,
    MissingStart,
    IoError(String),
}

impl std::fmt::Display for WigParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            WigParseError::EmptyLine => write!(f, "Empty line"),
            WigParseError::InvalidFormat(msg) => write!(f, "Invalid format: {}", msg),
            WigParseError::InvalidNumber(msg) => write!(f, "Invalid number: {}", msg),
            WigParseError::MissingChrom => write!(f, "Missing chrom parameter"),
            WigParseError::MissingSpan => write!(f, "Missing span parameter"),
            WigParseError::MissingStart => write!(f, "Missing start parameter"),
            WigParseError::IoError(msg) => write!(f, "IO error: {}", msg),
        }
    }
}

impl std::error::Error for WigParseError {}

/// Wiggle format type
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WigFormat {
    VariableStep,
    FixedStep,
}

/// Wiggle declaration line parameters
#[derive(Debug, Clone)]
pub struct WigDeclaration {
    pub format: WigFormat,
    pub chrom: String,
    pub span: u64,
    pub start: Option<u64>,  // Only for fixedStep
    pub step: Option<u64>,   // Only for fixedStep
}

impl WigDeclaration {
    /// Parse a declaration line (variableStep or fixedStep)
    pub fn parse(line: &str) -> Result<Self, WigParseError> {
        let line = line.trim();
        
        let (format, rest) = if line.starts_with("variableStep") {
            (WigFormat::VariableStep, &line[12..])
        } else if line.starts_with("fixedStep") {
            (WigFormat::FixedStep, &line[9..])
        } else {
            return Err(WigParseError::InvalidFormat(
                "Expected variableStep or fixedStep".to_string()
            ));
        };
        
        // Parse key=value pairs
        let mut chrom = None;
        let mut span = 1u64;
        let mut start = None;
        let mut step = None;
        
        for part in rest.split_whitespace() {
            if let Some((key, value)) = part.split_once('=') {
                match key {
                    "chrom" => chrom = Some(value.to_string()),
                    "span" => {
                        span = value.parse().map_err(|_| {
                            WigParseError::InvalidNumber(format!("span: {}", value))
                        })?;
                    }
                    "start" => {
                        start = Some(value.parse().map_err(|_| {
                            WigParseError::InvalidNumber(format!("start: {}", value))
                        })?);
                    }
                    "step" => {
                        step = Some(value.parse().map_err(|_| {
                            WigParseError::InvalidNumber(format!("step: {}", value))
                        })?);
                    }
                    _ => {} // Ignore unknown parameters
                }
            }
        }
        
        let chrom = chrom.ok_or(WigParseError::MissingChrom)?;
        
        // fixedStep requires start
        if format == WigFormat::FixedStep && start.is_none() {
            return Err(WigParseError::MissingStart);
        }
        
        Ok(Self {
            format,
            chrom,
            span,
            start,
            step,
        })
    }
}

/// A single Wiggle data point
#[derive(Debug, Clone)]
pub struct WigDataPoint {
    pub chrom: String,
    pub start: u64,  // 0-based
    pub end: u64,    // 0-based, exclusive
    pub value: f64,
}

/// bedGraph record for output
#[derive(Debug, Clone)]
pub struct BedGraphRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub value: f64,
}

impl BedGraphRecord {
    /// Format as bedGraph line
    pub fn to_line(&self) -> String {
        format!("{}\t{}\t{}\t{}", self.chrom, self.start, self.end, self.value)
    }
}

/// Conversion statistics
#[derive(Debug, Clone, Default)]
pub struct ConversionStats {
    pub total: usize,
    pub success: usize,
    pub failed: usize,
    pub merged: usize,
}

/// Parse a Wiggle file and yield data points
pub struct WigReader<R: BufRead> {
    reader: R,
    current_decl: Option<WigDeclaration>,
    current_pos: u64,  // For fixedStep
    line_buffer: String,
}

impl<R: BufRead> WigReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            current_decl: None,
            current_pos: 0,
            line_buffer: String::with_capacity(256),
        }
    }
}

impl<R: BufRead> Iterator for WigReader<R> {
    type Item = Result<WigDataPoint, WigParseError>;
    
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line_buffer.clear();
            match self.reader.read_line(&mut self.line_buffer) {
                Ok(0) => return None, // EOF
                Ok(_) => {}
                Err(e) => return Some(Err(WigParseError::IoError(e.to_string()))),
            }
            
            let line = self.line_buffer.trim();
            
            // Skip empty lines and comments
            if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
                continue;
            }
            
            // Check for declaration line
            if line.starts_with("variableStep") || line.starts_with("fixedStep") {
                match WigDeclaration::parse(line) {
                    Ok(decl) => {
                        if decl.format == WigFormat::FixedStep {
                            self.current_pos = decl.start.unwrap_or(1) - 1; // Convert to 0-based
                        }
                        self.current_decl = Some(decl);
                        continue;
                    }
                    Err(e) => return Some(Err(e)),
                }
            }
            
            // Parse data line
            // Check if it's a bedGraph line (4 columns: chrom start end value)
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 4 {
                // Try to parse as bedGraph format
                if let (Ok(start), Ok(end), Ok(value)) = (
                    parts[1].parse::<u64>(),
                    parts[2].parse::<u64>(),
                    parts[3].parse::<f64>(),
                ) {
                    return Some(Ok(WigDataPoint {
                        chrom: parts[0].to_string(),
                        start,
                        end,
                        value,
                    }));
                }
            }
            
            let decl = match &self.current_decl {
                Some(d) => d,
                None => {
                    return Some(Err(WigParseError::InvalidFormat(
                        "Data line before declaration".to_string()
                    )));
                }
            };
            
            match decl.format {
                WigFormat::VariableStep => {
                    // Format: position value
                    let parts: Vec<&str> = line.split_whitespace().collect();
                    if parts.len() < 2 {
                        return Some(Err(WigParseError::InvalidFormat(
                            format!("Expected position and value: {}", line)
                        )));
                    }
                    
                    let pos: u64 = match parts[0].parse() {
                        Ok(p) => p,
                        Err(_) => return Some(Err(WigParseError::InvalidNumber(parts[0].to_string()))),
                    };
                    let value: f64 = match parts[1].parse() {
                        Ok(v) => v,
                        Err(_) => return Some(Err(WigParseError::InvalidNumber(parts[1].to_string()))),
                    };
                    
                    // Wiggle uses 1-based coordinates, convert to 0-based
                    let start = pos - 1;
                    let end = start + decl.span;
                    
                    return Some(Ok(WigDataPoint {
                        chrom: decl.chrom.clone(),
                        start,
                        end,
                        value,
                    }));
                }
                WigFormat::FixedStep => {
                    // Format: value only
                    let value: f64 = match line.trim().parse() {
                        Ok(v) => v,
                        Err(_) => return Some(Err(WigParseError::InvalidNumber(line.to_string()))),
                    };
                    
                    let start = self.current_pos;
                    let end = start + decl.span;
                    
                    // Advance position for next data point
                    self.current_pos += decl.step.unwrap_or(decl.span);
                    
                    return Some(Ok(WigDataPoint {
                        chrom: decl.chrom.clone(),
                        start,
                        end,
                        value,
                    }));
                }
            }
        }
    }
}

/// Merge overlapping bedGraph records by splitting at boundaries and summing values.
/// After liftover, originally non-overlapping intervals can map to overlapping regions.
/// BigWig format forbids overlapping intervals, so we resolve them here.
fn merge_bedgraph_records(records: Vec<BedGraphRecord>) -> Vec<BedGraphRecord> {
    if records.is_empty() {
        return records;
    }

    let mut by_chrom: BTreeMap<String, Vec<BedGraphRecord>> = BTreeMap::new();
    for rec in records {
        by_chrom.entry(rec.chrom.clone()).or_default().push(rec);
    }

    let mut result = Vec::new();

    for (chrom, recs) in by_chrom {
        // Build sweep events: +value at start, -value at end
        let mut events: Vec<(u64, f64)> = Vec::with_capacity(recs.len() * 2);
        for rec in &recs {
            events.push((rec.start, rec.value));
            events.push((rec.end, -rec.value));
        }
        events.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)));

        // Sweep through events to build non-overlapping segments.
        //
        // Accumulating deltas in f64 does not cancel exactly: opening and closing
        // the same value can leave a residual around 1e-14. Left alone, that
        // residue leaks out as a phantom segment in the gap between two real
        // records, or as a tiny offset on the next real value. Snapping anything
        // below `EPS` back to zero keeps the sweep from drifting.
        const EPS: f64 = 1e-9;

        let mut merged_segments: Vec<(u64, u64, f64)> = Vec::new();
        let mut current_val = 0.0f64;
        let mut prev_pos = events[0].0;

        for (pos, delta) in &events {
            if *pos > prev_pos && current_val.abs() > EPS {
                merged_segments.push((prev_pos, *pos, current_val));
            }
            if *pos > prev_pos {
                prev_pos = *pos;
            }
            current_val += delta;
            if current_val.abs() < EPS {
                current_val = 0.0;
            }
        }

        // Coalesce adjacent segments with the same value
        for (seg_start, seg_end, val) in merged_segments {
            if let Some(last) = result.last_mut() {
                let last_rec: &mut BedGraphRecord = last;
                if last_rec.chrom == chrom
                    && last_rec.end == seg_start
                    && (last_rec.value - val).abs() < 1e-10
                {
                    last_rec.end = seg_end;
                    continue;
                }
            }
            result.push(BedGraphRecord {
                chrom: chrom.clone(),
                start: seg_start,
                end: seg_end,
                value: val,
            });
        }
    }

    result
}

/// Convert a single Wiggle data point
fn convert_wig_point(
    point: &WigDataPoint,
    mapper: &CoordinateMapper,
) -> Option<BedGraphRecord> {
    // Map coordinates
    let segments = mapper.map(&point.chrom, point.start, point.end, Strand::Plus)?;
    
    // Require single mapping
    if segments.len() != 1 {
        return None;
    }
    
    let seg = &segments[0];
    
    Some(BedGraphRecord {
        chrom: seg.target.chrom.clone(),
        start: seg.target.start,
        end: seg.target.end,
        value: point.value,
    })
}

/// Maximum size in bytes of one spooled record: `(start, end, value)`.
const SPOOL_RECORD_SIZE: u64 = 24;

/// A disk-backed spool of `(start, end, value)` triples, grouped by chromosome.
///
/// Wiggle liftover cannot write output as it reads. Both our output and
/// CrossMap's are grouped by chromosome in lexicographic order, while the input
/// may name chromosomes in any order and may revisit one more than once. Nor can
/// the merge run on a partial chromosome: splitting one chromosome across two
/// passes would emit two `variableStep` blocks where the output calls for one.
///
/// So each record is appended to a scratch file in a fixed-width binary form,
/// with a small in-memory index mapping chromosome name to the byte ranges it
/// occupies. Once the input is exhausted, chromosomes are replayed one at a
/// time: read back, merged, written out, dropped.
struct PointSpool {
    path: PathBuf,
    // `Option` so `Drop` can close the handle before unlinking (required on Windows).
    writer: Option<BufWriter<std::fs::File>>,
    offset: u64,
    // chrom -> [(byte offset, record count), ...]. More than one entry when the
    // input revisits a chromosome after emitting data for a different one.
    index: BTreeMap<String, Vec<(u64, u64)>>,
    count: usize,
}

impl PointSpool {
    fn create(path: &str) -> std::io::Result<Self> {
        let file = std::fs::File::create(path)?;
        Ok(Self {
            path: PathBuf::from(path),
            writer: Some(BufWriter::with_capacity(128 * 1024, file)),
            offset: 0,
            index: BTreeMap::new(),
            count: 0,
        })
    }

    fn push(&mut self, chrom: &str, start: u64, end: u64, value: f64) -> std::io::Result<()> {
        let mut buf = [0u8; SPOOL_RECORD_SIZE as usize];
        buf[0..8].copy_from_slice(&start.to_le_bytes());
        buf[8..16].copy_from_slice(&end.to_le_bytes());
        buf[16..24].copy_from_slice(&value.to_le_bytes());
        if let Some(writer) = self.writer.as_mut() {
            writer.write_all(&buf)?;
        }

        match self.index.get_mut(chrom) {
            Some(segments) => match segments.last_mut() {
                // Contiguous with this chromosome's previous run: extend it.
                Some(last) if last.0 + last.1 * SPOOL_RECORD_SIZE == self.offset => last.1 += 1,
                _ => segments.push((self.offset, 1)),
            },
            None => {
                self.index.insert(chrom.to_string(), vec![(self.offset, 1)]);
            }
        }

        self.offset += SPOOL_RECORD_SIZE;
        self.count += 1;
        Ok(())
    }

    fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// Chromosome names in lexicographic order — the order both tools write.
    fn chroms(&self) -> impl Iterator<Item = &str> {
        self.index.keys().map(|s| s.as_str())
    }

    /// Read back every record spooled for one chromosome, in input order.
    fn read_chrom(&mut self, chrom: &str) -> std::io::Result<Vec<BedGraphRecord>> {
        let segments = match self.index.get(chrom) {
            Some(segments) => segments.clone(),
            None => return Ok(Vec::new()),
        };
        if let Some(writer) = self.writer.as_mut() {
            writer.flush()?;
        }

        let total: u64 = segments.iter().map(|(_, count)| count).sum();
        let mut records = Vec::with_capacity(total as usize);
        // Decode straight from a buffered reader: reading the whole run into a
        // scratch buffer first would double this function's footprint.
        let mut file = BufReader::with_capacity(128 * 1024, std::fs::File::open(&self.path)?);
        let mut buf = [0u8; SPOOL_RECORD_SIZE as usize];
        for (offset, count) in segments {
            file.seek(std::io::SeekFrom::Start(offset))?;
            for _ in 0..count {
                file.read_exact(&mut buf)?;
                records.push(BedGraphRecord {
                    chrom: chrom.to_string(),
                    start: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
                    end: u64::from_le_bytes(buf[8..16].try_into().unwrap()),
                    value: f64::from_le_bytes(buf[16..24].try_into().unwrap()),
                });
            }
        }
        Ok(records)
    }
}

impl Drop for PointSpool {
    fn drop(&mut self) {
        // Close the write handle first: Windows refuses to unlink an open file.
        self.writer.take();
        let _ = std::fs::remove_file(&self.path);
    }
}

/// Map one batch of Wiggle/BigWig points and spool them, preserving input order.
///
/// The mapping (`convert_wig_point`, the chain lookup) runs across the pool, but
/// the two spool appends below run here, in the batch's original order — so the
/// scratch file's contents do not depend on the thread count. `pool` is `None`
/// for a single-threaded run and the mapping then runs inline; the result is the
/// same bytes either way.
fn flush_wig_batch(
    batch: &mut Vec<WigDataPoint>,
    mapper: &CoordinateMapper,
    pool: Option<&rayon::ThreadPool>,
    spool: &mut PointSpool,
    unmap_spool: &mut PointSpool,
    mapped_count: &mut usize,
    stats: &mut ConversionStats,
) -> Result<(), std::io::Error> {
    use rayon::prelude::*;

    if batch.is_empty() {
        return Ok(());
    }

    let converted: Vec<Option<BedGraphRecord>> = match pool {
        Some(pool) => pool.install(|| {
            batch
                .par_iter()
                .map(|point| convert_wig_point(point, mapper))
                .collect()
        }),
        None => batch.iter().map(|point| convert_wig_point(point, mapper)).collect(),
    };

    for (point, mapped) in batch.drain(..).zip(converted) {
        match mapped {
            Some(converted) => {
                spool.push(
                    &converted.chrom,
                    converted.start,
                    converted.end,
                    converted.value,
                )?;
                *mapped_count += 1;
                stats.success += 1;
            }
            None => {
                unmap_spool.push(&point.chrom, point.start, point.end, point.value)?;
                stats.failed += 1;
            }
        }
    }

    Ok(())
}

/// Convert a Wiggle file to Wiggle format (variableStep)
///
/// Records are spooled to a scratch file as they are read, then replayed one
/// chromosome at a time, so peak memory tracks the largest single chromosome
/// rather than the whole input (Requirement: bounded memory).
///
/// # Arguments
/// * `input` - Input Wiggle file path
/// * `output_prefix` - Output file prefix (will create .wig file)
/// * `mapper` - Coordinate mapper
/// * `threads` - Number of threads for parallel mapping (1 = sequential)
///
/// # Returns
/// Conversion statistics
pub fn convert_wig<P: AsRef<Path>>(
    input: P,
    output_prefix: P,
    mapper: &CoordinateMapper,
    threads: usize,
) -> Result<ConversionStats, std::io::Error> {
    let input_file = std::fs::File::open(input.as_ref())?;
    let reader = BufReader::with_capacity(128 * 1024, input_file);

    // Output files - use .wig extension for Wiggle format
    let output_path = format!("{}.wig", output_prefix.as_ref().display());
    let unmap_path = format!("{}.unmap.wig", output_prefix.as_ref().display());

    // Scratch files sit beside the output so they stay on one filesystem.
    let mapped_spool_path = format!("{}.tmp.mapped", output_prefix.as_ref().display());
    let unmapped_spool_path = format!("{}.tmp.unmapped", output_prefix.as_ref().display());

    let mut stats = ConversionStats::default();
    let mut mapped_count = 0usize;

    let mut spool = PointSpool::create(&mapped_spool_path)?;
    let mut unmap_spool = PointSpool::create(&unmapped_spool_path)?;

    // Phase 1: parse and map, spooling each record under its output chromosome.
    //
    // Only the per-record mapping runs across the pool. Parsing and the spool
    // append stay on this thread and in input order, so the spool holds exactly
    // the bytes it would single-threaded and phase 2 writes an identical file.
    let pool = crate::core::pipeline::build_pool(threads);
    let batch_size = crate::core::pipeline::batch_size_wig();
    let wig_reader = WigReader::new(reader);
    let mut batch: Vec<WigDataPoint> = Vec::with_capacity(batch_size);

    for result in wig_reader {
        match result {
            Ok(point) => {
                stats.total += 1;
                batch.push(point);
                if batch.len() >= batch_size {
                    flush_wig_batch(
                        &mut batch,
                        mapper,
                        pool.as_ref(),
                        &mut spool,
                        &mut unmap_spool,
                        &mut mapped_count,
                        &mut stats,
                    )?;
                }
            }
            Err(e) => {
                eprintln!("Warning: {}", e);
                stats.failed += 1;
            }
        }
    }
    flush_wig_batch(
        &mut batch,
        mapper,
        pool.as_ref(),
        &mut spool,
        &mut unmap_spool,
        &mut mapped_count,
        &mut stats,
    )?;

    // Phase 2: replay one chromosome at a time — merge, write, release.
    let mut merged_total = 0usize;
    {
        let mut output_file =
            BufWriter::with_capacity(128 * 1024, std::fs::File::create(&output_path)?);
        let chroms: Vec<String> = spool.chroms().map(|s| s.to_string()).collect();
        for chrom in chroms {
            let records = spool.read_chrom(&chrom)?;
            let merged = merge_bedgraph_records(records);
            merged_total += merged.len();
            write_wiggle_block(&mut output_file, &chrom, &merged)?;
        }
    }
    stats.merged = mapped_count.abs_diff(merged_total);

    // Unmapped points keep their source chromosome and, within a chromosome,
    // their input order.
    if !unmap_spool.is_empty() {
        let mut unmap_file =
            BufWriter::with_capacity(128 * 1024, std::fs::File::create(&unmap_path)?);
        let chroms: Vec<String> = unmap_spool.chroms().map(|s| s.to_string()).collect();
        for chrom in chroms {
            let records = unmap_spool.read_chrom(&chrom)?;
            write_wiggle_block(&mut unmap_file, &chrom, &records)?;
        }
    }

    Ok(stats)
}

/// Render a score the way Python's `repr` does.
///
/// CrossMap is a Python program, so its scores reach disk through `str(float)`:
/// shortest round-trip form, but always recognisable as a float (`44.0`, not
/// `44`). Rust's `Display` produces the same digits, minus the trailing `.0` on
/// whole numbers, so we add it back.
fn format_value(value: f64) -> String {
    let mut s = value.to_string();
    if !s.contains('.') && !s.contains('e') && !s.contains('E') && !s.contains("inf") && !s.contains("NaN") {
        s.push_str(".0");
    }
    s
}

/// Write one chromosome's records as a single `variableStep` block.
fn write_wiggle_block<W: Write>(
    output: &mut W,
    chrom: &str,
    records: &[BedGraphRecord],
) -> Result<(), std::io::Error> {
    // Determine span (use the most common span, default to 1)
    let span = if !records.is_empty() {
        records[0].end - records[0].start
    } else {
        1
    };

    // Write variableStep declaration
    if span > 1 {
        writeln!(output, "variableStep chrom={} span={}", chrom, span)?;
    } else {
        writeln!(output, "variableStep chrom={}", chrom)?;
    }

    // Write data points (convert 0-based to 1-based)
    for rec in records {
        writeln!(output, "{}\t{}", rec.start + 1, format_value(rec.value))?;
    }

    Ok(())
}

/// BigWig support module
pub mod bigwig {
    use super::*;
    use bigtools::BigWigRead;
    use std::collections::HashMap;
    
    /// Visit every interval in a BigWig file, chromosome by chromosome.
    ///
    /// Chromosomes are visited in the file's own order and intervals within a
    /// chromosome in coordinate order. The callback runs per interval, so callers
    /// can stream straight to disk instead of collecting the whole file — that
    /// `Vec<WigDataPoint>` was the single largest allocation on the BigWig path.
    pub fn for_each_bigwig_interval<P, F>(path: P, mut visit: F) -> Result<(), WigParseError>
    where
        P: AsRef<Path>,
        F: FnMut(&str, u64, u64, f64) -> Result<(), WigParseError>,
    {
        let mut reader = BigWigRead::open_file(path.as_ref().to_str().unwrap())
            .map_err(|e| WigParseError::IoError(e.to_string()))?;

        let chroms = reader.chroms().to_vec();

        for chrom_info in chroms {
            let chrom_name = chrom_info.name.clone();
            let chrom_len = chrom_info.length;

            // Read all intervals for this chromosome
            let intervals = reader
                .get_interval(&chrom_name, 0, chrom_len)
                .map_err(|e| WigParseError::IoError(e.to_string()))?;

            for interval in intervals {
                let interval = interval.map_err(|e| WigParseError::IoError(e.to_string()))?;
                visit(
                    &chrom_name,
                    interval.start as u64,
                    interval.end as u64,
                    interval.value as f64,
                )?;
            }
        }

        Ok(())
    }

    /// Read intervals from a BigWig file
    pub fn read_bigwig_intervals<P: AsRef<Path>>(
        path: P,
    ) -> Result<Vec<WigDataPoint>, WigParseError> {
        let mut points = Vec::new();
        for_each_bigwig_interval(path, |chrom, start, end, value| {
            points.push(WigDataPoint {
                chrom: chrom.to_string(),
                start,
                end,
                value,
            });
            Ok(())
        })?;
        Ok(points)
    }
    

    /// Write bedGraph records directly to a BigWig file using bigtools
    ///
    /// Records are spooled to a scratch file first so the streaming writer below
    /// can be reused without holding them all in memory.
    pub fn write_bigwig_direct<P: AsRef<Path>>(
        records: &[BedGraphRecord],
        output_path: P,
        chrom_sizes: &HashMap<String, u64>,
    ) -> Result<(), WigParseError> {
        let spool_path = format!("{}.records.tmp", output_path.as_ref().display());
        let mut spool = PointSpool::create(&spool_path)
            .map_err(|e| WigParseError::IoError(e.to_string()))?;
        for rec in records {
            spool
                .push(&rec.chrom, rec.start, rec.end, rec.value)
                .map_err(|e| WigParseError::IoError(e.to_string()))?;
        }

        let mut merged_total = 0usize;
        write_bigwig_from_spool(&mut spool, output_path.as_ref(), chrom_sizes, &mut merged_total)
    }

    /// Merge and write a spooled set of records out as a BigWig.
    ///
    /// Target sequences are visited in lexicographic order. Every sequence named
    /// in `chrom_sizes` gets an entry, with a placeholder `chrom\t0\t0\t0` line
    /// for those carrying no data — the same shape CrossMap writes. The number of
    /// merged records actually written is returned via `merged_total`.
    fn write_bigwig_from_spool(
        spool: &mut PointSpool,
        output_path: &Path,
        chrom_sizes: &HashMap<String, u64>,
        merged_total: &mut usize,
    ) -> Result<(), WigParseError> {
        use bigtools::BigWigWrite;
        use bigtools::beddata::BedParserStreamingIterator;

        // Create chrom sizes map for bigtools (u32)
        let chrom_map: HashMap<String, u32> = chrom_sizes
            .iter()
            .map(|(k, v)| (k.clone(), *v as u32))
            .collect();

        let mut targets: Vec<&str> = chrom_sizes.keys().map(|k| k.as_str()).collect();
        targets.sort_unstable();

        // Write to a temporary bedGraph file first, one chromosome at a time.
        let temp_bgr_path = format!("{}.temp.bedGraph", output_path.display());
        {
            let mut file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(&temp_bgr_path)
                .map_err(|e| WigParseError::IoError(e.to_string()))?);
            for chrom in targets {
                let records = spool
                    .read_chrom(chrom)
                    .map_err(|e| WigParseError::IoError(e.to_string()))?;
                if records.is_empty() {
                    writeln!(file, "{}\t0\t0\t0", chrom)
                        .map_err(|e| WigParseError::IoError(e.to_string()))?;
                    continue;
                }
                let merged = merge_bedgraph_records(records);
                *merged_total += merged.len();
                for rec in &merged {
                    writeln!(file, "{}\t{}\t{}\t{}", rec.chrom, rec.start, rec.end, rec.value)
                        .map_err(|e| WigParseError::IoError(e.to_string()))?;
                }
            }
        }

        // Use bigtools to convert bedGraph to BigWig
        let bedgraph_file = std::fs::File::open(&temp_bgr_path)
            .map_err(|e| WigParseError::IoError(e.to_string()))?;
        let vals = BedParserStreamingIterator::from_bedgraph_file(bedgraph_file, false);

        // Create tokio runtime for bigtools
        let runtime = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(2)
            .build()
            .map_err(|e| WigParseError::IoError(e.to_string()))?;

        // Create BigWig writer and write
        let writer = BigWigWrite::create_file(output_path, chrom_map)
            .map_err(|e| WigParseError::IoError(e.to_string()))?;

        writer.write(vals, runtime)
            .map_err(|e| WigParseError::IoError(format!("{:?}", e)))?;

        // Clean up temp file
        let _ = std::fs::remove_file(&temp_bgr_path);

        Ok(())
    }
    
    /// Map one batch of BigWig intervals and spool them, preserving input order.
    ///
    /// The mapping runs across the pool; the spool append and the `.unmap.bedGraph`
    /// write below stay on this thread and in the batch's original order, so both
    /// outputs are identical regardless of thread count. The unmap file is created
    /// lazily, on the first unmapped record, exactly as the single-threaded path
    /// does. `pool` is `None` for a single-threaded run.
    #[allow(clippy::too_many_arguments)]
    fn flush_bigwig_batch(
        batch: &mut Vec<WigDataPoint>,
        mapper: &CoordinateMapper,
        pool: Option<&rayon::ThreadPool>,
        spool: &mut PointSpool,
        unmap_file: &mut Option<BufWriter<std::fs::File>>,
        unmap_path: &str,
        mapped_count: &mut usize,
        stats: &mut ConversionStats,
    ) -> Result<(), std::io::Error> {
        use rayon::prelude::*;

        if batch.is_empty() {
            return Ok(());
        }

        let converted: Vec<Option<BedGraphRecord>> = match pool {
            Some(pool) => pool.install(|| {
                batch
                    .par_iter()
                    .map(|point| convert_wig_point(point, mapper))
                    .collect()
            }),
            None => batch.iter().map(|point| convert_wig_point(point, mapper)).collect(),
        };

        for (point, mapped) in batch.drain(..).zip(converted) {
            match mapped {
                Some(converted) => {
                    spool.push(
                        &converted.chrom,
                        converted.start,
                        converted.end,
                        converted.value,
                    )?;
                    *mapped_count += 1;
                    stats.success += 1;
                }
                None => {
                    if unmap_file.is_none() {
                        let file = std::fs::File::create(unmap_path)?;
                        *unmap_file = Some(BufWriter::with_capacity(64 * 1024, file));
                    }
                    if let Some(writer) = unmap_file.as_mut() {
                        writeln!(
                            writer,
                            "{}\t{}\t{}\t{}",
                            point.chrom, point.start, point.end, point.value
                        )?;
                    }
                    stats.failed += 1;
                }
            }
        }

        Ok(())
    }

    /// Convert a BigWig file to BigWig format
    ///
    /// Intervals are spooled as they are read and merged one chromosome at a
    /// time, so peak memory tracks the largest single chromosome rather than the
    /// whole input.
    ///
    /// # Arguments
    /// * `input` - Input BigWig file path
    /// * `output_prefix` - Output file prefix (will create .bw file)
    /// * `mapper` - Coordinate mapper
    /// * `threads` - Number of threads for parallel mapping (1 = sequential)
    ///
    /// # Returns
    /// Conversion statistics
    pub fn convert_bigwig<P: AsRef<Path>>(
        input: P,
        output_prefix: P,
        mapper: &CoordinateMapper,
        threads: usize,
    ) -> Result<ConversionStats, std::io::Error> {
        let unmap_path = format!("{}.unmap.bedGraph", output_prefix.as_ref().display());
        let spool_path = format!("{}.tmp.mapped", output_prefix.as_ref().display());

        let mut stats = ConversionStats::default();
        let mut mapped_count = 0usize;
        let mut spool = PointSpool::create(&spool_path)?;

        // Unmapped intervals are written straight through in input order, so this
        // stream costs nothing extra. (The Wiggle path instead groups `.unmap.wig`
        // by chromosome and therefore spools it — different output, different order.)
        let mut unmap_file: Option<BufWriter<std::fs::File>> = None;

        // As in convert_wig, only the mapping runs across the pool; the interval
        // reader and both output streams stay on this thread in input order, so
        // the spool and the `.unmap.bedGraph` are byte-for-byte thread-count
        // independent.
        let pool = crate::core::pipeline::build_pool(threads);
        let batch_size = crate::core::pipeline::batch_size_wig();
        let mut batch: Vec<WigDataPoint> = Vec::with_capacity(batch_size);

        for_each_bigwig_interval(&input, |chrom, start, end, value| {
            stats.total += 1;
            batch.push(WigDataPoint {
                chrom: chrom.to_string(),
                start,
                end,
                value,
            });
            if batch.len() >= batch_size {
                flush_bigwig_batch(
                    &mut batch,
                    mapper,
                    pool.as_ref(),
                    &mut spool,
                    &mut unmap_file,
                    &unmap_path,
                    &mut mapped_count,
                    &mut stats,
                )
                .map_err(|e| WigParseError::IoError(e.to_string()))?;
            }
            Ok(())
        })
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;

        flush_bigwig_batch(
            &mut batch,
            mapper,
            pool.as_ref(),
            &mut spool,
            &mut unmap_file,
            &unmap_path,
            &mut mapped_count,
            &mut stats,
        )
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;

        if let Some(mut writer) = unmap_file {
            writer.flush()?;
        }

        // Merge and write per chromosome, then hand off to bigtools.
        let bw_path = format!("{}.bw", output_prefix.as_ref().display());
        let mut merged_total = 0usize;
        if !spool.is_empty() {
            let target_sizes = mapper.target_sizes();
            write_bigwig_from_spool(&mut spool, Path::new(&bw_path), target_sizes, &mut merged_total)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        }
        stats.merged = mapped_count.abs_diff(merged_total);

        Ok(stats)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// Parse a BigWig's B+ chrom tree into `(name, id, length)` triples.
    fn read_chrom_tree(path: &std::path::Path) -> Vec<(String, u32, u32)> {
        let d = std::fs::read(path).unwrap();
        let rd_u32 = |o: usize| u32::from_le_bytes(d[o..o + 4].try_into().unwrap());
        let rd_u64 = |o: usize| u64::from_le_bytes(d[o..o + 8].try_into().unwrap());
        let ct = rd_u64(8) as usize;
        let key_size = rd_u32(ct + 8) as usize;
        let count = u16::from_le_bytes(d[ct + 34..ct + 36].try_into().unwrap()) as usize;
        let mut out = Vec::new();
        let mut o = ct + 36;
        for _ in 0..count {
            let name = String::from_utf8_lossy(&d[o..o + key_size])
                .trim_end_matches('\0')
                .to_string();
            let id = rd_u32(o + key_size);
            let len = rd_u32(o + key_size + 4);
            out.push((name, id, len));
            o += key_size + 8;
        }
        out
    }

    /// The dictionary must list every target sequence, not just the ones with data,
    /// with ids assigned in name-sorted order — the way CrossMap writes it.
    #[test]
    fn test_bigwig_chrom_tree_lists_all_targets() {
        let mut chrom_sizes = std::collections::HashMap::new();
        chrom_sizes.insert("chr1".to_string(), 1_000_000u64);
        chrom_sizes.insert("chr10".to_string(), 2_000_000u64);
        chrom_sizes.insert("chr2".to_string(), 3_000_000u64);

        // Only chr2 and chr10 carry data; chr1 must still appear in the dictionary.
        let records = vec![
            BedGraphRecord { chrom: "chr10".to_string(), start: 100, end: 200, value: 1.0 },
            BedGraphRecord { chrom: "chr2".to_string(), start: 300, end: 400, value: 2.0 },
        ];

        let path = std::env::temp_dir().join(format!(
            "fcm_bigwig_dict_test_{}.bw",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));

        bigwig::write_bigwig_direct(&records, &path, &chrom_sizes).unwrap();
        let tree = read_chrom_tree(&path);
        let _ = std::fs::remove_file(&path);

        assert_eq!(
            tree,
            vec![
                ("chr1".to_string(), 0, 1_000_000),
                ("chr10".to_string(), 1, 2_000_000),
                ("chr2".to_string(), 2, 3_000_000),
            ]
        );
    }

    #[test]
    fn test_variable_step_declaration() {
        let line = "variableStep chrom=chr1 span=10";
        let decl = WigDeclaration::parse(line).unwrap();
        
        assert_eq!(decl.format, WigFormat::VariableStep);
        assert_eq!(decl.chrom, "chr1");
        assert_eq!(decl.span, 10);
        assert!(decl.start.is_none());
        assert!(decl.step.is_none());
    }

    #[test]
    fn test_fixed_step_declaration() {
        let line = "fixedStep chrom=chr2 start=1000 step=100 span=50";
        let decl = WigDeclaration::parse(line).unwrap();
        
        assert_eq!(decl.format, WigFormat::FixedStep);
        assert_eq!(decl.chrom, "chr2");
        assert_eq!(decl.span, 50);
        assert_eq!(decl.start, Some(1000));
        assert_eq!(decl.step, Some(100));
    }

    #[test]
    fn test_variable_step_default_span() {
        let line = "variableStep chrom=chr1";
        let decl = WigDeclaration::parse(line).unwrap();
        
        assert_eq!(decl.span, 1); // Default span is 1
    }

    #[test]
    fn test_missing_chrom_error() {
        let line = "variableStep span=10";
        let result = WigDeclaration::parse(line);
        assert!(matches!(result, Err(WigParseError::MissingChrom)));
    }

    #[test]
    fn test_fixed_step_missing_start_error() {
        let line = "fixedStep chrom=chr1 step=100";
        let result = WigDeclaration::parse(line);
        assert!(matches!(result, Err(WigParseError::MissingStart)));
    }

    #[test]
    fn test_wig_reader_variable_step() {
        let wig_content = "\
variableStep chrom=chr1 span=10
1000 1.5
2000 2.5
3000 3.5
";
        let cursor = Cursor::new(wig_content.as_bytes());
        let reader = WigReader::new(std::io::BufReader::new(cursor));
        let points: Vec<_> = reader.collect();
        
        assert_eq!(points.len(), 3);
        
        let p0 = points[0].as_ref().unwrap();
        assert_eq!(p0.chrom, "chr1");
        assert_eq!(p0.start, 999); // 1000 - 1 (0-based)
        assert_eq!(p0.end, 1009);  // 999 + 10
        assert!((p0.value - 1.5).abs() < 1e-10);
        
        let p1 = points[1].as_ref().unwrap();
        assert_eq!(p1.start, 1999);
        assert_eq!(p1.end, 2009);
        
        let p2 = points[2].as_ref().unwrap();
        assert_eq!(p2.start, 2999);
        assert_eq!(p2.end, 3009);
    }

    #[test]
    fn test_wig_reader_fixed_step() {
        let wig_content = "\
fixedStep chrom=chr1 start=1000 step=100 span=50
1.0
2.0
3.0
";
        let cursor = Cursor::new(wig_content.as_bytes());
        let reader = WigReader::new(std::io::BufReader::new(cursor));
        let points: Vec<_> = reader.collect();
        
        assert_eq!(points.len(), 3);
        
        // First point: start=1000 (1-based) -> 999 (0-based), span=50
        let p0 = points[0].as_ref().unwrap();
        assert_eq!(p0.chrom, "chr1");
        assert_eq!(p0.start, 999);
        assert_eq!(p0.end, 1049);
        assert!((p0.value - 1.0).abs() < 1e-10);
        
        // Second point: 999 + 100 = 1099
        let p1 = points[1].as_ref().unwrap();
        assert_eq!(p1.start, 1099);
        assert_eq!(p1.end, 1149);
        
        // Third point: 1099 + 100 = 1199
        let p2 = points[2].as_ref().unwrap();
        assert_eq!(p2.start, 1199);
        assert_eq!(p2.end, 1249);
    }

    #[test]
    fn test_wig_reader_skip_comments() {
        let wig_content = "\
# This is a comment
track type=wiggle_0 name=\"test\"
browser position chr1:1000-2000
variableStep chrom=chr1 span=10
1000 1.5
";
        let cursor = Cursor::new(wig_content.as_bytes());
        let reader = WigReader::new(std::io::BufReader::new(cursor));
        let points: Vec<_> = reader.collect();
        
        assert_eq!(points.len(), 1);
        let p0 = points[0].as_ref().unwrap();
        assert_eq!(p0.chrom, "chr1");
    }

    #[test]
    fn test_wig_reader_multiple_chroms() {
        let wig_content = "\
variableStep chrom=chr1 span=10
1000 1.0
variableStep chrom=chr2 span=20
2000 2.0
";
        let cursor = Cursor::new(wig_content.as_bytes());
        let reader = WigReader::new(std::io::BufReader::new(cursor));
        let points: Vec<_> = reader.collect();
        
        assert_eq!(points.len(), 2);
        
        let p0 = points[0].as_ref().unwrap();
        assert_eq!(p0.chrom, "chr1");
        assert_eq!(p0.end - p0.start, 10);
        
        let p1 = points[1].as_ref().unwrap();
        assert_eq!(p1.chrom, "chr2");
        assert_eq!(p1.end - p1.start, 20);
    }

    #[test]
    fn test_bedgraph_record_to_line() {
        let rec = BedGraphRecord {
            chrom: "chr1".to_string(),
            start: 100,
            end: 200,
            value: 1.5,
        };
        
        assert_eq!(rec.to_line(), "chr1\t100\t200\t1.5");
    }

    #[test]
    fn test_merge_bedgraph_adjacent_same_value() {
        let records = vec![
            BedGraphRecord { chrom: "chr1".to_string(), start: 0, end: 100, value: 1.0 },
            BedGraphRecord { chrom: "chr1".to_string(), start: 100, end: 200, value: 1.0 },
            BedGraphRecord { chrom: "chr1".to_string(), start: 200, end: 300, value: 1.0 },
        ];

        let merged = merge_bedgraph_records(records);

        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 0);
        assert_eq!(merged[0].end, 300);
        assert!((merged[0].value - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_merge_bedgraph_different_values() {
        let records = vec![
            BedGraphRecord { chrom: "chr1".to_string(), start: 0, end: 100, value: 1.0 },
            BedGraphRecord { chrom: "chr1".to_string(), start: 100, end: 200, value: 2.0 },
        ];

        let merged = merge_bedgraph_records(records);

        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn test_merge_bedgraph_different_chroms() {
        let records = vec![
            BedGraphRecord { chrom: "chr1".to_string(), start: 0, end: 100, value: 1.0 },
            BedGraphRecord { chrom: "chr2".to_string(), start: 0, end: 100, value: 1.0 },
        ];

        let merged = merge_bedgraph_records(records);

        assert_eq!(merged.len(), 2);
    }

    #[test]
    fn test_merge_bedgraph_overlapping() {
        let records = vec![
            BedGraphRecord { chrom: "chr1".to_string(), start: 0, end: 150, value: 1.0 },
            BedGraphRecord { chrom: "chr1".to_string(), start: 100, end: 200, value: 1.0 },
        ];

        let merged = merge_bedgraph_records(records);

        // 0-100: val 1.0, 100-150: val 2.0 (summed), 150-200: val 1.0
        assert_eq!(merged.len(), 3);
        assert_eq!(merged[0].start, 0);
        assert_eq!(merged[0].end, 100);
        assert!((merged[0].value - 1.0).abs() < 1e-10);
        assert_eq!(merged[1].start, 100);
        assert_eq!(merged[1].end, 150);
        assert!((merged[1].value - 2.0).abs() < 1e-10);
        assert_eq!(merged[2].start, 150);
        assert_eq!(merged[2].end, 200);
        assert!((merged[2].value - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_merge_bedgraph_overlapping_different_values() {
        let records = vec![
            BedGraphRecord { chrom: "chr1".to_string(), start: 10, end: 50, value: 3.0 },
            BedGraphRecord { chrom: "chr1".to_string(), start: 30, end: 80, value: 5.0 },
        ];

        let merged = merge_bedgraph_records(records);

        // 10-30: 3.0, 30-50: 8.0, 50-80: 5.0
        assert_eq!(merged.len(), 3);
        assert_eq!((merged[0].start, merged[0].end), (10, 30));
        assert!((merged[0].value - 3.0).abs() < 1e-10);
        assert_eq!((merged[1].start, merged[1].end), (30, 50));
        assert!((merged[1].value - 8.0).abs() < 1e-10);
        assert_eq!((merged[2].start, merged[2].end), (50, 80));
        assert!((merged[2].value - 5.0).abs() < 1e-10);
    }

    /// Overlapping values that do not cancel exactly in f64 must not leave a
    /// phantom record in the gap before the next real one.
    #[test]
    fn test_merge_bedgraph_no_phantom_after_inexact_cancellation() {
        // 193.49 - 100.0 - 93.49 leaves a 1.4e-14 residue in f64.
        let records = vec![
            BedGraphRecord { chrom: "chr1".to_string(), start: 0, end: 10, value: 100.0 },
            BedGraphRecord { chrom: "chr1".to_string(), start: 5, end: 10, value: 93.49 },
            BedGraphRecord { chrom: "chr1".to_string(), start: 30, end: 40, value: 5.0 },
        ];

        let merged = merge_bedgraph_records(records);

        // 0-5: 100.0, 5-10: 193.49, then nothing until the real 30-40 record.
        assert_eq!(merged.len(), 3, "phantom record in the 10..30 gap: {:?}", merged);
        assert_eq!((merged[0].start, merged[0].end), (0, 5));
        assert!((merged[0].value - 100.0).abs() < 1e-10);
        assert_eq!((merged[1].start, merged[1].end), (5, 10));
        assert!((merged[1].value - 193.49).abs() < 1e-10);
        assert_eq!((merged[2].start, merged[2].end), (30, 40));
        assert!((merged[2].value - 5.0).abs() < 1e-10);
    }

    /// Scores must reach disk the way Python's `repr` writes them: 44.0, not 44.
    #[test]
    fn test_format_value_matches_python_repr() {
        assert_eq!(format_value(44.0), "44.0");
        assert_eq!(format_value(5.0), "5.0");
        assert_eq!(format_value(-0.0), "-0.0");
        assert_eq!(format_value(83.15), "83.15");
        assert_eq!(format_value(95.8), "95.8");
        assert_eq!(format_value(0.0), "0.0");
    }
}
