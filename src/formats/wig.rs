//! Wiggle/BigWig format adapter
//!
//! Handles Wiggle (variableStep, fixedStep) and BigWig format conversion.
//! Outputs bedGraph (.bgr) and optionally BigWig (.bw) files.
//!
//! **Validates: Requirements 9.1, 9.2, 9.3, 9.4, 9.5, 9.6**

use crate::core::{CoordinateMapper, Strand};
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

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

/// Merge overlapping bedGraph records with same value
fn merge_bedgraph_records(records: Vec<BedGraphRecord>) -> Vec<BedGraphRecord> {
    if records.is_empty() {
        return records;
    }
    
    // Group by chromosome
    let mut by_chrom: BTreeMap<String, Vec<BedGraphRecord>> = BTreeMap::new();
    for rec in records {
        by_chrom.entry(rec.chrom.clone()).or_default().push(rec);
    }
    
    let mut merged = Vec::new();
    
    for (chrom, mut recs) in by_chrom {
        // Sort by start position
        recs.sort_by_key(|r| r.start);
        
        let mut current: Option<BedGraphRecord> = None;
        
        for rec in recs {
            match current.take() {
                None => {
                    current = Some(rec);
                }
                Some(mut curr) => {
                    // Check if can merge (adjacent or overlapping with same value)
                    if rec.chrom == curr.chrom 
                        && rec.start <= curr.end 
                        && (rec.value - curr.value).abs() < 1e-10 
                    {
                        // Merge: extend end
                        curr.end = curr.end.max(rec.end);
                        current = Some(curr);
                    } else {
                        // Cannot merge, output current and start new
                        merged.push(BedGraphRecord {
                            chrom: chrom.clone(),
                            ..curr
                        });
                        current = Some(rec);
                    }
                }
            }
        }
        
        if let Some(curr) = current {
            merged.push(BedGraphRecord {
                chrom,
                ..curr
            });
        }
    }
    
    merged
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

/// Convert a Wiggle file to bedGraph format
///
/// # Arguments
/// * `input` - Input Wiggle file path
/// * `output_prefix` - Output file prefix (will create .bgr file)
/// * `mapper` - Coordinate mapper
///
/// # Returns
/// Conversion statistics
pub fn convert_wig<P: AsRef<Path>>(
    input: P,
    output_prefix: P,
    mapper: &CoordinateMapper,
) -> Result<ConversionStats, std::io::Error> {
    let input_file = std::fs::File::open(input.as_ref())?;
    let reader = BufReader::with_capacity(128 * 1024, input_file);
    
    // Output files
    let output_path = format!("{}.bgr", output_prefix.as_ref().display());
    let unmap_path = format!("{}.unmap.bgr", output_prefix.as_ref().display());
    
    let mut stats = ConversionStats::default();
    let mut converted_records = Vec::new();
    let mut unmapped_records = Vec::new();
    
    // Parse and convert
    let wig_reader = WigReader::new(reader);
    
    for result in wig_reader {
        match result {
            Ok(point) => {
                stats.total += 1;
                
                if let Some(converted) = convert_wig_point(&point, mapper) {
                    converted_records.push(converted);
                    stats.success += 1;
                } else {
                    unmapped_records.push(BedGraphRecord {
                        chrom: point.chrom,
                        start: point.start,
                        end: point.end,
                        value: point.value,
                    });
                    stats.failed += 1;
                }
            }
            Err(e) => {
                eprintln!("Warning: {}", e);
                stats.failed += 1;
            }
        }
    }
    
    // Merge overlapping records
    let original_count = converted_records.len();
    let merged_records = merge_bedgraph_records(converted_records);
    stats.merged = original_count - merged_records.len();
    
    // Write output with BufWriter for performance
    let mut output_file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(&output_path)?);
    for rec in &merged_records {
        writeln!(output_file, "{}", rec.to_line())?;
    }
    
    // Write unmapped
    if !unmapped_records.is_empty() {
        let mut unmap_file = BufWriter::with_capacity(64 * 1024, std::fs::File::create(&unmap_path)?);
        for rec in &unmapped_records {
            writeln!(unmap_file, "{}", rec.to_line())?;
        }
    }
    
    Ok(stats)
}

/// BigWig support module
pub mod bigwig {
    use super::*;
    use bigtools::BigWigRead;
    use std::collections::HashMap;
    
    /// Read intervals from a BigWig file
    pub fn read_bigwig_intervals<P: AsRef<Path>>(
        path: P,
    ) -> Result<Vec<WigDataPoint>, WigParseError> {
        let mut reader = BigWigRead::open_file(path.as_ref().to_str().unwrap())
            .map_err(|e| WigParseError::IoError(e.to_string()))?;
        
        let chroms = reader.chroms().to_vec();
        let mut points = Vec::new();
        
        for chrom_info in chroms {
            let chrom_name = chrom_info.name.clone();
            let chrom_len = chrom_info.length;
            
            // Read all intervals for this chromosome
            let intervals = reader
                .get_interval(&chrom_name, 0, chrom_len)
                .map_err(|e| WigParseError::IoError(e.to_string()))?;
            
            for interval in intervals {
                let interval = interval.map_err(|e| WigParseError::IoError(e.to_string()))?;
                points.push(WigDataPoint {
                    chrom: chrom_name.clone(),
                    start: interval.start as u64,
                    end: interval.end as u64,
                    value: interval.value as f64,
                });
            }
        }
        
        Ok(points)
    }
    
    /// Write bedGraph records to a BigWig file
    /// 
    /// Note: This is a simplified implementation that writes bedGraph first,
    /// then uses external tools (bedGraphToBigWig) for conversion.
    /// For production use, consider using bigtools CLI or bedGraphToBigWig.
    #[allow(dead_code)]
    pub fn write_bigwig_via_bedgraph<P: AsRef<Path>>(
        records: &[BedGraphRecord],
        output_path: P,
        chrom_sizes: &HashMap<String, u64>,
    ) -> Result<(), WigParseError> {
        // Write bedGraph file with BufWriter for performance
        let bgr_path = format!("{}.bgr", output_path.as_ref().display());
        {
            let mut file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(&bgr_path)
                .map_err(|e| WigParseError::IoError(e.to_string()))?);
            for rec in records {
                writeln!(file, "{}\t{}\t{}\t{}", rec.chrom, rec.start, rec.end, rec.value)
                    .map_err(|e| WigParseError::IoError(e.to_string()))?;
            }
        }
        
        // Write chrom.sizes file
        let sizes_path = format!("{}.chrom.sizes", output_path.as_ref().display());
        {
            let mut file = BufWriter::new(std::fs::File::create(&sizes_path)
                .map_err(|e| WigParseError::IoError(e.to_string()))?);
            for (chrom, size) in chrom_sizes {
                writeln!(file, "{}\t{}", chrom, size)
                    .map_err(|e| WigParseError::IoError(e.to_string()))?;
            }
        }
        
        // Try to use bedGraphToBigWig if available
        let result = std::process::Command::new("bedGraphToBigWig")
            .args(&[&bgr_path, &sizes_path, output_path.as_ref().to_str().unwrap()])
            .output();
        
        // Clean up temp files
        let _ = std::fs::remove_file(&sizes_path);
        
        match result {
            Ok(output) if output.status.success() => {
                let _ = std::fs::remove_file(&bgr_path);
                Ok(())
            }
            Ok(output) => {
                Err(WigParseError::IoError(format!(
                    "bedGraphToBigWig failed: {}",
                    String::from_utf8_lossy(&output.stderr)
                )))
            }
            Err(_) => {
                // bedGraphToBigWig not available, keep bedGraph file
                eprintln!("Warning: bedGraphToBigWig not found, keeping bedGraph output");
                Ok(())
            }
        }
    }
    
    /// Convert a BigWig file
    ///
    /// # Arguments
    /// * `input` - Input BigWig file path
    /// * `output_prefix` - Output file prefix (will create .bgr and optionally .bw files)
    /// * `mapper` - Coordinate mapper
    /// * `output_bigwig` - Whether to also write BigWig output (requires bedGraphToBigWig)
    ///
    /// # Returns
    /// Conversion statistics
    pub fn convert_bigwig<P: AsRef<Path>>(
        input: P,
        output_prefix: P,
        mapper: &CoordinateMapper,
        output_bigwig: bool,
    ) -> Result<ConversionStats, std::io::Error> {
        // Read BigWig intervals
        let points = read_bigwig_intervals(&input)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
        
        let mut stats = ConversionStats::default();
        let mut converted_records = Vec::new();
        let mut unmapped_records = Vec::new();
        
        // Convert each interval
        for point in points {
            stats.total += 1;
            
            if let Some(converted) = convert_wig_point(&point, mapper) {
                converted_records.push(converted);
                stats.success += 1;
            } else {
                unmapped_records.push(BedGraphRecord {
                    chrom: point.chrom,
                    start: point.start,
                    end: point.end,
                    value: point.value,
                });
                stats.failed += 1;
            }
        }
        
        // Merge overlapping records
        let original_count = converted_records.len();
        let merged_records = merge_bedgraph_records(converted_records);
        stats.merged = original_count - merged_records.len();
        
        // Write bedGraph output with BufWriter for performance
        let bgr_path = format!("{}.bgr", output_prefix.as_ref().display());
        let mut output_file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(&bgr_path)?);
        for rec in &merged_records {
            writeln!(output_file, "{}", rec.to_line())?;
        }
        
        // Write unmapped
        let unmap_path = format!("{}.unmap.bgr", output_prefix.as_ref().display());
        if !unmapped_records.is_empty() {
            let mut unmap_file = BufWriter::with_capacity(64 * 1024, std::fs::File::create(&unmap_path)?);
            for rec in &unmapped_records {
                writeln!(unmap_file, "{}", rec.to_line())?;
            }
        }
        
        // Optionally write BigWig output
        if output_bigwig && !merged_records.is_empty() {
            let bw_path = format!("{}.bw", output_prefix.as_ref().display());
            let target_sizes = mapper.target_sizes();
            
            write_bigwig_via_bedgraph(&merged_records, &bw_path, target_sizes)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))?;
        }
        
        Ok(stats)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

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
        
        assert_eq!(merged.len(), 1);
        assert_eq!(merged[0].start, 0);
        assert_eq!(merged[0].end, 200);
    }
}
