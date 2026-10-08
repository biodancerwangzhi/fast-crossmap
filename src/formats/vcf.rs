//! VCF format adapter
//!
//! Handles VCF format conversion with zero-copy parsing.
//!
//! **Validates: Requirements 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7**

use crate::core::{chr_template_from_contig_line, dna, update_chrom_id_by_template, CoordinateMapper, Strand};
use std::cell::{Cell, RefCell};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// File name of the reference genome, as CrossMap reports it in `assembly=`.
fn ref_basename<P: AsRef<Path>>(ref_genome: &Option<P>) -> Option<String> {
    ref_genome
        .as_ref()
        .and_then(|p| p.as_ref().file_name())
        .map(|n| n.to_string_lossy().into_owned())
}

/// Write the target `##contig` block that CrossMap emits in place of `#CHROM`.
///
/// The contigs come from the *reference* FASTA (not the chain), are sorted by
/// name, rewritten into the input's chromosome style, and tagged with the
/// reference file name — mirroring `cmmodule/mapvcf.py`.
fn write_contig_header(
    out: &mut impl Write,
    ref_reader: Option<&crate::core::fasta::FastaReader>,
    chr_template: &str,
    ref_basename: Option<&str>,
) -> std::io::Result<()> {
    let Some(rdr) = ref_reader else {
        return Ok(());
    };
    let mut contigs: Vec<(String, usize)> =
        rdr.references().iter().cloned().zip(rdr.lengths()).collect();
    contigs.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    for (chrom, len) in contigs {
        let id = update_chrom_id_by_template(&chrom, chr_template);
        match ref_basename {
            Some(assembly) => {
                writeln!(out, "##contig=<ID={},length={},assembly={}>", id, len, assembly)?
            }
            None => writeln!(out, "##contig=<ID={},length={}>", id, len)?,
        }
    }
    Ok(())
}

/// VCF record representation for output
#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub chrom: String,
    pub pos: u64,
    pub id: String,
    pub ref_allele: String,
    pub alt_alleles: Vec<String>,
    pub qual: String,
    pub filter: String,
    pub info: String,
    pub format: Option<String>,
    pub samples: Vec<String>,
}

/// Zero-copy VCF record view for parsing.
///
/// Only CHROM and POS are decoded eagerly; every other column is addressed
/// lazily by scanning the original line. Field bounds for the leading
/// (non-sample) columns are stored inline, so parsing a record performs
/// **no heap allocation** — which matters because VCF lines routinely carry
/// thousands of sample columns.
pub struct VcfRecordView<'a> {
    /// Original line (guaranteed valid UTF-8 by the reader)
    line: &'a str,
    /// Chromosome name (field 0)
    pub chrom: &'a str,
    /// Position, 1-based (field 1)
    pub pos: u64,
    /// Byte bounds of the first 9 fields (CHROM..FORMAT).
    /// Only the first `lead_count` entries are populated.
    lead: [(usize, usize); 9],
    /// Number of populated entries in `lead`
    lead_count: usize,
    /// Byte offset where sample columns begin (immediately after FORMAT),
    /// or `None` when the record has no sample columns.
    sample_start: Option<usize>,
    /// Cached INFO parsing
    info_parsed: Cell<bool>,
    info_cache: RefCell<Option<HashMap<String, String>>>,
}

impl<'a> VcfRecordView<'a> {
    /// Parse a VCF line without heap allocation.
    ///
    /// CHROM and POS are decoded immediately; all other columns are addressed
    /// lazily via [`field`](Self::field).
    pub fn parse(line: &'a [u8]) -> Result<Self, VcfParseError> {
        if line.is_empty() {
            return Err(VcfParseError::EmptyLine);
        }
        let s = std::str::from_utf8(line).map_err(|_| VcfParseError::InvalidUtf8("line"))?;

        // Locate the leading fields (and the start of the sample block) in a
        // single pass, storing bounds inline instead of in a heap Vec.
        let mut lead = [(0usize, 0usize); 9];
        let mut lead_count = 0usize;
        let mut sample_start = None;

        let mut field_start = 0usize;
        let mut idx = 0usize;
        for (i, &b) in line.iter().enumerate() {
            if b == b'\t' {
                if idx < 9 {
                    lead[idx] = (field_start, i);
                    lead_count = idx + 1;
                }
                field_start = i + 1;
                idx += 1;
                if idx == 9 {
                    // Closed the FORMAT column (field 8): sample columns follow.
                    sample_start = Some(i + 1);
                    break;
                }
            }
        }
        // Final field, when the line does not end with a tab.
        if sample_start.is_none() && idx < 9 {
            lead[idx] = (field_start, line.len());
            lead_count = idx + 1;
        }

        // VCF requires at least 8 fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
        if lead_count < 8 {
            return Err(VcfParseError::TooFewFields {
                expected: 8,
                found: lead_count,
            });
        }

        // Parse CHROM (field 0)
        let (cs, ce) = lead[0];
        let chrom = &s[cs..ce];

        // Parse POS (field 1)
        let (ps, pe) = lead[1];
        let pos_str = &s[ps..pe];
        let pos: u64 = pos_str
            .parse()
            .map_err(|_| VcfParseError::InvalidNumber("POS", pos_str.to_string()))?;

        Ok(Self {
            line: s,
            chrom,
            pos,
            lead,
            lead_count,
            sample_start,
            info_parsed: Cell::new(false),
            info_cache: RefCell::new(None),
        })
    }

    /// Get the number of fields in the record.
    pub fn field_count(&self) -> usize {
        match self.sample_start {
            None => self.lead_count,
            Some(s) => {
                let rest = &self.line[s..];
                if rest.is_empty() {
                    self.lead_count
                } else {
                    self.lead_count + rest.matches('\t').count() + 1
                }
            }
        }
    }

    /// Get field as string slice (lazy access).
    pub fn field(&self, index: usize) -> Option<&'a str> {
        let line = self.line;
        if index < self.lead_count {
            let (s, e) = self.lead[index];
            return Some(&line[s..e]);
        }
        // Walk the sample columns to reach `index`.
        let mut cur = self.sample_start?;
        let mut n = self.lead_count;
        while n < index {
            cur += line[cur..].find('\t')? + 1;
            n += 1;
        }
        if cur > line.len() {
            return None;
        }
        let end = line[cur..].find('\t').map(|t| cur + t).unwrap_or(line.len());
        Some(&line[cur..end])
    }

    /// Get ID field (field 2)
    pub fn id(&self) -> Option<&'a str> {
        self.field(2)
    }

    /// Get REF field (field 3)
    pub fn ref_allele(&self) -> Option<&'a str> {
        self.field(3)
    }

    /// Get ALT field (field 4)
    pub fn alt_alleles(&self) -> Option<&'a str> {
        self.field(4)
    }

    /// Get QUAL field (field 5)
    pub fn qual(&self) -> Option<&'a str> {
        self.field(5)
    }

    /// Get FILTER field (field 6)
    pub fn filter(&self) -> Option<&'a str> {
        self.field(6)
    }

    /// Get INFO field (field 7)
    pub fn info(&self) -> Option<&'a str> {
        self.field(7)
    }

    /// Get FORMAT field (field 8) if present
    pub fn format(&self) -> Option<&'a str> {
        self.field(8)
    }

    /// Visit each sample column (fields 9+) in order, without allocating.
    pub fn for_each_sample<F: FnMut(&'a str)>(&self, mut f: F) {
        let line = self.line;
        let mut cur = match self.sample_start {
            Some(s) if s < line.len() => s,
            _ => return,
        };
        loop {
            match line[cur..].find('\t') {
                Some(t) => {
                    f(&line[cur..cur + t]);
                    cur += t + 1;
                }
                None => {
                    f(&line[cur..]);
                    break;
                }
            }
        }
    }

    /// Collect the sample columns (fields 9+) into a vector.
    pub fn samples(&self) -> Vec<&'a str> {
        let mut v = Vec::new();
        self.for_each_sample(|s| v.push(s));
        v
    }

    /// The FORMAT column and all sample columns as a single contiguous slice
    /// (everything from FORMAT to the end of the line), or `None` when the
    /// record has fewer than 9 fields.
    ///
    /// The sample columns are never rewritten during liftover, so a whole-tail
    /// slice lets the output path copy them with one `push_str` instead of one
    /// call per sample. Equivalent to `format()` followed by
    /// [`for_each_sample`](Self::for_each_sample), including the trailing-tab
    /// and no-sample cases.
    pub fn format_and_samples(&self) -> Option<&'a str> {
        if self.lead_count < 9 {
            return None;
        }
        let (fs, fe) = self.lead[8];
        match self.sample_start {
            // `sample_start` is only set when at least one byte follows FORMAT's tab.
            Some(s) if s < self.line.len() => Some(&self.line[fs..]),
            _ => Some(&self.line[fs..fe]),
        }
    }

    /// Parse INFO field lazily (only when needed)
    pub fn parse_info(&self) -> HashMap<String, String> {
        if !self.info_parsed.get() {
            let info_str = self.info().unwrap_or(".");
            let mut map = HashMap::new();
            
            if info_str != "." {
                for item in info_str.split(';') {
                    if let Some(eq_pos) = item.find('=') {
                        let key = item[..eq_pos].to_string();
                        let value = item[eq_pos + 1..].to_string();
                        map.insert(key, value);
                    } else {
                        // Flag without value
                        map.insert(item.to_string(), String::new());
                    }
                }
            }
            
            *self.info_cache.borrow_mut() = Some(map.clone());
            self.info_parsed.set(true);
            map
        } else {
            self.info_cache.borrow().clone().unwrap_or_default()
        }
    }
    
    /// Get variant type based on REF and ALT lengths
    pub fn variant_type(&self) -> VariantType {
        let ref_len = self.ref_allele().map(|s| s.len()).unwrap_or(0);
        let alt = self.alt_alleles().unwrap_or(".");
        
        // Get first ALT allele for type determination
        let first_alt = alt.split(',').next().unwrap_or(".");
        let alt_len = first_alt.len();
        
        if ref_len == alt_len {
            VariantType::Substitution
        } else if alt_len > ref_len {
            VariantType::Insertion
        } else {
            VariantType::Deletion
        }
    }
}

/// Variant type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantType {
    Substitution,
    Insertion,
    Deletion,
}

/// VCF parsing error
#[derive(Debug, thiserror::Error)]
pub enum VcfParseError {
    #[error("Empty line")]
    EmptyLine,
    
    #[error("Too few fields: expected at least {expected}, found {found}")]
    TooFewFields { expected: usize, found: usize },
    
    #[error("Invalid UTF-8 in field: {0}")]
    InvalidUtf8(&'static str),
    
    #[error("Invalid number in field {0}: {1}")]
    InvalidNumber(&'static str, String),
    
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}

/// Conversion statistics
#[derive(Debug, Default, Clone)]
pub struct ConversionStats {
    pub total: usize,
    pub success: usize,
    pub failed: usize,
}

/// Convert a single VCF record.
///
/// On success the formatted output line is appended to `out` and `Ok(())` is
/// returned; otherwise `Err(reason)` carries the failure reason (the caller
/// supplies the original line when writing the unmap file).
fn convert_vcf_record(
    view: &VcfRecordView,
    mapper: &CoordinateMapper,
    ref_genome: Option<&crate::core::fasta::FastaReader>,
    no_comp_allele: bool,
    out: &mut String,
) -> Result<(), String> {
    // Map the first position of REF allele (VCF is 1-based)
    let start = view.pos - 1; // Convert to 0-based
    let end = start + 1; // Map only the first position

    let result = mapper.map(view.chrom, start, end, Strand::Plus);

    match result {
        Some(segments) if segments.len() == 1 => {
            let seg = &segments[0];
            let target_chrom = &seg.target.chrom;
            let target_start = seg.target.start;
            let target_end = seg.target.end;
            let target_strand = seg.target.strand;

            // Get original fields
            let ref_allele = view.ref_allele().unwrap_or("N");
            let alt_alleles_str = view.alt_alleles().unwrap_or(".");

            // Calculate new REF position based on strand and variant type
            let (new_pos, new_ref) = if let Some(ref_reader) = ref_genome {
                // Get REF from target reference genome
                let ref_start = target_start;
                let ref_end = ref_start + 1;

                match ref_reader.fetch(target_chrom, ref_start, ref_end) {
                    // `fetch` already returns an uppercased sequence.
                    Some(seq) => (target_start + 1, seq),
                    None => return Err("Fail(KeyError)".to_string()),
                }
            } else {
                // No reference genome provided, keep original REF
                (target_start + 1, ref_allele.to_string())
            };

            if new_ref.is_empty() {
                return Err("Fail(KeyError)".to_string());
            }

            // Process ALT alleles (CrossMap logic)
            let mut alt_alleles_updated: Vec<String> = Vec::new();
            for alt_allele in alt_alleles_str.split(',') {
                if dna::is_dna(alt_allele) {
                    let updated = if ref_allele.len() != alt_allele.len() {
                        // Indel: replace the first nucleotide with the new REF,
                        // keeping (and, on the minus strand, reverse-complementing)
                        // the remainder of the allele.
                        let first_char = new_ref.chars().next().unwrap_or('N');
                        if alt_allele.len() > 1 {
                            if target_strand == Strand::Minus {
                                format!("{}{}", first_char, dna::revcomp(&alt_allele[1..]))
                            } else {
                                format!("{}{}", first_char, &alt_allele[1..])
                            }
                        } else {
                            first_char.to_string()
                        }
                    } else if target_strand == Strand::Minus {
                        // Substitution on a flipped block
                        dna::revcomp(alt_allele)
                    } else {
                        // Forward-strand substitution
                        alt_allele.to_string()
                    };

                    // Add to list (will filter REF==ALT later, matching CrossMap)
                    alt_alleles_updated.push(updated);
                } else {
                    // Non-DNA allele (e.g., <DEL>, <INS>), keep as-is
                    alt_alleles_updated.push(alt_allele.to_string());
                }
            }

            // Filter out ALT alleles that equal REF
            // (CrossMap: alt_alleles_updated = [i for i in alt_alleles_updated if i != ref_allele])
            alt_alleles_updated.retain(|alt| alt != &new_ref);

            // CrossMap behavior: when alt_alleles_updated is empty after filtering,
            // it sets fields[4] = "" (empty string), then checks if fields[3] != fields[4].
            // Since REF != "", the record is output with empty ALT.
            let alt_joined = alt_alleles_updated.join(",");
            if !no_comp_allele && alt_joined == new_ref {
                return Err("Fail(REF==ALT)".to_string());
            }

            // Format the output line in place, appending to the caller's
            // reusable buffer so a record costs no allocation of its own.
            write_output_line(
                out,
                view,
                target_chrom,
                new_pos,
                &new_ref,
                &alt_alleles_updated,
                target_end,
            );

            Ok(())
        }
        // Multiple mappings
        Some(_) => Err("Fail(Multiple_hits)".to_string()),
        // No mapping found
        None => Err("Fail(Unmap)".to_string()),
    }
}

/// Append `info` to `out`, replacing every `END=<digits>` value with `new_end`.
///
/// Mirrors CrossMap: `re.sub(r'END=\d+', 'END=' + str(target_end), fields[7])`.
fn write_info_with_end(out: &mut String, info: &str, new_end: u64) {
    use std::fmt::Write as _;

    let bytes = info.as_bytes();
    let mut last = 0usize;
    let mut cursor = 0usize;
    while cursor + 4 < bytes.len() {
        // Match "END=" followed by at least one digit.
        if &bytes[cursor..cursor + 4] == b"END=" && bytes[cursor + 4].is_ascii_digit() {
            out.push_str(&info[last..cursor]);
            out.push_str("END=");
            cursor += 4;
            while cursor < bytes.len() && bytes[cursor].is_ascii_digit() {
                cursor += 1;
            }
            let _ = write!(out, "{}", new_end);
            last = cursor;
        } else {
            cursor += 1;
        }
    }
    // `last` is always a char boundary (it follows ASCII "END=" or digits).
    out.push_str(&info[last..]);
}

/// Append a successfully mapped VCF record to `out`: every column, from CHROM
/// through the sample block.
///
/// Appends rather than returning a fresh `String` so a caller converting a whole
/// file keeps one buffer and reuses its capacity from record to record.
fn write_output_line(
    out: &mut String,
    view: &VcfRecordView,
    chrom: &str,
    pos: u64,
    ref_allele: &str,
    alt_alleles: &[String],
    target_end: u64,
) {
    use std::fmt::Write as _;

    // CHROM
    out.push_str(chrom);
    out.push('\t');
    // POS
    let _ = write!(out, "{}", pos);
    out.push('\t');
    // ID
    out.push_str(view.id().unwrap_or("."));
    out.push('\t');
    // REF
    out.push_str(ref_allele);
    out.push('\t');
    // ALT
    for (i, alt) in alt_alleles.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(alt);
    }
    out.push('\t');
    // QUAL
    out.push_str(view.qual().unwrap_or("."));
    out.push('\t');
    // FILTER
    out.push_str(view.filter().unwrap_or("."));
    out.push('\t');
    // INFO - update END if present (CrossMap behavior)
    write_info_with_end(out, view.info().unwrap_or("."), target_end);
    // FORMAT and sample columns, carried through verbatim as one slice.
    if let Some(tail) = view.format_and_samples() {
        out.push('\t');
        out.push_str(tail);
    }
}

/// Convert a VCF file using the coordinate mapper
///
/// # Arguments
/// * `input` - Input VCF file path
/// * `output` - Output VCF file path for successfully mapped records
/// * `unmap` - Output file path for unmapped records (will be output.unmap)
/// * `mapper` - Coordinate mapper with loaded chain index
/// * `ref_genome` - Optional path to target reference genome FASTA
/// * `no_comp_allele` - If true, keep variants where REF==ALT
/// * `threads` - Number of threads for parallel processing (1 = sequential)
/// 
/// # Returns
/// Conversion statistics
pub fn convert_vcf<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    ref_genome: Option<P>,
    no_comp_allele: bool,
    threads: usize,
) -> Result<ConversionStats, VcfParseError> {
    if threads > 1 {
        convert_vcf_parallel(input, output, mapper, ref_genome, no_comp_allele, threads)
    } else {
        convert_vcf_sequential(input, output, mapper, ref_genome, no_comp_allele)
    }
}

/// Sequential VCF conversion (single-threaded, line-by-line)
fn convert_vcf_sequential<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    ref_genome: Option<P>,
    no_comp_allele: bool,
) -> Result<ConversionStats, VcfParseError> {
    let input_file = std::fs::File::open(input.as_ref())?;
    let mut reader: Box<dyn BufRead> = if input.as_ref().extension().and_then(|e| e.to_str()) == Some("gz") {
        Box::new(BufReader::with_capacity(
            128 * 1024,
            flate2::read::MultiGzDecoder::new(input_file),
        ))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, input_file))
    };

    let output_path = output.as_ref();
    let unmap_path = output_path.with_extension("vcf.unmap");

    let mut output_file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(output_path)?);
    let mut unmap_file = BufWriter::with_capacity(64 * 1024, std::fs::File::create(&unmap_path)?);

    let assembly = ref_basename(&ref_genome);
    let ref_reader = ref_genome
        .map(|p| crate::core::fasta::FastaReader::open(p.as_ref()))
        .transpose()?;

    // Chromosome style of the input's `##contig` lines; CrossMap applies it to
    // every target contig it writes.
    let mut chr_template = "chr1";

    let mut stats = ConversionStats::default();
    let mut line_buf = String::with_capacity(4096);
    // Reusable scratch buffer for building output lines without per-line allocation.
    let mut scratch = String::with_capacity(4096);

    loop {
        line_buf.clear();
        let bytes_read = reader.read_line(&mut line_buf)?;
        if bytes_read == 0 {
            break;
        }

        let line = line_buf.trim_end();

        if line.is_empty() {
            continue;
        }

        if line.starts_with('#') {
            if line.starts_with("##fileformat")
                || line.starts_with("##INFO")
                || line.starts_with("##FILTER")
                || line.starts_with("##FORMAT")
                || line.starts_with("##ALT")
                || line.starts_with("##SAMPLE")
                || line.starts_with("##PEDIGREE")
            {
                writeln!(output_file, "{}", line)?;
                writeln!(unmap_file, "{}", line)?;
            } else if line.starts_with("##assembly") || line.starts_with("##contig") {
                if line.starts_with("##contig") {
                    chr_template = chr_template_from_contig_line(line);
                }
                writeln!(unmap_file, "{}", line)?;
            } else if line.starts_with("#CHROM") {
                write_contig_header(
                    &mut output_file,
                    ref_reader.as_ref(),
                    chr_template,
                    assembly.as_deref(),
                )?;
                writeln!(output_file, "##liftOverProgram=FastCrossMap")?;
                writeln!(output_file, "{}", line)?;
                writeln!(unmap_file, "{}", line)?;
            } else {
                writeln!(output_file, "{}", line)?;
            }
            continue;
        }

        stats.total += 1;

        match VcfRecordView::parse(line.as_bytes()) {
            Ok(view) => {
                scratch.clear();
                match convert_vcf_record(
                    &view,
                    mapper,
                    ref_reader.as_ref(),
                    no_comp_allele,
                    &mut scratch,
                ) {
                    Ok(()) => {
                        output_file.write_all(scratch.as_bytes())?;
                        output_file.write_all(b"\n")?;
                        stats.success += 1;
                    }
                    Err(reason) => {
                        write!(unmap_file, "{}\t{}\n", line, reason)?;
                        stats.failed += 1;
                    }
                }
            }
            Err(_) => {
                write!(unmap_file, "{}\tFail(ParseError)\n", line)?;
                stats.failed += 1;
            }
        }
    }

    Ok(stats)
}

/// Strip trailing ASCII whitespace from a line, matching the `str::trim_end()`
/// the sequential path applies to each line it reads.
#[inline]
fn trim_line(buf: &[u8]) -> &[u8] {
    let mut end = buf.len();
    while end > 0 && buf[end - 1].is_ascii_whitespace() {
        end -= 1;
    }
    &buf[..end]
}

/// Upper bound on the bytes of raw input held in one batch.
///
/// Record length varies by three orders of magnitude across inputs (a
/// sites-only VCF line is ~50 bytes; a 2504-sample line is ~10 KB), so the
/// batch is capped by *bytes* as well as by record count. Without this a
/// many-sample VCF would hold tens of megabytes per batch in flight.
const MAX_BATCH_BYTES: usize = 4 << 20;

/// Parallel VCF conversion.
///
/// Reading, converting and writing run as one overlapped pipeline (see
/// [`crate::core::pipeline`]): a reader thread parses and decompresses lines
/// into batches, the rayon pool converts a whole batch at a time, and the
/// calling thread writes results out in input order. The bounded queues keep
/// peak memory proportional to the batch size rather than to the input size.
///
/// The conversion stage reuses one [`RecordSink`] per rayon partition rather
/// than allocating a `VcfOut` per record: a mapped record's line is formatted
/// into the sink's scratch (cleared, not freed, between records) and appended to
/// the sink's payload in place, which is what keeps the parallel path's CPU cost
/// near the sequential path's instead of well above it.
fn convert_vcf_parallel<P: AsRef<Path>>(
    input: P,
    output: P,
    mapper: &CoordinateMapper,
    ref_genome: Option<P>,
    no_comp_allele: bool,
    threads: usize,
) -> Result<ConversionStats, VcfParseError> {
    use crate::core::pipeline::{
        batch_size, run_ordered_pipeline_chunked, Batch, LineBatch, OutBatch, RecordSink,
    };

    let input_file = std::fs::File::open(input.as_ref())?;
    let mut reader: Box<dyn BufRead + Send> =
        if input.as_ref().extension().and_then(|e| e.to_str()) == Some("gz") {
            Box::new(BufReader::with_capacity(
                128 * 1024,
                flate2::read::MultiGzDecoder::new(input_file),
            ))
        } else {
            Box::new(BufReader::with_capacity(128 * 1024, input_file))
        };

    let assembly = ref_basename(&ref_genome);
    let ref_reader = ref_genome
        .map(|p| crate::core::fasta::FastaReader::open(p.as_ref()))
        .transpose()?;

    let output_path = output.as_ref();
    let unmap_path = output_path.with_extension("vcf.unmap");
    let mut output_file = BufWriter::with_capacity(128 * 1024, std::fs::File::create(output_path)?);
    let mut unmap_file = BufWriter::with_capacity(64 * 1024, std::fs::File::create(&unmap_path)?);

    // Chromosome style of the input's `##contig` lines; CrossMap applies it to
    // every target contig it writes.
    let mut chr_template = "chr1";

    // Phase 1: write the headers, stopping at the first data line.
    let mut line_buf: Vec<u8> = Vec::with_capacity(4096);
    let mut pending: Option<Vec<u8>> = None;
    loop {
        line_buf.clear();
        if reader.read_until(b'\n', &mut line_buf)? == 0 {
            break;
        }
        let line = trim_line(&line_buf);
        if line.is_empty() {
            continue;
        }
        if line[0] != b'#' {
            pending = Some(line.to_vec());
            break;
        }

        let text = std::str::from_utf8(line)
            .map_err(|_| VcfParseError::InvalidUtf8("line"))?;
        if text.starts_with("##fileformat")
            || text.starts_with("##INFO")
            || text.starts_with("##FILTER")
            || text.starts_with("##FORMAT")
            || text.starts_with("##ALT")
            || text.starts_with("##SAMPLE")
            || text.starts_with("##PEDIGREE")
        {
            writeln!(output_file, "{}", text)?;
            writeln!(unmap_file, "{}", text)?;
        } else if text.starts_with("##assembly") || text.starts_with("##contig") {
            if text.starts_with("##contig") {
                chr_template = chr_template_from_contig_line(text);
            }
            writeln!(unmap_file, "{}", text)?;
        } else if text.starts_with("#CHROM") {
            write_contig_header(
                &mut output_file,
                ref_reader.as_ref(),
                chr_template,
                assembly.as_deref(),
            )?;
            writeln!(output_file, "##liftOverProgram=FastCrossMap")?;
            writeln!(output_file, "{}", text)?;
            writeln!(unmap_file, "{}", text)?;
        } else {
            writeln!(output_file, "{}", text)?;
        }
    }

    // Phase 2: stream the data lines through the pipeline.
    let mut stats = ConversionStats::default();
    let mut first = pending;

    run_ordered_pipeline_chunked(
        threads,
        |line: &[u8], sink: &mut RecordSink| {
            match VcfRecordView::parse(line) {
                Ok(view) => match convert_vcf_record(
                    &view,
                    mapper,
                    ref_reader.as_ref(),
                    no_comp_allele,
                    sink.text_mut(),
                ) {
                    Ok(()) => {
                        sink.flush_text();
                        true
                    }
                    Err(reason) => {
                        sink.write(reason.as_bytes());
                        false
                    }
                },
                Err(_) => {
                    sink.write(b"Fail(ParseError)");
                    false
                }
            }
        },
        // Fill one batch. The line is appended to the batch buffer directly —
        // no per-record allocation and no second copy of the line.
        |batch: &mut LineBatch| {
            if let Some(line) = first.take() {
                batch.push_line(&line);
            }
            while batch.len() < batch_size() && batch.buffer.len() < MAX_BATCH_BYTES {
                line_buf.clear();
                if reader.read_until(b'\n', &mut line_buf)? == 0 {
                    break;
                }
                let line = trim_line(&line_buf);
                if line.is_empty() {
                    continue;
                }
                batch.push_line(line);
            }
            Ok(batch.len())
        },
        |batch: &LineBatch, results: &OutBatch| {
            // A short result batch would silently drop records through `zip`,
            // so check the count rather than trusting it.
            if results.len() != batch.len() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!(
                        "internal error: {} converted records for {} input lines",
                        results.len(),
                        batch.len()
                    ),
                ));
            }
            for ((start, end), (payload, mapped)) in batch.lines.iter().zip(results.iter()) {
                if mapped {
                    output_file.write_all(payload)?;
                    output_file.write_all(b"\n")?;
                    stats.success += 1;
                } else {
                    // For a mapped-but-rejected record this matches the
                    // sequential path, which records the reason it reconstructed;
                    // a parse failure writes the same `Fail(ParseError)` marker.
                    unmap_file.write_all(&batch.buffer[*start..*end])?;
                    unmap_file.write_all(b"\t")?;
                    unmap_file.write_all(payload)?;
                    unmap_file.write_all(b"\n")?;
                    stats.failed += 1;
                }
            }
            stats.total += batch.len();
            Ok(())
        },
    )?;

    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_vcf_record_view_basic() {
        let line = b"chr1\t12345\trs123\tA\tG\t30\tPASS\tDP=100";
        let view = VcfRecordView::parse(line).unwrap();
        
        assert_eq!(view.chrom, "chr1");
        assert_eq!(view.pos, 12345);
        assert_eq!(view.id(), Some("rs123"));
        assert_eq!(view.ref_allele(), Some("A"));
        assert_eq!(view.alt_alleles(), Some("G"));
        assert_eq!(view.qual(), Some("30"));
        assert_eq!(view.filter(), Some("PASS"));
        assert_eq!(view.info(), Some("DP=100"));
    }
    
    #[test]
    fn test_vcf_record_view_with_samples() {
        let line = b"chr1\t12345\t.\tA\tG\t.\t.\t.\tGT:DP\t0/1:30\t1/1:25";
        let view = VcfRecordView::parse(line).unwrap();
        
        assert_eq!(view.chrom, "chr1");
        assert_eq!(view.pos, 12345);
        assert_eq!(view.format(), Some("GT:DP"));
        assert_eq!(view.samples(), vec!["0/1:30", "1/1:25"]);
    }
    
    #[test]
    fn test_vcf_record_view_too_few_fields() {
        let line = b"chr1\t12345\trs123";
        let result = VcfRecordView::parse(line);
        assert!(matches!(result, Err(VcfParseError::TooFewFields { .. })));
    }
    
    #[test]
    fn test_vcf_record_view_empty_line() {
        let line = b"";
        let result = VcfRecordView::parse(line);
        assert!(matches!(result, Err(VcfParseError::EmptyLine)));
    }
    
    #[test]
    fn test_variant_type_detection() {
        // Substitution
        let line = b"chr1\t100\t.\tA\tG\t.\t.\t.";
        let view = VcfRecordView::parse(line).unwrap();
        assert_eq!(view.variant_type(), VariantType::Substitution);
        
        // Insertion
        let line = b"chr1\t100\t.\tA\tAG\t.\t.\t.";
        let view = VcfRecordView::parse(line).unwrap();
        assert_eq!(view.variant_type(), VariantType::Insertion);
        
        // Deletion
        let line = b"chr1\t100\t.\tAG\tA\t.\t.\t.";
        let view = VcfRecordView::parse(line).unwrap();
        assert_eq!(view.variant_type(), VariantType::Deletion);
    }
    
    #[test]
    fn test_info_parsing() {
        let line = b"chr1\t100\t.\tA\tG\t.\t.\tDP=100;AF=0.5;DB";
        let view = VcfRecordView::parse(line).unwrap();
        let info = view.parse_info();
        
        assert_eq!(info.get("DP"), Some(&"100".to_string()));
        assert_eq!(info.get("AF"), Some(&"0.5".to_string()));
        assert_eq!(info.get("DB"), Some(&"".to_string())); // Flag
    }
    
    #[test]
    fn test_multi_allelic() {
        let line = b"chr1\t100\t.\tA\tG,T,C\t.\t.\t.";
        let view = VcfRecordView::parse(line).unwrap();
        
        assert_eq!(view.alt_alleles(), Some("G,T,C"));
    }
}
