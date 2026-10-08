//! Ordered conversion pipeline
//!
//! A three-stage pipeline — **read → convert → write** — with all three stages
//! running concurrently:
//!
//! ```text
//! thread A (reader)     read + decompress, fill a batch   ──tx1──▶
//! thread B (converter)  convert one batch across the pool ──tx2──▶
//! calling thread        write batches in order
//! ```
//!
//! Output order always matches input order. The conversion thread is the only
//! producer of the channel the writer consumes and it handles one batch at a
//! time (farming that batch out across the rayon pool), so batches arrive in
//! submission order and no reordering buffer is needed.
//!
//! The bounded channels apply backpressure, so peak memory stays proportional
//! to `batch_size × threads` rather than to the size of the input.

use std::io;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::thread;
use std::time::{Duration, Instant};

use rayon::prelude::*;

/// Default number of records converted as one unit.
pub const DEFAULT_BATCH_SIZE: usize = 4096;

/// Number of records converted as one unit.
///
/// `FCM_BATCH` overrides the format's own default, for measuring how batch size
/// interacts with the worker count: a batch that is too small spends its time in
/// the channel hand-off and rayon's split/collect rather than in conversion.
pub fn batch_size_or(default: usize) -> usize {
    std::env::var("FCM_BATCH")
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|&n| n > 0)
        .unwrap_or(default)
}

/// Batch size for the line-oriented formats.
pub fn batch_size() -> usize {
    batch_size_or(DEFAULT_BATCH_SIZE)
}

/// How many batches may be queued between two stages.
///
/// Each stage is allowed to run this far ahead of the next one, which is what
/// lets reading, converting and writing overlap without letting memory grow
/// with the input.
///
/// `FCM_QUEUE_DEPTH` overrides it: a deeper queue lets a fast stage run further
/// ahead of a slow one, at the cost of holding more batches in memory.
fn queue_depth(threads: usize) -> usize {
    std::env::var("FCM_QUEUE_DEPTH")
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|&n| n >= 2)
        .unwrap_or((threads * 2).max(2))
}

/// The result of converting one record.
#[derive(Debug, Clone)]
pub struct Conversion<Out> {
    /// The converted record.
    pub output: Out,
    /// Whether the record was mapped. Unmapped records are written to the
    /// unmap file by the caller, using the original input.
    pub mapped: bool,
}

/// A group of input records converted together.
///
/// Implementors own their storage. The only thing that varies between formats
/// is how the records are laid out and therefore how they are iterated in
/// parallel, which is what [`convert_parallel`](Self::convert_parallel) exposes.
pub trait Batch: Send + Default {
    /// The record type the conversion function sees.
    type Item: Sync + ?Sized;

    /// Number of records in the batch.
    fn len(&self) -> usize;

    /// Whether the batch holds no records.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Convert every record in the batch, returning one result per record in
    /// input order.
    fn convert_parallel<Out, F>(&self, convert: F, pool: &rayon::ThreadPool) -> Vec<Conversion<Out>>
    where
        Out: Default + Send,
        F: Fn(&Self::Item, &mut Out) -> bool + Sync;
}

/// Convert a batch of owned items, allocating one result per record.
pub struct ItemBatch<T> {
    pub items: Vec<T>,
}

impl<T> Default for ItemBatch<T> {
    fn default() -> Self {
        Self { items: Vec::new() }
    }
}

impl<T: Send + Sync> Batch for ItemBatch<T> {
    type Item = T;

    fn len(&self) -> usize {
        self.items.len()
    }

    fn convert_parallel<Out, F>(&self, convert: F, pool: &rayon::ThreadPool) -> Vec<Conversion<Out>>
    where
        Out: Default + Send,
        F: Fn(&T, &mut Out) -> bool + Sync,
    {
        pool.install(|| {
            self.items
                .par_iter()
                .map(|item| {
                    let mut output = Out::default();
                    let mapped = convert(item, &mut output);
                    Conversion { output, mapped }
                })
                .collect()
        })
    }
}

/// A batch of lines read from a text input.
///
/// Lines are stored back to back in `buffer`; `lines[i]` is the byte range of
/// line `i`. Keeping one allocation per batch instead of one per line is what
/// makes the line-oriented formats fast — a batch of 4096 ten-kilobyte VCF
/// records costs two allocations, not 4096.
pub struct LineBatch {
    /// Concatenated line bytes.
    pub buffer: Vec<u8>,
    /// `(start, end)` byte range of each line within `buffer`.
    pub lines: Vec<(usize, usize)>,
}

impl LineBatch {
    /// Number of bytes of line storage reserved up front for a new batch.
    ///
    /// Worth a large-ish reservation because VCF/GFF records are long; the
    /// vector grows if the batch's lines need more, and stays bounded by
    /// `batch_size × longest line`.
    const INITIAL_BUFFER_CAPACITY: usize = 1 << 20;

    /// Append one line to the batch, recording its byte range.
    ///
    /// The line is taken verbatim, including any trailing newline; conversion
    /// functions receive it as-is, matching what a `read_until` loop yields.
    pub fn push_line(&mut self, line: &[u8]) {
        let start = self.buffer.len();
        self.buffer.extend_from_slice(line);
        self.lines.push((start, self.buffer.len()));
    }

    /// Convert every line, reusing one output buffer per rayon partition.
    ///
    /// [`Batch::convert_parallel`] hands back an owned `Out` per record, and for
    /// the line formats `Out` owns a `String` that starts empty — so every record
    /// pays a fresh growth sequence (8, 16, 32 … bytes) plus the copies that go
    /// with it, and then a free. On a 1.1-million-record VCF that is the whole
    /// reason the parallel path costs ~2× the CPU of the sequential one rather
    /// than the same.
    ///
    /// Here each partition keeps one [`RecordSink`] for the whole batch: its text
    /// scratch is cleared (not freed) between records and its byte buffer is
    /// appended to in place, so after the first few records a record costs no
    /// allocation at all.
    pub fn convert_parallel_chunked<F>(&self, convert: F, pool: &rayon::ThreadPool) -> OutBatch
    where
        F: Fn(&[u8], &mut RecordSink) -> bool + Sync,
    {
        let buffer = &self.buffer;
        let parts: Vec<RecordSink> = pool.install(|| {
            self.lines
                .par_iter()
                .fold(RecordSink::default, |mut sink, &(start, end)| {
                    sink.begin();
                    let mapped = convert(&buffer[start..end], &mut sink);
                    sink.end(mapped);
                    sink
                })
                .collect()
        });
        OutBatch { parts }
    }
}

impl Default for LineBatch {
    fn default() -> Self {
        Self {
            buffer: Vec::with_capacity(Self::INITIAL_BUFFER_CAPACITY),
            lines: Vec::new(),
        }
    }
}

impl Batch for LineBatch {
    type Item = [u8];

    fn len(&self) -> usize {
        self.lines.len()
    }

    fn convert_parallel<Out, F>(&self, convert: F, pool: &rayon::ThreadPool) -> Vec<Conversion<Out>>
    where
        Out: Default + Send,
        F: Fn(&[u8], &mut Out) -> bool + Sync,
    {
        let buffer = &self.buffer;
        pool.install(|| {
            self.lines
                .par_iter()
                .map(|&(start, end)| {
                    let mut output = Out::default();
                    let mapped = convert(&buffer[start..end], &mut output);
                    Conversion { output, mapped }
                })
                .collect()
        })
    }
}

/// Per-partition sink for the records of one batch.
///
/// Two ways to fill a record's payload, matching the two shapes the line formats
/// need: `text_mut` for the ones that format a line (the expensive case, and the
/// one this type exists to make cheap) and `write` for the ones that pass bytes
/// through unchanged. The text scratch is appended to the payload by
/// [`flush_text`](Self::flush_text), so the two may be combined in whatever order
/// the format needs.
///
/// A record's payload is one flat, append-only stretch of `buffer` with no
/// per-record allocation anywhere: the scratch holds its capacity from record to
/// record, and the payload grows into the same `Vec` for the whole batch.
#[derive(Default)]
pub struct RecordSink {
    /// Concatenated payload of the records committed so far.
    buffer: Vec<u8>,
    /// Scratch for the record being built; cleared — not freed — by `begin`, so
    /// its capacity survives from record to record.
    text: String,
    /// Where the record being built starts within `buffer`.
    start: usize,
    /// `(start, end)` of every committed payload within `buffer`.
    ranges: Vec<(usize, usize)>,
    /// Whether each committed record mapped, parallel to `ranges`.
    mapped: Vec<bool>,
}

impl RecordSink {
    /// Begin a new record's payload.
    #[inline]
    pub fn begin(&mut self) {
        self.start = self.buffer.len();
        self.text.clear();
    }

    /// The reusable text scratch for the record being built.
    ///
    /// Its contents are appended to the payload by
    /// [`flush_text`](Self::flush_text), which must be called before
    /// [`end`](Self::end) if the scratch was used.
    #[inline]
    pub fn text_mut(&mut self) -> &mut String {
        &mut self.text
    }

    /// Append the text scratch to the current record's payload.
    #[inline]
    pub fn flush_text(&mut self) {
        self.buffer.extend_from_slice(self.text.as_bytes());
    }

    /// Append bytes to the current record's payload.
    #[inline]
    pub fn write(&mut self, bytes: &[u8]) {
        self.buffer.extend_from_slice(bytes);
    }

    /// Append one byte to the current record's payload.
    #[inline]
    pub fn push(&mut self, byte: u8) {
        self.buffer.push(byte);
    }

    /// Finish the current record, recording whether it mapped.
    #[inline]
    pub fn end(&mut self, mapped: bool) {
        self.ranges.push((self.start, self.buffer.len()));
        self.mapped.push(mapped);
    }
}

/// One batch's converted records, laid out the way the input is: a handful of
/// flat buffers plus per-record byte ranges, instead of an owned value per
/// record.
///
/// The partitions are kept as they came back rather than concatenated into one
/// buffer, which saves a second copy of the whole batch's output. They are in
/// input order, so [`iter`](Self::iter) yields the records in input order.
pub struct OutBatch {
    parts: Vec<RecordSink>,
}

impl OutBatch {
    /// Number of records in the batch.
    pub fn len(&self) -> usize {
        self.parts.iter().map(|p| p.ranges.len()).sum()
    }

    /// Whether the batch holds no records.
    pub fn is_empty(&self) -> bool {
        self.parts.iter().all(|p| p.ranges.is_empty())
    }

    /// Every record's payload in input order, paired with whether it mapped.
    pub fn iter(&self) -> impl Iterator<Item = (&[u8], bool)> + '_ {
        self.parts.iter().flat_map(|part| {
            part.ranges
                .iter()
                .zip(part.mapped.iter())
                .map(move |(&(start, end), &mapped)| (&part.buffer[start..end], mapped))
        })
    }
}

/// Per-stage wall-clock totals, for diagnosing where a pipeline spends its time.
///
/// Enabled by setting `FCM_PIPELINE_STATS`; when off, each stage's contribution
/// is a single relaxed load. The three totals overlap in real time (the stages
/// run concurrently), so compare them against each other and against the total
/// runtime, not as a breakdown that sums to it.
#[derive(Default)]
struct StageTiming {
    enabled: bool,
    read: AtomicU64,
    convert: AtomicU64,
    write: AtomicU64,
    batches: AtomicU64,
}

impl StageTiming {
    fn new() -> Self {
        Self {
            enabled: std::env::var_os("FCM_PIPELINE_STATS").is_some(),
            ..Default::default()
        }
    }

    /// Record `elapsed` against one stage. A no-op unless stats are enabled.
    fn add(&self, stage: &AtomicU64, elapsed: Duration) {
        if self.enabled {
            stage.fetch_add(elapsed.as_nanos() as u64, Ordering::Relaxed);
        }
    }

    fn report(&self, threads: usize) {
        if !self.enabled {
            return;
        }
        let ms = |a: &AtomicU64| a.load(Ordering::Relaxed) as f64 / 1e6;
        let total = ms(&self.read) + ms(&self.convert) + ms(&self.write);
        eprintln!(
            "[pipeline t={}] batches={} read={:.0}ms convert={:.0}ms write={:.0}ms",
            threads,
            self.batches.load(Ordering::Relaxed),
            ms(&self.read),
            ms(&self.convert),
            ms(&self.write),
        );
        if total > 0.0 {
            eprintln!(
                "[pipeline t={}] share: read={:.0}% convert={:.0}% write={:.0}%",
                threads,
                100.0 * ms(&self.read) / total,
                100.0 * ms(&self.convert) / total,
                100.0 * ms(&self.write) / total,
            );
        }
    }
}

/// Number of threads in the conversion pool.
///
/// `FCM_CONV_CAP` overrides it, for measuring how much the conversion stage
/// actually benefits from extra workers.
fn convert_workers(threads: usize) -> usize {
    match std::env::var("FCM_CONV_CAP").ok().and_then(|v| v.parse().ok()) {
        Some(cap) => threads.min(cap).max(1),
        None => threads,
    }
}

/// Run a conversion as an overlapped, order-preserving pipeline.
///
/// * `threads` — worker threads for the conversion stage (must be > 1; callers
///   are expected to keep a separate sequential path for a single thread).
/// * `convert` — converts one record. Returns whether it was mapped.
/// * `read` — appends up to one batch of records to the given (initially empty)
///   batch and returns how many it appended. Returning 0 means end of input.
/// * `write` — writes one batch and its per-record results, in input order.
///
/// A `read` or `write` error is propagated; the pipeline shuts down cleanly
/// either way.
///
/// The stages run on scoped threads, so `convert` may borrow (for instance the
/// mapper and the reference reader) instead of having to own `'static` clones.
pub fn run_ordered_pipeline<B, Out, F, R, W>(
    threads: usize,
    convert: F,
    mut read: R,
    mut write: W,
) -> io::Result<()>
where
    B: Batch + Send,
    Out: Default + Send,
    F: Fn(&B::Item, &mut Out) -> bool + Sync + Send,
    R: FnMut(&mut B) -> io::Result<usize> + Send,
    W: FnMut(&B, &[Conversion<Out>]) -> io::Result<()>,
{
    let depth = queue_depth(threads);
    // Borrowed (not moved) so both stages can share it.
    let timing = &StageTiming::new();

    let (batch_tx, batch_rx): (SyncSender<B>, Receiver<B>) = sync_channel(depth);
    let (result_tx, result_rx) = sync_channel(depth);

    thread::scope(|scope| -> io::Result<()> {
        // Stage 1: read and decompress.
        let reader = scope.spawn(move || -> io::Result<()> {
            loop {
                let mut batch = B::default();
                let t0 = Instant::now();
                let n = read(&mut batch)?;
                timing.add(&timing.read, t0.elapsed());
                if n == 0 {
                    break;
                }
                // The converter has stopped (writer failed); nothing left to do.
                if batch_tx.send(batch).is_err() {
                    break;
                }
            }
            Ok(())
        });

        // Stage 2: convert each batch across the worker pool.
        let converter = scope.spawn(move || {
            let pool = match rayon::ThreadPoolBuilder::new()
                .num_threads(convert_workers(threads))
                .build()
            {
                Ok(pool) => pool,
                Err(_) => return,
            };
            while let Ok(batch) = batch_rx.recv() {
                let t0 = Instant::now();
                let results = batch.convert_parallel(&convert, &pool);
                timing.add(&timing.convert, t0.elapsed());
                timing.batches.fetch_add(1, Ordering::Relaxed);
                if result_tx.send((batch, results)).is_err() {
                    break;
                }
            }
        });

        // Stage 3: write, on the calling thread.
        let mut write_result = Ok(());
        while let Ok((batch, results)) = result_rx.recv() {
            let t0 = Instant::now();
            let r = write(&batch, &results);
            timing.add(&timing.write, t0.elapsed());
            if let Err(e) = r {
                write_result = Err(e);
                break;
            }
        }
        // Dropping the receiver unblocks the converter, which in turn unblocks
        // the reader, so both threads are guaranteed to finish.
        drop(result_rx);

        let reader_result = reader.join().unwrap_or_else(|_| {
            Err(io::Error::new(io::ErrorKind::Other, "reader thread panicked"))
        });
        let _ = converter.join();

        timing.report(threads);

        write_result?;
        reader_result
    })
}

/// Run a line-oriented conversion as the same overlapped pipeline, but with a
/// reusable [`RecordSink`] per rayon partition instead of an owned `Out` per
/// record.
///
/// Same stages, same ordering guarantee and same shutdown as
/// [`run_ordered_pipeline`]; only the conversion stage's allocation behaviour
/// differs.
///
/// Takes `&[u8]` lines and hands the sink's payload back as bytes, so the writer
/// needs no knowledge of how the format built the record.
pub fn run_ordered_pipeline_chunked<F, R, W>(
    threads: usize,
    convert: F,
    mut read: R,
    mut write: W,
) -> io::Result<()>
where
    F: Fn(&[u8], &mut RecordSink) -> bool + Sync + Send,
    R: FnMut(&mut LineBatch) -> io::Result<usize> + Send,
    W: FnMut(&LineBatch, &OutBatch) -> io::Result<()>,
{
    let depth = queue_depth(threads);
    let timing = &StageTiming::new();

    let (batch_tx, batch_rx): (SyncSender<LineBatch>, Receiver<LineBatch>) = sync_channel(depth);
    let (result_tx, result_rx) = sync_channel(depth);

    thread::scope(|scope| -> io::Result<()> {
        let reader = scope.spawn(move || -> io::Result<()> {
            loop {
                let mut batch = LineBatch::default();
                let t0 = Instant::now();
                let n = read(&mut batch)?;
                timing.add(&timing.read, t0.elapsed());
                if n == 0 {
                    break;
                }
                if batch_tx.send(batch).is_err() {
                    break;
                }
            }
            Ok(())
        });

        let converter = scope.spawn(move || {
            let pool = match rayon::ThreadPoolBuilder::new()
                .num_threads(convert_workers(threads))
                .build()
            {
                Ok(pool) => pool,
                Err(_) => return,
            };
            while let Ok(batch) = batch_rx.recv() {
                let t0 = Instant::now();
                let results = batch.convert_parallel_chunked(&convert, &pool);
                timing.add(&timing.convert, t0.elapsed());
                timing.batches.fetch_add(1, Ordering::Relaxed);
                if result_tx.send((batch, results)).is_err() {
                    break;
                }
            }
        });

        let mut write_result = Ok(());
        while let Ok((batch, results)) = result_rx.recv() {
            let t0 = Instant::now();
            let r = write(&batch, &results);
            timing.add(&timing.write, t0.elapsed());
            if let Err(e) = r {
                write_result = Err(e);
                break;
            }
        }
        drop(result_rx);

        let reader_result = reader.join().unwrap_or_else(|_| {
            Err(io::Error::new(io::ErrorKind::Other, "reader thread panicked"))
        });
        let _ = converter.join();

        timing.report(threads);

        write_result?;
        reader_result
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fmt::Write as _;
    use std::sync::Mutex;

    /// Records the order in which results were written.
    #[derive(Default)]
    struct Recorder {
        written: Mutex<Vec<usize>>,
    }

    impl Recorder {
        fn order(&self) -> Vec<usize> {
            self.written.lock().unwrap().clone()
        }
    }

    #[test]
    fn line_batch_parallel_conversion_preserves_order() {
        let mut batch = LineBatch::default();
        for i in 0..1000usize {
            batch.push_line(format!("{}\n", i).as_bytes());
        }
        assert_eq!(batch.len(), 1000);

        let pool = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        let results =
            batch.convert_parallel(|line, out: &mut String| {
                let text = std::str::from_utf8(line).unwrap().trim_end();
                out.push_str(text);
                text.parse::<usize>().unwrap() % 2 == 0
            }, &pool);

        assert_eq!(results.len(), 1000);
        for (i, r) in results.iter().enumerate() {
            assert_eq!(r.output, i.to_string());
            assert_eq!(r.mapped, i % 2 == 0);
        }
    }

    #[test]
    fn item_batch_parallel_conversion_preserves_order() {
        let batch = ItemBatch { items: (0..500u64).collect() };
        let pool = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        let results = batch.convert_parallel(|n: &u64, out: &mut String| {
            write!(out, "{}", n * 2).unwrap();
            *n < 250
        }, &pool);

        for (i, r) in results.iter().enumerate() {
            assert_eq!(r.output, (i as u64 * 2).to_string());
            assert_eq!(r.mapped, i < 250);
        }
    }

    #[test]
    fn pipeline_writes_batches_in_input_order() {
        const TOTAL: usize = 20_000;
        const BATCH: usize = 512;

        let recorder = Recorder::default();
        let mut next = 0usize;

        run_ordered_pipeline::<LineBatch, String, _, _, _>(
            4,
            |line, out: &mut String| {
                let text = std::str::from_utf8(line).unwrap().trim_end();
                out.push_str(text);
                true
            },
            move |batch: &mut LineBatch| {
                while batch.len() < BATCH && next < TOTAL {
                    batch.push_line(format!("{}\n", next).as_bytes());
                    next += 1;
                }
                Ok(batch.len())
            },
            |batch, results| {
                let mut written = recorder.written.lock().unwrap();
                assert_eq!(batch.len(), results.len());
                for r in results {
                    written.push(r.output.parse().unwrap());
                }
                Ok(())
            },
        )
        .unwrap();

        let order = recorder.order();
        assert_eq!(order.len(), TOTAL);
        assert_eq!(order, (0..TOTAL).collect::<Vec<_>>());
    }

    #[test]
    fn chunked_pipeline_writes_batches_in_input_order_and_keeps_payloads() {
        const TOTAL: usize = 20_000;
        const BATCH: usize = 512;

        let recorder = Recorder::default();
        let mut next = 0usize;

        run_ordered_pipeline_chunked(
            4,
            |line, sink: &mut RecordSink| {
                // Each record writes through both routes: text scratch plus a
                // pass-through byte, so `flush_text` ordering is exercised too.
                let text = std::str::from_utf8(line).unwrap().trim_end();
                sink.text_mut().push_str(text);
                sink.flush_text();
                sink.push(b':');
                // Deliberately leave `text` dirty: `begin` must clear it.
                sink.text_mut().push_str("LEFTOVER");
                let n: usize = text.parse().unwrap();
                n % 2 == 0
            },
            move |batch: &mut LineBatch| {
                while batch.len() < BATCH && next < TOTAL {
                    batch.push_line(format!("{}\n", next).as_bytes());
                    next += 1;
                }
                Ok(batch.len())
            },
            |batch, results| {
                let mut written = recorder.written.lock().unwrap();
                assert_eq!(batch.len(), results.len());
                for (payload, mapped) in results.iter() {
                    let text = std::str::from_utf8(payload).unwrap();
                    let (number, rest) = text.split_once(':').unwrap();
                    // The record's payload is exactly the scratch that was
                    // flushed: the "LEFTOVER" left in `text` afterwards must not
                    // appear, because `begin` clears the scratch each record.
                    assert_eq!(rest, "");
                    assert!(!text.contains("LEFTOVER"), "text scratch leaked between records");
                    written.push(number.parse().unwrap());
                    assert_eq!(mapped, number.parse::<usize>().unwrap() % 2 == 0);
                }
                Ok(())
            },
        )
        .unwrap();

        let order = recorder.order();
        assert_eq!(order.len(), TOTAL);
        assert_eq!(order, (0..TOTAL).collect::<Vec<_>>());
    }

    #[test]
    fn pipeline_propagates_write_error_and_shuts_down() {
        let err = run_ordered_pipeline::<LineBatch, String, _, _, _>(
            4,
            |line, out: &mut String| {
                out.push_str(std::str::from_utf8(line).unwrap().trim_end());
                true
            },
            |batch: &mut LineBatch| {
                while batch.len() < 128 {
                    batch.push_line(b"x\n");
                }
                Ok(batch.len())
            },
            |_, _| Err(io::Error::new(io::ErrorKind::Other, "boom")),
        )
        .unwrap_err();
        assert_eq!(err.to_string(), "boom");
    }
}
