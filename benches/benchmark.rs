//! Performance benchmarks for FastCrossMap
//!
//! Run with: cargo bench
//!
//! **Validates: Requirements 12.1, 12.2, 12.3, 12.4, 12.5**

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use fast_crossmap::core::{ChainIndex, CoordinateMapper, ChromStyle, Strand};
use std::path::PathBuf;

/// Chain file path for benchmarks
const CHAIN_FILE: &str = "ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz";

/// Check if chain file exists
fn chain_file_exists() -> bool {
    PathBuf::from(CHAIN_FILE).exists()
}

/// Benchmark chain file loading
fn bench_chain_loading(c: &mut Criterion) {
    if !chain_file_exists() {
        eprintln!("Skipping chain loading benchmark: chain file not found");
        return;
    }
    
    c.bench_function("chain_load_gz", |b| {
        b.iter(|| {
            let index = ChainIndex::from_chain_file(black_box(CHAIN_FILE)).unwrap();
            black_box(index)
        })
    });
}

/// Benchmark single coordinate mapping
fn bench_single_mapping(c: &mut Criterion) {
    if !chain_file_exists() {
        eprintln!("Skipping single mapping benchmark: chain file not found");
        return;
    }
    
    let index = ChainIndex::from_chain_file(CHAIN_FILE).unwrap();
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    c.bench_function("map_single_coord", |b| {
        b.iter(|| {
            let result = mapper.map(
                black_box("chr1"),
                black_box(1000000),
                black_box(1000100),
                black_box(Strand::Plus),
            );
            black_box(result)
        })
    });
}

/// Benchmark batch coordinate mapping
fn bench_batch_mapping(c: &mut Criterion) {
    if !chain_file_exists() {
        eprintln!("Skipping batch mapping benchmark: chain file not found");
        return;
    }
    
    let index = ChainIndex::from_chain_file(CHAIN_FILE).unwrap();
    let mapper = CoordinateMapper::new(index, ChromStyle::AsIs);
    
    // Generate test coordinates
    let coords: Vec<(u64, u64)> = (0..1000)
        .map(|i| (10000 + i * 1000, 10100 + i * 1000))
        .collect();
    
    let mut group = c.benchmark_group("batch_mapping");
    
    for size in [100, 500, 1000].iter() {
        group.throughput(Throughput::Elements(*size as u64));
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            b.iter(|| {
                for (start, end) in coords.iter().take(size) {
                    let result = mapper.map("chr1", *start, *end, Strand::Plus);
                    black_box(result);
                }
            })
        });
    }
    
    group.finish();
}

/// Benchmark interval query
fn bench_interval_query(c: &mut Criterion) {
    if !chain_file_exists() {
        eprintln!("Skipping interval query benchmark: chain file not found");
        return;
    }
    
    let index = ChainIndex::from_chain_file(CHAIN_FILE).unwrap();
    
    c.bench_function("interval_query", |b| {
        b.iter(|| {
            let result = index.query_intervals(
                black_box("chr1"),
                black_box(1000000),
                black_box(2000000),
            );
            black_box(result)
        })
    });
}

/// Benchmark chromosome name normalization
fn bench_chrom_normalization(c: &mut Criterion) {
    use fast_crossmap::core::normalize_chrom;
    
    let chroms = ["chr1", "1", "CHR1", "Chr1", "chrX", "X", "chrM", "MT"];
    
    c.bench_function("chrom_normalize", |b| {
        b.iter(|| {
            for chrom in &chroms {
                let result = normalize_chrom(black_box(chrom));
                black_box(result);
            }
        })
    });
}

/// Benchmark DNA reverse complement
fn bench_revcomp(c: &mut Criterion) {
    use fast_crossmap::core::dna::revcomp;
    
    let sequences = [
        "ACGT",
        "ACGTACGTACGT",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
    ];
    
    let mut group = c.benchmark_group("revcomp");
    
    for seq in &sequences {
        group.throughput(Throughput::Bytes(seq.len() as u64));
        group.bench_with_input(BenchmarkId::from_parameter(seq.len()), seq, |b, seq| {
            b.iter(|| {
                let result = revcomp(black_box(seq));
                black_box(result)
            })
        });
    }
    
    group.finish();
}

/// Benchmark BED line parsing
fn bench_bed_parsing(c: &mut Criterion) {
    use fast_crossmap::formats::BedRecordView;
    
    let lines = [
        b"chr1\t1000\t2000".as_slice(),
        b"chr1\t1000\t2000\tgene1\t500\t+".as_slice(),
        b"chr1\t1000\t2000\tgene1\t500\t+\t1100\t1900\t0,0,0\t2\t100,100\t0,900".as_slice(),
    ];
    
    let mut group = c.benchmark_group("bed_parsing");
    
    for (i, line) in lines.iter().enumerate() {
        let name = match i {
            0 => "BED3",
            1 => "BED6",
            2 => "BED12",
            _ => "unknown",
        };
        group.bench_with_input(BenchmarkId::from_parameter(name), line, |b, line| {
            b.iter(|| {
                let result = BedRecordView::parse(black_box(line));
                black_box(result)
            })
        });
    }
    
    group.finish();
}

/// Benchmark VCF line parsing
fn bench_vcf_parsing(c: &mut Criterion) {
    use fast_crossmap::formats::VcfRecordView;
    
    let line = b"chr1\t12345\trs123\tA\tG\t30\tPASS\tDP=100;AF=0.5\tGT:DP\t0/1:30";
    
    c.bench_function("vcf_parsing", |b| {
        b.iter(|| {
            let result = VcfRecordView::parse(black_box(line.as_slice()));
            black_box(result)
        })
    });
}

criterion_group!(
    benches,
    bench_chain_loading,
    bench_single_mapping,
    bench_batch_mapping,
    bench_interval_query,
    bench_chrom_normalization,
    bench_revcomp,
    bench_bed_parsing,
    bench_vcf_parsing,
);

criterion_main!(benches);
