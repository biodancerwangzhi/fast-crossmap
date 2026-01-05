# FastCrossMap

[![Release](https://img.shields.io/github/v/release/example/fast-crossmap)](https://github.com/example/fast-crossmap/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/example/fast-crossmap/workflows/Release/badge.svg)](https://github.com/example/fast-crossmap/actions)

**FastCrossMap** is a high-performance genome coordinate liftover tool written in Rust. It achieves **10-20x speedup** over Python CrossMap while maintaining **100% output compatibility**.

## Features

- 🚀 **10-20x faster** than CrossMap (single-threaded)
- 🧵 **Multi-threading support** with near-linear scalability
- 💾 **64x less memory** usage (16MB vs 1GB for BAM processing)
- ✅ **100% compatible** with CrossMap output (strict mode)
- 📦 **8 file formats** supported: BED, BAM/SAM/CRAM, VCF, GVCF, GFF/GTF, Wiggle, BigWig, MAF
- 🗜️ **Compressed file support**: .gz, .bz2 for both chain files and input files
- 🖥️ **Cross-platform**: Linux, macOS, Windows

## Performance

| Format | FastCrossMap (1T) | CrossMap | Speedup |
|--------|-------------------|----------|---------|
| BED (297K records) | 0.28s | 2.85s | **10x** |
| BAM (1.3GB) | 86s | 679s | **8x** |

| Tool | Peak Memory (BAM) |
|------|-------------------|
| FastCrossMap | 16 MB |
| CrossMap | 1,022 MB |

## Installation

### Pre-built Binaries

Download from [Releases](https://github.com/example/fast-crossmap/releases):

```bash
# Linux
wget https://github.com/example/fast-crossmap/releases/download/v0.1.0/fast-crossmap-linux-x64.tar.gz
tar -xzf fast-crossmap-linux-x64.tar.gz
./fast-crossmap --help

# macOS (Apple Silicon)
wget https://github.com/example/fast-crossmap/releases/download/v0.1.0/fast-crossmap-macos-arm64.tar.gz

# Windows
# Download fast-crossmap-windows-x64.zip from Releases
```

### Build from Source

```bash
# Requires Rust 1.70+
cargo build --release
./target/release/fast-crossmap --help
```

## Usage

### Basic Usage

```bash
# BED format
fast-crossmap bed hg19ToHg38.chain.gz input.bed output.bed

# BAM format
fast-crossmap bam hg19ToHg38.chain.gz input.bam output.bam

# VCF format
fast-crossmap vcf hg19ToHg38.chain.gz input.vcf output.vcf

# GFF/GTF format
fast-crossmap gff hg19ToHg38.chain.gz input.gff output.gff
```

### Multi-threading

```bash
# Use 4 threads
fast-crossmap bed -t 4 hg19ToHg38.chain.gz input.bed output.bed

# Use all available cores
fast-crossmap bed -t 0 hg19ToHg38.chain.gz input.bed output.bed
```

### Compatibility Mode

```bash
# Strict mode: 100% compatible with CrossMap output
fast-crossmap --compat-mode strict bed chain.gz input.bed output.bed

# Improved mode (default): optimized logic
fast-crossmap --compat-mode improved bed chain.gz input.bed output.bed
```

### Compressed Files

```bash
# Supports .gz and .bz2 for both chain and input files
fast-crossmap bed hg19ToHg38.chain.gz input.bed.gz output.bed
```

## Supported Formats

| Format | Description | Multi-threading |
|--------|-------------|-----------------|
| BED | BED3/BED6/BED12 | ✅ |
| BAM/SAM/CRAM | Alignment files | ✅ |
| VCF | Variant Call Format | ✅ |
| GVCF | Genomic VCF | ✅ |
| GFF/GTF | Gene annotations | ✅ |
| Wiggle | Coverage tracks | ✅ |
| BigWig | Binary Wiggle | ✅ |
| MAF | Multiple Alignment Format | ✅ |

## Comparison with Other Tools

| Feature | FastCrossMap | CrossMap | liftOver | FastRemap |
|---------|:------------:|:--------:|:--------:|:---------:|
| BED | ✅ | ✅ | ✅ | ✅ |
| BAM/SAM | ✅ | ✅ | ❌ | ✅ |
| VCF | ✅ | ✅ | ❌ | ❌ |
| GFF/GTF | ✅ | ✅ | ✅ | ❌ |
| Multi-threading | ✅ | ❌ | ❌ | ⚠️* |
| Compressed chain | ✅ | ✅ | ✅ | ❌ |
| Windows support | ✅ | ⚠️ | ❌ | ❌ |

*FastRemap uses internal 4 threads (not user-controllable)

## Algorithm

FastCrossMap uses the **BITS (Binary Indexed Tree Search)** algorithm via `rust-lapper` for O(log n + k) interval queries, combined with:

- **Zero-copy parsing**: Minimal memory allocation during file parsing
- **Parallel pipeline**: Rayon-based work-stealing scheduler
- **Memory-mapped I/O**: Efficient handling of large files
- **Cache-friendly data structures**: Contiguous memory layout for CPU cache efficiency

## Citation

If you use FastCrossMap in your research, please cite:

```bibtex
@software{fastcrossmap,
  title = {FastCrossMap: High-Performance Genome Coordinate Liftover},
  author = {FastCrossMap Contributors},
  year = {2026},
  url = {https://github.com/example/fast-crossmap}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- [CrossMap](https://crossmap.readthedocs.io/) - Original Python implementation
- [rust-lapper](https://docs.rs/rust-lapper/) - Interval tree library
- [Rayon](https://docs.rs/rayon/) - Parallel processing library
