# Changelog

All notable changes to FastCrossMap will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-01-06

### Added
- Initial release of FastCrossMap
- Support for 8 file formats: BED, BAM/SAM/CRAM, VCF, GVCF, GFF/GTF, Wiggle, BigWig, MAF
- Multi-threading support with `-t` option
- Compressed file support (.gz, .bz2) for both chain files and input files
- Two compatibility modes: `strict` (100% CrossMap compatible) and `improved` (optimized)
- Cross-platform support: Linux, macOS, Windows (Windows without BAM support)
- Property-based testing suite
- Benchmark scripts for performance comparison

### Performance
- 10-20x faster than CrossMap (single-threaded)
- 64x less memory usage for BAM processing
- Near-linear multi-threading scalability

[Unreleased]: https://github.com/biodancerwangzhi/fast-crossmap/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.1.0
