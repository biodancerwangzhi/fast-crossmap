# Changelog

All notable changes to FastCrossMap will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.0] - 2026-10-08

### Added
- VCF/GVCF: gzip-compressed input (`.vcf.gz` / `.gvcf.gz`) is now read directly

### Changed
- **BAM/SAM backend switched from rust-htslib to noodles**: BAM/SAM is now pure Rust,
  so users no longer need a C toolchain (gcc / cmake / zlib-dev) to build
- `-t/--threads` help text now states the actual default; it has always been 1
- VCF/GVCF parallel path rewritten as an ordered pipeline (read / convert / write
  overlap, output order preserved)
- Windows binaries now include BAM/SAM: the pure-Rust backend dropped the htslib
  dependency that used to make them unbuildable there
- Documentation no longer claims CRAM support; only BAM and SAM are implemented

### Fixed
- `##contig` header order now follows the reference `.fai` index instead of being
  dependent on hash iteration order, so output is byte-stable across runs
- Windows release build produced no executable: the `cli` feature became required by
  the binary, and the Windows job disabled it along with the default features

## [0.4.0] - 2026-01-11

### Added
- Wiki documentation (Installation, QuickStart, Advanced, FAQ, Contributing)
- VCF: target reference genome is now a required positional argument and REF alleles
  are updated at the destination coordinate (`--no-comp-allele` to keep REF==ALT variants)
- BAM/SAM/CRAM support
- `tokio` runtime for BigWig writing

### Changed
- GVCF, MAF, Wiggle and VCF conversion paths reworked
- CLI rework in `src/main.rs` (positional reference genome, dropped `-c/--compress`)

## [0.3.0] - 2026-01-08

### Added
- BED: interval handling reworked
- Wiggle output
- `LICENSE`

### Removed
- `.devcontainer/cross-compile/` (cross-compilation is handled by the release workflow)

## [0.2.0] - 2026-01-06

### Changed
- README: expanded installation and usage documentation

## [0.1.1] - 2026-01-05

### Fixed
- Windows build and release workflow

### Changed
- BAM/SAM/CRAM support moved behind a feature flag, so Windows builds exclude it

## [0.1.0] - 2026-01-05

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

[Unreleased]: https://github.com/biodancerwangzhi/fast-crossmap/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.5.0
[0.4.0]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.4.0
[0.3.0]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.3.0
[0.2.0]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.2.0
[0.1.1]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.1.1
[0.1.0]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.1.0
