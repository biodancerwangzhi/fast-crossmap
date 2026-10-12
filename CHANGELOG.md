# Changelog

All notable changes to FastCrossMap will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.7.0] - 2026-10-12

### Added
- `-t` / `--threads` for the `maf` and `bigwig` subcommands, so every format
  can now be parallelised from the CLI
- `-t 0` now means "use every available core" on all subcommands, instead of
  being rejected as an invalid thread count

### Changed
- MAF: records are parsed and mapped in parallel batches and then flushed in
  input order, so the bytes written are identical to the single-threaded run
- Wiggle / BigWig: the parallel batch size is raised to 262,144 records, which
  removes the slowdown at `-t 2` where a worker pool was woken for every small
  batch; the pool construction is now shared with the other format pipelines
## [0.6.0] - 2026-10-11

### Added
- BAM/SAM: paired-end liftover. Mate coordinates (`RNEXT` / `PNEXT` / `TLEN`)
  are rewritten against the destination, the paired FLAG bits are rebuilt, and
  the nine paired-end categories CrossMap reports are tallied and printed, so
  the counts can be reconciled against CrossMap's own run for the same input

### Changed
- Chain chromosome names are interned into a shared `Arc<str>` instead of one
  `String` per block: a block carries two pointers rather than two allocations,
  and a third copy that used to be cloned into the index is gone. Cost drops
  from **241 to 145 bytes per chain block** (flat across 54k-432k blocks); a BED
  run's peak RSS falls from 16.7 MB to 11.7 MB
- Wiggle / BigWig: records are now spooled to a scratch file and merged one
  chromosome at a time, so peak memory tracks the largest single chromosome
  instead of the whole input, which is what put these two formats above CrossMap
  on memory. Output is byte-for-byte unchanged

### Fixed
- BAM/SAM: `@HD` is now written as `VN:1.0 SO:coordinate`, as CrossMap does,
  instead of noodles' default `VN:1.6` with no sort order

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
- `##contig` lines are now emitted in sorted contig-name order, as CrossMap does,
  instead of in hash iteration order, so output is byte-stable across runs
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

[Unreleased]: https://github.com/biodancerwangzhi/fast-crossmap/compare/v0.6.0...HEAD
[0.6.0]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.6.0
[0.5.0]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.5.0
[0.4.0]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.4.0
[0.3.0]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.3.0
[0.2.0]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.2.0
[0.1.1]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.1.1
[0.1.0]: https://github.com/biodancerwangzhi/fast-crossmap/releases/tag/v0.1.0
