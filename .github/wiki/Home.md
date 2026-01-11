# FastCrossMap Wiki

Welcome to **FastCrossMap** - a high-performance genome coordinate liftover tool written in Rust.

## 🚀 Key Features

- **10-20x faster** than CrossMap (single-threaded)
- **Multi-threading support** with near-linear scalability
- **64x less memory** usage (16MB vs 1GB for BAM processing)
- **100% compatible** with CrossMap output (strict mode)
- **8 file formats**: BED, BAM/SAM/CRAM, VCF, GVCF, GFF/GTF, Wiggle, BigWig, MAF
- **Compressed file support**: .gz, .bz2 for both chain files and input files
- **Cross-platform**: Linux, macOS, Windows

## 📚 Documentation

| Page | Description |
|------|-------------|
| [Installation](Installation) | Install FastCrossMap + download chain files |
| [Quick Start](QuickStart) | Get started in 5 minutes |
| [Advanced Usage](Advanced) | Detailed format guides, CLI reference, multi-threading, performance |
| [FAQ](FAQ) | Common questions and troubleshooting |
| [Contributing](Contributing) | How to contribute + development setup |

## 📊 Performance at a Glance

| Format | FastCrossMap (1T) | CrossMap | Speedup |
|--------|-------------------|----------|---------|
| BED (297K records) | 0.28s | 2.85s | **10x** |
| BAM (1.3GB) | 86s | 679s | **8x** |

| Tool | Peak Memory (BAM) |
|------|-------------------|
| FastCrossMap | 16 MB |
| CrossMap | 1,022 MB |

## 🔗 Links

- **GitHub**: https://github.com/biodancerwangzhi/fast-crossmap
- **Issues**: https://github.com/biodancerwangzhi/fast-crossmap/issues
- **Releases**: https://github.com/biodancerwangzhi/fast-crossmap/releases
- **Contact**: szxszx@foxmail.com

## 📖 Citation

```bibtex
@software{fastcrossmap,
  title = {FastCrossMap: High-Performance Genome Coordinate Liftover},
  author = {Wang, Zhi},
  year = {2026},
  url = {https://github.com/biodancerwangzhi/fast-crossmap}
}
```

---
**Version**: v0.4.0 | **Author**: Wang Zhi | **License**: MIT
