# FAQ & Troubleshooting

## General Questions

**Q: What is FastCrossMap?**  
A: A high-performance genome coordinate liftover tool written in Rust.

**Q: How is it different from CrossMap?**  
A: 10-20x faster, 64x less memory, multi-threading support, same output compatibility.

**Q: Which formats are supported?**  
A: BED, BAM/SAM/CRAM, VCF, GVCF, GFF/GTF, Wiggle, BigWig, MAF (8 formats total).

---

## Installation

**Q: Do I need Rust installed?**  
A: No, if using pre-built binaries. Yes, if building from source.

**Q: Does it work on Windows?**  
A: Yes, but BAM/SAM/CRAM support is not available. Use Linux/macOS or WSL for BAM files.

**Q: "Command not found" error?**  
A: Add to PATH: `export PATH="/path/to/fast-crossmap:$PATH"` or use full path.

**Q: "Permission denied" error?**  
A: Make executable: `chmod +x fast-crossmap`

---

## Usage

**Q: Where do I get chain files?**  
A: Download from UCSC: https://hgdownload.cse.ucsc.edu/downloads.html

**Q: Why does VCF/GVCF/MAF require a reference genome?**  
A: After coordinate conversion, the REF allele at the new position may differ. The reference genome is used to fetch the correct REF allele.

**Q: What is the `.unmap` file?**  
A: Records that couldn't be converted, with failure reasons.

**Q: Can I use compressed files?**  
A: Yes! Supports .gz and .bz2 for both chain files and input files.

**Q: How many threads should I use?**  
A: BED: 1-2 threads (I/O-bound), BAM: 4-8 threads (CPU-bound), VCF: 4 threads.

---

## Performance

**Q: Why is FastCrossMap so fast?**  
A: BITS algorithm (rust-lapper), zero-copy parsing, Rust's performance, streaming architecture.

**Q: Does multi-threading always help?**  
A: No. BED is I/O-bound (limited benefit). BAM is CPU-bound (significant benefit with 4-8 threads).

---

## Troubleshooting

### Many unmapped records

**Cause:** Wrong chain file direction or assembly mismatch.

**Solution:** Verify you're using the correct chain file:
- Converting hg19 → hg38? Use `hg19ToHg38.over.chain.gz`
- Converting hg38 → hg19? Use `hg38ToHg19.over.chain.gz`

### "Invalid BED format"

**Cause:** BED file has fewer than 3 columns.

**Solution:** Ensure BED has at least 3 columns (chr, start, end).

### Program hangs (VCF/GVCF/MAF)

**Cause:** Reference genome loading issue.

**Solution:** Ensure reference genome FASTA is indexed:
```bash
samtools faidx hg38.fa
```

### Many variants marked as failed (VCF)

**Cause:** Target chromosome not in reference genome.

**Solution:** Check if your reference genome contains all chromosomes. Alt contigs (e.g., `chr1_random`) may not be present in primary assembly files.

### Slow performance

**Solution:** Use multi-threading:
```bash
fast-crossmap bam -t 8 chain.gz input.bam output.bam
```

### Windows: BAM not supported

**Solution:** Use Linux, macOS, or WSL (Windows Subsystem for Linux).

---

## Getting Help

1. Check this FAQ
2. Search [GitHub Issues](https://github.com/biodancerwangzhi/fast-crossmap/issues)
3. Open a new issue with:
   - FastCrossMap version: `fast-crossmap --version`
   - Operating system
   - Command that caused the issue
   - Error message
4. Email: szxszx@foxmail.com

---

## Next Steps

- [Advanced Usage](Advanced) - Detailed CLI reference
- [Contributing](Contributing) - How to contribute
