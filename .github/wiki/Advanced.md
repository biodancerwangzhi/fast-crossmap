# Advanced Usage

Complete reference for FastCrossMap commands, options, and performance tuning.

---

## 1. Command-Line Reference

### Global Options

```bash
fast-crossmap [OPTIONS] <SUBCOMMAND>
```

| Option | Description |
|--------|-------------|
| `--compat-mode <MODE>` | `strict` (CrossMap-identical) or `improved` (default) |
| `-h, --help` | Print help |
| `-V, --version` | Print version |

### Common Options (all formats)

| Option | Description |
|--------|-------------|
| `-t, --threads <N>` | Thread count (0 = all cores) |
| `--chromid <MODE>` | Chromosome ID style: `a` (as-is), `s` (short: `1`), `l` (long: `chr1`) |

---

## 2. Format-Specific Usage

### BED Format

```bash
fast-crossmap bed [OPTIONS] <CHAIN> <INPUT> [OUTPUT]
```

**Examples:**
```bash
fast-crossmap bed hg19ToHg38.chain.gz input.bed output.bed
fast-crossmap bed -t 4 hg19ToHg38.chain.gz input.bed.gz output.bed.gz
```

**Supported variants:** BED3, BED6, BED12

---

### BAM/SAM/CRAM Format

```bash
fast-crossmap bam [OPTIONS] <CHAIN> <INPUT> <OUTPUT>
```

**Examples:**
```bash
fast-crossmap bam hg19ToHg38.chain.gz input.bam output.bam
fast-crossmap bam -t 8 hg19ToHg38.chain.gz input.bam output.bam
```

**Note:** BAM support not available on Windows. Use Linux/macOS or WSL.

---

### VCF Format

```bash
fast-crossmap vcf [OPTIONS] <CHAIN> <INPUT> <REFGENOME> [OUTPUT]
```

**Parameters:**
- `<REFGENOME>`: Target reference genome FASTA (required)

**Options:**
- `--no-comp-allele`: Do not compare REF allele with reference genome

**Examples:**
```bash
fast-crossmap vcf hg19ToHg38.chain.gz input.vcf hg38.fa output.vcf
fast-crossmap vcf -t 4 --no-comp-allele hg19ToHg38.chain.gz input.vcf hg38.fa output.vcf
```

---

### GVCF Format

```bash
fast-crossmap gvcf [OPTIONS] <CHAIN> <INPUT> <REFGENOME> [OUTPUT]
```

**Examples:**
```bash
fast-crossmap gvcf hg19ToHg38.chain.gz input.g.vcf hg38.fa output.g.vcf
fast-crossmap gvcf -t 4 hg19ToHg38.chain.gz input.g.vcf.gz hg38.fa output.g.vcf.gz
```

**Special handling:**
- Non-variant blocks (with `END=`): Keep original REF allele
- Variant records: Fetch new REF from target reference genome

---

### GFF/GTF Format

```bash
fast-crossmap gff [OPTIONS] <CHAIN> <INPUT> [OUTPUT]
```

**Examples:**
```bash
fast-crossmap gff hg19ToHg38.chain.gz genes.gff3 genes.hg38.gff3
fast-crossmap gff hg19ToHg38.chain.gz genes.gtf genes.hg38.gtf
```

---

### Wiggle Format

```bash
fast-crossmap wig [OPTIONS] <CHAIN> <INPUT> [OUTPUT]
```

**Example:**
```bash
fast-crossmap wig hg19ToHg38.chain.gz input.wig output.wig
```

---

### BigWig Format

```bash
fast-crossmap bigwig [OPTIONS] <CHAIN> <INPUT> [OUTPUT]
```

**Example:**
```bash
fast-crossmap bigwig hg19ToHg38.chain.gz input.bw output.bw
```

---

### MAF Format

```bash
fast-crossmap maf [OPTIONS] <CHAIN> <INPUT> <REFGENOME> -b <BUILD> [OUTPUT]
```

**Parameters:**
- `<REFGENOME>`: Target reference genome FASTA (required)
- `-b, --build <BUILD>`: Target genome build name (required, e.g., `hg38`, `GRCh38`)

**Example:**
```bash
fast-crossmap maf hg19ToHg38.chain.gz input.maf hg38.fa -b hg38 output.maf
```

---

## 3. Multi-threading

### Usage

```bash
fast-crossmap bed -t 4 chain.gz input.bed output.bed
fast-crossmap bam -t 0 chain.gz input.bam output.bam  # 0 = all cores
```

### Recommended Thread Counts

| File Size | Threads |
|-----------|---------|
| < 100 MB | 1-2 |
| 100 MB - 1 GB | 4 |
| 1 GB - 5 GB | 8 |
| > 5 GB | 8-16 |

### Performance by Format

| Format | Speedup Potential | Recommendation |
|--------|-------------------|----------------|
| BED | Limited (I/O-bound) | 1-2 threads |
| BAM | Significant (CPU-bound) | 4-8 threads |
| VCF | Moderate | 4 threads |

**Memory:** Multi-threading does NOT significantly increase memory (~16 MB → ~20 MB with 8 threads)

---

## 4. Compatibility Mode

### Strict Mode (CrossMap-identical)

```bash
fast-crossmap --compat-mode strict bed chain.gz input.bed output.bed
```

**Use when:**
- Need exact CrossMap compatibility
- Reproducing published analyses
- Comparing results with CrossMap

### Improved Mode (Default)

```bash
fast-crossmap bed chain.gz input.bed output.bed
```

**Use when:**
- Want optimal performance
- Starting new analyses

| Aspect | Strict | Improved |
|--------|--------|----------|
| CrossMap compatibility | 100% | ~99.9% |
| Performance | Fast | Faster |

---

## 5. Performance Comparison

### Speed (vs CrossMap baseline)

| Format | FastCrossMap (1T) | CrossMap | liftOver | FastRemap (4T) |
|--------|-------------------|----------|----------|----------------|
| BED (297K) | **0.28s (10x)** | 2.85s (1x) | 3.28s | 1.67s |
| BAM (1.3GB) | **86s (8x)** | 679s (1x) | N/A | 58s |

### Memory

| Tool | Peak Memory (BAM) |
|------|-------------------|
| **FastCrossMap** | **16 MB** |
| CrossMap | 1,022 MB |
| FastRemap | 70 MB |

FastCrossMap uses **64x less memory** than CrossMap.

### Feature Comparison

| Feature | FastCrossMap | CrossMap | liftOver | FastRemap |
|---------|:------------:|:--------:|:--------:|:---------:|
| BED | ✅ | ✅ | ✅ | ✅ |
| BAM/SAM/CRAM | ✅ | ✅ | ❌ | ✅ |
| VCF/GVCF | ✅ | ✅ | ❌ | ❌ |
| GFF/GTF | ✅ | ✅ | ✅ | ❌ |
| Wiggle/BigWig | ✅ | ✅ | ❌ | ❌ |
| MAF | ✅ | ✅ | ❌ | ❌ |
| Multi-threading | ✅ | ❌ | ❌ | ⚠️ (fixed 4T) |
| Compressed chain | ✅ | ✅ | ✅ | ❌ |

---

## 6. Why Reference Genome is Required (VCF/GVCF/MAF)

After coordinate conversion, the genomic position changes. The REF allele at the new position may differ:

1. **Original**: chr1:1000 REF=A (in hg19)
2. **After liftover**: chr1:2000 (in hg38)
3. **New REF**: Must fetch from hg38 reference at position 2000

This ensures the output file is valid against the target reference genome.

---

## Next Steps

- [FAQ](FAQ) - Common questions and troubleshooting
- [Contributing](Contributing) - How to contribute
