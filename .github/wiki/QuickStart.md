# Quick Start

Get started with FastCrossMap in 5 minutes!

## Step 1: Download Required Files

```bash
# Chain file: Human hg19 → hg38
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

# Reference genome (required for VCF/GVCF/MAF)
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

## Step 2: Convert Your File

### BED Format
```bash
fast-crossmap bed hg19ToHg38.over.chain.gz input.bed output.bed
```

### BAM Format
```bash
fast-crossmap bam hg19ToHg38.over.chain.gz input.bam output.bam
```

### VCF Format (requires reference genome)
```bash
fast-crossmap vcf hg19ToHg38.over.chain.gz input.vcf hg38.fa output.vcf
```

### GVCF Format (requires reference genome)
```bash
fast-crossmap gvcf hg19ToHg38.over.chain.gz input.g.vcf hg38.fa output.g.vcf
```

### GFF/GTF Format
```bash
fast-crossmap gff hg19ToHg38.over.chain.gz input.gff output.gff
fast-crossmap gff hg19ToHg38.over.chain.gz input.gtf output.gtf
```

### Wiggle Format
```bash
fast-crossmap wig hg19ToHg38.over.chain.gz input.wig output.wig
```

### BigWig Format
```bash
fast-crossmap bigwig hg19ToHg38.over.chain.gz input.bw output.bw
```

### MAF Format (requires reference genome and build name)
```bash
fast-crossmap maf hg19ToHg38.over.chain.gz input.maf hg38.fa -b hg38 output.maf
```

## Step 3: Use Multi-threading

```bash
# Use 4 threads for faster processing
fast-crossmap bed -t 4 hg19ToHg38.over.chain.gz input.bed output.bed
fast-crossmap bam -t 8 hg19ToHg38.over.chain.gz input.bam output.bam
fast-crossmap vcf -t 4 hg19ToHg38.over.chain.gz input.vcf hg38.fa output.vcf
```

## Output Files

FastCrossMap generates:
- **Main output**: Successfully converted records
- **`.unmap` file**: Records that couldn't be converted

## Next Steps

- [Advanced Usage](Advanced) - Detailed CLI reference, format-specific options, performance tuning
- [FAQ](FAQ) - Common questions and troubleshooting
