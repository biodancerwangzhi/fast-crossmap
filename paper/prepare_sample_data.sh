#!/bin/bash
# prepare_sample_data.sh - Generate minimal sample data for tutorial/verification
#
# Creates small test files (100 records each) from existing benchmark data,
# then runs FastCrossMap to generate expected output.
#
# Usage: bash paper/prepare_sample_data.sh

set -e
cd "$(dirname "$0")/.."

SAMPLE_DIR="paper/sample_data"
INPUT_DIR="$SAMPLE_DIR/input"
OUTPUT_DIR="$SAMPLE_DIR/output"
DATA_DIR="paper/data"
CHAIN="$DATA_DIR/hg19ToHg38.over.chain.gz"
FCM="./target/release/fast-crossmap"

mkdir -p "$INPUT_DIR" "$OUTPUT_DIR"

echo "============================================================"
echo " Preparing sample data (8 formats)"
echo "============================================================"

# --- 1. BED ---
echo "[1/8] BED"
zcat "$DATA_DIR/bed_k562_dnase.bed.gz" | grep -v "^#" | head -100 > "$INPUT_DIR/sample.bed"
echo "  $(wc -l < "$INPUT_DIR/sample.bed") records"

# --- 2. VCF ---
echo "[2/8] VCF"
# Keep header + first 100 variant lines
zcat "$DATA_DIR/vcf_1kg_chr22.vcf.gz" | head -500 | grep "^#" > "$INPUT_DIR/sample.vcf"
zcat "$DATA_DIR/vcf_1kg_chr22.vcf.gz" | grep -v "^#" | head -100 >> "$INPUT_DIR/sample.vcf"
echo "  $(grep -cv '^#' "$INPUT_DIR/sample.vcf") variants"

# --- 3. GFF ---
echo "[3/8] GFF"
zcat "$DATA_DIR/gff_gencode_v44_lncrna_grch37.gff3.gz" | head -10 | grep "^#" > "$INPUT_DIR/sample.gff3"
zcat "$DATA_DIR/gff_gencode_v44_lncrna_grch37.gff3.gz" | grep -v "^#" | head -100 >> "$INPUT_DIR/sample.gff3"
echo "  $(grep -cv '^#' "$INPUT_DIR/sample.gff3") features"

# --- 4. SAM ---
echo "[4/8] SAM"
# Keep header + first 100 alignment lines
grep "^@" "$DATA_DIR/sam_k562_ctcf.sam" > "$INPUT_DIR/sample.sam"
grep -v "^@" "$DATA_DIR/sam_k562_ctcf.sam" | head -100 >> "$INPUT_DIR/sample.sam"
echo "  $(grep -cv '^@' "$INPUT_DIR/sample.sam") reads"

# --- 5. BAM (from sample SAM) ---
echo "[5/8] BAM"
if command -v samtools &>/dev/null; then
    samtools view -bS "$INPUT_DIR/sample.sam" > "$INPUT_DIR/sample.bam"
    echo "  $(samtools view -c "$INPUT_DIR/sample.bam") reads"
else
    echo "  ⚠ samtools not found, skipping BAM"
fi

# --- 6. WIG ---
echo "[6/8] WIG"
# Take header line + first 100 data lines from a real WIG
head -1 "$DATA_DIR/wig_k562_h3k27ac.wig" > "$INPUT_DIR/sample.wig"
grep -v "^variableStep\|^fixedStep\|^track\|^#" "$DATA_DIR/wig_k562_h3k27ac.wig" | head -99 >> "$INPUT_DIR/sample.wig"
# Actually we need proper wig format with chrom headers
# Redo: take first ~110 lines which includes header + data
head -110 "$DATA_DIR/wig_k562_h3k27ac.wig" > "$INPUT_DIR/sample.wig"
DATA_LINES=$(grep -cv "^variableStep\|^fixedStep\|^track\|^#\|^$" "$INPUT_DIR/sample.wig" || echo 0)
echo "  $DATA_LINES data points"

# --- 7. BigWig ---
echo "[7/8] BigWig"
if command -v wigToBigWig &>/dev/null; then
    # Need chrom sizes for wigToBigWig
    echo "  ⚠ wigToBigWig conversion skipped (use sample.wig with FCM wig subcommand)"
else
    echo "  ⚠ wigToBigWig not available, copying a small slice is not feasible"
    echo "  Users can test BigWig with the full benchmark files"
fi

# --- 8. MAF ---
echo "[8/8] MAF"
head -1 "$DATA_DIR/maf_synthetic_small.maf" > "$INPUT_DIR/sample.maf"
head -2 "$DATA_DIR/maf_synthetic_small.maf" | tail -1 >> "$INPUT_DIR/sample.maf"
grep -v "^#\|^Hugo_Symbol" "$DATA_DIR/maf_synthetic_small.maf" | head -100 >> "$INPUT_DIR/sample.maf"
echo "  $(grep -cv '^#\|^Hugo_Symbol' "$INPUT_DIR/sample.maf") mutations"

# --- Generate expected output ---
echo ""
echo "============================================================"
echo " Generating expected output with FastCrossMap"
echo "============================================================"

if [ ! -f "$FCM" ]; then
    echo "ERROR: fast-crossmap binary not found. Run: cargo build --release"
    exit 1
fi

# BED
echo "[1/6] BED output"
$FCM bed "$CHAIN" "$INPUT_DIR/sample.bed" "$OUTPUT_DIR/sample_hg38.bed"

# VCF (needs reference genome)
echo "[2/6] VCF output"
if [ -f "$DATA_DIR/hg38.fa" ]; then
    $FCM vcf "$CHAIN" "$INPUT_DIR/sample.vcf" "$DATA_DIR/hg38.fa" "$OUTPUT_DIR/sample_hg38.vcf"
else
    echo "  ⚠ hg38.fa not found, skipping VCF output"
fi

# GFF
echo "[3/6] GFF output"
$FCM gff "$CHAIN" "$INPUT_DIR/sample.gff3" "$OUTPUT_DIR/sample_hg38.gff3"

# SAM/BAM (uses bam subcommand for both)
echo "[4/6] SAM output"
$FCM bam "$CHAIN" "$INPUT_DIR/sample.sam" "$OUTPUT_DIR/sample_hg38.sam"

# WIG
echo "[5/6] WIG output"
$FCM wig "$CHAIN" "$INPUT_DIR/sample.wig" "$OUTPUT_DIR/sample_hg38.wig"

# MAF (needs reference genome)
echo "[6/6] MAF output"
if [ -f "$DATA_DIR/hg38.fa" ]; then
    $FCM maf -b GRCh38 "$CHAIN" "$INPUT_DIR/sample.maf" "$DATA_DIR/hg38.fa" "$OUTPUT_DIR/sample_hg38.maf"
else
    echo "  ⚠ hg38.fa not found, skipping MAF output"
fi

# --- Summary ---
echo ""
echo "============================================================"
echo " Sample data summary"
echo "============================================================"
echo ""
echo "Input files:"
ls -lh "$INPUT_DIR"/ 2>/dev/null
echo ""
echo "Output files:"
ls -lh "$OUTPUT_DIR"/ 2>/dev/null
echo ""
echo "Done! Sample data ready in $SAMPLE_DIR/"
