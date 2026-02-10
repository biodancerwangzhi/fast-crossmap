#!/bin/bash
# =============================================================================
# 01_download_data.sh - Download public data for paper reproduction
# =============================================================================
# 
# Data sources:
#   - Chain file: UCSC Genome Browser
#   - BED file: ENCODE Project (DNase-seq peaks)
#   - BAM file: ENCODE Project (ChIP-seq alignments)
#   - VCF file: 1000 Genomes Project
#
# Usage: bash paper/01_download_data.sh
# =============================================================================

set -e  # Exit immediately on error

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Create directory structure
echo -e "${GREEN}[1/6] Creating directory structure...${NC}"
mkdir -p paper/data
mkdir -p paper/results
mkdir -p paper/figures

DATA_DIR="paper/data"

# =============================================================================
# 1. Download Chain file (UCSC)
# =============================================================================
echo -e "${GREEN}[2/6] Downloading Chain file (hg19 -> hg38)...${NC}"
CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
CHAIN_FILE="${DATA_DIR}/hg19ToHg38.over.chain.gz"

if [ -f "$CHAIN_FILE" ]; then
    echo -e "${YELLOW}  Chain file already exists, skipping download${NC}"
else
    wget -q --show-progress -O "$CHAIN_FILE" "$CHAIN_URL"
    echo -e "${GREEN}  ✓ Chain file download complete${NC}"
fi

# =============================================================================
# 2. Download BED file (ENCODE DNase-seq peaks)
# =============================================================================
echo -e "${GREEN}[3/6] Downloading BED file (ENCODE DNase-seq peaks)...${NC}"

# ENCODE DNase-seq peaks for K562 cell line (hg19)
# Source: ENCODE Project
# Records: ~300K (representative of typical usage)
BED_URL="https://www.encodeproject.org/files/ENCFF001WBV/@@download/ENCFF001WBV.bed.gz"
BED_FILE="${DATA_DIR}/encode_dnase_peaks.bed.gz"

# Large cCRE dataset (~2.35M records, for stress testing)
# BED_URL="https://downloads.wenglab.org/Registry-V4/cCRE.hg19.bed"
# BED_FILE="${DATA_DIR}/encode_ccre.bed"

if [ -f "$BED_FILE" ]; then
    echo -e "${YELLOW}  BED file already exists, skipping download${NC}"
else
    wget -q --show-progress -O "$BED_FILE" "$BED_URL"
    echo -e "${GREEN}  ✓ BED file download complete${NC}"
fi

# =============================================================================
# 3. Download BAM file (ENCODE ChIP-seq)
# =============================================================================
echo -e "${GREEN}[4/6] Downloading BAM file (ENCODE ChIP-seq)...${NC}"

# ENCODE ChIP-seq alignments for K562 cell line (hg19)
# File: ENCFF000PED - K562 CTCF ChIP-seq alignments (~1.5GB)
# Selected a medium-sized BAM file for testing
BAM_URL="https://www.encodeproject.org/files/ENCFF000PED/@@download/ENCFF000PED.bam"
BAM_FILE="${DATA_DIR}/encode_chipseq.bam"

if [ -f "$BAM_FILE" ]; then
    echo -e "${YELLOW}  BAM file already exists, skipping download${NC}"
else
    echo -e "${YELLOW}  Note: BAM file is large (~1.5GB), download may take a few minutes...${NC}"
    wget -q --show-progress -O "$BAM_FILE" "$BAM_URL"
    echo -e "${GREEN}  ✓ BAM file download complete${NC}"
fi

# =============================================================================
# 4. Download VCF file (1000 Genomes Phase 3)
# =============================================================================
echo -e "${GREEN}[5/6] Downloading VCF file (1000 Genomes chr22)...${NC}"

# 1000 Genomes Phase 3 - chromosome 22 variants
VCF_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
VCF_FILE="${DATA_DIR}/1000g_chr22.vcf.gz"

if [ -f "$VCF_FILE" ]; then
    echo -e "${YELLOW}  VCF file already exists, skipping download${NC}"
else
    echo -e "${YELLOW}  Note: VCF file is large (~200MB), download may take a few minutes...${NC}"
    wget -q --show-progress -O "$VCF_FILE" "$VCF_URL" || {
        echo -e "${YELLOW}  FTP download failed, creating test VCF file...${NC}"
        # Create a small test VCF
        cat > "${DATA_DIR}/test.vcf" << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr22,length=50818468>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr22	16050075	rs587697622	A	G	100	PASS	.
chr22	16050115	rs587755077	G	A	100	PASS	.
chr22	16050213	rs587654921	C	T	100	PASS	.
EOF
        gzip -k "${DATA_DIR}/test.vcf"
        mv "${DATA_DIR}/test.vcf.gz" "$VCF_FILE"
    }
    echo -e "${GREEN}  ✓ VCF file download complete${NC}"
fi

# =============================================================================
# 5. Verify file integrity
# =============================================================================
echo -e "${GREEN}[6/6] Verifying file integrity...${NC}"

echo ""
echo "=========================================="
echo "Downloaded files:"
echo "=========================================="
ls -lh ${DATA_DIR}/

echo ""
echo "=========================================="
echo "File line counts:"
echo "=========================================="

# Chain file
if [ -f "$CHAIN_FILE" ]; then
    CHAIN_LINES=$(zcat "$CHAIN_FILE" 2>/dev/null | wc -l || echo "N/A")
    echo "Chain file: $CHAIN_LINES lines"
fi

# BED file
if [ -f "$BED_FILE" ]; then
    BED_LINES=$(zcat "$BED_FILE" 2>/dev/null | wc -l || echo "N/A")
    echo "BED file: $BED_LINES lines"
fi

# VCF file
if [ -f "$VCF_FILE" ]; then
    VCF_LINES=$(zcat "$VCF_FILE" 2>/dev/null | grep -v "^#" | wc -l || echo "N/A")
    echo "VCF file: $VCF_LINES variants"
fi

# BAM file
if [ -f "$BAM_FILE" ]; then
    BAM_SIZE=$(ls -lh "$BAM_FILE" | awk '{print $5}')
    echo "BAM file: $BAM_SIZE"
fi

echo ""
echo -e "${GREEN}=========================================="
echo "✓ Data download complete!"
echo "==========================================${NC}"
echo ""
echo "Next step: run python paper/02_benchmark_bed.py"
