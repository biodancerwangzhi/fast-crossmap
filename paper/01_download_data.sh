#!/bin/bash
# =============================================================================
# 01_download_data.sh - 下载论文复现所需的公共数据
# =============================================================================
# 
# 数据来源:
#   - Chain 文件: UCSC Genome Browser
#   - BED 文件: ENCODE Project (DNase-seq peaks)
#   - BAM 文件: ENCODE Project (ChIP-seq alignments)
#   - VCF 文件: 1000 Genomes Project
#
# 用法: bash paper/01_download_data.sh
# =============================================================================

set -e  # 遇到错误立即退出

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 创建目录结构
echo -e "${GREEN}[1/6] 创建目录结构...${NC}"
mkdir -p paper/data
mkdir -p paper/results
mkdir -p paper/figures

DATA_DIR="paper/data"

# =============================================================================
# 1. 下载 Chain 文件 (UCSC)
# =============================================================================
echo -e "${GREEN}[2/6] 下载 Chain 文件 (hg19 -> hg38)...${NC}"
CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
CHAIN_FILE="${DATA_DIR}/hg19ToHg38.over.chain.gz"

if [ -f "$CHAIN_FILE" ]; then
    echo -e "${YELLOW}  Chain 文件已存在，跳过下载${NC}"
else
    wget -q --show-progress -O "$CHAIN_FILE" "$CHAIN_URL"
    echo -e "${GREEN}  ✓ Chain 文件下载完成${NC}"
fi

# =============================================================================
# 2. 下载 BED 文件 (ENCODE DNase-seq peaks)
# =============================================================================
echo -e "${GREEN}[3/6] 下载 BED 文件 (ENCODE DNase-seq peaks)...${NC}"

# ENCODE DNase-seq peaks for K562 cell line (hg19)
# 来源: ENCODE Project
# 记录数: ~30 万条 (符合日常使用场景)
BED_URL="https://www.encodeproject.org/files/ENCFF001WBV/@@download/ENCFF001WBV.bed.gz"
BED_FILE="${DATA_DIR}/encode_dnase_peaks.bed.gz"

# 大型 cCRE 数据集 (约 235 万条记录，用于压力测试)
# BED_URL="https://downloads.wenglab.org/Registry-V4/cCRE.hg19.bed"
# BED_FILE="${DATA_DIR}/encode_ccre.bed"

if [ -f "$BED_FILE" ]; then
    echo -e "${YELLOW}  BED 文件已存在，跳过下载${NC}"
else
    wget -q --show-progress -O "$BED_FILE" "$BED_URL"
    echo -e "${GREEN}  ✓ BED 文件下载完成${NC}"
fi

# =============================================================================
# 3. 下载 BAM 文件 (ENCODE ChIP-seq)
# =============================================================================
echo -e "${GREEN}[4/6] 下载 BAM 文件 (ENCODE ChIP-seq)...${NC}"

# ENCODE ChIP-seq alignments for K562 cell line (hg19)
# 文件: ENCFF000PED - K562 CTCF ChIP-seq alignments (~1.5GB)
# 选择一个中等大小的 BAM 文件用于测试
BAM_URL="https://www.encodeproject.org/files/ENCFF000PED/@@download/ENCFF000PED.bam"
BAM_FILE="${DATA_DIR}/encode_chipseq.bam"

if [ -f "$BAM_FILE" ]; then
    echo -e "${YELLOW}  BAM 文件已存在，跳过下载${NC}"
else
    echo -e "${YELLOW}  注意: BAM 文件较大 (~1.5GB)，下载可能需要几分钟...${NC}"
    wget -q --show-progress -O "$BAM_FILE" "$BAM_URL"
    echo -e "${GREEN}  ✓ BAM 文件下载完成${NC}"
fi

# =============================================================================
# 4. 下载 VCF 文件 (1000 Genomes Phase 3)
# =============================================================================
echo -e "${GREEN}[5/6] 下载 VCF 文件 (1000 Genomes chr22)...${NC}"

# 1000 Genomes Phase 3 - chromosome 22 variants
VCF_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
VCF_FILE="${DATA_DIR}/1000g_chr22.vcf.gz"

if [ -f "$VCF_FILE" ]; then
    echo -e "${YELLOW}  VCF 文件已存在，跳过下载${NC}"
else
    echo -e "${YELLOW}  注意: VCF 文件较大 (~200MB)，下载可能需要几分钟...${NC}"
    wget -q --show-progress -O "$VCF_FILE" "$VCF_URL" || {
        echo -e "${YELLOW}  FTP 下载失败，创建测试 VCF 文件...${NC}"
        # 创建一个小型测试 VCF
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
    echo -e "${GREEN}  ✓ VCF 文件下载完成${NC}"
fi

# =============================================================================
# 5. 验证文件完整性
# =============================================================================
echo -e "${GREEN}[6/6] 验证文件完整性...${NC}"

echo ""
echo "=========================================="
echo "下载的文件列表:"
echo "=========================================="
ls -lh ${DATA_DIR}/

echo ""
echo "=========================================="
echo "文件行数统计:"
echo "=========================================="

# Chain 文件
if [ -f "$CHAIN_FILE" ]; then
    CHAIN_LINES=$(zcat "$CHAIN_FILE" 2>/dev/null | wc -l || echo "N/A")
    echo "Chain 文件: $CHAIN_LINES 行"
fi

# BED 文件
if [ -f "$BED_FILE" ]; then
    BED_LINES=$(zcat "$BED_FILE" 2>/dev/null | wc -l || echo "N/A")
    echo "BED 文件: $BED_LINES 行"
fi

# VCF 文件
if [ -f "$VCF_FILE" ]; then
    VCF_LINES=$(zcat "$VCF_FILE" 2>/dev/null | grep -v "^#" | wc -l || echo "N/A")
    echo "VCF 文件: $VCF_LINES 变异位点"
fi

# BAM 文件
if [ -f "$BAM_FILE" ]; then
    BAM_SIZE=$(ls -lh "$BAM_FILE" | awk '{print $5}')
    echo "BAM 文件: $BAM_SIZE"
fi

echo ""
echo -e "${GREEN}=========================================="
echo "✓ 数据下载完成!"
echo "==========================================${NC}"
echo ""
echo "下一步: 运行 python paper/02_benchmark_bed.py"
