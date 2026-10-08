#!/bin/bash
# =============================================================================
# 01_download_data.sh - Download public data for paper reproduction
# =============================================================================
#
# Data sources:
#   - Chain file: UCSC Genome Browser (hg19 -> hg38)
#   - BED files: ENCODE Project (DNase-seq peaks, multiple cell lines)
#   - BAM files: ENCODE Project (ChIP-seq alignments, multiple experiments)
#   - VCF files: 1000 Genomes (chr22, chr1, chr11, chr20, chr2)
#   - GFF files: GENCODE (GRCh37), RefSeq (GRCh37)
#   - Reference genome: hg38 (UCSC, for VCF liftover)
#
# Usage: bash paper/01_download_data.sh [--skip-large] [--bed-only] [--vcf-only] ...
# =============================================================================

set -e

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

DATA_DIR="paper/data"
mkdir -p "$DATA_DIR"
mkdir -p paper/results
mkdir -p paper/figures

# Parse arguments
SKIP_LARGE=false
SECTIONS=""
for arg in "$@"; do
    case $arg in
        --skip-large) SKIP_LARGE=true ;;
        --bed-only)   SECTIONS="bed" ;;
        --bam-only)   SECTIONS="bam" ;;
        --vcf-only)   SECTIONS="vcf" ;;
        --gff-only)   SECTIONS="gff" ;;
        --bigwig-only) SECTIONS="bigwig" ;;
        --wig-only)   SECTIONS="wig" ;;
        --maf-only)   SECTIONS="maf" ;;
        --sam-only)   SECTIONS="sam" ;;
        --ref-only)   SECTIONS="ref" ;;
        --all)        SECTIONS="" ;;
    esac
done

should_run() {
    [ -z "$SECTIONS" ] || echo "$SECTIONS" | grep -q "$1"
}

download_file() {
    local url="$1"
    local output="$2"
    local desc="$3"

    if [ -f "$output" ]; then
        echo -e "${YELLOW}  [skip] $desc already exists${NC}"
        return 0
    fi

    echo -e "${BLUE}  Downloading $desc ...${NC}"
    if wget -q --show-progress -c -O "${output}" "$url" 2>/dev/null; then
        echo -e "${GREEN}  ✓ $desc${NC}"
    elif curl -fSL -C - -o "${output}" "$url" 2>/dev/null; then
        echo -e "${GREEN}  ✓ $desc (via curl)${NC}"
    else
        rm -f "${output}"
        echo -e "${RED}  ✗ Failed to download $desc${NC}"
        echo -e "${RED}    URL: $url${NC}"
        return 1
    fi
}

# =============================================================================
# 1. Chain file (required by all benchmarks)
# =============================================================================
echo -e "${GREEN}[1/11] Chain file (hg19 -> hg38)${NC}"
download_file \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" \
    "${DATA_DIR}/hg19ToHg38.over.chain.gz" \
    "hg19ToHg38 chain file"

# =============================================================================
# 2. BED files (ENCODE DNase-seq peaks, hg19, multiple cell lines)
# =============================================================================
if should_run "bed"; then
echo -e "${GREEN}[2/11] BED files (ENCODE DNase-seq peaks, hg19)${NC}"

# Dataset 1: K562 DNase-seq peaks (~297K records) — original dataset
download_file \
    "https://www.encodeproject.org/files/ENCFF001WBV/@@download/ENCFF001WBV.bed.gz" \
    "${DATA_DIR}/bed_k562_dnase.bed.gz" \
    "BED: K562 DNase-seq peaks (ENCFF001WBV)"

# Dataset 2: GM12878 DNase-seq peaks (hg19)
download_file \
    "https://www.encodeproject.org/files/ENCFF001WFH/@@download/ENCFF001WFH.bed.gz" \
    "${DATA_DIR}/bed_gm12878_dnase.bed.gz" \
    "BED: GM12878 DNase-seq peaks (ENCFF001WFH)"

# Dataset 3: H1-hESC DNase-seq peaks (hg19)
download_file \
    "https://www.encodeproject.org/files/ENCFF001WDU/@@download/ENCFF001WDU.bed.gz" \
    "${DATA_DIR}/bed_h1hesc_dnase.bed.gz" \
    "BED: H1-hESC DNase-seq peaks (ENCFF001WDU)"

# Dataset 4: HepG2 DNase-seq narrowPeak (hg19, UCSC ENCODE)
download_file \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseHepg2Pk.narrowPeak.gz" \
    "${DATA_DIR}/bed_hepg2_dnase.narrowPeak.gz" \
    "BED: HepG2 DNase-seq narrowPeak (UCSC ENCODE)"

# Dataset 5: A549 DNase-seq narrowPeak (hg19, UCSC ENCODE)
download_file \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeOpenChromDnase/wgEncodeOpenChromDnaseA549Pk.narrowPeak.gz" \
    "${DATA_DIR}/bed_a549_dnase.narrowPeak.gz" \
    "BED: A549 DNase-seq narrowPeak (UCSC ENCODE)"

# Bonus: ENCODE cCRE (GRCh38, not used in main benchmark — different coordinate system)
if [ "$SKIP_LARGE" = false ]; then
download_file \
    "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed" \
    "${DATA_DIR}/bed_ccre_grch38.bed" \
    "BED: ENCODE cCRE registry (GRCh38, ~2.35M records, bonus)"
fi

fi

# =============================================================================
# 3. BAM files (ENCODE + 1000 Genomes, hg19, multiple experiments & sizes)
# =============================================================================
if should_run "bam"; then
echo -e "${GREEN}[3/11] BAM files (hg19)${NC}"

# Dataset 1: K562 CTCF ChIP-seq (~1.3GB) — original dataset
download_file \
    "https://www.encodeproject.org/files/ENCFF000PED/@@download/ENCFF000PED.bam" \
    "${DATA_DIR}/bam_k562_ctcf.bam" \
    "BAM: K562 CTCF ChIP-seq (ENCFF000PED, ~1.3GB)"

# Dataset 2: K562 ChIP-seq (ENCODE, different experiment)
download_file \
    "https://www.encodeproject.org/files/ENCFF000PEE/@@download/ENCFF000PEE.bam" \
    "${DATA_DIR}/bam_k562_chipseq2.bam" \
    "BAM: K562 ChIP-seq (ENCFF000PEE)"

# Dataset 3: 1000 Genomes NA12878 chr11 (~672MB, different source)
download_file \
    "https://1000genomes.s3.amazonaws.com/phase3/data/NA12878/alignment/NA12878.chrom11.ILLUMINA.bwa.CEU.low_coverage.20121211.bam" \
    "${DATA_DIR}/bam_na12878_chr11.bam" \
    "BAM: 1000G NA12878 chr11 (~672MB)"

# Dataset 4 (optional small): 1000 Genomes NA12878 chr20 (~297MB)
download_file \
    "https://1000genomes.s3.amazonaws.com/phase3/data/NA12878/alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam" \
    "${DATA_DIR}/bam_na12878_chr20.bam" \
    "BAM: 1000G NA12878 chr20 (~297MB)"

# Dataset 5: 1000 Genomes HG00096 chr11 (~661MB, different individual)
download_file \
    "https://1000genomes.s3.amazonaws.com/phase3/data/HG00096/alignment/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam" \
    "${DATA_DIR}/bam_hg00096_chr11.bam" \
    "BAM: 1000G HG00096 chr11 (~661MB)"

fi

# =============================================================================
# 4. VCF files (1000 Genomes, ClinVar — GRCh37/hg19)
# =============================================================================
if should_run "vcf"; then
echo -e "${GREEN}[4/11] VCF files (GRCh37/hg19)${NC}"

# Dataset 1: 1000 Genomes Phase 3 chr22 (~200MB, ~1M variants)
download_file \
    "https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" \
    "${DATA_DIR}/vcf_1kg_chr22.vcf.gz" \
    "VCF: 1000 Genomes chr22 (~1M variants)"

# Dataset 2: 1000 Genomes Phase 3 chr1 (~1.2GB, ~6M variants, large)
if [ "$SKIP_LARGE" = false ]; then
download_file \
    "https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" \
    "${DATA_DIR}/vcf_1kg_chr1.vcf.gz" \
    "VCF: 1000 Genomes chr1 (~6M variants, large)"
fi

# Dataset 3: 1000 Genomes Phase 3 chr11 (~731MB, ~3M variants)
download_file \
    "https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" \
    "${DATA_DIR}/vcf_1kg_chr11.vcf.gz" \
    "VCF: 1000 Genomes chr11 (~3M variants)"

# Dataset 4: 1000 Genomes Phase 3 chr20 (~312MB, medium)
download_file \
    "https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" \
    "${DATA_DIR}/vcf_1kg_chr20.vcf.gz" \
    "VCF: 1000 Genomes chr20 (~1.8M variants)"

# Dataset 5: 1000 Genomes Phase 3 chr2 (~1.2GB, large)
if [ "$SKIP_LARGE" = false ]; then
download_file \
    "https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" \
    "${DATA_DIR}/vcf_1kg_chr2.vcf.gz" \
    "VCF: 1000 Genomes chr2 (~3.6M variants, large)"
fi

fi

# =============================================================================
# 5. BigWig files (UCSC/ENCODE signal tracks, hg19)
# =============================================================================
if should_run "bigwig"; then
echo -e "${GREEN}[5/11] BigWig files (hg19)${NC}"

# Dataset 1: K562 H3K27ac signal (~292MB)
download_file \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig" \
    "${DATA_DIR}/bigwig_k562_h3k27ac.bigWig" \
    "BigWig: K562 H3K27ac signal (~292MB)"

# Dataset 2: GM12878 H3K4me3 signal (~434MB)
download_file \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneGm12878H3k4me3StdSig.bigWig" \
    "${DATA_DIR}/bigwig_gm12878_h3k4me3.bigWig" \
    "BigWig: GM12878 H3K4me3 signal (~434MB)"

# Dataset 3: K562 H3K4me1 signal (~355MB)
download_file \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H3k4me1StdSig.bigWig" \
    "${DATA_DIR}/bigwig_k562_h3k4me1.bigWig" \
    "BigWig: K562 H3K4me1 signal (~355MB)"

# Dataset 4: K562 H3K36me3 signal (~376MB)
download_file \
    "https://genome-browser.s3.amazonaws.com/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H3k36me3StdSig.bigWig" \
    "${DATA_DIR}/bigwig_k562_h3k36me3.bigWig" \
    "BigWig: K562 H3K36me3 signal (~376MB)"

# Dataset 5: GM12878 H3K27ac signal (~250MB)
download_file \
    "https://genome-browser.s3.amazonaws.com/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneGm12878H3k27acStdSig.bigWig" \
    "${DATA_DIR}/bigwig_gm12878_h3k27ac.bigWig" \
    "BigWig: GM12878 H3K27ac signal (~250MB)"

fi

# =============================================================================
# 6. Wiggle files (generated from BigWig or synthetic)
# =============================================================================
if should_run "wig"; then
echo -e "${GREEN}[6/11] Wiggle files${NC}"

# Generate Wiggle from BigWig if bigWigToWig is available,
# otherwise generate synthetic Wiggle data
if command -v bigWigToWig &>/dev/null; then
    for bw in "${DATA_DIR}"/bigwig_*.bigWig; do
        [ -f "$bw" ] || continue
        base=$(basename "$bw" .bigWig)
        wig="${DATA_DIR}/wig_${base#bigwig_}.wig"
        if [ -f "$wig" ]; then
            echo -e "${YELLOW}  [skip] $(basename $wig) already exists${NC}"
        else
            echo -e "${BLUE}  Converting $(basename $bw) -> $(basename $wig) ...${NC}"
            bigWigToWig "$bw" "$wig"
            echo -e "${GREEN}  ✓ $(basename $wig)${NC}"
        fi
    done
else
    echo -e "${YELLOW}  bigWigToWig not found, generating synthetic Wiggle data...${NC}"
    python3 -c "
import random
random.seed(42)
chroms = {'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430,
          'chr7': 159138663, 'chr11': 135006516, 'chr20': 63025520, 'chr22': 51304566}
sizes = [200000, 150000, 100000]
names = ['large', 'medium', 'small']
for size, name in zip(sizes, names):
    fname = '${DATA_DIR}/wig_synthetic_' + name + '.wig'
    with open(fname, 'w') as f:
        per_chrom = size // len(chroms)
        for chrom, length in chroms.items():
            f.write(f'variableStep chrom={chrom} span=10\n')
            positions = sorted(random.sample(range(1, length - 100), per_chrom))
            for pos in positions:
                f.write(f'{pos}\t{random.uniform(0, 100):.2f}\n')
    print(f'  Generated {fname} ({size} records)')
"
fi

fi

# =============================================================================
# 7. MAF files (synthetic realistic data for benchmark)
# =============================================================================
if should_run "maf"; then
echo -e "${GREEN}[7/11] MAF files${NC}"

echo -e "${BLUE}  Generating synthetic MAF datasets...${NC}"
python3 -c "
import random
random.seed(42)
chroms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X']
chrom_lengths = {
    '1':249250621,'2':243199373,'3':198022430,'4':191154276,'5':180915260,
    '6':171115067,'7':159138663,'8':146364022,'9':141213431,'10':135534747,
    '11':135006516,'12':133851895,'13':115169878,'14':107349540,'15':102531392,
    '16':90354753,'17':81195210,'18':78077248,'19':59128983,'20':63025520,
    '21':48129895,'22':51304566,'X':155270560
}
var_classes = ['Missense_Mutation','Silent','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site','In_Frame_Del','In_Frame_Ins']
var_types = ['SNP','DEL','INS']
bases = 'ACGT'
header = '#version 2.4\nHugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\n'

sizes = [1000000, 500000, 200000, 100000, 50000]
names = ['xlarge', 'large', 'medium', 'small', 'tiny']
for size, name in zip(sizes, names):
    fname = '${DATA_DIR}/maf_synthetic_' + name + '.maf'
    with open(fname, 'w') as f:
        f.write(header)
        for i in range(size):
            ch = random.choice(chroms)
            pos = random.randint(1, chrom_lengths[ch] - 20)
            span = random.randint(1, 10)
            ref = ''.join(random.choice(bases) for _ in range(random.randint(1,3)))
            alt1 = ''.join(random.choice(bases) for _ in range(random.randint(1,3)))
            alt2 = ''.join(random.choice(bases) for _ in range(random.randint(1,3)))
            f.write(f'GENE{i:06d}\t{random.randint(1000,99999)}\tTestCenter\tGRCh37\t{ch}\t{pos}\t{pos+span}\t{random.choice([\"+\",\"-\"])}\t{random.choice(var_classes)}\t{random.choice(var_types)}\t{ref}\t{alt1}\t{alt2}\tSAMPLE_{i%1000:04d}\n')
    print(f'  Generated {fname} ({size} records)')
"

fi

# =============================================================================
# 8. SAM files (converted from BAM via samtools)
# =============================================================================
if should_run "sam"; then
echo -e "${GREEN}[8/11] SAM files (converted from BAM)${NC}"

if ! command -v samtools &>/dev/null; then
    echo -e "${RED}  samtools not found, cannot convert BAM to SAM${NC}"
else
    for bam in "${DATA_DIR}"/bam_*.bam; do
        [ -f "$bam" ] || continue
        base=$(basename "$bam" .bam)
        sam="${DATA_DIR}/sam_${base#bam_}.sam"
        if [ -f "$sam" ]; then
            echo -e "${YELLOW}  [skip] $(basename $sam) already exists${NC}"
        else
            echo -e "${BLUE}  Converting $(basename $bam) -> $(basename $sam) ...${NC}"
            samtools view -h "$bam" > "$sam"
            echo -e "${GREEN}  ✓ $(basename $sam) ($(ls -lh $sam | awk '{print $5}'))${NC}"
        fi
    done
fi

fi

# =============================================================================
# 9. GFF files (GENCODE — GRCh37/hg19)
# =============================================================================
if should_run "gff"; then
echo -e "${GREEN}[9/11] GFF files (GRCh37/hg19)${NC}"

# Dataset 1: GENCODE v44 GRCh37 mapped (comprehensive gene annotation)
download_file \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.annotation.gff3.gz" \
    "${DATA_DIR}/gff_gencode_v44_grch37.gff3.gz" \
    "GFF: GENCODE v44 GRCh37 annotation"

# Dataset 2: GENCODE v44 GRCh37 basic annotation (smaller subset)
download_file \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.basic.annotation.gff3.gz" \
    "${DATA_DIR}/gff_gencode_v44_basic_grch37.gff3.gz" \
    "GFF: GENCODE v44 basic GRCh37 annotation"

# Dataset 3: GENCODE v44 long non-coding RNA (different gene set)
download_file \
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.long_noncoding_RNAs.gff3.gz" \
    "${DATA_DIR}/gff_gencode_v44_lncrna_grch37.gff3.gz" \
    "GFF: GENCODE v44 lncRNA GRCh37 annotation"

# Dataset 4: RefSeq GRCh37 latest (different annotation source)
download_file \
    "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz" \
    "${DATA_DIR}/gff_refseq_grch37.gff.gz" \
    "GFF: RefSeq GRCh37 latest annotation"

# Dataset 5: Ensembl GRCh37 release 87 (different annotation source)
download_file \
    "https://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz" \
    "${DATA_DIR}/gff_ensembl_grch37.gff3.gz" \
    "GFF: Ensembl GRCh37 release 87 annotation"

fi

# =============================================================================
# 10. Reference genomes (required for VCF/GVCF/MAF liftover and re-alignment benchmark)
# =============================================================================
if should_run "ref"; then
echo -e "${GREEN}[10/12] Reference genomes${NC}"

# hg38 (required for VCF/GVCF/MAF liftover)
HG38_FA="${DATA_DIR}/hg38.fa"
HG38_GZ="${DATA_DIR}/hg38.fa.gz"

if [ -f "$HG38_FA" ]; then
    echo -e "${YELLOW}  [skip] hg38.fa already exists${NC}"
elif [ -f "$HG38_GZ" ]; then
    echo -e "${BLUE}  Decompressing hg38.fa.gz ...${NC}"
    gunzip -k "$HG38_GZ"
    echo -e "${GREEN}  ✓ hg38.fa decompressed${NC}"
else
    download_file \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" \
        "$HG38_GZ" \
        "hg38 reference genome (~900MB compressed)"
    if [ -f "$HG38_GZ" ]; then
        echo -e "${BLUE}  Decompressing hg38.fa.gz ...${NC}"
        gunzip -k "$HG38_GZ"
        echo -e "${GREEN}  ✓ hg38.fa decompressed${NC}"
    fi
fi

# hg19 (required for re-alignment benchmark: wgsim + bwa-mem2)
HG19_FA="${DATA_DIR}/hg19.fa"
HG19_GZ="${DATA_DIR}/hg19.fa.gz"

if [ -f "$HG19_FA" ]; then
    echo -e "${YELLOW}  [skip] hg19.fa already exists${NC}"
elif [ -f "$HG19_GZ" ]; then
    echo -e "${BLUE}  Decompressing hg19.fa.gz ...${NC}"
    gunzip -k "$HG19_GZ"
    echo -e "${GREEN}  ✓ hg19.fa decompressed${NC}"
else
    download_file \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz" \
        "$HG19_GZ" \
        "hg19 reference genome (~900MB compressed)"
    if [ -f "$HG19_GZ" ]; then
        echo -e "${BLUE}  Decompressing hg19.fa.gz ...${NC}"
        gunzip -k "$HG19_GZ"
        echo -e "${GREEN}  ✓ hg19.fa decompressed${NC}"
    fi
fi

fi

# =============================================================================
# 11. Verify all files
# =============================================================================
echo -e "${GREEN}[11/11] Verifying downloaded files${NC}"
echo ""
echo "=========================================="
echo "Downloaded files:"
echo "=========================================="

for f in "${DATA_DIR}"/*; do
    [ -f "$f" ] || continue
    size=$(ls -lh "$f" | awk '{print $5}')
    name=$(basename "$f")
    printf "  %-45s %8s\n" "$name" "$size"
done

echo ""
echo "=========================================="
echo "Record counts:"
echo "=========================================="

# BED files
for f in "${DATA_DIR}"/bed_*.bed.gz "${DATA_DIR}"/bed_*.bed; do
    [ -f "$f" ] || continue
    name=$(basename "$f")
    if [[ "$f" == *.gz ]]; then
        count=$(zcat "$f" 2>/dev/null | grep -cv "^#" || echo "N/A")
    else
        count=$(grep -cv "^#" "$f" 2>/dev/null || echo "N/A")
    fi
    printf "  %-40s %s records\n" "$name" "$count"
done

# VCF files
for f in "${DATA_DIR}"/vcf_*.vcf.gz; do
    [ -f "$f" ] || continue
    name=$(basename "$f")
    count=$(zcat "$f" 2>/dev/null | grep -cv "^#" || echo "N/A")
    printf "  %-40s %s variants\n" "$name" "$count"
done

# GFF files
for f in "${DATA_DIR}"/gff_*.gff3.gz; do
    [ -f "$f" ] || continue
    name=$(basename "$f")
    count=$(zcat "$f" 2>/dev/null | grep -cv "^#" || echo "N/A")
    printf "  %-40s %s features\n" "$name" "$count"
done

# BAM files (size + read count if samtools available)
for f in "${DATA_DIR}"/bam_*.bam; do
    [ -f "$f" ] || continue
    name=$(basename "$f")
    size=$(ls -lh "$f" | awk '{print $5}')
    if command -v samtools &>/dev/null; then
        reads=$(samtools view -c "$f" 2>/dev/null || echo "N/A")
        printf "  %-40s %s  (%s reads)\n" "$name" "$size" "$reads"
    else
        printf "  %-40s %s\n" "$name" "$size"
    fi
done

# SAM files
for f in "${DATA_DIR}"/sam_*.sam; do
    [ -f "$f" ] || continue
    name=$(basename "$f")
    count=$(grep -cv "^@" "$f" 2>/dev/null || echo "N/A")
    size=$(ls -lh "$f" | awk '{print $5}')
    printf "  %-40s %s  (%s reads)\n" "$name" "$size" "$count"
done

# BigWig files
for f in "${DATA_DIR}"/bigwig_*.bigWig; do
    [ -f "$f" ] || continue
    name=$(basename "$f")
    size=$(ls -lh "$f" | awk '{print $5}')
    printf "  %-40s %s\n" "$name" "$size"
done

# Wiggle files
for f in "${DATA_DIR}"/wig_*.wig; do
    [ -f "$f" ] || continue
    name=$(basename "$f")
    count=$(grep -cv "^variableStep\|^fixedStep\|^track\|^#\|^$" "$f" 2>/dev/null || echo "N/A")
    printf "  %-40s %s data points\n" "$name" "$count"
done

# MAF files
for f in "${DATA_DIR}"/maf_*.maf; do
    [ -f "$f" ] || continue
    name=$(basename "$f")
    count=$(grep -cv "^#\|^Hugo_Symbol" "$f" 2>/dev/null || echo "N/A")
    printf "  %-40s %s mutations\n" "$name" "$count"
done

echo ""
echo -e "${GREEN}=========================================="
echo "✓ Data download complete!"
echo "==========================================${NC}"
echo ""
echo "Next steps:"
echo "  python paper/02_benchmark_bed.py     # BED benchmark"
echo "  python paper/02c_benchmark_vcf.py    # VCF benchmark"
echo "  python paper/02d_benchmark_gff.py    # GFF benchmark"
echo "  python paper/02e_benchmark_wig.py    # Wiggle benchmark"
echo "  python paper/02f_benchmark_bigwig.py # BigWig benchmark"
echo "  python paper/02g_benchmark_maf.py    # MAF benchmark"
echo "  python paper/02h_benchmark_sam.py    # SAM benchmark"
echo "  python paper/03_benchmark_bam.py     # BAM benchmark"
