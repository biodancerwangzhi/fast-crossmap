#!/bin/bash
# run_all_benchmarks.sh - Run all 8 format benchmarks sequentially
#
# Prerequisites:
#   - cargo build --release (fast-crossmap binary)
#   - pip install CrossMap
#   - conda install -c bioconda fastremap-bio (FastRemap)
#   - liftOver binary in PATH
#   - samtools in PATH (for BAM stats)
#   - All data files in paper/data/ (run: bash paper/01_download_data.sh)
#
# Usage: bash paper/run_all_benchmarks.sh [format...]
#   No args = run all 8 formats
#   With args = run only specified formats, e.g.: bash paper/run_all_benchmarks.sh bed vcf

set -e
cd "$(dirname "$0")/.."

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

FORMATS=(bed gff maf bigwig wig vcf bam sam)

if [ $# -gt 0 ]; then
    FORMATS=("$@")
fi

echo "============================================================"
echo " FastCrossMap Benchmark Suite"
echo " Formats: ${FORMATS[*]}"
echo " Started: $(date)"
echo "============================================================"
echo ""

# Check prerequisites
echo -e "${YELLOW}Checking prerequisites...${NC}"
MISSING=0
if [ ! -f target/release/fast-crossmap ]; then
    echo -e "${RED}  ✗ fast-crossmap binary not found. Run: cargo build --release${NC}"
    MISSING=1
fi
command -v CrossMap >/dev/null 2>&1 || { echo -e "${RED}  ✗ CrossMap not found. Run: pip install CrossMap${NC}"; MISSING=1; }
command -v liftOver >/dev/null 2>&1 || { echo -e "${YELLOW}  ⚠ liftOver not found (only needed for BED/GFF)${NC}"; }
command -v FastRemap >/dev/null 2>&1 || { echo -e "${YELLOW}  ⚠ FastRemap not found. Run: conda install -c bioconda fastremap-bio${NC}"; }
command -v samtools >/dev/null 2>&1 || { echo -e "${YELLOW}  ⚠ samtools not found (needed for BAM benchmark)${NC}"; }
if [ ! -f paper/data/hg19ToHg38.over.chain.gz ]; then
    echo -e "${RED}  ✗ Chain file not found. Run: bash paper/01_download_data.sh${NC}"
    MISSING=1
fi
if [ $MISSING -eq 1 ]; then
    echo -e "${RED}Missing required prerequisites. Aborting.${NC}"
    exit 1
fi
echo -e "${GREEN}  ✓ Prerequisites OK${NC}"
echo ""

SCRIPT_MAP=(
    "bed:paper/02_benchmark_bed.py"
    "vcf:paper/02c_benchmark_vcf.py"
    "gff:paper/02d_benchmark_gff.py"
    "wig:paper/02e_benchmark_wig.py"
    "bigwig:paper/02f_benchmark_bigwig.py"
    "maf:paper/02g_benchmark_maf.py"
    "sam:paper/02h_benchmark_sam.py"
    "bam:paper/03_benchmark_bam.py"
)

PASSED=0
FAILED=0
SKIPPED=0

for fmt in "${FORMATS[@]}"; do
    SCRIPT=""
    for entry in "${SCRIPT_MAP[@]}"; do
        key="${entry%%:*}"
        val="${entry#*:}"
        if [ "$key" = "$fmt" ]; then
            SCRIPT="$val"
            break
        fi
    done

    if [ -z "$SCRIPT" ]; then
        echo -e "${RED}Unknown format: $fmt${NC}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    if [ ! -f "$SCRIPT" ]; then
        echo -e "${RED}Script not found: $SCRIPT${NC}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    echo ""
    echo -e "${GREEN}============================================================${NC}"
    echo -e "${GREEN} [$fmt] Starting: $SCRIPT${NC}"
    echo -e "${GREEN} Time: $(date)${NC}"
    echo -e "${GREEN}============================================================${NC}"

    if python3 "$SCRIPT" 2>&1; then
        echo -e "${GREEN}  ✓ [$fmt] DONE${NC}"
        PASSED=$((PASSED + 1))
    else
        echo -e "${RED}  ✗ [$fmt] FAILED${NC}"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "============================================================"
echo " Benchmark Suite Complete"
echo " Finished: $(date)"
echo " Passed: $PASSED  Failed: $FAILED  Skipped: $SKIPPED"
echo "============================================================"
echo ""
echo "Results in paper/results/:"
ls -lh paper/results/benchmark_*.json 2>/dev/null
