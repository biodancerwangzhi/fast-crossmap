#!/bin/bash
# setup_benchmark_env.sh - Set up benchmark environment on a fresh server
#
# Prerequisites: conda (miniconda/mambaforge), Rust toolchain (rustup)
#
# Usage:
#   bash paper/setup_benchmark_env.sh          # full setup
#   bash paper/setup_benchmark_env.sh --skip-data  # skip data download

set -e
cd "$(dirname "$0")/.."
ROOT=$(pwd)

SKIP_DATA=false
if [ "$1" = "--skip-data" ]; then
    SKIP_DATA=true
fi

echo "============================================================"
echo " FastCrossMap Benchmark Environment Setup"
echo " Root: $ROOT"
echo "============================================================"
echo ""

# --- 1. Build FastCrossMap ---
echo "[1/5] Building FastCrossMap..."
if command -v cargo &>/dev/null; then
    cargo build --release
    echo "  ✓ fast-crossmap binary: $ROOT/target/release/fast-crossmap"
else
    echo "  ✗ Rust/cargo not found. Install: curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh"
    exit 1
fi
echo ""

# --- 2. Install Python tools ---
echo "[2/5] Installing comparison tools..."

# CrossMap
if command -v CrossMap &>/dev/null; then
    echo "  ✓ CrossMap already installed"
else
    echo "  Installing CrossMap..."
    pip3 install CrossMap -i https://pypi.org/simple/ || pip install CrossMap
    echo "  ✓ CrossMap installed"
fi

# FastRemap
if command -v FastRemap &>/dev/null; then
    echo "  ✓ FastRemap already installed"
else
    echo "  Installing FastRemap (conda)..."
    conda install -y -c bioconda fastremap-bio 2>/dev/null || echo "  ⚠ FastRemap install failed (optional, BED-only)"
fi

# liftOver
if command -v liftOver &>/dev/null; then
    echo "  ✓ liftOver already installed"
else
    echo "  Downloading liftOver..."
    wget -q "https://genome-browser.s3.amazonaws.com/admin/exe/linux.x86_64/liftOver" -O liftOver && chmod +x liftOver
    sudo mv liftOver /usr/local/bin/ 2>/dev/null || mv liftOver "$HOME/.local/bin/" 2>/dev/null || {
        mkdir -p "$ROOT/bin" && mv liftOver "$ROOT/bin/"
        export PATH="$ROOT/bin:$PATH"
    }
    echo "  ✓ liftOver installed"
fi

# samtools
if command -v samtools &>/dev/null; then
    echo "  ✓ samtools already installed"
else
    echo "  Installing samtools (conda)..."
    conda install -y -c bioconda samtools 2>/dev/null || echo "  ⚠ samtools install failed"
fi

echo ""

# --- 3. Download data ---
if [ "$SKIP_DATA" = false ]; then
    echo "[3/5] Downloading benchmark data..."
    bash paper/01_download_data.sh
    echo ""
else
    echo "[3/5] Skipping data download (--skip-data)"
    echo ""
fi

# --- 4. Verify ---
echo "[4/5] Verifying installation..."
echo ""
echo "  Tools:"
echo "  fast-crossmap: $(./target/release/fast-crossmap --version 2>&1 || echo 'NOT FOUND')"
echo "  CrossMap:      $(CrossMap --version 2>&1 | head -1 || echo 'NOT FOUND')"
echo "  liftOver:      $(command -v liftOver 2>/dev/null && echo 'OK' || echo 'NOT FOUND')"
echo "  FastRemap:     $(command -v FastRemap 2>/dev/null && echo 'OK' || echo 'NOT FOUND (optional)')"
echo "  samtools:      $(samtools --version 2>&1 | head -1 || echo 'NOT FOUND')"
echo ""
echo "  Data files:"
if [ -d paper/data ]; then
    FILE_COUNT=$(ls paper/data/ 2>/dev/null | wc -l)
    TOTAL_SIZE=$(du -sh paper/data/ 2>/dev/null | cut -f1)
    echo "  $FILE_COUNT files, $TOTAL_SIZE total"
    echo ""
    echo "  Chain file: $(ls -lh paper/data/hg19ToHg38.over.chain.gz 2>/dev/null | awk '{print $5}' || echo 'MISSING')"
    echo "  hg38.fa:    $(ls -lh paper/data/hg38.fa 2>/dev/null | awk '{print $5}' || echo 'MISSING (needed for VCF/MAF)')"
else
    echo "  paper/data/ not found"
fi
echo ""

# --- 5. Instructions ---
echo "[5/5] Ready!"
echo ""
echo "============================================================"
echo " To run all benchmarks:"
echo "   bash paper/run_all_benchmarks.sh"
echo ""
echo " To run specific formats:"
echo "   bash paper/run_all_benchmarks.sh bed vcf gff"
echo ""
echo " Results will be saved to paper/results/"
echo "============================================================"
