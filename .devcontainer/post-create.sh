#!/bin/bash
set -e

echo "=== Setting up FastCrossMap development environment ==="

cd /workspace

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to verify a tool
verify_tool() {
    local tool_name=$1
    local verify_cmd=$2
    local description=$3
    
    echo -n "  Checking $description... "
    if eval "$verify_cmd" > /dev/null 2>&1; then
        echo -e "${GREEN}✓${NC}"
        return 0
    else
        echo -e "${RED}✗ FAILED${NC}"
        echo -e "${RED}    Error: $tool_name is not accessible or not working properly${NC}"
        return 1
    fi
}

FAILED_TOOLS=()

# Verify installations
echo ""
echo "[1/6] Verifying Rust..."
if ! verify_tool "Rust" "rustc --version" "rustc"; then
    FAILED_TOOLS+=("Rust (rustc)")
fi
if ! verify_tool "Cargo" "cargo --version" "cargo"; then
    FAILED_TOOLS+=("Rust (cargo)")
fi

echo ""
echo "[2/6] Verifying CrossMap (Python)..."
if ! verify_tool "CrossMap" "CrossMap --version 2>/dev/null || python -c 'import CrossMap'" "CrossMap"; then
    FAILED_TOOLS+=("CrossMap")
fi

echo ""
echo "[3/6] Verifying samtools..."
if ! verify_tool "samtools" "samtools --version" "samtools"; then
    FAILED_TOOLS+=("samtools")
fi

echo ""
echo "[4/6] Verifying FastRemap..."
if ! verify_tool "FastRemap" "FastRemap --version 2>&1 | head -1" "FastRemap"; then
    # Try alternative command
    if ! verify_tool "FastRemap" "which FastRemap" "FastRemap (binary)"; then
        FAILED_TOOLS+=("FastRemap")
    fi
fi

echo ""
echo "[5/6] Verifying UCSC liftOver..."
if ! verify_tool "liftOver" "liftOver 2>&1 | head -1" "liftOver"; then
    # liftOver prints usage to stderr when called without args
    if ! verify_tool "liftOver" "which liftOver" "liftOver (binary)"; then
        FAILED_TOOLS+=("UCSC liftOver")
    fi
fi

echo ""
echo "[6/6] Building and verifying FastCrossMap..."
echo "  Building release binary..."
if cargo build --release --quiet 2>/dev/null; then
    echo -e "  Build: ${GREEN}✓${NC}"
    if ! verify_tool "fast-crossmap" "./target/release/fast-crossmap --version" "fast-crossmap binary"; then
        FAILED_TOOLS+=("FastCrossMap")
    fi
else
    echo -e "  Build: ${RED}✗ FAILED${NC}"
    FAILED_TOOLS+=("FastCrossMap (build failed)")
fi

# Create directories
mkdir -p /workspace/test_data
mkdir -p /workspace/benchmark_results
mkdir -p /workspace/results

# Make scripts executable
chmod +x scripts/*.sh 2>/dev/null || true
chmod +x scripts/*.py 2>/dev/null || true

echo ""
echo "=========================================="

# Report results
if [ ${#FAILED_TOOLS[@]} -eq 0 ]; then
    echo -e "  ${GREEN}FastCrossMap Dev Environment Ready!${NC}"
    echo "=========================================="
    echo ""
    echo "Available tools for benchmarking:"
    echo "  - fast-crossmap : FastCrossMap (Rust)"
    echo "  - CrossMap      : CrossMap (Python)"
    echo "  - liftOver      : UCSC liftOver"
    echo "  - FastRemap     : FastRemap"
    echo ""
    echo "Other tools:"
    echo "  - cargo/rustc   : Rust compiler"
    echo "  - samtools      : BAM/SAM utilities"
    echo "  - valgrind      : Memory debugging"
    echo "  - heaptrack     : Heap profiler"
    echo "  - hyperfine     : Benchmarking"
    echo ""
    echo "Quick start:"
    echo "  cargo build --release"
    echo "  cargo test"
    echo ""
else
    echo -e "  ${RED}Environment Setup INCOMPLETE${NC}"
    echo "=========================================="
    echo ""
    echo -e "${RED}The following tools failed verification:${NC}"
    for tool in "${FAILED_TOOLS[@]}"; do
        echo -e "  ${RED}✗ $tool${NC}"
    done
    echo ""
    echo "Please check the installation and try again."
    echo "You may need to rebuild the container or install missing packages."
    exit 1
fi
