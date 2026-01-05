#!/bin/bash
# Test strict mode compatibility with CrossMap
# This script compares FastCrossMap --compat-mode=strict output with CrossMap output

set -e

CHAIN_FILE="${1:-ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz}"
TEST_FILE="${2:-test_data/compat_test.bed}"
FORMAT="${3:-bed}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=== Strict Mode Compatibility Test ==="
echo "Chain file: $CHAIN_FILE"
echo "Test file: $TEST_FILE"
echo "Format: $FORMAT"
echo ""

# Check prerequisites
if [ ! -f "$CHAIN_FILE" ]; then
    echo -e "${RED}Error: Chain file not found: $CHAIN_FILE${NC}"
    exit 1
fi

if [ ! -f "$TEST_FILE" ]; then
    echo -e "${YELLOW}Test file not found. Generating test data...${NC}"
    mkdir -p test_data
    python scripts/generate_test_data.py --format "$FORMAT" --count 1000 --output "$TEST_FILE"
fi

# Check if CrossMap is available
if ! command -v CrossMap &> /dev/null; then
    echo -e "${RED}Error: CrossMap not installed${NC}"
    echo "Install with: pip install CrossMap"
    exit 1
fi

# Create temp directory
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

echo "Running CrossMap..."
CrossMap "$FORMAT" "$CHAIN_FILE" "$TEST_FILE" "$TEMP_DIR/crossmap_output.$FORMAT" 2>/dev/null || true

echo "Running FastCrossMap (strict mode)..."
cargo run --release -- --compat-mode=strict "$FORMAT" "$CHAIN_FILE" "$TEST_FILE" "$TEMP_DIR/fast_strict_output.$FORMAT" 2>/dev/null

echo "Running FastCrossMap (improved mode)..."
cargo run --release -- --compat-mode=improved "$FORMAT" "$CHAIN_FILE" "$TEST_FILE" "$TEMP_DIR/fast_improved_output.$FORMAT" 2>/dev/null

echo ""
echo "=== Comparison Results ==="

# Compare strict mode with CrossMap
echo ""
echo "--- Strict Mode vs CrossMap ---"
python scripts/compare_outputs.py "$TEMP_DIR/crossmap_output.$FORMAT" "$TEMP_DIR/fast_strict_output.$FORMAT" --format "$FORMAT"

# Compare improved mode with CrossMap
echo ""
echo "--- Improved Mode vs CrossMap ---"
python scripts/compare_outputs.py "$TEMP_DIR/crossmap_output.$FORMAT" "$TEMP_DIR/fast_improved_output.$FORMAT" --format "$FORMAT"

# Byte-level comparison (sorted)
echo ""
echo "=== Byte-Level Comparison (sorted) ==="

# Sort and compare
sort "$TEMP_DIR/crossmap_output.$FORMAT" > "$TEMP_DIR/crossmap_sorted.$FORMAT" 2>/dev/null || true
sort "$TEMP_DIR/fast_strict_output.$FORMAT" > "$TEMP_DIR/fast_strict_sorted.$FORMAT" 2>/dev/null || true

if diff -q "$TEMP_DIR/crossmap_sorted.$FORMAT" "$TEMP_DIR/fast_strict_sorted.$FORMAT" > /dev/null 2>&1; then
    echo -e "${GREEN}✓ Strict mode output is BYTE-IDENTICAL to CrossMap (after sorting)${NC}"
else
    echo -e "${YELLOW}✗ Strict mode output differs from CrossMap${NC}"
    echo ""
    echo "First 10 differences:"
    diff "$TEMP_DIR/crossmap_sorted.$FORMAT" "$TEMP_DIR/fast_strict_sorted.$FORMAT" | head -20 || true
fi

echo ""
echo "=== Test Complete ==="
