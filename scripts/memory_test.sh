#!/bin/bash
# Memory profiling for FastCrossMap
#
# Usage: ./scripts/memory_test.sh <chain_file> <input_file> [format] [threads]
#
# Requirements:
#   - heaptrack (recommended) or valgrind
#   - /usr/bin/time with -v option (GNU time)

set -e

CHAIN_FILE=${1:?"Usage: $0 <chain_file> <input_file> [format] [threads]"}
INPUT_FILE=${2:?"Usage: $0 <chain_file> <input_file> [format] [threads]"}
FORMAT=${3:-bed}
THREADS=${4:-1}

OUTPUT_DIR="benchmark_results"
mkdir -p "$OUTPUT_DIR"

echo "=== FastCrossMap Memory Profiling ==="
echo "Chain file: $CHAIN_FILE"
echo "Input file: $INPUT_FILE"
echo "Format:     $FORMAT"
echo "Threads:    $THREADS"
echo ""

# Build release binary
echo "Building release binary..."
cargo build --release --quiet

# Count input lines
INPUT_LINES=$(grep -v "^#" "$INPUT_FILE" | wc -l)
echo "Input records: $INPUT_LINES"
echo ""

# Method 1: GNU time (basic memory stats)
echo "=== Method 1: GNU time ==="
/usr/bin/time -v ./target/release/fast-crossmap "$FORMAT" "$CHAIN_FILE" "$INPUT_FILE" \
    "$OUTPUT_DIR/memory_test_out.$FORMAT" -t "$THREADS" 2>&1 | \
    grep -E "(Maximum resident set size|User time|System time|Elapsed)"

# Method 2: heaptrack (detailed heap analysis)
if command -v heaptrack &> /dev/null; then
    echo ""
    echo "=== Method 2: heaptrack ==="
    heaptrack -o "$OUTPUT_DIR/heaptrack_fastcrossmap" \
        ./target/release/fast-crossmap "$FORMAT" "$CHAIN_FILE" "$INPUT_FILE" \
        "$OUTPUT_DIR/memory_test_out.$FORMAT" -t "$THREADS"
    
    echo "Heaptrack output saved to $OUTPUT_DIR/heaptrack_fastcrossmap.gz"
    echo "View with: heaptrack_gui $OUTPUT_DIR/heaptrack_fastcrossmap.gz"
    
    # Print summary if heaptrack_print is available
    if command -v heaptrack_print &> /dev/null; then
        echo ""
        echo "=== Heaptrack Summary ==="
        heaptrack_print "$OUTPUT_DIR/heaptrack_fastcrossmap.gz" 2>/dev/null | head -30
    fi
fi

# Method 3: valgrind massif (alternative)
if command -v valgrind &> /dev/null && [ ! -f "$OUTPUT_DIR/heaptrack_fastcrossmap.gz" ]; then
    echo ""
    echo "=== Method 3: Valgrind Massif ==="
    valgrind --tool=massif --massif-out-file="$OUTPUT_DIR/massif.out" \
        ./target/release/fast-crossmap "$FORMAT" "$CHAIN_FILE" "$INPUT_FILE" \
        "$OUTPUT_DIR/memory_test_out.$FORMAT" -t "$THREADS" 2>&1 | tail -5
    
    echo "Massif output saved to $OUTPUT_DIR/massif.out"
    echo "View with: ms_print $OUTPUT_DIR/massif.out"
fi

# Compare with CrossMap memory usage
echo ""
echo "=== CrossMap Memory Usage (for comparison) ==="
/usr/bin/time -v CrossMap "$FORMAT" "$CHAIN_FILE" "$INPUT_FILE" \
    "$OUTPUT_DIR/crossmap_memory_test_out.$FORMAT" 2>&1 | \
    grep -E "(Maximum resident set size|User time|System time|Elapsed)"

echo ""
echo "=== Memory Test Complete ==="
echo "Results saved to $OUTPUT_DIR/"
