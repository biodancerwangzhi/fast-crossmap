#!/bin/bash
# Benchmark FastCrossMap vs Python CrossMap
#
# Usage: ./scripts/benchmark.sh <chain_file> <input_file> [format] [threads]
#
# Examples:
#   ./scripts/benchmark.sh chain.gz input.bed bed 4
#   ./scripts/benchmark.sh chain.gz input.vcf vcf 8

set -e

CHAIN_FILE=${1:?"Usage: $0 <chain_file> <input_file> [format] [threads]"}
INPUT_FILE=${2:?"Usage: $0 <chain_file> <input_file> [format] [threads]"}
FORMAT=${3:-bed}
THREADS=${4:-4}

OUTPUT_DIR="benchmark_results"
mkdir -p "$OUTPUT_DIR"

echo "=== FastCrossMap vs CrossMap Benchmark ==="
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

# Warm up (load chain file into cache)
echo "Warming up..."
CrossMap "$FORMAT" "$CHAIN_FILE" "$INPUT_FILE" /dev/null 2>/dev/null || true

# Run benchmark with hyperfine
echo "Running benchmark..."
hyperfine \
    --warmup 2 \
    --min-runs 5 \
    --export-json "$OUTPUT_DIR/benchmark_${FORMAT}.json" \
    --export-markdown "$OUTPUT_DIR/benchmark_${FORMAT}.md" \
    --command-name "CrossMap" \
    "CrossMap $FORMAT $CHAIN_FILE $INPUT_FILE $OUTPUT_DIR/crossmap_out.$FORMAT" \
    --command-name "FastCrossMap (1 thread)" \
    "./target/release/fast-crossmap $FORMAT $CHAIN_FILE $INPUT_FILE $OUTPUT_DIR/fast_1t_out.$FORMAT -t 1" \
    --command-name "FastCrossMap ($THREADS threads)" \
    "./target/release/fast-crossmap $FORMAT $CHAIN_FILE $INPUT_FILE $OUTPUT_DIR/fast_${THREADS}t_out.$FORMAT -t $THREADS"

# Compare outputs
echo ""
echo "=== Output Comparison ==="
python scripts/compare_outputs.py "$OUTPUT_DIR/crossmap_out.$FORMAT" "$OUTPUT_DIR/fast_${THREADS}t_out.$FORMAT" --format "$FORMAT"

# Calculate speedup
echo ""
echo "=== Performance Summary ==="
if [ -f "$OUTPUT_DIR/benchmark_${FORMAT}.json" ]; then
    python3 -c "
import json
with open('$OUTPUT_DIR/benchmark_${FORMAT}.json') as f:
    data = json.load(f)
    results = data['results']
    crossmap_time = results[0]['mean']
    fast_1t_time = results[1]['mean']
    fast_mt_time = results[2]['mean']
    print(f'CrossMap:              {crossmap_time:.3f}s')
    print(f'FastCrossMap (1 thread): {fast_1t_time:.3f}s ({crossmap_time/fast_1t_time:.1f}x speedup)')
    print(f'FastCrossMap ($THREADS threads): {fast_mt_time:.3f}s ({crossmap_time/fast_mt_time:.1f}x speedup)')
    print(f'Parallel efficiency:   {fast_1t_time/fast_mt_time/$THREADS*100:.1f}%')
"
fi

echo ""
echo "Results saved to $OUTPUT_DIR/"
