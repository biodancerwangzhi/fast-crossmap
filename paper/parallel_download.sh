#!/bin/bash
# parallel_download.sh - Multi-threaded file download using curl range requests
# Retries failed chunks up to 3 times with no timeout on individual chunks.
# Usage: bash parallel_download.sh <URL> <OUTPUT_FILE> [THREADS]

URL="$1"
OUTPUT="$2"
THREADS="${3:-16}"
MAX_RETRIES=3

if [ -z "$URL" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: $0 <URL> <OUTPUT_FILE> [THREADS]"
    exit 1
fi

if [ -f "$OUTPUT" ]; then
    EXISTING_SIZE=$(stat -c%s "$OUTPUT")
    if [ "$EXISTING_SIZE" -gt 1000 ]; then
        echo "[skip] $OUTPUT already exists ($(( EXISTING_SIZE / 1024 / 1024 ))MB)"
        exit 0
    fi
    rm -f "$OUTPUT"
fi

TOTAL_SIZE=$(curl -sI --max-time 30 "$URL" | grep -i '^content-length:' | tail -1 | tr -d '\r' | awk '{print $2}')

if [ -z "$TOTAL_SIZE" ] || ! [[ "$TOTAL_SIZE" =~ ^[0-9]+$ ]] || [ "$TOTAL_SIZE" -eq 0 ]; then
    echo "[warn] Cannot get file size, falling back to single-thread download"
    curl -fSL --retry 3 -o "$OUTPUT" "$URL"
    exit $?
fi

CHUNK_SIZE=$(( (TOTAL_SIZE + THREADS - 1) / THREADS ))
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

SIZE_MB=$(( TOTAL_SIZE / 1024 / 1024 ))
echo "[download] $OUTPUT (${SIZE_MB}MB, ${THREADS} threads)"

download_chunk() {
    local idx=$1 start=$2 end=$3 outfile=$4
    local attempt
    for attempt in $(seq 1 $MAX_RETRIES); do
        if curl -s --retry 2 --retry-delay 5 -r "${start}-${end}" -o "$outfile" "$URL"; then
            local got=$(stat -c%s "$outfile" 2>/dev/null || echo 0)
            local expected=$(( end - start + 1 ))
            if [ "$got" -eq "$expected" ]; then
                return 0
            fi
        fi
        echo "[retry] chunk $idx attempt $attempt failed, retrying..." >&2
    done
    return 1
}

PIDS=""
for i in $(seq 0 $((THREADS - 1))); do
    START=$((i * CHUNK_SIZE))
    END=$(( (i + 1) * CHUNK_SIZE - 1 ))
    if [ $END -ge $TOTAL_SIZE ]; then
        END=$((TOTAL_SIZE - 1))
    fi
    download_chunk "$i" "$START" "$END" "${TMPDIR}/part_$(printf '%03d' $i)" &
    PIDS="$PIDS $!"
done

FAILED=0
for pid in $PIDS; do
    if ! wait $pid; then
        FAILED=$((FAILED + 1))
    fi
done

if [ $FAILED -gt 0 ]; then
    echo "[error] $FAILED chunks failed after $MAX_RETRIES retries for $OUTPUT"
    exit 1
fi

cat ${TMPDIR}/part_* > "$OUTPUT"
ACTUAL_SIZE=$(stat -c%s "$OUTPUT")

if [ "$ACTUAL_SIZE" -ne "$TOTAL_SIZE" ]; then
    echo "[error] Size mismatch: expected $TOTAL_SIZE, got $ACTUAL_SIZE"
    rm -f "$OUTPUT"
    exit 1
fi

echo "[done] $OUTPUT (${SIZE_MB}MB)"
