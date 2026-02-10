#!/usr/bin/env python3
"""
02b_benchmark_bed_multithread.py - BED format multi-thread scalability test

Test FastCrossMap performance with different thread counts
Used to generate data for Figure 1(b)

Usage: python paper/02b_benchmark_bed_multithread.py
Output: paper/results/benchmark_bed_multithread.json
"""

import subprocess
import time
import json
from pathlib import Path
from datetime import datetime

# =============================================================================
# 配置
# =============================================================================
DATA_DIR = Path("paper/data")
RESULTS_DIR = Path("paper/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Test files
CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
BED_FILE = DATA_DIR / "encode_dnase_peaks.bed.gz"

# Thread counts to test
THREAD_COUNTS = [1, 2, 4, 8, 16]

# Number of runs per configuration
NUM_RUNS = 5


def count_bed_records(bed_file):
    """Count BED file records (supports .gz compression)"""
    import gzip
    
    count = 0
    bed_path = Path(bed_file)
    
    if str(bed_path).endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    with opener(bed_path, mode) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                count += 1
    return count


def run_fastcrossmap(chain_file, input_file, output_file, threads=1):
    """Run FastCrossMap and return execution time"""
    cmd = [
        "./fast-crossmap-linux-x64/fast-crossmap", "bed",
        "-t", str(threads),
        str(chain_file),
        str(input_file),
        str(output_file)
    ]
    
    start = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - start
    
    return {
        "success": result.returncode == 0,
        "time": elapsed,
        "stderr": result.stderr
    }


def main():
    print("=" * 60)
    print("FastCrossMap BED Multi-Thread Scalability Test")
    print("=" * 60)
    
    # Check files
    if not CHAIN_FILE.exists():
        print(f"Error: Chain file not found: {CHAIN_FILE}")
        print("Please run first: bash paper/01_download_data.sh")
        return
    
    if not BED_FILE.exists():
        print(f"Error: BED file not found: {BED_FILE}")
        print("Please run first: bash paper/01_download_data.sh")
        return
    
    # If BED file is .gz format, decompress first
    import gzip
    import shutil
    
    bed_file_to_use = BED_FILE
    if str(BED_FILE).endswith('.gz'):
        bed_file_unzipped = Path(str(BED_FILE)[:-3])  # Remove .gz suffix
        if not bed_file_unzipped.exists():
            print(f"Decompressing BED file: {BED_FILE} -> {bed_file_unzipped}")
            with gzip.open(BED_FILE, 'rb') as f_in:
                with open(bed_file_unzipped, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        bed_file_to_use = bed_file_unzipped
    
    # Count records
    num_records = count_bed_records(BED_FILE)
    print(f"Input file: {bed_file_to_use}")
    print(f"Records: {num_records:,}")
    print(f"Thread counts: {THREAD_COUNTS}")
    print(f"Runs per configuration: {NUM_RUNS}")
    print()
    
    results = []
    
    for threads in THREAD_COUNTS:
        print(f"\nTesting {threads} threads...")
        output_file = RESULTS_DIR / f"fastcrossmap_mt{threads}_output.bed"
        
        times = []
        for run in range(NUM_RUNS):
            result = run_fastcrossmap(CHAIN_FILE, bed_file_to_use, output_file, threads)
            if result["success"]:
                times.append(result["time"])
                print(f"  Run {run+1}: {result['time']:.3f}s")
            else:
                print(f"  Run {run+1}: FAILED - {result['stderr'][:100]}")
        
        if times:
            avg_time = sum(times) / len(times)
            min_time = min(times)
            max_time = max(times)
            throughput = num_records / avg_time
            
            results.append({
                "threads": threads,
                "execution_time_sec": avg_time,
                "min_time_sec": min_time,
                "max_time_sec": max_time,
                "all_times": times,
                "throughput_rec_per_sec": throughput,
                "success": True
            })
            
            print(f"  Average: {avg_time:.3f}s (min: {min_time:.3f}s, max: {max_time:.3f}s)")
            print(f"  Throughput: {throughput:,.0f} records/sec")
        else:
            results.append({
                "threads": threads,
                "success": False,
                "error": "All runs failed"
            })
    
    # Calculate speedup
    if results and results[0]["success"]:
        baseline = results[0]["execution_time_sec"]
        print("\n" + "=" * 60)
        print("Scalability Analysis")
        print("=" * 60)
        for r in results:
            if r["success"]:
                speedup = baseline / r["execution_time_sec"]
                efficiency = speedup / r["threads"] * 100
                print(f"{r['threads']}T: {r['execution_time_sec']:.3f}s, "
                      f"Speedup: {speedup:.2f}x, Efficiency: {efficiency:.1f}%")
    
    # Save results
    output_data = {
        "timestamp": datetime.now().isoformat(),
        "format": "BED",
        "input_file": str(BED_FILE),
        "input_records": num_records,
        "chain_file": str(CHAIN_FILE),
        "num_runs": NUM_RUNS,
        "results": results
    }
    
    output_file = RESULTS_DIR / "benchmark_bed_multithread.json"
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nResults saved to: {output_file}")
    print("\nNext step: python paper/04_plot_performance.py")


if __name__ == "__main__":
    main()
