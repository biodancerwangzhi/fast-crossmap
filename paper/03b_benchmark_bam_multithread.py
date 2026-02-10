#!/usr/bin/env python3
"""
03b_benchmark_bam_multithread.py - BAM format multi-thread scalability test

Test FastCrossMap performance with different thread counts
Used to generate data for Figure 1(d)

Usage: python paper/03b_benchmark_bam_multithread.py
Output: paper/results/benchmark_bam_multithread.json
"""

import subprocess
import time
import json
import os
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
BAM_FILE = DATA_DIR / "encode_chipseq.bam"

# Thread counts to test
THREAD_COUNTS = [1, 2, 4, 8, 16]

# Number of runs per configuration
NUM_RUNS = 5


def get_file_size_mb(filepath):
    """Get file size (MB)"""
    return os.path.getsize(filepath) / (1024 * 1024)


def run_fastcrossmap_bam(chain_file, input_file, output_file, threads=1):
    """Run FastCrossMap BAM conversion and return execution time"""
    cmd = [
        "./fast-crossmap-linux-x64/fast-crossmap", "bam",
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
    print("FastCrossMap BAM Multi-Thread Scalability Test")
    print("=" * 60)
    
    # Check files
    if not CHAIN_FILE.exists():
        print(f"Error: Chain file not found: {CHAIN_FILE}")
        print("Please run first: bash paper/01_download_data.sh")
        return
    
    if not BAM_FILE.exists():
        print(f"Error: BAM file not found: {BAM_FILE}")
        print("Please run first: bash paper/01_download_data.sh")
        return
    
    # Get file size
    file_size_mb = get_file_size_mb(BAM_FILE)
    print(f"Input file: {BAM_FILE}")
    print(f"File size: {file_size_mb:.2f} MB")
    print(f"Thread counts: {THREAD_COUNTS}")
    print(f"Runs per configuration: {NUM_RUNS}")
    print()
    
    results = []
    
    for threads in THREAD_COUNTS:
        print(f"\nTesting {threads} threads...")
        output_file = RESULTS_DIR / f"fastcrossmap_bam_mt{threads}_output.bam"
        
        times = []
        for run in range(NUM_RUNS):
            result = run_fastcrossmap_bam(CHAIN_FILE, BAM_FILE, output_file, threads)
            if result["success"]:
                times.append(result["time"])
                print(f"  Run {run+1}: {result['time']:.2f}s")
            else:
                print(f"  Run {run+1}: FAILED - {result['stderr'][:100]}")
        
        if times:
            avg_time = sum(times) / len(times)
            min_time = min(times)
            max_time = max(times)
            throughput = file_size_mb / avg_time
            
            results.append({
                "threads": threads,
                "execution_time_sec": avg_time,
                "min_time_sec": min_time,
                "max_time_sec": max_time,
                "all_times": times,
                "throughput_mb_per_sec": throughput,
                "success": True
            })
            
            print(f"  Average: {avg_time:.2f}s (min: {min_time:.2f}s, max: {max_time:.2f}s)")
            print(f"  Throughput: {throughput:.2f} MB/sec")
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
                print(f"{r['threads']}T: {r['execution_time_sec']:.2f}s, "
                      f"Speedup: {speedup:.2f}x, Efficiency: {efficiency:.1f}%")
    
    # Save results
    output_data = {
        "timestamp": datetime.now().isoformat(),
        "format": "BAM",
        "input_file": str(BAM_FILE),
        "input_size_mb": file_size_mb,
        "chain_file": str(CHAIN_FILE),
        "num_runs": NUM_RUNS,
        "results": results
    }
    
    output_file = RESULTS_DIR / "benchmark_bam_multithread.json"
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nResults saved to: {output_file}")
    print("\nNext step: python paper/04_plot_performance.py")


if __name__ == "__main__":
    main()
