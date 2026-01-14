#!/usr/bin/env python3
"""
benchmark_quick.py - Quick benchmark using sample data

This script runs a quick benchmark (~1 minute) using the small sample datasets
included in the repository. For full benchmark reproduction, see README.md.

Usage:
    python paper/benchmark_quick.py
"""

import os
import sys
import time
import subprocess
import platform
from pathlib import Path

# Configuration
SAMPLE_DIR = Path("paper/sample_data")
RESULTS_DIR = Path("paper/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Detect platform
IS_WINDOWS = platform.system() == "Windows"

# Tool paths
if IS_WINDOWS:
    FASTCROSSMAP = Path("target/release/fast-crossmap.exe")
else:
    FASTCROSSMAP = Path("target/release/fast-crossmap")


def check_prerequisites():
    """Check if all required files exist"""
    print("Checking prerequisites...")
    
    # Check FastCrossMap binary
    if not FASTCROSSMAP.exists():
        print(f"  ✗ FastCrossMap not found: {FASTCROSSMAP}")
        print("    Please build first: cargo build --release")
        return False
    print(f"  ✓ FastCrossMap: {FASTCROSSMAP}")
    
    # Check sample data
    chain_file = SAMPLE_DIR / "hg19ToHg38.over.chain.gz"
    bed_file = SAMPLE_DIR / "sample.bed"
    sam_file = SAMPLE_DIR / "sample.sam"
    
    missing = []
    for f in [chain_file, bed_file, sam_file]:
        if f.exists():
            print(f"  ✓ {f.name}")
        else:
            print(f"  ✗ {f.name} not found")
            missing.append(f)
    
    if missing:
        print("\n  Missing sample data. Run: python paper/generate_sample_data.py")
        return False
    
    return True


def run_benchmark(name, cmd, input_file):
    """Run a single benchmark and return timing"""
    print(f"\n  Running {name}...")
    
    # Get input file size
    input_size = input_file.stat().st_size
    
    # Run benchmark
    start = time.perf_counter()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )
        elapsed = time.perf_counter() - start
        
        if result.returncode != 0:
            print(f"    ✗ Failed: {result.stderr[:200]}")
            return None
        
        return {
            'name': name,
            'time': elapsed,
            'input_size': input_size,
            'throughput': input_size / elapsed / 1024 / 1024  # MB/s
        }
    except subprocess.TimeoutExpired:
        print(f"    ✗ Timeout (>5 min)")
        return None
    except Exception as e:
        print(f"    ✗ Error: {e}")
        return None


def benchmark_bed():
    """Benchmark BED file conversion"""
    print("\n" + "=" * 50)
    print("BED Benchmark (10,000 records)")
    print("=" * 50)
    
    chain_file = SAMPLE_DIR / "hg19ToHg38.over.chain.gz"
    bed_file = SAMPLE_DIR / "sample.bed"
    output_file = RESULTS_DIR / "quick_bed_output.bed"
    
    results = []
    
    # FastCrossMap (1 thread)
    cmd = [str(FASTCROSSMAP), "bed", str(chain_file), str(bed_file), str(output_file), "-t", "1"]
    result = run_benchmark("FastCrossMap (1T)", cmd, bed_file)
    if result:
        results.append(result)
        print(f"    ✓ Time: {result['time']:.3f}s")
    
    # FastCrossMap (4 threads)
    cmd = [str(FASTCROSSMAP), "bed", str(chain_file), str(bed_file), str(output_file), "-t", "4"]
    result = run_benchmark("FastCrossMap (4T)", cmd, bed_file)
    if result:
        results.append(result)
        print(f"    ✓ Time: {result['time']:.3f}s")
    
    # Count output records
    if output_file.exists():
        with open(output_file) as f:
            output_count = sum(1 for _ in f)
        print(f"\n  Output records: {output_count}")
    
    return results


def benchmark_sam():
    """Benchmark SAM file conversion"""
    # Skip SAM benchmark on Windows
    if IS_WINDOWS:
        print("\n" + "=" * 50)
        print("SAM Benchmark - SKIPPED (Windows)")
        print("=" * 50)
        print("  Note: BAM/SAM conversion is not supported on Windows.")
        print("  Please use Linux or macOS for BAM/SAM benchmarks.")
        return []
    
    print("\n" + "=" * 50)
    print("SAM Benchmark (1,000 reads)")
    print("=" * 50)
    
    chain_file = SAMPLE_DIR / "hg19ToHg38.over.chain.gz"
    sam_file = SAMPLE_DIR / "sample.sam"
    output_file = RESULTS_DIR / "quick_sam_output.sam"
    
    results = []
    
    # FastCrossMap (1 thread)
    cmd = [str(FASTCROSSMAP), "bam", str(chain_file), str(sam_file), str(output_file), "-t", "1"]
    result = run_benchmark("FastCrossMap (1T)", cmd, sam_file)
    if result:
        results.append(result)
        print(f"    ✓ Time: {result['time']:.3f}s")
    
    return results


def print_summary(bed_results, sam_results):
    """Print benchmark summary"""
    print("\n" + "=" * 50)
    print("Summary")
    print("=" * 50)
    
    print("\nBED Benchmark:")
    print(f"  {'Tool':<25} {'Time (s)':<12} {'Records/s':<15}")
    print("  " + "-" * 50)
    for r in bed_results:
        # Estimate records (10,000 for sample)
        records_per_sec = 10000 / r['time']
        print(f"  {r['name']:<25} {r['time']:<12.3f} {records_per_sec:<15,.0f}")
    
    if sam_results:
        print("\nSAM Benchmark:")
        print(f"  {'Tool':<25} {'Time (s)':<12}")
        print("  " + "-" * 50)
        for r in sam_results:
            print(f"  {r['name']:<25} {r['time']:<12.3f}")
    
    print("\n" + "=" * 50)
    print("Quick benchmark completed!")
    print("For full benchmark with ENCODE data, see paper/README.md")
    print("=" * 50)


def main():
    print("=" * 50)
    print("FastCrossMap Quick Benchmark")
    print("=" * 50)
    print()
    
    if not check_prerequisites():
        sys.exit(1)
    
    bed_results = benchmark_bed()
    sam_results = benchmark_sam()
    
    print_summary(bed_results, sam_results)


if __name__ == "__main__":
    main()
