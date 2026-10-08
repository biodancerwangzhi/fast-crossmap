#!/usr/bin/env python3
"""
03b_benchmark_bam_multithread.py - BAM format multi-thread scalability test

Test FastCrossMap performance with different thread counts
Used to generate data for Figure 1(d)

The input is read into page cache before each timed run: the thread counts are
tested in order and the first run would otherwise pay for the disk read on any
input larger than free RAM -- a distortion that would land on the t=1 baseline
(see `warm_page_cache`).

Usage: python paper/03b_benchmark_bam_multithread.py
Output: paper/results/benchmark_bam_multithread.json
"""

import subprocess
import shutil
import time
import json
import os
from pathlib import Path
from typing import Optional
from datetime import datetime

# =============================================================================
# 配置
# =============================================================================
DATA_DIR = Path("paper/data")
RESULTS_DIR = Path("paper/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def cleanup_outputs(*targets) -> None:
    """Delete a dataset's tool outputs once its numbers are already recorded.

    Every dataset writes each tool's mapped output at full input size -- the VCF
    ones are ~11 GB apiece -- so keeping them costs tens of GB that the results
    JSON already summarises. Removing the directory as soon as that dataset's
    results are in keeps peak disk use to one dataset instead of all of them,
    which is what lets a full re-run fit on a smaller volume. Set
    FCM_KEEP_OUTPUTS=1 in the environment to keep the files for inspection.

    Only paths inside RESULTS_DIR are ever removed, so a mistyped argument
    cannot delete an input or anything else outside the results tree.
    """
    if os.environ.get("FCM_KEEP_OUTPUTS"):
        print("    FCM_KEEP_OUTPUTS set: keeping outputs")
        return
    root = RESULTS_DIR.resolve()
    for target in targets:
        if target is None:
            continue
        target = Path(target)
        if not target.exists():
            continue
        try:
            target.resolve().relative_to(root)
        except ValueError:
            print(f"    refusing to clean {target}: outside {RESULTS_DIR}")
            continue
        freed = 0
        if target.is_dir():
            for p in target.rglob("*"):
                if p.is_file():
                    try:
                        freed += p.stat().st_size
                    except OSError:
                        pass
            shutil.rmtree(target, ignore_errors=True)
        else:
            try:
                freed = target.stat().st_size
            except OSError:
                freed = 0
            try:
                target.unlink()
            except OSError:
                pass
        print(f"    cleaned {target}  ({freed / 1e9:.2f} GB freed)".rstrip())


# FastCrossMap binary (override with FCM_BIN=/path/to/fast-crossmap)
FCM_BIN = os.environ.get("FCM_BIN", "./target/release/fast-crossmap")

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


def warm_page_cache(path: Optional[Path], chunk_size: int = 8 << 20) -> None:
    """Read `path` once and throw it away, so the timed run starts from RAM.

    Here the same command is repeated NUM_RUNS times per thread count and the
    thread counts are tested in order (1, 2, 4, 8, 16), so whoever runs first
    pays for the disk read whenever the input is larger than the machine's free
    RAM. On a 1.3 GB BAM that is not a small effect: it lands on the t=1 runs,
    which are exactly the baseline the whole speedup table is divided by, so a
    cold read there would understate every other thread count. The 2026-10-07
    GFF run shows the size of the distortion -- FastCrossMap recorded 80 s on a
    1.87 GB input against 5.20 s for the same file once it was in cache.

    Warming before *every* timed run rather than once per file keeps the
    comparison fair even when the input does not fit in RAM: each run then
    starts from whatever the cache holds after the same warm-up read, instead of
    from whatever the previous run happened to leave behind. The read is outside
    the timed region, so it cannot inflate a measurement.
    """
    if path is None or not path.exists():
        return
    with open(path, 'rb') as f:
        while f.read(chunk_size):
            pass


def run_fastcrossmap_bam(chain_file, input_file, output_file, threads=1):
    """Run FastCrossMap BAM conversion and return execution time"""
    cmd = [
        FCM_BIN, "bam",
        "-t", str(threads),
        str(chain_file),
        str(input_file),
        str(output_file)
    ]

    # Bring the input into page cache first, so the elapsed time below measures
    # the conversion and not who happened to read the file first.
    warm_page_cache(input_file)
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
    cleanup_outputs(*RESULTS_DIR.glob('fastcrossmap_bam_mt*_output.*'))
    print("\nNext step: python paper/04_plot_performance.py")


if __name__ == "__main__":
    main()
