#!/usr/bin/env python3
"""
03_benchmark_bam.py - BAM format benchmark

Run FastCrossMap vs CrossMap vs FastRemap benchmark
(liftOver does not support BAM format)

Usage: python paper/03_benchmark_bam.py
Output: paper/results/benchmark_bam.json
"""

import json
import os
import subprocess
import time
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional

# =============================================================================
# 配置
# =============================================================================
DATA_DIR = Path("paper/data")
RESULTS_DIR = Path("paper/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Input files
CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
BAM_FILE = DATA_DIR / "encode_chipseq.bam"

# Number of runs
NUM_RUNS = 5


@dataclass
class BenchmarkResult:
    """Benchmark result"""
    tool: str
    format: str
    input_file: str
    input_size_mb: float
    execution_time_sec: float
    throughput_mb_per_sec: float
    peak_memory_mb: float
    success: bool
    error_message: Optional[str] = None


def get_file_size_mb(file_path: Path) -> float:
    """Get file size (MB)"""
    return file_path.stat().st_size / (1024 * 1024)


def run_with_time(cmd: list[str]) -> tuple[float, float, bool, str]:
    """
    Run command and measure time and memory.
    Returns: (execution_time_sec, peak_memory_mb, success, error_message)
    """
    import os
    
    start_time = time.time()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )
        elapsed = time.time() - start_time
        
        # Try to get memory via /usr/bin/time
        peak_memory_mb = 0
        try:
            time_result = subprocess.run(
                ["/usr/bin/time", "-v"] + cmd,
                capture_output=True,
                text=True,
                timeout=3600
            )
            for line in time_result.stderr.split('\n'):
                if 'Maximum resident set size' in line:
                    peak_memory_kb = int(line.split(':')[1].strip())
                    peak_memory_mb = peak_memory_kb / 1024
                    break
        except:
            try:
                import psutil
                peak_memory_mb = psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)
            except:
                peak_memory_mb = 0
        
        if result.returncode != 0:
            return elapsed, peak_memory_mb, False, result.stderr[:500]
        
        return elapsed, peak_memory_mb, True, ""
        
    except subprocess.TimeoutExpired:
        return 3600, 0, False, "Timeout after 3600 seconds"
    except Exception as e:
        return 0, 0, False, str(e)


def benchmark_fastcrossmap(bam_file: Path, chain_file: Path, output_dir: Path) -> BenchmarkResult:
    """Benchmark FastCrossMap"""
    print("  Running FastCrossMap...")
    output_file = output_dir / "fastcrossmap_output.bam"
    
    # FastCrossMap uses positional args
    cmd = [
        "./fast-crossmap-linux-x64/fast-crossmap",
        "bam",
        str(chain_file),
        str(bam_file),
        str(output_file)
    ]
    
    times = []
    memories = []
    success = False
    error_msg = ""
    
    for i in range(NUM_RUNS):
        print(f"    Run {i+1}/{NUM_RUNS}...")
        elapsed, memory, ok, err = run_with_time(cmd)
        if ok:
            times.append(elapsed)
            memories.append(memory)
            success = True
        else:
            error_msg = err
    
    input_size = get_file_size_mb(bam_file)
    
    if not times:
        return BenchmarkResult(
            tool="FastCrossMap",
            format="BAM",
            input_file=str(bam_file),
            input_size_mb=input_size,
            execution_time_sec=0,
            throughput_mb_per_sec=0,
            peak_memory_mb=0,
            success=False,
            error_message=error_msg
        )
    
    avg_time = sum(times) / len(times)
    avg_memory = sum(memories) / len(memories)
    
    return BenchmarkResult(
        tool="FastCrossMap",
        format="BAM",
        input_file=str(bam_file),
        input_size_mb=round(input_size, 2),
        execution_time_sec=round(avg_time, 2),
        throughput_mb_per_sec=round(input_size / avg_time, 2),
        peak_memory_mb=round(avg_memory, 2),
        success=success
    )


def benchmark_crossmap(bam_file: Path, chain_file: Path, output_dir: Path) -> BenchmarkResult:
    """Benchmark CrossMap"""
    print("  Running CrossMap...")
    output_file = output_dir / "crossmap_output.bam"
    
    cmd = [
        "CrossMap", "bam",
        "-a",  # Output all reads
        str(chain_file),
        str(bam_file),
        str(output_file)
    ]
    
    times = []
    memories = []
    success = False
    error_msg = ""
    
    for i in range(NUM_RUNS):
        print(f"    Run {i+1}/{NUM_RUNS}...")
        elapsed, memory, ok, err = run_with_time(cmd)
        if ok:
            times.append(elapsed)
            memories.append(memory)
            success = True
        else:
            error_msg = err
    
    input_size = get_file_size_mb(bam_file)
    
    if not times:
        return BenchmarkResult(
            tool="CrossMap",
            format="BAM",
            input_file=str(bam_file),
            input_size_mb=input_size,
            execution_time_sec=0,
            throughput_mb_per_sec=0,
            peak_memory_mb=0,
            success=False,
            error_message=error_msg
        )
    
    avg_time = sum(times) / len(times)
    avg_memory = sum(memories) / len(memories)
    
    return BenchmarkResult(
        tool="CrossMap",
        format="BAM",
        input_file=str(bam_file),
        input_size_mb=round(input_size, 2),
        execution_time_sec=round(avg_time, 2),
        throughput_mb_per_sec=round(input_size / avg_time, 2),
        peak_memory_mb=round(avg_memory, 2),
        success=success
    )


def benchmark_fastremap(bam_file: Path, chain_file: Path, output_dir: Path) -> BenchmarkResult:
    """Benchmark FastRemap"""
    print("  Running FastRemap...")
    output_file = output_dir / "fastremap_output.bam"
    unmap_file = output_dir / "fastremap_output.bam.unmap"
    
    # FastRemap does not support .gz, needs uncompressed chain file
    chain_unzipped = chain_file.parent / "hg19ToHg38.over.chain"
    if not chain_unzipped.exists():
        print("    Decompressing chain file for FastRemap...")
        import subprocess
        subprocess.run(["gunzip", "-k", str(chain_file)], check=True)
    
    cmd = [
        "FastRemap",
        "-f", "bam",
        "-c", str(chain_unzipped),
        "-i", str(bam_file),
        "-o", str(output_file),
        "-u", str(unmap_file)
    ]
    
    times = []
    memories = []
    success = False
    error_msg = ""
    
    for i in range(NUM_RUNS):
        print(f"    Run {i+1}/{NUM_RUNS}...")
        elapsed, memory, ok, err = run_with_time(cmd)
        if ok:
            times.append(elapsed)
            memories.append(memory)
            success = True
        else:
            error_msg = err
    
    input_size = get_file_size_mb(bam_file)
    
    if not times:
        return BenchmarkResult(
            tool="FastRemap",
            format="BAM",
            input_file=str(bam_file),
            input_size_mb=input_size,
            execution_time_sec=0,
            throughput_mb_per_sec=0,
            peak_memory_mb=0,
            success=False,
            error_message=error_msg
        )
    
    avg_time = sum(times) / len(times)
    avg_memory = sum(memories) / len(memories)
    
    return BenchmarkResult(
        tool="FastRemap",
        format="BAM",
        input_file=str(bam_file),
        input_size_mb=round(input_size, 2),
        execution_time_sec=round(avg_time, 2),
        throughput_mb_per_sec=round(input_size / avg_time, 2),
        peak_memory_mb=round(avg_memory, 2),
        success=success
    )


def main():
    print("=" * 60)
    print("BAM Format Benchmark")
    print("=" * 60)
    
    # Check input files
    if not BAM_FILE.exists():
        print(f"Error: BAM file not found: {BAM_FILE}")
        print("Please run first: bash paper/01_download_data.sh")
        return
    
    if not CHAIN_FILE.exists():
        print(f"Error: Chain file not found: {CHAIN_FILE}")
        print("Please run first: bash paper/01_download_data.sh")
        return
    
    # Get file size
    input_size = get_file_size_mb(BAM_FILE)
    print(f"\nInput file: {BAM_FILE}")
    print(f"File size: {input_size:.2f} MB")
    print(f"Number of runs: {NUM_RUNS}")
    print()
    print("Note: liftOver does not support BAM format")
    print()
    
    # Create output directory
    output_dir = RESULTS_DIR / "bam_benchmark"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run benchmarks
    results = []
    
    print("[1/3] FastCrossMap")
    results.append(benchmark_fastcrossmap(BAM_FILE, CHAIN_FILE, output_dir))
    
    print("[2/3] CrossMap")
    results.append(benchmark_crossmap(BAM_FILE, CHAIN_FILE, output_dir))
    
    print("[3/3] FastRemap")
    results.append(benchmark_fastremap(BAM_FILE, CHAIN_FILE, output_dir))
    
    # Save results
    output_json = RESULTS_DIR / "benchmark_bam.json"
    with open(output_json, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "input_file": str(BAM_FILE),
            "input_size_mb": input_size,
            "num_runs": NUM_RUNS,
            "results": [asdict(r) for r in results]
        }, f, indent=2)
    
    print(f"\nResults saved to: {output_json}")
    
    # Print summary
    print("\n" + "=" * 60)
    print("Benchmark Results Summary")
    print("=" * 60)
    print(f"{'Tool':<15} {'Time(s)':<12} {'Throughput(MB/s)':<15} {'Memory(MB)':<12} {'Status':<10}")
    print("-" * 60)
    
    for r in results:
        status = "✓" if r.success else "✗"
        print(f"{r.tool:<15} {r.execution_time_sec:<12.2f} {r.throughput_mb_per_sec:<15.2f} {r.peak_memory_mb:<12.1f} {status:<10}")
    
    # Calculate speedup
    fc_result = next((r for r in results if r.tool == "FastCrossMap" and r.success), None)
    cm_result = next((r for r in results if r.tool == "CrossMap" and r.success), None)
    fr_result = next((r for r in results if r.tool == "FastRemap" and r.success), None)
    
    if fc_result and cm_result:
        speedup = cm_result.execution_time_sec / fc_result.execution_time_sec
        print(f"\nFastCrossMap vs CrossMap: {speedup:.1f}x speedup")
    
    if fc_result and fr_result:
        speedup = fr_result.execution_time_sec / fc_result.execution_time_sec
        print(f"FastCrossMap vs FastRemap: {speedup:.1f}x speedup")
    
    print("\nNext step: python paper/04_plot_performance.py")


if __name__ == "__main__":
    main()
