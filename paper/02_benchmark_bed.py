#!/usr/bin/env python3
"""
02_benchmark_bed.py - BED format benchmark

Run 4-way benchmark: FastCrossMap, CrossMap, liftOver, FastRemap
Measure execution time, throughput, memory usage

Usage: python paper/02_benchmark_bed.py
Output: paper/results/benchmark_bed.json
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
BED_FILE = DATA_DIR / "encode_dnase_peaks.bed.gz"

# Uncompressed chain file (required by FastRemap)
CHAIN_FILE_UNZIPPED = DATA_DIR / "hg19ToHg38.over.chain"

# Number of runs
NUM_RUNS = 5


@dataclass
class BenchmarkResult:
    """Benchmark result"""
    tool: str
    format: str
    input_file: str
    input_records: int
    execution_time_sec: float
    throughput_rec_per_sec: float
    peak_memory_mb: float
    mapped_records: int
    unmapped_records: int
    all_times: list  # 所有运行的时间
    success: bool
    error_message: Optional[str] = None


def count_bed_records(bed_file: Path) -> int:
    """Count BED file records (supports .gz compression)"""
    import gzip
    
    count = 0
    if str(bed_file).endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    with opener(bed_file, mode) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                count += 1
    return count


def count_output_records(output_file: Path) -> tuple[int, int]:
    """Count mapped and unmapped records in output file"""
    mapped = 0
    unmapped = 0
    
    if output_file.exists():
        with open(output_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    mapped += 1
    
    unmap_file = Path(str(output_file) + ".unmap")
    if unmap_file.exists():
        with open(unmap_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    unmapped += 1
    
    return mapped, unmapped


def run_with_time(cmd: list[str], output_file: Path) -> tuple[float, float, bool, str]:
    """
    Run command and measure time and memory.
    Returns: (execution_time_sec, peak_memory_mb, success, error_message)
    """
    import resource
    import os
    
    start_time = time.time()
    try:
        # Run command directly, measure with Python time
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        elapsed = time.time() - start_time
        
        # Try to get memory usage via /usr/bin/time, fallback to 0
        peak_memory_mb = 0
        try:
            # Try running with time -v to get memory info
            time_result = subprocess.run(
                ["/usr/bin/time", "-v"] + cmd,
                capture_output=True,
                text=True,
                timeout=600
            )
            for line in time_result.stderr.split('\n'):
                if 'Maximum resident set size' in line:
                    peak_memory_kb = int(line.split(':')[1].strip())
                    peak_memory_mb = peak_memory_kb / 1024
                    break
        except:
            # If /usr/bin/time is unavailable, try psutil or set to 0
            try:
                import psutil
                # Rough estimate: use current process memory as reference
                peak_memory_mb = psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)
            except:
                peak_memory_mb = 0
        
        if result.returncode != 0:
            return elapsed, peak_memory_mb, False, result.stderr[:500]
        
        return elapsed, peak_memory_mb, True, ""
        
    except subprocess.TimeoutExpired:
        return 600, 0, False, "Timeout after 600 seconds"
    except Exception as e:
        return 0, 0, False, str(e)


def benchmark_fastcrossmap(bed_file: Path, chain_file: Path, output_dir: Path) -> BenchmarkResult:
    """Benchmark FastCrossMap"""
    print("  Running FastCrossMap...")
    output_file = output_dir / "fastcrossmap_output.bed"
    
    # FastCrossMap uses positional args: fast-crossmap bed <CHAIN> <INPUT> [OUTPUT]
    cmd = [
        "./fast-crossmap-linux-x64/fast-crossmap",
        "bed",
        str(chain_file),
        str(bed_file),
        str(output_file)
    ]
    
    times = []
    memories = []
    success = False
    error_msg = ""
    
    for i in range(NUM_RUNS):
        elapsed, memory, ok, err = run_with_time(cmd, output_file)
        if ok:
            times.append(elapsed)
            memories.append(memory)
            success = True
        else:
            error_msg = err
    
    if not times:
        return BenchmarkResult(
            tool="FastCrossMap",
            format="BED",
            input_file=str(bed_file),
            input_records=count_bed_records(bed_file),
            execution_time_sec=0,
            throughput_rec_per_sec=0,
            peak_memory_mb=0,
            mapped_records=0,
            unmapped_records=0,
            all_times=[],
            success=False,
            error_message=error_msg
        )
    
    avg_time = sum(times) / len(times)
    avg_memory = sum(memories) / len(memories)
    input_records = count_bed_records(bed_file)
    mapped, unmapped = count_output_records(output_file)
    
    return BenchmarkResult(
        tool="FastCrossMap",
        format="BED",
        input_file=str(bed_file),
        input_records=input_records,
        execution_time_sec=round(avg_time, 3),
        throughput_rec_per_sec=round(input_records / avg_time, 0),
        peak_memory_mb=round(avg_memory, 2),
        mapped_records=mapped,
        unmapped_records=unmapped,
        all_times=times,
        success=success
    )


def benchmark_crossmap(bed_file: Path, chain_file: Path, output_dir: Path) -> BenchmarkResult:
    """Benchmark CrossMap"""
    print("  Running CrossMap...")
    output_file = output_dir / "crossmap_output.bed"
    
    cmd = [
        "CrossMap", "bed",
        str(chain_file),
        str(bed_file),
        str(output_file)
    ]
    
    times = []
    memories = []
    success = False
    error_msg = ""
    
    for i in range(NUM_RUNS):
        elapsed, memory, ok, err = run_with_time(cmd, output_file)
        if ok:
            times.append(elapsed)
            memories.append(memory)
            success = True
        else:
            error_msg = err
    
    if not times:
        return BenchmarkResult(
            tool="CrossMap",
            format="BED",
            input_file=str(bed_file),
            input_records=count_bed_records(bed_file),
            execution_time_sec=0,
            throughput_rec_per_sec=0,
            peak_memory_mb=0,
            mapped_records=0,
            unmapped_records=0,
            all_times=[],
            success=False,
            error_message=error_msg
        )
    
    avg_time = sum(times) / len(times)
    avg_memory = sum(memories) / len(memories)
    input_records = count_bed_records(bed_file)
    mapped, unmapped = count_output_records(output_file)
    
    return BenchmarkResult(
        tool="CrossMap",
        format="BED",
        input_file=str(bed_file),
        input_records=input_records,
        execution_time_sec=round(avg_time, 3),
        throughput_rec_per_sec=round(input_records / avg_time, 0),
        peak_memory_mb=round(avg_memory, 2),
        mapped_records=mapped,
        unmapped_records=unmapped,
        all_times=times,
        success=success
    )


def benchmark_liftover(bed_file: Path, chain_file: Path, output_dir: Path) -> BenchmarkResult:
    """Benchmark UCSC liftOver"""
    print("  Running liftOver...")
    output_file = output_dir / "liftover_output.bed"
    unmap_file = output_dir / "liftover_output.bed.unmap"
    
    # liftOver has strict BED format requirements; narrowPeak format column 7 float values cause errors
    # Need to convert to standard BED6 format first
    import gzip
    
    bed6_file = output_dir / "input_bed6.bed"
    
    # Choose open method based on file extension
    if str(bed_file).endswith('.gz'):
        fin = gzip.open(bed_file, 'rt')
    else:
        fin = open(bed_file, 'r')
    
    try:
        with open(bed6_file, 'w') as fout:
            for line in fin:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 6:
                    # Take only first 6 columns: chrom, start, end, name, score, strand
                    # Convert float score to integer
                    try:
                        score = int(float(fields[4])) if len(fields) > 4 else 0
                    except:
                        score = 0
                    strand = fields[5] if len(fields) > 5 else '.'
                    fout.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\t{score}\t{strand}\n")
                elif len(fields) >= 3:
                    # Minimum 3 columns
                    fout.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t.\t0\t.\n")
    finally:
        fin.close()
    
    cmd = [
        "liftOver",
        str(bed6_file),
        str(chain_file),
        str(output_file),
        str(unmap_file)
    ]
    
    times = []
    memories = []
    success = False
    error_msg = ""
    
    for i in range(NUM_RUNS):
        elapsed, memory, ok, err = run_with_time(cmd, output_file)
        if ok:
            times.append(elapsed)
            memories.append(memory)
            success = True
        else:
            error_msg = err
    
    if not times:
        return BenchmarkResult(
            tool="liftOver",
            format="BED",
            input_file=str(bed_file),
            input_records=count_bed_records(bed_file),
            execution_time_sec=0,
            throughput_rec_per_sec=0,
            peak_memory_mb=0,
            mapped_records=0,
            unmapped_records=0,
            all_times=[],
            success=False,
            error_message=error_msg
        )
    
    avg_time = sum(times) / len(times)
    avg_memory = sum(memories) / len(memories)
    input_records = count_bed_records(bed_file)
    mapped, unmapped = count_output_records(output_file)
    
    return BenchmarkResult(
        tool="liftOver",
        format="BED",
        input_file=str(bed_file),
        input_records=input_records,
        execution_time_sec=round(avg_time, 3),
        throughput_rec_per_sec=round(input_records / avg_time, 0),
        peak_memory_mb=round(avg_memory, 2),
        mapped_records=mapped,
        unmapped_records=unmapped,
        all_times=times,
        success=success
    )


def benchmark_fastremap(bed_file: Path, chain_file: Path, output_dir: Path) -> BenchmarkResult:
    """Benchmark FastRemap"""
    print("  Running FastRemap...")
    output_file = output_dir / "fastremap_output.bed"
    unmap_file = output_dir / "fastremap_output.bed.unmap"
    
    # FastRemap does not support .gz, needs uncompressed chain file
    chain_unzipped = CHAIN_FILE_UNZIPPED
    if not chain_unzipped.exists():
        print("    Decompressing chain file for FastRemap...")
        subprocess.run(["gunzip", "-k", str(chain_file)], check=True)
    
    cmd = [
        "FastRemap",
        "-f", "bed",
        "-c", str(chain_unzipped),
        "-i", str(bed_file),
        "-o", str(output_file),
        "-u", str(unmap_file)
    ]
    
    times = []
    memories = []
    success = False
    error_msg = ""
    
    for i in range(NUM_RUNS):
        elapsed, memory, ok, err = run_with_time(cmd, output_file)
        if ok:
            times.append(elapsed)
            memories.append(memory)
            success = True
        else:
            error_msg = err
    
    if not times:
        return BenchmarkResult(
            tool="FastRemap",
            format="BED",
            input_file=str(bed_file),
            input_records=count_bed_records(bed_file),
            execution_time_sec=0,
            throughput_rec_per_sec=0,
            peak_memory_mb=0,
            mapped_records=0,
            unmapped_records=0,
            all_times=[],
            success=False,
            error_message=error_msg
        )
    
    avg_time = sum(times) / len(times)
    avg_memory = sum(memories) / len(memories)
    input_records = count_bed_records(bed_file)
    mapped, unmapped = count_output_records(output_file)
    
    return BenchmarkResult(
        tool="FastRemap",
        format="BED",
        input_file=str(bed_file),
        input_records=input_records,
        execution_time_sec=round(avg_time, 3),
        throughput_rec_per_sec=round(input_records / avg_time, 0),
        peak_memory_mb=round(avg_memory, 2),
        mapped_records=mapped,
        unmapped_records=unmapped,
        all_times=times,
        success=success
    )


def main():
    print("=" * 60)
    print("BED Format Benchmark")
    print("=" * 60)
    
    # Check input files
    if not BED_FILE.exists():
        print(f"Error: BED file not found: {BED_FILE}")
        print("Please run first: bash paper/01_download_data.sh")
        return
    
    if not CHAIN_FILE.exists():
        print(f"Error: Chain file not found: {CHAIN_FILE}")
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
    
    # Count input records
    input_records = count_bed_records(BED_FILE)
    print(f"\nInput file: {bed_file_to_use}")
    print(f"Records: {input_records:,}")
    print(f"Number of runs: {NUM_RUNS}")
    print()
    
    # Create output directory
    output_dir = RESULTS_DIR / "bed_benchmark"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run benchmarks (using decompressed file)
    results = []
    
    print("[1/4] FastCrossMap")
    results.append(benchmark_fastcrossmap(bed_file_to_use, CHAIN_FILE, output_dir))
    
    print("[2/4] CrossMap")
    results.append(benchmark_crossmap(bed_file_to_use, CHAIN_FILE, output_dir))
    
    print("[3/4] liftOver")
    results.append(benchmark_liftover(bed_file_to_use, CHAIN_FILE, output_dir))
    
    print("[4/4] FastRemap")
    results.append(benchmark_fastremap(bed_file_to_use, CHAIN_FILE, output_dir))
    
    # Save results
    output_json = RESULTS_DIR / "benchmark_bed.json"
    with open(output_json, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "input_file": str(BED_FILE),
            "input_records": input_records,
            "num_runs": NUM_RUNS,
            "results": [asdict(r) for r in results]
        }, f, indent=2)
    
    print(f"\nResults saved to: {output_json}")
    
    # Print summary
    print("\n" + "=" * 60)
    print("Benchmark Results Summary")
    print("=" * 60)
    print(f"{'Tool':<15} {'Time(s)':<10} {'Throughput(rec/s)':<15} {'Memory(MB)':<10} {'Status':<10}")
    print("-" * 60)
    
    for r in results:
        status = "✓" if r.success else "✗"
        print(f"{r.tool:<15} {r.execution_time_sec:<10.3f} {r.throughput_rec_per_sec:<15,.0f} {r.peak_memory_mb:<10.1f} {status:<10}")
    
    print("\nNext step: python paper/03_benchmark_bam.py")


if __name__ == "__main__":
    main()
