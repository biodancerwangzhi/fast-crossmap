#!/usr/bin/env python3
"""
05b_memory_scalability.py - Test FastCrossMap memory scalability

Purpose: Demonstrate the advantage of streaming architecture - memory usage is independent of file size

Test procedure:
1. Generate subsets of different sizes from the original BAM file
2. Run FastCrossMap on each file and sample memory
3. Save results to paper/results/memory_scalability.json

Usage: python paper/05b_memory_scalability.py
Output: paper/results/memory_scalability.json
"""

import json
import subprocess
import time
import psutil
import threading
from pathlib import Path
from datetime import datetime

# =============================================================================
# 配置
# =============================================================================
RESULTS_DIR = Path("paper/results")
DATA_DIR = Path("paper/data")
TEMP_DIR = Path("paper/data/temp_scalability")

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)
TEMP_DIR.mkdir(parents=True, exist_ok=True)

# FastCrossMap executable path
FASTCROSSMAP_BIN = "target/release/fast-crossmap"

# Chain file
CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"

# Source BAM file (needs to be large enough, at least 1GB+)
SOURCE_BAM = DATA_DIR / "encode_chipseq.bam"

# Test file sizes (MB)
# Test from small to large, proving memory usage doesn't grow with file size
TEST_SIZES_MB = [50, 100, 200, 500, 1000, 2000]

# Memory sampling interval (seconds)
SAMPLE_INTERVAL = 0.5


def check_dependencies():
    """Check dependency tools"""
    print("Checking dependencies...")
    
    # Check FastCrossMap
    if not Path(FASTCROSSMAP_BIN).exists():
        print(f"Error: FastCrossMap not found: {FASTCROSSMAP_BIN}")
        print("Please build first: cargo build --release")
        return False
    
    # Check samtools
    try:
        subprocess.run(["samtools", "--version"], 
                      capture_output=True, check=True)
        print("  ✓ samtools")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: samtools not installed")
        print("Please install: conda install -c bioconda samtools")
        return False
    
    # Check Chain file
    if not CHAIN_FILE.exists():
        print(f"Error: Chain file not found: {CHAIN_FILE}")
        print("Please run first: python paper/01_download_data.sh")
        return False
    
    print("  ✓ Chain file")
    
    # Check source BAM file
    if not SOURCE_BAM.exists():
        print(f"Error: Source BAM file not found: {SOURCE_BAM}")
        print("Please run first: python paper/01_download_data.sh")
        return False
    
    print("  ✓ Source BAM file")
    
    return True


def get_file_size_mb(filepath):
    """Get file size (MB)"""
    return filepath.stat().st_size / (1024 * 1024)


def create_bam_subset(source_bam, output_bam, target_size_mb):
    """
    Create a subset of specified size from source BAM file.
    
    Args:
        source_bam: Source BAM file path
        output_bam: Output BAM file path
        target_size_mb: Target file size (MB)
    
    Returns:
        Actual file size (MB)
    """
    print(f"\nGenerating {target_size_mb} MB BAM subset...")
    
    # Get source file size
    source_size_mb = get_file_size_mb(source_bam)
    print(f"  Source file size: {source_size_mb:.2f} MB")
    
    if target_size_mb >= source_size_mb:
        print(f"  Target size >= source file, copying directly")
        subprocess.run(["cp", str(source_bam), str(output_bam)], check=True)
        return source_size_mb
    
    # Calculate extraction ratio
    ratio = target_size_mb / source_size_mb
    print(f"  Extraction ratio: {ratio:.2%}")
    
    # Use samtools view to extract subset
    # -s parameter specifies sampling ratio (needs random seed)
    seed = 42  # Fixed seed for reproducibility
    subsample_fraction = f"{seed}.{int(ratio * 100)}"
    
    cmd = [
        "samtools", "view",
        "-b",  # Output BAM format
        "-s", subsample_fraction,  # Sampling ratio
        "-o", str(output_bam),  # Output file
        str(source_bam)
    ]
    
    print(f"  Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    # Index BAM file
    print(f"  Indexing BAM file...")
    subprocess.run(["samtools", "index", str(output_bam)], check=True)
    
    actual_size_mb = get_file_size_mb(output_bam)
    print(f"  ✓ Generation complete: {actual_size_mb:.2f} MB")
    
    return actual_size_mb


def memory_sampler(process, sample_interval, results):
    """
    Memory sampling thread.
    
    Args:
        process: psutil.Process object
        sample_interval: Sampling interval (seconds)
        results: Result dictionary (for storing sampling data)
    """
    results["memory_samples"] = []
    results["sample_times"] = []
    start_time = time.time()
    
    try:
        while process.is_running() and process.status() != psutil.STATUS_ZOMBIE:
            try:
                # Get memory usage (RSS)
                mem_info = process.memory_info()
                memory_mb = mem_info.rss / (1024 * 1024)
                
                elapsed = time.time() - start_time
                results["memory_samples"].append(memory_mb)
                results["sample_times"].append(elapsed)
                
                time.sleep(sample_interval)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
    except Exception as e:
        print(f"    Memory sampling error: {e}")


def run_fastcrossmap_with_memory_profiling(input_bam, output_bam, chain_file):
    """
    Run FastCrossMap and sample memory.
    
    Args:
        input_bam: Input BAM file
        output_bam: Output BAM file
        chain_file: Chain file
    
    Returns:
        {
            "execution_time_sec": float,
            "peak_memory_mb": float,
            "memory_samples": [float],
            "sample_times": [float]
        }
    """
    print(f"  Running FastCrossMap...")
    
    # Build command
    cmd = [
        FASTCROSSMAP_BIN,
        "bam",
        str(chain_file),
        str(input_bam),
        str(output_bam)
    ]
    
    # Start process
    start_time = time.time()
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # Get psutil Process object
    ps_process = psutil.Process(process.pid)
    
    # Start memory sampling thread
    results = {}
    sampler_thread = threading.Thread(
        target=memory_sampler,
        args=(ps_process, SAMPLE_INTERVAL, results)
    )
    sampler_thread.start()
    
    # Wait for process to complete
    stdout, stderr = process.communicate()
    execution_time = time.time() - start_time
    
    # Wait for sampling thread to complete
    sampler_thread.join()
    
    # Check if successful
    if process.returncode != 0:
        print(f"    Error: FastCrossMap failed")
        print(f"    stderr: {stderr.decode()}")
        return None
    
    # Calculate peak memory
    peak_memory_mb = max(results["memory_samples"]) if results["memory_samples"] else 0
    
    print(f"  ✓ Complete: {execution_time:.2f}s, Peak memory: {peak_memory_mb:.2f} MB")
    
    return {
        "execution_time_sec": execution_time,
        "peak_memory_mb": peak_memory_mb,
        "memory_samples": results["memory_samples"],
        "sample_times": results["sample_times"]
    }


def main():
    print("=" * 60)
    print("FastCrossMap Memory Scalability Test")
    print("=" * 60)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Check dependencies
    if not check_dependencies():
        return
    
    # Test results
    test_results = []
    
    # Test for each target size
    for target_size_mb in TEST_SIZES_MB:
        print(f"\n{'=' * 60}")
        print(f"Test file size: {target_size_mb} MB")
        print(f"{'=' * 60}")
        
        # Generate BAM subset
        subset_bam = TEMP_DIR / f"subset_{target_size_mb}mb.bam"
        output_bam = TEMP_DIR / f"output_{target_size_mb}mb.bam"
        
        try:
            actual_size_mb = create_bam_subset(SOURCE_BAM, subset_bam, target_size_mb)
            
            # Run FastCrossMap and sample memory
            result = run_fastcrossmap_with_memory_profiling(
                subset_bam, output_bam, CHAIN_FILE
            )
            
            if result:
                test_results.append({
                    "target_size_mb": target_size_mb,
                    "actual_size_mb": actual_size_mb,
                    "execution_time_sec": result["execution_time_sec"],
                    "peak_memory_mb": result["peak_memory_mb"],
                    "memory_samples": result["memory_samples"],
                    "sample_times": result["sample_times"]
                })
            
            # Clean up temp files (keep subset files, delete output files)
            if output_bam.exists():
                output_bam.unlink()
            if (output_bam.parent / f"{output_bam.name}.unmap").exists():
                (output_bam.parent / f"{output_bam.name}.unmap").unlink()
        
        except Exception as e:
            print(f"  Error: {e}")
            continue
    
    # Save results
    output_file = RESULTS_DIR / "memory_scalability.json"
    
    output_data = {
        "tool": "FastCrossMap",
        "test_date": datetime.now().isoformat(),
        "source_bam": str(SOURCE_BAM),
        "chain_file": str(CHAIN_FILE),
        "sample_interval_sec": SAMPLE_INTERVAL,
        "test_results": test_results
    }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\n{'=' * 60}")
    print("Test Complete")
    print(f"{'=' * 60}")
    print(f"Results saved to: {output_file}")
    
    # Print summary
    print(f"\n{'=' * 60}")
    print("Test Summary")
    print(f"{'=' * 60}")
    print(f"{'File Size (MB)':<15} {'Exec Time (s)':<15} {'Peak Memory (MB)':<15}")
    print("-" * 45)
    
    for result in test_results:
        print(f"{result['actual_size_mb']:<15.2f} "
              f"{result['execution_time_sec']:<15.2f} "
              f"{result['peak_memory_mb']:<15.2f}")
    
    # Calculate memory usage variation
    if len(test_results) >= 2:
        min_mem = min(r["peak_memory_mb"] for r in test_results)
        max_mem = max(r["peak_memory_mb"] for r in test_results)
        mem_variation = max_mem - min_mem
        
        print(f"\nMemory usage variation:")
        print(f"  Min: {min_mem:.2f} MB")
        print(f"  Max: {max_mem:.2f} MB")
        print(f"  Variation: {mem_variation:.2f} MB ({mem_variation/min_mem*100:.1f}%)")
        
        if mem_variation < 10:
            print(f"  ✓ Memory usage is nearly constant, confirming streaming architecture effectiveness")
    
    print(f"\nNext step: python paper/06b_plot_memory_scalability.py")


if __name__ == "__main__":
    main()
