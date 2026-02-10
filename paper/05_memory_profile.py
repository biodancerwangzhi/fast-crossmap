#!/usr/bin/env python3
"""
05_memory_profile.py - Memory efficiency test

Profile memory usage for FastCrossMap, CrossMap, FastRemap
Sample memory usage (RSS) once per second, recording complete memory curves

Usage: python paper/05_memory_profile.py
Output: paper/results/memory_profile.json
"""

import json
import os
import subprocess
import time
import threading
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import List, Optional
import psutil

# =============================================================================
# 配置
# =============================================================================
DATA_DIR = Path("paper/data")
RESULTS_DIR = Path("paper/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Input files - use BAM file for memory testing
CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
BAM_FILE = DATA_DIR / "encode_chipseq.bam"

# Uncompressed chain file (required by FastRemap)
CHAIN_FILE_UNZIPPED = DATA_DIR / "hg19ToHg38.over.chain"

# Sampling interval (seconds)
SAMPLE_INTERVAL = 1.0


@dataclass
class MemoryProfile:
    """Memory profiling result"""
    tool: str
    format: str
    input_file: str
    execution_time_sec: float
    peak_memory_mb: float
    memory_samples: List[float]  # 每秒采样的内存值 (MB)
    sample_times: List[float]    # 采样时间点 (相对于开始时间的秒数)
    success: bool
    error_message: Optional[str] = None


class MemorySampler:
    """Memory sampler - periodically samples process memory in a background thread"""
    
    def __init__(self, pid: int, interval: float = 1.0):
        self.pid = pid
        self.interval = interval
        self.samples = []
        self.sample_times = []
        self.start_time = None
        self.stop_flag = threading.Event()
        self.thread = None
        
    def _sample_loop(self):
        """Sampling loop"""
        try:
            process = psutil.Process(self.pid)
            while not self.stop_flag.is_set():
                try:
                    # Get RSS (Resident Set Size) memory
                    mem_info = process.memory_info()
                    rss_mb = mem_info.rss / (1024 * 1024)
                    
                    elapsed = time.time() - self.start_time
                    self.samples.append(rss_mb)
                    self.sample_times.append(elapsed)
                    
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    # Process has ended or access denied
                    break
                
                # Wait for next sample
                self.stop_flag.wait(self.interval)
                
        except Exception as e:
            print(f"    Sampling thread error: {e}")
    
    def start(self):
        """Start sampling"""
        self.start_time = time.time()
        self.stop_flag.clear()
        self.thread = threading.Thread(target=self._sample_loop, daemon=True)
        self.thread.start()
    
    def stop(self):
        """Stop sampling"""
        self.stop_flag.set()
        if self.thread:
            self.thread.join(timeout=2.0)
    
    def get_results(self):
        """Get sampling results"""
        return self.samples, self.sample_times


def run_with_memory_profile(cmd: List[str], output_file: Path) -> tuple[float, List[float], List[float], bool, str]:
    """
    Run command and sample memory.
    Returns: (execution_time_sec, memory_samples_mb, sample_times, success, error_message)
    """
    try:
        # Start process
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Start memory sampler
        sampler = MemorySampler(process.pid, interval=SAMPLE_INTERVAL)
        sampler.start()
        
        start_time = time.time()
        
        # Wait for process to complete
        try:
            stdout, stderr = process.communicate(timeout=1800)  # 30 minute timeout
            elapsed = time.time() - start_time
            
            # Stop sampler
            sampler.stop()
            
            # Get sampling results
            memory_samples, sample_times = sampler.get_results()
            
            if process.returncode != 0:
                return elapsed, memory_samples, sample_times, False, stderr[:500]
            
            return elapsed, memory_samples, sample_times, True, ""
            
        except subprocess.TimeoutExpired:
            process.kill()
            sampler.stop()
            return 1800, [], [], False, "Timeout after 1800 seconds"
            
    except Exception as e:
        return 0, [], [], False, str(e)


def profile_fastcrossmap(bam_file: Path, chain_file: Path, output_dir: Path) -> MemoryProfile:
    """Profile FastCrossMap memory usage"""
    print("  Profiling FastCrossMap...")
    output_file = output_dir / "fastcrossmap_memory.bam"
    
    cmd = [
        "./target/release/fast-crossmap",
        "bam",
        str(chain_file),
        str(bam_file),
        str(output_file)
    ]
    
    elapsed, mem_samples, sample_times, success, error_msg = run_with_memory_profile(cmd, output_file)
    
    if not success or not mem_samples:
        return MemoryProfile(
            tool="FastCrossMap",
            format="BAM",
            input_file=str(bam_file),
            execution_time_sec=elapsed,
            peak_memory_mb=0,
            memory_samples=[],
            sample_times=[],
            success=False,
            error_message=error_msg if not success else "No memory samples collected"
        )
    
    peak_memory = max(mem_samples)
    
    return MemoryProfile(
        tool="FastCrossMap",
        format="BAM",
        input_file=str(bam_file),
        execution_time_sec=round(elapsed, 2),
        peak_memory_mb=round(peak_memory, 2),
        memory_samples=[round(m, 2) for m in mem_samples],
        sample_times=[round(t, 2) for t in sample_times],
        success=True
    )


def profile_crossmap(bam_file: Path, chain_file: Path, output_dir: Path) -> MemoryProfile:
    """Profile CrossMap memory usage"""
    print("  Profiling CrossMap...")
    output_file = output_dir / "crossmap_memory.bam"
    
    cmd = [
        "CrossMap", "bam",
        str(chain_file),
        str(bam_file),
        str(output_file)
    ]
    
    elapsed, mem_samples, sample_times, success, error_msg = run_with_memory_profile(cmd, output_file)
    
    if not success or not mem_samples:
        return MemoryProfile(
            tool="CrossMap",
            format="BAM",
            input_file=str(bam_file),
            execution_time_sec=elapsed,
            peak_memory_mb=0,
            memory_samples=[],
            sample_times=[],
            success=False,
            error_message=error_msg if not success else "No memory samples collected"
        )
    
    peak_memory = max(mem_samples)
    
    return MemoryProfile(
        tool="CrossMap",
        format="BAM",
        input_file=str(bam_file),
        execution_time_sec=round(elapsed, 2),
        peak_memory_mb=round(peak_memory, 2),
        memory_samples=[round(m, 2) for m in mem_samples],
        sample_times=[round(t, 2) for t in sample_times],
        success=True
    )


def profile_fastremap(bam_file: Path, chain_file: Path, output_dir: Path) -> MemoryProfile:
    """Profile FastRemap memory usage"""
    print("  Profiling FastRemap...")
    output_file = output_dir / "fastremap_memory.bam"
    unmap_file = output_dir / "fastremap_memory.bam.unmap"
    
    # FastRemap does not support .gz, needs uncompressed chain file
    chain_unzipped = CHAIN_FILE_UNZIPPED
    if not chain_unzipped.exists():
        print("    Decompressing chain file for FastRemap...")
        subprocess.run(["gunzip", "-k", str(chain_file)], check=True)
    
    cmd = [
        "FastRemap",
        "-f", "bam",
        "-c", str(chain_unzipped),
        "-i", str(bam_file),
        "-o", str(output_file),
        "-u", str(unmap_file)
    ]
    
    elapsed, mem_samples, sample_times, success, error_msg = run_with_memory_profile(cmd, output_file)
    
    if not success or not mem_samples:
        return MemoryProfile(
            tool="FastRemap",
            format="BAM",
            input_file=str(bam_file),
            execution_time_sec=elapsed,
            peak_memory_mb=0,
            memory_samples=[],
            sample_times=[],
            success=False,
            error_message=error_msg if not success else "No memory samples collected"
        )
    
    peak_memory = max(mem_samples)
    
    return MemoryProfile(
        tool="FastRemap",
        format="BAM",
        input_file=str(bam_file),
        execution_time_sec=round(elapsed, 2),
        peak_memory_mb=round(peak_memory, 2),
        memory_samples=[round(m, 2) for m in mem_samples],
        sample_times=[round(t, 2) for t in sample_times],
        success=True
    )


def main():
    print("=" * 60)
    print("Memory Efficiency Test")
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
    
    # Check if psutil is installed
    try:
        import psutil
    except ImportError:
        print("Error: psutil library is required")
        print("Please run: pip install psutil")
        return
    
    # Get file size
    file_size_mb = BAM_FILE.stat().st_size / (1024 * 1024)
    print(f"\nInput file: {BAM_FILE}")
    print(f"File size: {file_size_mb:.2f} MB")
    print(f"Sampling interval: {SAMPLE_INTERVAL} seconds")
    print()
    
    # Create output directory
    output_dir = RESULTS_DIR / "memory_profile"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run memory profiling
    results = []
    
    print("[1/3] FastCrossMap")
    results.append(profile_fastcrossmap(BAM_FILE, CHAIN_FILE, output_dir))
    
    print("[2/3] CrossMap")
    results.append(profile_crossmap(BAM_FILE, CHAIN_FILE, output_dir))
    
    print("[3/3] FastRemap")
    results.append(profile_fastremap(BAM_FILE, CHAIN_FILE, output_dir))
    
    # Save results
    output_json = RESULTS_DIR / "memory_profile.json"
    with open(output_json, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "input_file": str(BAM_FILE),
            "input_size_mb": round(file_size_mb, 2),
            "sample_interval_sec": SAMPLE_INTERVAL,
            "results": [asdict(r) for r in results]
        }, f, indent=2)
    
    print(f"\nResults saved to: {output_json}")
    
    # Print summary
    print("\n" + "=" * 60)
    print("Memory Profiling Results Summary")
    print("=" * 60)
    print(f"{'Tool':<15} {'Time(s)':<10} {'Peak Memory(MB)':<15} {'Samples':<10} {'Status':<10}")
    print("-" * 60)
    
    for r in results:
        status = "✓" if r.success else "✗"
        sample_count = len(r.memory_samples)
        print(f"{r.tool:<15} {r.execution_time_sec:<10.2f} {r.peak_memory_mb:<15.2f} {sample_count:<10} {status:<10}")
    
    # Print memory efficiency comparison
    print("\n" + "=" * 60)
    print("Memory Efficiency Comparison")
    print("=" * 60)
    
    successful_results = [r for r in results if r.success]
    if len(successful_results) >= 2:
        # 找到 FastCrossMap 和 CrossMap
        fc = next((r for r in successful_results if r.tool == "FastCrossMap"), None)
        cm = next((r for r in successful_results if r.tool == "CrossMap"), None)
        
        if fc and cm:
            mem_ratio = cm.peak_memory_mb / fc.peak_memory_mb
            print(f"FastCrossMap peak memory: {fc.peak_memory_mb:.2f} MB")
            print(f"CrossMap peak memory: {cm.peak_memory_mb:.2f} MB")
            print(f"Memory savings: {mem_ratio:.2f}x (CrossMap uses {mem_ratio:.1f}x memory)")
    
    print("\nNext step: python paper/06_plot_memory.py")


if __name__ == "__main__":
    main()
