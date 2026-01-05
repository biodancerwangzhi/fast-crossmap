#!/usr/bin/env python3
"""
Comprehensive 4-way benchmark comparison for genome coordinate conversion tools.

Compares: FastCrossMap, CrossMap (Python), UCSC liftOver, FastRemap

Usage:
    python scripts/run_benchmarks.py --chain <chain_file> --input <input_file> --format <bed|bam>
    python scripts/run_benchmarks.py --generate-test-data  # Generate test data first
    python scripts/run_benchmarks.py --quick-test  # Run with small test data

Features:
    - Cold start vs warm start timing
    - Peak RSS memory measurement
    - Throughput calculation (records/sec)
    - JSON and CSV output
"""

import argparse
import csv
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Tool-format support matrix
TOOL_FORMAT_SUPPORT = {
    "fastcrossmap": ["bed", "vcf", "gff", "bam", "sam", "wig", "bigwig", "maf", "gvcf", "region"],
    "crossmap": ["bed", "vcf", "gff", "bam", "sam", "wig", "bigwig", "maf", "gvcf", "region"],
    "liftover": ["bed"],  # UCSC liftOver only supports BED
    "fastremap": ["bed", "bam"],  # FastRemap supports BED and BAM
}

# Command templates for each tool
TOOL_COMMANDS = {
    "fastcrossmap": "./target/release/fast-crossmap {format} {chain} {input} {output} -t {threads}",
    "crossmap": "CrossMap {format} {chain} {input} {output}",
    "liftover": "liftOver {input} {chain} {output} {unmapped}",
    "fastremap": "FastRemap -f {format} -c {chain} -i {input} -u {unmapped} -o {output}",
}

@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""
    tool: str
    format: str
    input_records: int
    execution_time_sec: float
    cold_start_time_sec: float
    warm_start_time_sec: float
    peak_rss_mb: float
    throughput_records_per_sec: float
    exit_code: int
    supported: bool
    error_message: str = ""
    output_records: int = 0
    unmapped_records: int = 0

@dataclass
class BenchmarkMetadata:
    """Metadata for benchmark run."""
    timestamp: str
    system: Dict[str, str]
    tool_versions: Dict[str, str]
    input_file: str
    input_file_size_mb: float
    chain_file: str
    format: str
    threads: int

@dataclass
class BenchmarkReport:
    """Complete benchmark report."""
    metadata: BenchmarkMetadata
    results: List[BenchmarkResult]


class BenchmarkRunner:
    """Runs benchmarks for genome coordinate conversion tools."""
    
    def __init__(self, chain_file: Path, threads: int = 4, output_dir: Path = Path("results")):
        self.chain_file = chain_file
        self.threads = threads
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # FastRemap doesn't support .gz files, prepare uncompressed version
        self.chain_file_uncompressed = self._prepare_uncompressed_chain()
    
    def _prepare_uncompressed_chain(self) -> Path:
        """Prepare uncompressed chain file for tools that don't support .gz (e.g., FastRemap)."""
        if not str(self.chain_file).endswith('.gz'):
            return self.chain_file
        
        # Check if uncompressed version exists
        uncompressed = Path(str(self.chain_file)[:-3])  # Remove .gz
        if uncompressed.exists():
            return uncompressed
        
        # Try to decompress
        print(f"Note: Decompressing chain file for FastRemap (doesn't support .gz)...")
        try:
            import gzip
            import shutil
            with gzip.open(self.chain_file, 'rb') as f_in:
                with open(uncompressed, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"  Created: {uncompressed}")
            return uncompressed
        except Exception as e:
            print(f"  Warning: Could not decompress chain file: {e}")
            return self.chain_file
        
    def get_tool_version(self, tool: str) -> str:
        """Get version string for a tool."""
        try:
            if tool == "fastcrossmap":
                result = subprocess.run(
                    ["./target/release/fast-crossmap", "--version"],
                    capture_output=True, text=True, timeout=10
                )
                return result.stdout.strip() or result.stderr.strip() or "unknown"
            elif tool == "crossmap":
                result = subprocess.run(
                    ["CrossMap", "--version"],
                    capture_output=True, text=True, timeout=10
                )
                return result.stdout.strip() or result.stderr.strip() or "unknown"
            elif tool == "liftover":
                result = subprocess.run(
                    ["liftOver"],
                    capture_output=True, text=True, timeout=10
                )
                # liftOver prints version in usage message
                for line in result.stderr.split('\n'):
                    if 'liftOver' in line.lower():
                        return line.strip()
                return "unknown"
            elif tool == "fastremap":
                result = subprocess.run(
                    ["FastRemap", "--version"],
                    capture_output=True, text=True, timeout=10
                )
                return result.stdout.strip() or result.stderr.strip() or "unknown"
        except Exception as e:
            return f"error: {e}"
        return "unknown"
    
    def get_system_info(self) -> Dict[str, str]:
        """Get system information."""
        return {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "python_version": platform.python_version(),
            "cpu_count": str(os.cpu_count()),
        }
    
    def count_records(self, file_path: Path, format: str) -> int:
        """Count records in input file."""
        if not file_path.exists():
            return 0
        count = 0
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if not line.startswith('#') and line.strip():
                        count += 1
        except:
            pass
        return count
    
    def run_with_time(self, cmd: str, timeout: int = 1800) -> Tuple[float, float, int]:
        """
        Run command and measure execution time and Peak RSS.
        Returns: (execution_time_sec, peak_rss_mb, exit_code)
        
        Uses Python's time module for timing and psutil for memory (if available).
        Falls back to /usr/bin/time on Linux if available.
        """
        # Try to use psutil for cross-platform memory measurement
        try:
            import psutil
            return self._run_with_psutil(cmd, timeout)
        except ImportError:
            pass
        
        # Try GNU time on Linux
        if shutil.which("time"):
            return self._run_with_gnu_time(cmd, timeout)
        
        # Fallback: just measure time, no memory
        return self._run_simple(cmd, timeout)
    
    def _run_with_psutil(self, cmd: str, timeout: int) -> Tuple[float, float, int]:
        """Run with psutil for memory tracking."""
        import psutil
        
        start_time = time.time()
        peak_rss_mb = 0.0
        
        try:
            process = subprocess.Popen(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            
            try:
                ps_process = psutil.Process(process.pid)
                while process.poll() is None:
                    try:
                        # Get memory of process and all children
                        mem_info = ps_process.memory_info()
                        current_rss = mem_info.rss / (1024 * 1024)
                        
                        # Also check children
                        for child in ps_process.children(recursive=True):
                            try:
                                child_mem = child.memory_info()
                                current_rss += child_mem.rss / (1024 * 1024)
                            except (psutil.NoSuchProcess, psutil.AccessDenied):
                                pass
                        
                        peak_rss_mb = max(peak_rss_mb, current_rss)
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                    time.sleep(0.01)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
            
            stdout, stderr = process.communicate(timeout=timeout)
            elapsed = time.time() - start_time
            
            return elapsed, peak_rss_mb, process.returncode
            
        except subprocess.TimeoutExpired:
            process.kill()
            return timeout, peak_rss_mb, -1
        except Exception as e:
            print(f"    Error: {e}")
            return 0, 0, -1
    
    def _run_with_gnu_time(self, cmd: str, timeout: int) -> Tuple[float, float, int]:
        """Run with GNU time for memory measurement."""
        # Use 'time' command (might be shell builtin or GNU time)
        time_cmd = f"command time -v {cmd} 2>&1"
        
        start_time = time.time()
        try:
            result = subprocess.run(
                time_cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=timeout,
                executable='/bin/bash'
            )
            elapsed = time.time() - start_time
            
            # Parse Peak RSS from time output
            peak_rss_kb = 0
            output = result.stdout + result.stderr
            for line in output.split('\n'):
                if 'Maximum resident set size' in line:
                    match = re.search(r'(\d+)', line)
                    if match:
                        peak_rss_kb = int(match.group(1))
                        break
            
            peak_rss_mb = peak_rss_kb / 1024.0
            return elapsed, peak_rss_mb, result.returncode
            
        except subprocess.TimeoutExpired:
            return timeout, 0, -1
        except Exception as e:
            print(f"    Error: {e}")
            return 0, 0, -1
    
    def _run_simple(self, cmd: str, timeout: int) -> Tuple[float, float, int]:
        """Simple run without memory measurement."""
        start_time = time.time()
        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                timeout=timeout
            )
            elapsed = time.time() - start_time
            return elapsed, 0.0, result.returncode
        except subprocess.TimeoutExpired:
            return timeout, 0, -1
        except Exception as e:
            print(f"    Error: {e}")
            return 0, 0, -1
    
    def build_command(self, tool: str, format: str, input_file: Path, 
                      output_file: Path, unmapped_file: Path) -> str:
        """Build command string for a tool."""
        template = TOOL_COMMANDS.get(tool, "")
        
        # FastRemap doesn't support .gz chain files
        chain_file = self.chain_file_uncompressed if tool == "fastremap" else self.chain_file
        
        return template.format(
            format=format,
            chain=str(chain_file),
            input=str(input_file),
            output=str(output_file),
            unmapped=str(unmapped_file),
            threads=self.threads
        )
    
    def run_benchmark(self, tool: str, format: str, input_file: Path,
                      num_runs: int = 3) -> BenchmarkResult:
        """Run benchmark for a single tool."""
        
        # Check if tool supports this format
        if format not in TOOL_FORMAT_SUPPORT.get(tool, []):
            return BenchmarkResult(
                tool=tool,
                format=format,
                input_records=0,
                execution_time_sec=0,
                cold_start_time_sec=0,
                warm_start_time_sec=0,
                peak_rss_mb=0,
                throughput_records_per_sec=0,
                exit_code=-1,
                supported=False,
                error_message="Format not supported"
            )
        
        # Count input records
        input_records = self.count_records(input_file, format)
        
        # Setup output files
        output_file = self.output_dir / f"{tool}_output.{format}"
        unmapped_file = self.output_dir / f"{tool}_output.{format}.unmap"
        
        # Build command
        cmd = self.build_command(tool, format, input_file, output_file, unmapped_file)
        
        print(f"  Running {tool}...")
        print(f"    Command: {cmd}")
        
        # Cold start run (first run)
        cold_time, cold_rss, cold_exit = self.run_with_time(cmd)
        
        if cold_exit != 0:
            return BenchmarkResult(
                tool=tool,
                format=format,
                input_records=input_records,
                execution_time_sec=cold_time,
                cold_start_time_sec=cold_time,
                warm_start_time_sec=0,
                peak_rss_mb=cold_rss,
                throughput_records_per_sec=0,
                exit_code=cold_exit,
                supported=True,
                error_message=f"Command failed with exit code {cold_exit}"
            )
        
        # Warm start runs (subsequent runs)
        warm_times = []
        warm_rss_values = []
        
        for i in range(num_runs - 1):
            warm_time, warm_rss, warm_exit = self.run_with_time(cmd)
            if warm_exit == 0:
                warm_times.append(warm_time)
                warm_rss_values.append(warm_rss)
        
        # Calculate averages
        avg_warm_time = sum(warm_times) / len(warm_times) if warm_times else cold_time
        avg_rss = sum(warm_rss_values) / len(warm_rss_values) if warm_rss_values else cold_rss
        
        # Use warm start time as the representative execution time
        exec_time = avg_warm_time
        throughput = input_records / exec_time if exec_time > 0 else 0
        
        # Count output records
        output_records = self.count_records(output_file, format)
        unmapped_records = self.count_records(unmapped_file, format)
        
        return BenchmarkResult(
            tool=tool,
            format=format,
            input_records=input_records,
            execution_time_sec=exec_time,
            cold_start_time_sec=cold_time,
            warm_start_time_sec=avg_warm_time,
            peak_rss_mb=avg_rss,
            throughput_records_per_sec=throughput,
            exit_code=0,
            supported=True,
            output_records=output_records,
            unmapped_records=unmapped_records
        )
    
    def run_all_benchmarks(self, input_file: Path, format: str,
                           tools: List[str] = None) -> BenchmarkReport:
        """Run benchmarks for all tools."""
        if tools is None:
            tools = ["fastcrossmap", "crossmap", "liftover", "fastremap"]
        
        print(f"\n{'='*60}")
        print(f"Running 4-way benchmark comparison")
        print(f"{'='*60}")
        print(f"Chain file: {self.chain_file}")
        print(f"Input file: {input_file}")
        print(f"Format: {format}")
        print(f"Threads: {self.threads}")
        print(f"{'='*60}\n")
        
        # Collect metadata
        metadata = BenchmarkMetadata(
            timestamp=datetime.now().isoformat(),
            system=self.get_system_info(),
            tool_versions={tool: self.get_tool_version(tool) for tool in tools},
            input_file=str(input_file),
            input_file_size_mb=input_file.stat().st_size / (1024 * 1024) if input_file.exists() else 0,
            chain_file=str(self.chain_file),
            format=format,
            threads=self.threads
        )
        
        # Run benchmarks
        results = []
        for tool in tools:
            result = self.run_benchmark(tool, format, input_file)
            results.append(result)
            
            # Print result summary
            if result.supported:
                print(f"    Cold start: {result.cold_start_time_sec:.3f}s")
                print(f"    Warm start: {result.warm_start_time_sec:.3f}s")
                print(f"    Peak RSS: {result.peak_rss_mb:.1f} MB")
                print(f"    Throughput: {result.throughput_records_per_sec:.0f} records/sec")
            else:
                print(f"    Status: Not Supported")
            print()
        
        return BenchmarkReport(metadata=metadata, results=results)


def export_to_json(report: BenchmarkReport, output_path: Path):
    """Export benchmark report to JSON."""
    data = {
        "metadata": asdict(report.metadata),
        "results": [asdict(r) for r in report.results]
    }
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"JSON report saved to: {output_path}")


def export_to_csv(report: BenchmarkReport, output_path: Path):
    """Export benchmark results to CSV."""
    if not report.results:
        return
    
    fieldnames = list(asdict(report.results[0]).keys())
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for result in report.results:
            writer.writerow(asdict(result))
    print(f"CSV report saved to: {output_path}")


def generate_test_data(output_dir: Path, num_records: int = 10000) -> Tuple[Path, Path]:
    """Generate test BED data for benchmarking."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    bed_file = output_dir / "benchmark_test.bed"
    
    print(f"Generating {num_records} test BED records...")
    
    # Generate random BED records for human chromosomes
    import random
    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    
    with open(bed_file, 'w') as f:
        for i in range(num_records):
            chrom = random.choice(chromosomes)
            start = random.randint(1000000, 200000000)
            end = start + random.randint(100, 10000)
            name = f"region_{i}"
            score = random.randint(0, 1000)
            strand = random.choice(['+', '-'])
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
    
    print(f"Test data generated: {bed_file}")
    return bed_file


def print_summary(report: BenchmarkReport):
    """Print a summary table of benchmark results."""
    print(f"\n{'='*80}")
    print("BENCHMARK SUMMARY")
    print(f"{'='*80}")
    print(f"{'Tool':<15} {'Status':<12} {'Time (s)':<12} {'Memory (MB)':<12} {'Throughput':<15}")
    print(f"{'-'*80}")
    
    for r in report.results:
        status = "✓ OK" if r.supported and r.exit_code == 0 else "✗ N/A" if not r.supported else "✗ FAIL"
        time_str = f"{r.execution_time_sec:.3f}" if r.supported else "-"
        mem_str = f"{r.peak_rss_mb:.1f}" if r.supported else "-"
        tp_str = f"{r.throughput_records_per_sec:.0f} rec/s" if r.supported else "-"
        print(f"{r.tool:<15} {status:<12} {time_str:<12} {mem_str:<12} {tp_str:<15}")
    
    print(f"{'='*80}")
    
    # Calculate speedups relative to CrossMap
    crossmap_result = next((r for r in report.results if r.tool == "crossmap" and r.supported), None)
    if crossmap_result and crossmap_result.execution_time_sec > 0:
        print("\nSpeedup vs CrossMap:")
        for r in report.results:
            if r.supported and r.execution_time_sec > 0 and r.tool != "crossmap":
                speedup = crossmap_result.execution_time_sec / r.execution_time_sec
                print(f"  {r.tool}: {speedup:.1f}x")


def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive 4-way benchmark for genome coordinate conversion tools"
    )
    parser.add_argument("--chain", type=Path, 
                        default=Path("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz"),
                        help="Path to chain file")
    parser.add_argument("--input", type=Path, help="Path to input file")
    parser.add_argument("--format", choices=["bed", "bam", "vcf", "gff"], default="bed",
                        help="Input file format")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--output-dir", type=Path, default=Path("results"),
                        help="Output directory for results")
    parser.add_argument("--generate-test-data", action="store_true",
                        help="Generate test data before running")
    parser.add_argument("--test-records", type=int, default=10000,
                        help="Number of test records to generate")
    parser.add_argument("--quick-test", action="store_true",
                        help="Run quick test with small dataset")
    parser.add_argument("--tools", nargs="+", 
                        choices=["fastcrossmap", "crossmap", "liftover", "fastremap"],
                        help="Specific tools to benchmark")
    
    args = parser.parse_args()
    
    # Quick test mode
    if args.quick_test:
        args.generate_test_data = True
        args.test_records = 1000
    
    # Generate test data if requested
    if args.generate_test_data:
        test_data_dir = Path("test_data")
        args.input = generate_test_data(test_data_dir, args.test_records)
    
    # Validate inputs
    if not args.input or not args.input.exists():
        print("Error: Input file not specified or does not exist")
        print("Use --generate-test-data to create test data, or specify --input")
        sys.exit(1)
    
    if not args.chain.exists():
        print(f"Error: Chain file not found: {args.chain}")
        sys.exit(1)
    
    # Build FastCrossMap if needed
    if not Path("./target/release/fast-crossmap").exists():
        print("Building FastCrossMap...")
        subprocess.run(["cargo", "build", "--release"], check=True)
    
    # Run benchmarks
    runner = BenchmarkRunner(
        chain_file=args.chain,
        threads=args.threads,
        output_dir=args.output_dir
    )
    
    report = runner.run_all_benchmarks(
        input_file=args.input,
        format=args.format,
        tools=args.tools
    )
    
    # Export results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    json_path = args.output_dir / f"benchmark_{args.format}_{timestamp}.json"
    csv_path = args.output_dir / f"benchmark_{args.format}_{timestamp}.csv"
    
    export_to_json(report, json_path)
    export_to_csv(report, csv_path)
    
    # Also save as latest
    export_to_json(report, args.output_dir / f"benchmark_{args.format}_latest.json")
    export_to_csv(report, args.output_dir / f"benchmark_{args.format}_latest.csv")
    
    # Print summary
    print_summary(report)
    
    print(f"\nResults saved to: {args.output_dir}/")


if __name__ == "__main__":
    main()
