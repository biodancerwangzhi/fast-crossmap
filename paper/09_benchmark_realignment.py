#!/usr/bin/env python3
"""
09_benchmark_realignment.py - BAM liftover vs re-alignment benchmark

Compares FastCrossMap BAM liftover against a full BWA-MEM re-alignment
pipeline at multiple data scales and thread counts. Uses wgsim to generate
simulated reads from hg19, then benchmarks both approaches for converting
coordinates to hg38.

Usage: python3 paper/09_benchmark_realignment.py [--scales 1M,5M] [--threads 1,4,8]
Output: paper/results/benchmark_realignment.json
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, asdict, field
from datetime import datetime
from pathlib import Path
from typing import Optional

# =============================================================================
# Configuration
# =============================================================================
DATA_DIR = Path("paper/data")
SIM_DIR = DATA_DIR / "simulated"
RESULTS_DIR = Path("paper/results")

CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
HG19_FA = DATA_DIR / "hg19.fa"
HG38_FA = DATA_DIR / "hg38.fa"
FCM_BIN = Path("./target/release/fast-crossmap")

NUM_RUNS = 3

SCALES = {
    "1M":  1_000_000,
    "5M":  5_000_000,
    "10M": 10_000_000,
}

DEFAULT_THREADS = [1, 4, 8]


# =============================================================================
# Data classes
# =============================================================================
@dataclass
class RunResult:
    wall_time_sec: float
    peak_memory_mb: float
    success: bool
    error: str = ""


@dataclass
class BenchmarkEntry:
    method: str
    scale: str
    num_reads: int
    threads: int
    mean_time_sec: float
    std_time_sec: float
    min_time_sec: float
    max_time_sec: float
    peak_memory_mb: float
    all_times: list = field(default_factory=list)
    success: bool = True
    error: str = ""


# =============================================================================
# Utility functions
# =============================================================================
def check_tool(name: str, cmd: list[str]) -> bool:
    try:
        subprocess.run(cmd, capture_output=True, timeout=10)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def check_prerequisites() -> bool:
    ok = True
    tools = [
        ("bwa", ["bwa"]),
        ("samtools", ["samtools", "version"]),
        ("wgsim",    ["wgsim"]),  # wgsim prints usage to stderr on no args
        ("fast-crossmap", [str(FCM_BIN), "--version"]),
    ]
    for name, cmd in tools:
        if name in ("wgsim", "bwa"):
            # these return non-zero with no args but that means they exist
            found = shutil.which(name) is not None
        else:
            found = check_tool(name, cmd)
        status = "OK" if found else "MISSING"
        print(f"  {name}: {status}")
        if not found:
            ok = False

    for ref, label in [(HG19_FA, "hg19.fa"), (HG38_FA, "hg38.fa")]:
        exists = ref.exists()
        print(f"  {label}: {'OK' if exists else 'MISSING'} ({ref})")
        if not exists:
            ok = False

    if not CHAIN_FILE.exists():
        print(f"  chain file: MISSING ({CHAIN_FILE})")
        ok = False
    else:
        print(f"  chain file: OK")

    return ok


# Per-run cap for one re-alignment run (FASTQ conversion + BWA-MEM + sort +
# markdup over the whole file). It scales with the read count and the machine
# is often busy, so keep this generous -- a 10M-read run can take hours, and
# when the cap fires the run is scored as a failure rather than retried.
DEFAULT_TIMEOUT = 43200  # 12 h


def run_timed(cmd: list[str], timeout: int = DEFAULT_TIMEOUT) -> RunResult:
    try:
        result = subprocess.run(
            ["/usr/bin/time", "-v"] + cmd,
            capture_output=True, text=True, timeout=timeout
        )
        wall_time = None
        peak_mem = 0
        for line in result.stderr.split('\n'):
            if 'Elapsed (wall clock)' in line:
                # Format: h:mm:ss or m:ss or m:ss.ss
                match = re.search(r'(\d+):(\d+)[.:](\d+)', line.split(': ', 1)[1].strip())
                if match:
                    parts = line.split(': ', 1)[1].strip().split(':')
                    if len(parts) == 3:
                        wall_time = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
                    elif len(parts) == 2:
                        wall_time = int(parts[0]) * 60 + float(parts[1])
            if 'Maximum resident set size' in line:
                peak_mem = int(line.split(':')[1].strip()) / 1024  # KB -> MB

        if wall_time is None:
            wall_time = 0
            # Fallback: if /usr/bin/time parsing failed, not reliable
            return RunResult(0, peak_mem, False, "Could not parse wall time")

        if result.returncode != 0:
            err = result.stderr[:500] if result.stderr else "non-zero exit"
            return RunResult(wall_time, peak_mem, False, err)

        return RunResult(wall_time, peak_mem, True)

    except subprocess.TimeoutExpired:
        return RunResult(timeout, 0, False, f"Timeout after {timeout}s")
    except Exception as e:
        return RunResult(0, 0, False, str(e))


def run_pipeline_timed(cmds_str: str, timeout: int = DEFAULT_TIMEOUT) -> RunResult:
    """Time a shell pipeline string (e.g. 'cmd1 | cmd2 | cmd3').

    Run with `pipefail`: without it, `samtools fastq ... | bwa mem ... | samtools
    sort ...` reports the exit status of `sort` alone, so a missing `bwa` (or any
    upstream crash) exits 0 and leaves behind a valid, empty BAM that scores as a
    successful run.
    """
    try:
        result = subprocess.run(
            ["/usr/bin/time", "-v", "bash", "-o", "pipefail", "-c", cmds_str],
            capture_output=True, text=True, timeout=timeout
        )
        wall_time = None
        peak_mem = 0
        for line in result.stderr.split('\n'):
            if 'Elapsed (wall clock)' in line:
                parts = line.split(': ', 1)[1].strip().split(':')
                if len(parts) == 3:
                    wall_time = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
                elif len(parts) == 2:
                    wall_time = int(parts[0]) * 60 + float(parts[1])
            if 'Maximum resident set size' in line:
                peak_mem = int(line.split(':')[1].strip()) / 1024

        if wall_time is None:
            return RunResult(0, peak_mem, False, "Could not parse wall time")
        if result.returncode != 0:
            return RunResult(wall_time, peak_mem, False, result.stderr[:500])
        return RunResult(wall_time, peak_mem, True)

    except subprocess.TimeoutExpired:
        return RunResult(timeout, 0, False, f"Timeout after {timeout}s")
    except Exception as e:
        return RunResult(0, 0, False, str(e))


# =============================================================================
# BWA index
# =============================================================================
def count_bam_reads(bam: Path) -> int:
    """Number of alignment records, or -1 if samtools is unavailable."""
    if shutil.which("samtools") is None:
        return -1
    try:
        r = subprocess.run(["samtools", "view", "-c", str(bam)],
                           capture_output=True, text=True, timeout=1800)
        return int(r.stdout.strip()) if r.returncode == 0 else -1
    except (ValueError, subprocess.SubprocessError):
        return -1


def ensure_bwa_index(ref_fa: Path):
    idx = Path(str(ref_fa) + ".bwt")
    if idx.exists():
        print(f"  BWA index exists for {ref_fa.name}")
        return
    print(f"  Building BWA index for {ref_fa.name} (this takes ~10 min)...")
    result = subprocess.run(
        ["bwa", "index", str(ref_fa)],
        capture_output=True, text=True, timeout=7200
    )
    if result.returncode != 0:
        print(f"  ERROR: bwa index failed: {result.stderr[:300]}")
        sys.exit(1)
    print(f"  Index built.")


# =============================================================================
# Simulate reads with wgsim
# =============================================================================
def simulate_reads(num_reads: int, scale_name: str) -> tuple[Path, Path]:
    SIM_DIR.mkdir(parents=True, exist_ok=True)
    fq1 = SIM_DIR / f"sim_{scale_name}_R1.fq"
    fq2 = SIM_DIR / f"sim_{scale_name}_R2.fq"

    if fq1.exists() and fq2.exists():
        print(f"  Simulated reads already exist for {scale_name}")
        return fq1, fq2

    print(f"  Generating {num_reads:,} paired-end reads with wgsim...")
    # Extract chr1 from hg19 for simulation
    chr1_fa = SIM_DIR / "hg19_chr1.fa"
    if not chr1_fa.exists():
        print(f"  Extracting chr1 from hg19...")
        subprocess.run(
            ["samtools", "faidx", str(HG19_FA), "chr1"],
            capture_output=True, text=True, check=True, timeout=300
        )
        result = subprocess.run(
            ["samtools", "faidx", str(HG19_FA), "chr1"],
            capture_output=True, text=True, check=True, timeout=300
        )
        chr1_fa.write_text(result.stdout)

    cmd = [
        "wgsim",
        "-N", str(num_reads),
        "-1", "150", "-2", "150",
        "-r", "0", "-R", "0", "-X", "0",
        str(chr1_fa),
        str(fq1), str(fq2),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    if result.returncode != 0:
        print(f"  ERROR: wgsim failed: {result.stderr[:300]}")
        sys.exit(1)

    print(f"  Generated: {fq1.name} ({fq1.stat().st_size / 1e6:.0f} MB), "
          f"{fq2.name} ({fq2.stat().st_size / 1e6:.0f} MB)")
    return fq1, fq2


# =============================================================================
# Generate hg19 BAM (align simulated reads to hg19)
# =============================================================================
def generate_hg19_bam(fq1: Path, fq2: Path, scale_name: str, threads: int = 4) -> Path:
    bam_out = SIM_DIR / f"sim_{scale_name}_hg19.bam"

    if bam_out.exists() and count_bam_reads(bam_out) != 0:
        print(f"  hg19 BAM already exists for {scale_name}")
        return bam_out
    if bam_out.exists():
        print(f"  {bam_out.name} exists but holds 0 reads -- regenerating")

    print(f"  Aligning to hg19 with BWA-MEM ({threads} threads)...")
    pipeline = (
        f"bwa mem -t {threads} {HG19_FA} {fq1} {fq2} "
        f"| samtools sort -@ {threads} -o {bam_out}"
    )
    # pipefail so a failed bwa is not masked by a successful sort
    result = subprocess.run(
        ["bash", "-o", "pipefail", "-c", pipeline],
        capture_output=True, text=True, timeout=DEFAULT_TIMEOUT
    )
    if result.returncode != 0:
        print(f"  ERROR: alignment failed: {result.stderr[:500]}")
        sys.exit(1)

    # An empty BAM is not a valid result. Without this check the whole benchmark
    # below runs happily on zero reads and reports meaningless numbers.
    n_reads = count_bam_reads(bam_out)
    if n_reads == 0:
        print(f"  ERROR: {bam_out.name} contains 0 reads -- alignment produced nothing.")
        print(f"         stderr: {result.stderr[:300]}")
        sys.exit(1)
    print(f"  {n_reads:,} reads aligned")

    subprocess.run(["samtools", "index", str(bam_out)], check=True, timeout=300)
    size_mb = bam_out.stat().st_size / 1e6
    print(f"  Generated: {bam_out.name} ({size_mb:.0f} MB)")
    return bam_out


# =============================================================================
# Benchmark: Re-alignment pipeline
# =============================================================================
def benchmark_realignment(bam_in: Path, scale_name: str, threads: int) -> list[RunResult]:
    """Full re-alignment: BAM→FASTQ→BWA-MEM(hg38)→sort→markdup."""
    results = []
    out_prefix = SIM_DIR / f"realign_{scale_name}_t{threads}"
    out_bam = Path(f"{out_prefix}_hg38.bam")

    pipeline = (
        f"samtools bam2fq {bam_in} "
        f"| bwa mem -t {threads} {HG38_FA} -p - "
        f"| samtools sort -@ {threads} -o {out_prefix}_sorted.bam && "
        f"samtools markdup -@ {threads} {out_prefix}_sorted.bam {out_bam} && "
        f"rm -f {out_prefix}_sorted.bam"
    )

    for i in range(NUM_RUNS):
        print(f"      Run {i+1}/{NUM_RUNS}...")
        # Clean up previous output. This also drops the partial file left behind
        # by a run killed at the timeout -- otherwise the next run's sort step
        # would read it back instead of the current run's alignments.
        for f in [out_bam, Path(f"{out_prefix}_sorted.bam")]:
            f.unlink(missing_ok=True)
        r = run_pipeline_timed(pipeline)
        results.append(r)
        if not r.success:
            print(f"      FAILED: {r.error[:100]}")

    # Clean up
    out_bam.unlink(missing_ok=True)
    return results


# =============================================================================
# Benchmark: FastCrossMap liftover
# =============================================================================
def benchmark_liftover(bam_in: Path, scale_name: str, threads: int) -> list[RunResult]:
    results = []
    out_bam = SIM_DIR / f"liftover_{scale_name}_t{threads}_hg38.bam"

    cmd = [
        str(FCM_BIN), "bam",
        "-t", str(threads),
        str(CHAIN_FILE), str(bam_in), str(out_bam),
    ]

    for i in range(NUM_RUNS):
        print(f"      Run {i+1}/{NUM_RUNS}...")
        out_bam.unlink(missing_ok=True)
        Path(str(out_bam) + ".unmap").unlink(missing_ok=True)
        r = run_timed(cmd)
        results.append(r)
        if not r.success:
            print(f"      FAILED: {r.error[:100]}")

    # Clean up
    out_bam.unlink(missing_ok=True)
    Path(str(out_bam) + ".unmap").unlink(missing_ok=True)
    return results


# =============================================================================
# Aggregate results
# =============================================================================
def aggregate(method: str, scale: str, num_reads: int, threads: int,
              runs: list[RunResult]) -> BenchmarkEntry:
    ok_times = [r.wall_time_sec for r in runs if r.success]
    ok_mems = [r.peak_memory_mb for r in runs if r.success]

    if not ok_times:
        return BenchmarkEntry(
            method=method, scale=scale, num_reads=num_reads, threads=threads,
            mean_time_sec=0, std_time_sec=0, min_time_sec=0, max_time_sec=0,
            peak_memory_mb=0, all_times=[], success=False,
            error=runs[0].error if runs else "no runs"
        )

    import statistics
    mean_t = statistics.mean(ok_times)
    std_t = statistics.stdev(ok_times) if len(ok_times) > 1 else 0
    peak_mem = max(ok_mems) if ok_mems else 0

    return BenchmarkEntry(
        method=method, scale=scale, num_reads=num_reads, threads=threads,
        mean_time_sec=round(mean_t, 2),
        std_time_sec=round(std_t, 2),
        min_time_sec=round(min(ok_times), 2),
        max_time_sec=round(max(ok_times), 2),
        peak_memory_mb=round(peak_mem, 1),
        all_times=[round(t, 2) for t in ok_times],
        success=True,
    )


# =============================================================================
# Main
# =============================================================================
def main():
    parser = argparse.ArgumentParser(description="BAM liftover vs re-alignment benchmark")
    parser.add_argument("--scales", type=str, default=None,
                        help="Comma-separated scales to test (e.g. 1M,5M,10M)")
    parser.add_argument("--threads", type=str, default=None,
                        help="Comma-separated thread counts (e.g. 1,4,8)")
    parser.add_argument("--prep-threads", type=int, default=4,
                        help="Threads for data preparation (default: 4)")
    args = parser.parse_args()

    scales = list(SCALES.keys())
    if args.scales:
        scales = [s.strip() for s in args.scales.split(",")]
        for s in scales:
            if s not in SCALES:
                print(f"ERROR: unknown scale '{s}'. Available: {list(SCALES.keys())}")
                sys.exit(1)

    thread_counts = DEFAULT_THREADS
    if args.threads:
        thread_counts = [int(t) for t in args.threads.split(",")]

    print("=" * 60)
    print(" BAM Liftover vs Re-alignment Benchmark")
    print("=" * 60)
    print(f"  Scales: {scales}")
    print(f"  Threads: {thread_counts}")
    print(f"  Runs per config: {NUM_RUNS}")
    print()

    # --- Prerequisites ---
    print("[1/5] Checking prerequisites...")
    if not check_prerequisites():
        print("\nERROR: Missing prerequisites. Install them and try again.")
        sys.exit(1)
    print()

    # --- BWA indexes ---
    print("[2/5] Ensuring BWA indexes...")
    ensure_bwa_index(HG19_FA)
    ensure_bwa_index(HG38_FA)
    print()

    all_entries = []

    for scale_name in scales:
        num_reads = SCALES[scale_name]
        print("=" * 60)
        print(f" Scale: {scale_name} ({num_reads:,} reads)")
        print("=" * 60)

        # --- Simulate reads ---
        print(f"[3/5] Simulating reads...")
        fq1, fq2 = simulate_reads(num_reads, scale_name)

        # --- Generate hg19 BAM ---
        print(f"[4/5] Generating hg19 BAM...")
        bam_in = generate_hg19_bam(fq1, fq2, scale_name, threads=args.prep_threads)
        print()

        # --- Benchmark ---
        print(f"[5/5] Running benchmarks...")
        for t in thread_counts:
            print(f"\n  --- {scale_name}, {t} thread(s) ---")

            print(f"    [A] Re-alignment pipeline (BWA-MEM)...")
            realign_runs = benchmark_realignment(bam_in, scale_name, t)
            entry_r = aggregate("realignment", scale_name, num_reads, t, realign_runs)
            all_entries.append(entry_r)

            # Do not burn hours on a doomed config: if every re-alignment run
            # failed, the cause is the environment (missing bwa, timed-out run,
            # no disk), not the thread count -- the liftover side is skipped and
            # the failure is reported by the summary at the end.
            if not entry_r.success:
                print(f"    [B] FastCrossMap liftover... SKIPPED "
                      f"(re-alignment failed: {entry_r.error[:60]})")
                liftover_runs = [
                    RunResult(0, 0, False, "skipped: re-alignment failed")
                    for _ in range(NUM_RUNS)
                ]
                entry_l = aggregate("liftover", scale_name, num_reads, t, liftover_runs)
                all_entries.append(entry_l)
                print()
                continue

            print(f"    [B] FastCrossMap liftover...")
            liftover_runs = benchmark_liftover(bam_in, scale_name, t)
            entry_l = aggregate("liftover", scale_name, num_reads, t, liftover_runs)
            all_entries.append(entry_l)

            # Summary line
            if entry_r.success and entry_l.success and entry_l.mean_time_sec > 0:
                speedup = entry_r.mean_time_sec / entry_l.mean_time_sec
                print(f"    => Realign: {entry_r.mean_time_sec:.1f}s, "
                      f"Liftover: {entry_l.mean_time_sec:.1f}s, "
                      f"Speedup: {speedup:.1f}x")
            print()

    # --- Save results ---
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output_file = RESULTS_DIR / "benchmark_realignment.json"
    output = {
        "benchmark": "realignment_vs_liftover",
        "date": datetime.now().isoformat(),
        "num_runs": NUM_RUNS,
        "scales": scales,
        "thread_counts": thread_counts,
        "results": [asdict(e) for e in all_entries],
    }
    output_file.write_text(json.dumps(output, indent=2))
    print(f"\nResults saved to {output_file}")

    # --- Print summary table ---
    print("\n" + "=" * 80)
    print(" SUMMARY")
    print("=" * 80)
    print(f"{'Scale':<8} {'Threads':<8} {'Realign (s)':<14} {'Liftover (s)':<14} {'Speedup':<10} {'Mem R (MB)':<12} {'Mem L (MB)':<12}")
    print("-" * 80)
    for i in range(0, len(all_entries), 2):
        r = all_entries[i]      # realignment
        l = all_entries[i + 1]  # liftover
        if r.success and l.success and l.mean_time_sec > 0:
            speedup = f"{r.mean_time_sec / l.mean_time_sec:.1f}x"
        else:
            speedup = "N/A"
        print(f"{r.scale:<8} {r.threads:<8} {r.mean_time_sec:<14.1f} {l.mean_time_sec:<14.1f} {speedup:<10} {r.peak_memory_mb:<12.0f} {l.peak_memory_mb:<12.0f}")
    print("=" * 80)

    # A run that failed everywhere wrote a JSON full of zeros. Say so loudly:
    # those entries must not end up quoted as measurements.
    if any(not e.success for e in all_entries):
        print()
        print("!! Some configurations FAILED and were recorded as zeros.")
        print("!! Check the errors above before using this JSON.")
        failed = sorted({(e.method, e.scale, e.threads, e.error)
                         for e in all_entries if not e.success})
        for method, scale, threads, err in failed[:10]:
            print(f"   {scale} t={threads} {method}: {err[:90]}")


if __name__ == "__main__":
    main()
