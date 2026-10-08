#!/usr/bin/env python3
"""
02g_benchmark_maf.py - MAF format benchmark (multi-dataset)

Run FastCrossMap vs CrossMap benchmark for MAF format.
(liftOver and FastRemap do not support MAF)

Usage: python paper/02g_benchmark_maf.py
Output: paper/results/benchmark_maf.json
"""

import json
import statistics
import os
import subprocess
import shutil
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Optional

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


CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
REF_GENOME = DATA_DIR / "hg38.fa"
NUM_RUNS = 5

MAF_DATASETS = [
    {"name": "Synthetic xlarge", "file": DATA_DIR / "maf_synthetic_xlarge.maf",
     "source": "Synthetic (1M mutations)"},
    {"name": "Synthetic large", "file": DATA_DIR / "maf_synthetic_large.maf",
     "source": "Synthetic (500K mutations)"},
    {"name": "Synthetic medium", "file": DATA_DIR / "maf_synthetic_medium.maf",
     "source": "Synthetic (200K mutations)"},
    {"name": "Synthetic small", "file": DATA_DIR / "maf_synthetic_small.maf",
     "source": "Synthetic (100K mutations)"},
    {"name": "Synthetic tiny", "file": DATA_DIR / "maf_synthetic_tiny.maf",
     "source": "Synthetic (50K mutations)"},
]


@dataclass
class BenchmarkResult:
    tool: str
    format: str
    dataset_name: str
    input_file: str
    input_mutations: int
    input_size_mb: float
    execution_time_sec: float
    throughput_mutations_per_sec: float
    peak_memory_mb: float
    all_times: list
    success: bool
    error_message: Optional[str] = None
    # `execution_time_sec` is the arithmetic mean. These formats run from under a
    # second to several minutes, and on the short rows a single scheduling hiccup
    # moves the mean by tens of percent (published spread: GFF 6.67-35.41 s on
    # one FastCrossMap row), so the robust statistics are recorded too.
    min_time_sec: Optional[float] = None
    median_time_sec: Optional[float] = None
    min_throughput: Optional[float] = None
    median_throughput: Optional[float] = None
    std_time_sec: Optional[float] = None
    # Per-run `File system inputs` (KB) and major page faults, in the
    # same order as `all_times`. A non-zero entry means that run read
    # its input from the device rather than from page cache.
    all_disk_read_kb: Optional[list] = None
    all_major_faults: Optional[list] = None


def count_maf_mutations(maf_file: Path) -> int:
    count = 0
    with open(maf_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('Hugo_Symbol'):
                continue
            if line.strip():
                count += 1
    return count


def get_file_size_mb(file_path: Path) -> float:
    return file_path.stat().st_size / (1024 * 1024)


def run_with_time(cmd: list[str]) -> tuple[float, float, bool, str]:
    # Wall clock and peak RSS both come from ONE `/usr/bin/time -v` execution.
    # The previous version ran the command a second time whenever the wall-clock
    # string failed to parse (and, in the other shape, always ran it a second
    # time just to read `Maximum resident set size`), so every data point was an
    # average over two executions of the same command.
    try:
        result = subprocess.run(
            ["/usr/bin/time", "-v"] + cmd,
            capture_output=True, text=True, timeout=3600
        )
    except subprocess.TimeoutExpired:
        return 3600, 0, False, "Timeout after 3600 seconds", 0, 0
    except Exception as e:
        return 0, 0, False, str(e), 0, 0

    elapsed = 0.0
    peak_memory_mb = 0.0
    disk_read_kb = 0
    major_faults = 0
    for line in result.stderr.split('\n'):
        if 'File system inputs' in line:
            # Bytes this run actually fetched from the device. This is the only
            # per-run number that separates "the input was in page cache" from
            # "the timed command paid for the disk read": peak RSS cannot show
            # it, and the wall clock only reveals it when the tool is faster
            # than the disk, which is why the GFF run looked like a slow tool.
            try:
                disk_read_kb = int(line.split(':')[1].strip())
            except ValueError:
                pass
        elif 'Major (requiring I/O) page faults' in line:
            try:
                major_faults = int(line.split(':')[1].strip())
            except ValueError:
                pass
        elif 'Maximum resident set size' in line:
            try:
                peak_memory_mb = int(line.split(':')[1].strip()) / 1024
            except ValueError:
                pass
        elif 'Elapsed (wall clock)' in line:
            # The line is "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.23";
            # the timestamp is the final whitespace-separated field. Splitting
            # the whole line on ':' instead lands in the "(h:mm:ss or m:ss)"
            # parenthetical and yields zeros.
            parts = line.split()[-1].split(':')
            try:
                if len(parts) == 3:
                    elapsed = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
                elif len(parts) == 2:
                    elapsed = int(parts[0]) * 60 + float(parts[1])
            except ValueError:
                elapsed = 0.0

    if result.returncode != 0:
        return elapsed, peak_memory_mb, False, result.stderr[:500], disk_read_kb, major_faults
    return elapsed, peak_memory_mb, True, "", disk_read_kb, major_faults


def _output_snapshot(outputs: list) -> dict:
    """Size and mtime of each expected output, taken before the timed runs."""
    snapshot = {}
    for path in outputs:
        try:
            st = path.stat()
            snapshot[str(path)] = (st.st_mtime_ns, st.st_size)
        except OSError:
            snapshot[str(path)] = None
    return snapshot


def _count_lines(path: Path, skip_hash: bool, skip_at: bool = False,
                 skip_fasta: bool = False) -> int:
    """Records in `path`, using the per-format rules in `count_output_records`."""
    count = 0
    fasta = False
    with open(path, 'rb') as f:
        for line in f:
            if fasta:
                if line.startswith(b'>'):
                    fasta = False
                continue
            if skip_fasta and line.startswith(b'>'):
                fasta = True
                continue
            if line in (b'\n', b''):
                continue
            if skip_hash and line.startswith(b'#'):
                continue
            if skip_at and line.startswith(b'@'):
                continue
            count += 1
    return count


def count_output_records(outputs: list) -> tuple[str, bool]:
    """Mapped-record count per output, for the run log, plus "did it map any".

    Counting bytes cannot tell a converted file from a header: the GFF run had a
    dataset whose seqids matched nothing in the chain file, where no tool mapped
    anything, and a 34 KB FastCrossMap output made of comment lines was recorded
    as a 27.51 s success. A count is what separates "this tool lost on this
    dataset" from "this tool was merely slower". What counts as a record differs
    per format:

      * BED: no header, every line is a record.
      * VCF: `#` lines are headers.
      * GFF: `#` headers, plus a `##FASTA` section and the sequence under it.
      * SAM: `@` and `#` lines are headers.
      * MAF: line 1 is `#version`, the rest are aligned blocks.
      * Wiggle: `fixedStep` / `variableStep` / `track` lines are headers.
      * BAM / BigWig: binary, so there is no line structure to count; the size
        is reported instead.
      * bgzip input (`.gz`): a member count would be a *block* count, not a
        record count, so it is labelled rather than counted.
    """
    parts = []
    any_records = False
    for path in outputs:
        name = path.name.lower()
        if not path.exists():
            parts.append(f"{path.name}: missing")
            continue
        if name.endswith('.gz'):
            parts.append(f"{path.name}: bgzf")
            any_records = any_records or path.stat().st_size > 0
            continue
        try:
            if name.endswith(('.vcf', '.gff', '.gff3')):
                n = _count_lines(path, skip_hash=True, skip_fasta=True)
            elif name.endswith('.sam'):
                n = _count_lines(path, skip_hash=True, skip_at=True)
            elif name.endswith('.maf'):
                n = _count_lines(path, skip_hash=False) - 1
            elif name.endswith('.wig'):
                n = _count_lines(path, skip_hash=True) - sum(
                    1 for line in open(path, 'rb')
                    if line.startswith((b'fixedStep', b'variableStep', b'track')))
            elif name.endswith(('.bam', '.bw', '.bigwig')):
                size = path.stat().st_size
                parts.append(f"{path.name}: {size:,} B (binary)")
                any_records = any_records or size > 0
                continue
            else:  # BED
                n = _count_lines(path, skip_hash=False)
            n = max(n, 0)
            parts.append(f"{path.name}: {n:,} records")
            any_records = any_records or n > 0
        except OSError as e:
            parts.append(f"{path.name}: unreadable ({e.strerror})")
    return "; ".join(parts), any_records


def _check_outputs(outputs: list, snapshot: dict) -> tuple[bool, str]:
    """True when the timed command left a non-empty, freshly written output.

    `run_with_time` only looks at the exit status, and every tool here can exit
    0 without producing usable output -- an unrecognised flag, a refused
    overwrite, a truncated stream. Requiring a new mtime as well as a non-zero
    size also rejects a stale file left behind by an earlier run.
    """
    for path in outputs:
        try:
            st = path.stat()
        except OSError:
            continue
        if st.st_size == 0:
            continue
        previous = snapshot.get(str(path))
        if previous is None or st.st_mtime_ns > previous[0]:
            return True, ""
    return False, "no output written: " + ", ".join(str(p) for p in outputs)


def warm_page_cache(path: Optional[Path], chunk_size: int = 8 << 20) -> None:
    """Read `path` once and throw it away, so the timed run starts from RAM.

    The tools for a dataset run in a fixed order, so whenever an input is larger
    than the machine's free RAM the first tool pays for the disk read and the
    ones after it read from cache. That cost is not a property of the tool being
    timed, and it is large: in the 2026-10-07 GFF results FastCrossMap recorded
    80 s on the 1.87 GB GENCODE file and 80 s again on the 1.11 GB basic file --
    the same number for 60% of the input, which is what a disk-bound run looks
    like, not a conversion. Rerunning that file by hand on the same server, with
    the file already in cache, gave 5.20 s against 80 s in the recorded run.

    Warming before *every* timed run rather than once per dataset keeps the
    comparison fair even when the input does not fit in RAM: each tool then
    starts from whatever the cache holds after the same warm-up read, instead of
    from whatever the previous tool happened to leave behind. The read is
    outside the timed region, so it cannot inflate a measurement.
    """
    if path is None or not path.exists():
        return
    with open(path, 'rb') as f:
        while f.read(chunk_size):
            pass


def run_tool_benchmark(tool_name: str, cmd: list[str], num_mutations: int,
                       size_mb: float, dataset_name: str, input_file: str,
                       outputs: list = (),
                       warm_file: Optional[Path] = None,
                       warm_ref: Optional[Path] = None) -> BenchmarkResult:
    times = []
    memories = []
    disk_reads = []
    major_faults_list = []
    success = False
    error_msg = ""
    snapshot = _output_snapshot(outputs)

    for i in range(NUM_RUNS):
        print(f"    Run {i+1}/{NUM_RUNS}...")
        # Bring the input into page cache first, so the wall clock below measures
        # the conversion and not who happened to read the file first.
        warm_page_cache(warm_file)
        # The reference FASTA is read too (scattered reads for the bases
        # under each mapped interval), so it needs the same treatment.
        warm_page_cache(warm_ref)
        elapsed, memory, ok, err, disk_kb, major = run_with_time(cmd)
        if ok and elapsed > 0:
            times.append(elapsed)
            memories.append(memory)
            disk_reads.append(disk_kb)
            major_faults_list.append(major)
            success = True
        elif ok:
            # /usr/bin/time printed no wall-clock value: the run is not a
            # measurement, and averaging it in would divide by zero below.
            error_msg = "no elapsed time recorded"
        else:
            error_msg = err

    counts, mapped_any = count_output_records(outputs)
    print(f"    Output: {counts}")
    if any(k > 1024 for k in disk_reads):
        print(f"    Disk reads (KB): {disk_reads}  <-- not all runs were cached")
    if not mapped_any and not error_msg:
        error_msg = f"no records mapped [{counts}]"

    written, output_err = _check_outputs(outputs, snapshot)
    if not times or not written or not mapped_any:
        print(f"    FAILED: {output_err or error_msg}")
        return BenchmarkResult(
            tool=tool_name, format="MAF", dataset_name=dataset_name,
            input_file=input_file, input_mutations=num_mutations,
            input_size_mb=size_mb, execution_time_sec=0,
            throughput_mutations_per_sec=0, peak_memory_mb=0,
            all_times=[], success=False, error_message=output_err or error_msg
        )

    avg_time = sum(times) / len(times)
    avg_memory = sum(memories) / len(memories)
    min_time = min(times)
    median_time = statistics.median(times)

    return BenchmarkResult(
        tool=tool_name, format="MAF", dataset_name=dataset_name,
        input_file=input_file, input_mutations=num_mutations,
        input_size_mb=round(size_mb, 2),
        execution_time_sec=round(avg_time, 2),
        throughput_mutations_per_sec=round(num_mutations / avg_time, 0) if num_mutations > 0 else 0,
        peak_memory_mb=round(avg_memory, 2),
        all_times=times, success=success,
        all_disk_read_kb=disk_reads, all_major_faults=major_faults_list,
        min_time_sec=round(min_time, 2),
        median_time_sec=round(median_time, 2),
        min_throughput=round(num_mutations / min_time, 0) if num_mutations > 0 else 0,
        median_throughput=round(num_mutations / median_time, 0) if num_mutations > 0 else 0,
        std_time_sec=round(statistics.pstdev(times), 2) if len(times) > 1 else 0.0,
    )


def benchmark_dataset(maf_file: Path, dataset_name: str, output_dir: Path) -> list[BenchmarkResult]:
    num_mutations = count_maf_mutations(maf_file)
    size_mb = get_file_size_mb(maf_file)
    print(f"  Mutations: {num_mutations:,} | Size: {size_mb:.1f} MB")

    ds_dir = output_dir / maf_file.stem
    ds_dir.mkdir(parents=True, exist_ok=True)
    results = []

    # 1. FastCrossMap
    print(f"  [1/2] FastCrossMap on {dataset_name}")
    out = ds_dir / "fastcrossmap_output.maf"
    cmd = ["./target/release/fast-crossmap", "maf", "-b", "GRCh38",
           str(CHAIN_FILE), str(maf_file), str(REF_GENOME), str(out)]
    results.append(run_tool_benchmark(
        "FastCrossMap", cmd, num_mutations, size_mb, dataset_name, str(maf_file), outputs=[out], warm_file=maf_file,
                       warm_ref=REF_GENOME))

    # 2. CrossMap
    print(f"  [2/2] CrossMap on {dataset_name}")
    out = ds_dir / "crossmap_output.maf"
    cmd = ["CrossMap", "maf", str(CHAIN_FILE), str(maf_file),
           str(REF_GENOME), "GRCh38", str(out)]
    results.append(run_tool_benchmark(
        "CrossMap", cmd, num_mutations, size_mb, dataset_name, str(maf_file), outputs=[out], warm_file=maf_file,
                       warm_ref=REF_GENOME))

    cleanup_outputs(ds_dir)
    return results


def main():
    print("=" * 60)
    print("MAF Format Benchmark (Multi-Dataset)")
    print("=" * 60)
    print("Note: liftOver and FastRemap do not support MAF\n")

    if not CHAIN_FILE.exists():
        print(f"Error: Chain file not found: {CHAIN_FILE}")
        print("Please run first: bash paper/01_download_data.sh")
        return

    if not REF_GENOME.exists():
        print(f"Error: Reference genome not found: {REF_GENOME}")
        print("Please run: bash paper/01_download_data.sh --ref-only")
        return

    datasets = []
    for ds in MAF_DATASETS:
        f = ds["file"]
        if f.exists():
            datasets.append((f, ds["name"], ds["source"]))
        else:
            print(f"Warning: {ds['name']} not found at {f}, skipping")

    if not datasets:
        print("Error: No MAF datasets found. Run bash paper/01_download_data.sh")
        return

    output_dir = RESULTS_DIR / "maf_benchmark"
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = []
    for maf_file, name, source in datasets:
        print(f"\n{'─'*60}")
        print(f"Dataset: {name} ({source})")
        print(f"File: {maf_file}")
        print(f"{'─'*60}")
        ds_results = benchmark_dataset(maf_file, name, output_dir)
        all_results.extend(ds_results)

    output_json = RESULTS_DIR / "benchmark_maf.json"
    with open(output_json, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "num_runs": NUM_RUNS,
            "num_datasets": len(datasets),
            "datasets": [
                {"name": n, "file": str(fp), "source": s}
                for fp, n, s in datasets
            ],
            "results": [asdict(r) for r in all_results]
        }, f, indent=2)

    print(f"\nResults saved to: {output_json}")

    # Summary
    print(f"\n{'='*70}")
    print("Benchmark Results Summary")
    print(f"{'='*70}")
    print(f"{'Dataset':<22} {'Tool':<15} {'mean(s)':<10} {'min(s)':<9} {'med(s)':<9} {'Mut/s':<14} {'Mem(MB)':<10} {'OK'}")
    print("-" * 93)
    for r in all_results:
        ok = "✓" if r.success else "✗"
        print(f"{r.dataset_name:<22} {r.tool:<15} {r.execution_time_sec:<10.2f} {(r.min_time_sec or 0):<9.2f} {(r.median_time_sec or 0):<9.2f} "
              f"{r.throughput_mutations_per_sec:<14,.0f} {r.peak_memory_mb:<10.1f} {ok}")

    # Speedup
    print(f"\n{'='*60}")
    print("Speedup Summary (CrossMap / FastCrossMap)")
    print(f"{'='*60}")
    print(f"{'Dataset':<22} {'mean':>8} {'min':>8} {'median':>8}")
    for _, name, _ in datasets:
        fc = next((r for r in all_results if r.dataset_name == name
                   and r.tool == "FastCrossMap" and r.success), None)
        cm = next((r for r in all_results if r.dataset_name == name
                   and r.tool == "CrossMap" and r.success), None)
        if not (fc and cm):
            continue
        stats = (fc.execution_time_sec, fc.min_time_sec or 0, fc.median_time_sec or 0,
                 cm.execution_time_sec, cm.min_time_sec or 0, cm.median_time_sec or 0)
        if min(stats) <= 0:
            continue
        by_mean = cm.execution_time_sec / fc.execution_time_sec
        by_min = cm.min_time_sec / fc.min_time_sec
        by_med = cm.median_time_sec / fc.median_time_sec
        print(f"  {name:<20} {by_mean:7.1f}x {by_min:7.1f}x {by_med:7.1f}x")

    print("\nNext step: python paper/02h_benchmark_sam.py")


if __name__ == "__main__":
    main()
