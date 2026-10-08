#!/usr/bin/env python3
"""
02_benchmark_bed.py - BED format benchmark (multi-dataset)

Run 4-way benchmark: FastCrossMap, CrossMap, liftOver, FastRemap
across multiple BED datasets for reproducibility validation.

Usage: python paper/02_benchmark_bed.py
Output: paper/results/benchmark_bed.json
"""

import gzip
import json
import os
import shutil
import statistics
import subprocess
from dataclasses import dataclass, asdict, field
from datetime import datetime
from pathlib import Path
from typing import Optional

# =============================================================================
# Configuration
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


CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
CHAIN_FILE_UNZIPPED = DATA_DIR / "hg19ToHg38.over.chain"

NUM_RUNS = 5

BED_DATASETS = [
    {
        "name": "K562 DNase-seq",
        "file": DATA_DIR / "bed_k562_dnase.bed.gz",
        "legacy_file": DATA_DIR / "encode_dnase_peaks.bed.gz",
        "source": "ENCODE ENCFF001WBV",
    },
    {
        "name": "GM12878 DNase-seq",
        "file": DATA_DIR / "bed_gm12878_dnase.bed.gz",
        "source": "ENCODE ENCFF001WFH",
    },
    {
        "name": "H1-hESC DNase-seq",
        "file": DATA_DIR / "bed_h1hesc_dnase.bed.gz",
        "source": "ENCODE ENCFF001WDU",
    },
    {
        "name": "HepG2 DNase-seq",
        "file": DATA_DIR / "bed_hepg2_dnase.narrowPeak.gz",
        "source": "UCSC ENCODE OpenChrom HepG2",
    },
    {
        "name": "A549 DNase-seq",
        "file": DATA_DIR / "bed_a549_dnase.narrowPeak.gz",
        "source": "UCSC ENCODE OpenChrom A549",
    },
]


@dataclass
class BenchmarkResult:
    tool: str
    format: str
    dataset_name: str
    input_file: str
    input_records: int
    execution_time_sec: float
    throughput_rec_per_sec: float
    peak_memory_mb: float
    mapped_records: int
    unmapped_records: int
    all_times: list
    success: bool
    error_message: Optional[str] = None
    # `execution_time_sec` is the arithmetic mean. FastCrossMap finishes these
    # datasets in well under a second, where one scheduling hiccup moves the mean
    # by tens of percent, so the robust statistics are recorded too.
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


def count_bed_records(bed_file: Path) -> int:
    count = 0
    if str(bed_file).endswith('.gz'):
        opener, mode = gzip.open, 'rt'
    else:
        opener, mode = open, 'r'
    with opener(bed_file, mode) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                count += 1
    return count


def count_mapped_unmapped(output_file: Path,
                          unmap_file: Optional[Path] = None) -> tuple[int, int]:
    """Record counts for `output_file` and its companion `.unmap`.

    These are the counts written into the JSON, so they must stay independent of
    `count_output_records`, which reports in a single string for the run log.
    `unmap_file` defaults to `output_file + ".unmap"`, which is where
    FastCrossMap, CrossMap and liftOver put it. FastRemap appends its own `-f`
    suffix and its unmap path cannot be derived, so its caller passes the real
    file.
    """
    if unmap_file is None:
        unmap_file = Path(str(output_file) + ".unmap")
    mapped = 0
    unmapped = 0
    if output_file.exists():
        with open(output_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    mapped += 1
    if unmap_file.exists():
        with open(unmap_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    unmapped += 1
    return mapped, unmapped


def run_with_time(cmd: list[str]) -> tuple[float, float, bool, str]:
    # Wall clock and peak RSS both come from ONE `/usr/bin/time -v` execution.
    # The previous version ran the command a second time whenever the wall-clock
    # string failed to parse (and, in the other shape, always ran it a second
    # time just to read `Maximum resident set size`), so every data point was an
    # average over two executions of the same command.
    try:
        result = subprocess.run(
            ["/usr/bin/time", "-v"] + cmd,
            capture_output=True, text=True, timeout=600
        )
    except subprocess.TimeoutExpired:
        return 600, 0, False, "Timeout after 600 seconds", 0, 0
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


def ensure_decompressed(bed_file: Path) -> Path:
    """Decompress `bed_file` and return a path with a `.bed` extension.

    FastCrossMap, CrossMap and liftOver dispatch on the contents, but FastRemap
    picks its parser from the *input file name*: any other extension aborts with
    `File Extension: <ext> not supported or unknown` on stdout and exit code 1,
    before doing any work. Two of the BED datasets ship as `*.narrowPeak.gz`, so
    their decompressed names are rewritten -- a symlink, leaving the original
    file untouched.
    """
    if not str(bed_file).endswith('.gz'):
        return bed_file
    decompressed = Path(str(bed_file)[:-3])
    if not decompressed.exists():
        print(f"    Decompressing {bed_file.name} ...")
        with gzip.open(bed_file, 'rb') as f_in, open(decompressed, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    if decompressed.suffix != '.bed':
        as_bed = decompressed.with_suffix('.bed')
        if not as_bed.exists():
            as_bed.symlink_to(decompressed.name)
        return as_bed
    return decompressed


def ensure_chain_unzipped() -> Path:
    if not CHAIN_FILE_UNZIPPED.exists():
        print("    Decompressing chain file for FastRemap...")
        subprocess.run(["gunzip", "-k", str(CHAIN_FILE)], check=True)
    return CHAIN_FILE_UNZIPPED


def prepare_bed6_for_liftover(bed_file: Path, output_dir: Path) -> Path:
    bed6_file = output_dir / "input_bed6.bed"
    with open(bed_file, 'r') as fin, open(bed6_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 6:
                try:
                    score = int(float(fields[4]))
                except Exception:
                    score = 0
                strand = fields[5] if len(fields) > 5 else '.'
                fout.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{fields[3]}\t{score}\t{strand}\n")
            elif len(fields) >= 3:
                fout.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t.\t0\t.\n")
    return bed6_file


# PLACEHOLDER_BENCHMARK_FUNCTIONS


def warm_page_cache(path: Optional[Path], chunk_size: int = 8 << 20) -> None:
    """Read `path` once and throw it away, so the timed run starts from RAM.

    The four tools for a dataset run in a fixed order, and only the record count
    after the timed loop reads the input -- so whichever tool goes first pays for
    the disk read and the ones after it read from cache. The BED inputs are small
    enough that this is a sub-second effect, but a sub-second effect is the whole
    measurement here: a first run of 0.81 s against 0.20 s for the next four
    (GM12878, 2026-10-07) is a 4x distortion of run 1, and run 1 is what a
    reader sees if they stop at the first line of `all_times`.

    Warming before *every* timed run rather than once per dataset keeps the
    comparison fair, and the read sits outside the timed region, so it cannot
    inflate a measurement.
    """
    if path is None or not path.exists():
        return
    with open(path, 'rb') as f:
        while f.read(chunk_size):
            pass


def run_tool_benchmark(tool_name: str, cmd: list[str], bed_file: Path,
                       dataset_name: str, output_file: Path,
                       unmap_file: Optional[Path] = None,
                       warm_file: Optional[Path] = None) -> BenchmarkResult:
    times = []
    memories = []
    disk_reads = []
    major_faults_list = []
    success = False
    error_msg = ""

    # Validate the artefact the tool was actually asked to write. FastRemap
    # appends the `-f` suffix itself, so its real output is not the path given
    # to `-o`; the caller passes the real one here.
    snapshot = _output_snapshot([output_file])

    for i in range(NUM_RUNS):
        print(f"    Run {i+1}/{NUM_RUNS}...")
        # Bring the input into page cache first, so the wall clock below
        # measures the conversion and not who happened to read the file first.
        warm_page_cache(warm_file)
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

    input_records = count_bed_records(bed_file)

    counts, mapped_any = count_output_records([output_file])
    print(f"    Output: {counts}")
    if any(k > 1024 for k in disk_reads):
        print(f"    Disk reads (KB): {disk_reads}  <-- not all runs were cached")
    if not mapped_any and not error_msg:
        error_msg = f"no records mapped [{counts}]"

    written, output_err = _check_outputs([output_file], snapshot)
    if not times or not written or not mapped_any:
        print(f"    FAILED: {output_err or error_msg}")
        return BenchmarkResult(
            tool=tool_name, format="BED", dataset_name=dataset_name,
            input_file=str(bed_file), input_records=input_records,
            execution_time_sec=0, throughput_rec_per_sec=0, peak_memory_mb=0,
            mapped_records=0, unmapped_records=0, all_times=[],
            success=False, error_message=output_err or error_msg
        )

    avg_time = sum(times) / len(times)
    avg_memory = sum(memories) / len(memories)
    min_time = min(times)
    median_time = statistics.median(times)
    mapped, unmapped = count_mapped_unmapped(output_file, unmap_file)

    return BenchmarkResult(
        tool=tool_name, format="BED", dataset_name=dataset_name,
        input_file=str(bed_file), input_records=input_records,
        execution_time_sec=round(avg_time, 3),
        throughput_rec_per_sec=round(input_records / avg_time, 0),
        peak_memory_mb=round(avg_memory, 2),
        mapped_records=mapped, unmapped_records=unmapped,
        all_times=times, success=success,
        all_disk_read_kb=disk_reads, all_major_faults=major_faults_list,
        min_time_sec=round(min_time, 3),
        median_time_sec=round(median_time, 3),
        min_throughput=round(input_records / min_time, 0),
        median_throughput=round(input_records / median_time, 0),
        std_time_sec=round(statistics.pstdev(times), 3) if len(times) > 1 else 0.0,
    )


def benchmark_dataset(bed_gz: Path, dataset_name: str, output_dir: Path) -> list[BenchmarkResult]:
    bed_file = ensure_decompressed(bed_gz)
    ds_dir = output_dir / bed_gz.stem.replace('.bed', '')
    ds_dir.mkdir(parents=True, exist_ok=True)

    results = []

    # 1. FastCrossMap
    print(f"  [1/4] FastCrossMap on {dataset_name}")
    out = ds_dir / "fastcrossmap_output.bed"
    cmd = ["./target/release/fast-crossmap", "bed",
           str(CHAIN_FILE), str(bed_file), str(out)]
    results.append(run_tool_benchmark("FastCrossMap", cmd, bed_gz, dataset_name, out,
                                      warm_file=bed_file))

    # 2. CrossMap
    print(f"  [2/4] CrossMap on {dataset_name}")
    out = ds_dir / "crossmap_output.bed"
    cmd = ["CrossMap", "bed", str(CHAIN_FILE), str(bed_file), str(out)]
    results.append(run_tool_benchmark("CrossMap", cmd, bed_gz, dataset_name, out,
                                      warm_file=bed_file))

    # 3. liftOver
    print(f"  [3/4] liftOver on {dataset_name}")
    out = ds_dir / "liftover_output.bed"
    unmap = ds_dir / "liftover_output.bed.unmap"
    bed6 = prepare_bed6_for_liftover(bed_file, ds_dir)
    cmd = ["liftOver", str(bed6), str(CHAIN_FILE), str(out), str(unmap)]
    results.append(run_tool_benchmark("liftOver", cmd, bed_gz, dataset_name, out,
                                      warm_file=bed_file))

    # 4. FastRemap
    print(f"  [4/4] FastRemap on {dataset_name}")
    out = ds_dir / "fastremap_output.bed"
    unmap = ds_dir / "fastremap_output.unmap"
    chain_plain = ensure_chain_unzipped()
    # FastRemap appends the `-f` suffix to `-o` itself, so the prefix is passed
    # without one -- otherwise it writes `fastremap_output.bed.bed`.
    cmd = ["FastRemap", "-f", "bed", "-c", str(chain_plain),
           "-i", str(bed_file), "-o", str(ds_dir / "fastremap_output"),
           "-u", str(unmap)]
    results.append(run_tool_benchmark("FastRemap", cmd, bed_gz, dataset_name,
                                      out, unmap, warm_file=bed_file))

    cleanup_outputs(ds_dir)
    return results


def main():
    print("=" * 60)
    print("BED Format Benchmark (Multi-Dataset)")
    print("=" * 60)

    if not CHAIN_FILE.exists():
        print(f"Error: Chain file not found: {CHAIN_FILE}")
        print("Please run first: bash paper/01_download_data.sh")
        return

    # Resolve datasets: support legacy file names
    datasets = []
    for ds in BED_DATASETS:
        f = ds["file"]
        if not f.exists() and "legacy_file" in ds:
            f = ds["legacy_file"]
        if f.exists():
            datasets.append((f, ds["name"], ds["source"]))
        else:
            print(f"Warning: {ds['name']} not found at {ds['file']}, skipping")

    if not datasets:
        print("Error: No BED datasets found. Run bash paper/01_download_data.sh")
        return

    output_dir = RESULTS_DIR / "bed_benchmark"
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = []
    for bed_file, name, source in datasets:
        records = count_bed_records(bed_file)
        print(f"\n{'─'*60}")
        print(f"Dataset: {name} ({source})")
        print(f"File: {bed_file} | Records: {records:,}")
        print(f"{'─'*60}")
        ds_results = benchmark_dataset(bed_file, name, output_dir)
        all_results.extend(ds_results)

    # Save results
    output_json = RESULTS_DIR / "benchmark_bed.json"
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
    print(f"\n{'='*60}")
    print("Benchmark Results Summary")
    print(f"{'='*60}")
    print(f"{'Dataset':<22} {'Tool':<15} {'mean(s)':<10} {'min(s)':<9} {'med(s)':<9} "
          f"{'Throughput':<14} {'Mem(MB)':<10} {'OK'}")
    print("-" * 95)
    for r in all_results:
        ok = "✓" if r.success else "✗"
        print(f"{r.dataset_name:<22} {r.tool:<15} {r.execution_time_sec:<10.3f} "
              f"{(r.min_time_sec if r.min_time_sec is not None else 0):<9.3f} "
              f"{(r.median_time_sec if r.median_time_sec is not None else 0):<9.3f} "
              f"{r.throughput_rec_per_sec:<14,.0f} {r.peak_memory_mb:<10.1f} {ok}")

    # Speedup summary per dataset. `execution_time_sec` is the mean, and on the
    # sub-second BED runs a single scheduling hiccup moves it by tens of percent,
    # so the min- and median-based speedups are reported alongside it.
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
        if min(fc.execution_time_sec, fc.min_time_sec or 0,
               fc.median_time_sec or 0, cm.execution_time_sec,
               cm.min_time_sec or 0, cm.median_time_sec or 0) <= 0:
            continue
        by_mean = cm.execution_time_sec / fc.execution_time_sec
        by_min = cm.min_time_sec / fc.min_time_sec
        by_med = cm.median_time_sec / fc.median_time_sec
        print(f"  {name:<20} {by_mean:7.1f}x {by_min:7.1f}x {by_med:7.1f}x")

    print("\nNext step: python paper/02b_benchmark_bed_multithread.py")


if __name__ == "__main__":
    main()
