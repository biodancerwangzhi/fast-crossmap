#!/usr/bin/env python3
"""
02d_benchmark_gff.py - GFF format benchmark (multi-dataset)

Run FastCrossMap vs CrossMap vs liftOver benchmark for GFF format.
(FastRemap does not support GFF)

The input is read into page cache before each timed run: the tools run in a fixed
order, so without that warm-up the first tool pays for the disk read on any input
larger than free RAM and the later ones do not (see `warm_page_cache`).

Usage: python paper/02d_benchmark_gff.py
Output: paper/results/benchmark_gff.json
"""

import gzip
import json
import statistics
import os
import subprocess
import shutil
from dataclasses import dataclass, asdict
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

NUM_RUNS = 5

GFF_DATASETS = [
    {
        "name": "GENCODE v44 full",
        "file": DATA_DIR / "gff_gencode_v44_grch37.gff3.gz",
        "source": "GENCODE v44 GRCh37 mapped",
    },
    {
        "name": "GENCODE v44 basic",
        "file": DATA_DIR / "gff_gencode_v44_basic_grch37.gff3.gz",
        "source": "GENCODE v44 basic GRCh37",
    },
    {
        "name": "GENCODE v44 lncRNA",
        "file": DATA_DIR / "gff_gencode_v44_lncrna_grch37.gff3.gz",
        "source": "GENCODE v44 lncRNA GRCh37",
    },
    {
        "name": "RefSeq GRCh37",
        "file": DATA_DIR / "gff_refseq_grch37.gff.gz",
        "source": "NCBI RefSeq GRCh37",
    },
    {
        "name": "Ensembl GRCh37",
        "file": DATA_DIR / "gff_ensembl_grch37.gff3.gz",
        "source": "Ensembl release 87 GRCh37",
    },
]


@dataclass
class BenchmarkResult:
    tool: str
    format: str
    dataset_name: str
    input_file: str
    input_features: int
    input_size_mb: float
    execution_time_sec: float
    throughput_features_per_sec: float
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
    # Records in the converted output, and how many of the output's records the
    # input had (`##FASTA` sequences and comment lines are not counted on either
    # side). A run that maps almost nothing is a different failure from a slow
    # run, and only a count can tell them apart.
    output_records: int = 0
    mapped_percent: Optional[float] = None


# UCSC long names, the spelling the hg19->hg38 chain file uses, which is also
# what the GENCODE datasets already ship.
_CHR_NAMES = ({str(n): f"chr{n}" for n in range(1, 23)}
              | {"X": "chrX", "Y": "chrY", "MT": "chrM"})
# RefSeq GRCh37 accessions -> UCSC name.
_REFSEQ_NAMES = {
    "NC_000001.10": "chr1", "NC_000002.11": "chr2", "NC_000003.11": "chr3",
    "NC_000004.11": "chr4", "NC_000005.10": "chr5", "NC_000006.11": "chr6",
    "NC_000007.13": "chr7", "NC_000008.10": "chr8", "NC_000009.11": "chr9",
    "NC_000010.10": "chr10", "NC_000011.9": "chr11", "NC_000012.11": "chr12",
    "NC_000013.10": "chr13", "NC_000014.8": "chr14", "NC_000015.9": "chr15",
    "NC_000016.9": "chr16", "NC_000017.10": "chr17", "NC_000018.9": "chr18",
    "NC_000019.9": "chr19", "NC_000020.10": "chr20", "NC_000021.8": "chr21",
    "NC_000022.10": "chr22", "NC_000023.10": "chrX", "NC_000024.9": "chrY",
    "NC_012920.1": "chrM",
}


def count_gff_features(gff_file: Path) -> int:
    count = 0
    if str(gff_file).endswith('.gz'):
        opener, mode = gzip.open, 'rt'
    else:
        opener, mode = open, 'r'
    with opener(gff_file, mode) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                count += 1
    return count


def get_file_size_mb(file_path: Path) -> float:
    return file_path.stat().st_size / (1024 * 1024)


def normalise_seqids(gff_file: Path) -> Path:
    """Rewrite the seqids in `gff_file` to UCSC long names, if they need it.

    Two of the five GFF datasets name their sequences in ways the hg19->hg38
    chain does not use, so *every* record fails to map:

      RefSeq  `NC_000001.10` / `NT_167207.1`   (RefSeq accessions)
      Ensembl `1` / `MT`                       (bare names)

    `--chromid` cannot fix this: it is an *output* style switch, so `--chromid l`
    rewrites `1` to `chr1` only after a mapping that has already failed. It is
    the *input* names that have to match the chain.

    The normalised copy is what every tool in the benchmark reads, so this is
    input preparation sitting next to the `.gz` decompression rather than a
    per-tool special case, and the comparison stays like-for-like. `NT_*`/`NW_*`
    unplaced scaffolds have no `chr*` equivalent and are left alone; those
    records stay unmapped, which is the honest outcome.

    Nothing is rewritten in place: the normalised file lands beside the original
    as `<stem>.renamed<ext>`, which also makes a second run a no-op.
    """
    if gff_file.suffix not in ('.gff', '.gff3'):
        return gff_file
    if gff_file.name.startswith('gff_refseq'):
        table = _REFSEQ_NAMES
    elif gff_file.name.startswith('gff_ensembl'):
        table = _CHR_NAMES
    else:
        return gff_file  # GENCODE already uses chr* names

    out = gff_file.with_name(gff_file.stem + '.renamed' + gff_file.suffix)
    if out.exists():
        return out

    renamed = kept = total = 0
    seen = set()
    with open(gff_file, 'r') as fin, open(out, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            total += 1
            tab = line.find('\t')
            seqid = line[:tab] if tab > 0 else ''
            seen.add(seqid)
            if seqid in table:
                fout.write(table[seqid] + line[tab:])
                renamed += 1
            else:
                fout.write(line)
                kept += 1

    print(f"    Seqid normalisation -> {out.name}: "
          f"{len({k for k in seen if k in table})} names renamed, "
          f"{renamed:,}/{total:,} records rewritten, "
          f"{len(seen - set(table))} seqid(s) left as-is")
    stray = sorted(k for k in seen - set(table)
                   if not k.startswith(('NT_', 'NW_')))
    if stray:
        print(f"      WARNING: seqids with no chr* equivalent: {stray[:8]}")
    return out


def prepare_input(gff_gz: Path) -> Path:
    """Decompressed, seqid-normalised copy of `gff_gz` -- what the tools read."""
    path = ensure_decompressed(gff_gz)
    renamed = normalise_seqids(path)
    if renamed != path:
        print(f"    Using seqid-normalised input: {renamed.name}")
    return renamed


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
    """Records in `path`, using the per-format rules named in the docstring of
    `count_output_records`."""
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

    Counting bytes cannot tell a lifted-over file from an empty header. Two of
    the five GFF datasets name their sequences in a way the chain file does not
    (`NC_000001.10`, and Ensembl's bare `1`), and *no* tool maps anything from
    them: CrossMap's RefSeq output was 0 bytes and FastCrossMap's was 34 KB of
    comment lines, yet the latter was still recorded as a 27.51 s success. A
    record count is what separates "this tool lost on this dataset" from "this
    tool was merely slower".

    What counts as a record differs per format:

      * BED: no header, every line is a record.
      * GFF / VCF: `#` lines are headers; in a GFF a `##FASTA` section and the
        sequence under it are not records either.
      * SAM: `@` lines are headers, as are `#` lines.
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
            # A bgzip member is a compressed block, not a record.
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
    """True when the timed command left a freshly written, non-empty output.

    `run_with_time` only looks at the exit status, and every tool here can exit
    0 without producing usable output -- an unrecognised flag, a refused
    overwrite, a truncated stream, or an input whose seqids match nothing in
    the chain file. Requiring a new mtime as well as a non-zero size also
    rejects a stale file left behind by an earlier run.
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


def ensure_decompressed(gff_file: Path) -> Path:
    if not str(gff_file).endswith('.gz'):
        return gff_file
    import shutil
    decompressed = Path(str(gff_file)[:-3])
    if not decompressed.exists():
        print(f"    Decompressing {gff_file.name} ...")
        with gzip.open(gff_file, 'rb') as f_in, open(decompressed, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return decompressed


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


def run_tool_benchmark(tool_name: str, cmd: list[str], num_features: int,
                       size_mb: float, dataset_name: str, input_file: str,
                       outputs: list = (), warm_file: Optional[Path] = None,
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

    # Logged for every row: "0 records" is the fingerprint of a run that could
    # not map anything, which is a different failure from a slow run and has to
    # be visible without reading the output files afterwards.
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
            tool=tool_name, format="GFF", dataset_name=dataset_name,
            input_file=input_file, input_features=num_features,
            input_size_mb=size_mb, execution_time_sec=0,
            throughput_features_per_sec=0, peak_memory_mb=0,
            all_times=[], success=False, error_message=output_err or error_msg
        )

    avg_time = sum(times) / len(times)
    avg_memory = sum(memories) / len(memories)
    min_time = min(times)
    median_time = statistics.median(times)
    # What the tool actually wrote, against the records the input offers: the
    # two together give the mapped fraction, which the pass/fail check alone
    # cannot express (a tool can write 99% of a file and still be wrong).
    output_records = sum(
        int(part.rsplit(': ', 1)[1].replace(' records', '').replace(',', ''))
        for part in counts.split('; ') if part.endswith(' records'))

    return BenchmarkResult(
        tool=tool_name, format="GFF", dataset_name=dataset_name,
        input_file=input_file, input_features=num_features,
        input_size_mb=round(size_mb, 2),
        execution_time_sec=round(avg_time, 2),
        throughput_features_per_sec=round(num_features / avg_time, 0) if num_features > 0 else 0,
        peak_memory_mb=round(avg_memory, 2),
        all_times=times, success=success,
        all_disk_read_kb=disk_reads, all_major_faults=major_faults_list,
        min_time_sec=round(min_time, 2),
        median_time_sec=round(median_time, 2),
        min_throughput=round(num_features / min_time, 0) if num_features > 0 else 0,
        median_throughput=round(num_features / median_time, 0) if num_features > 0 else 0,
        std_time_sec=round(statistics.pstdev(times), 2) if len(times) > 1 else 0.0,
        output_records=output_records,
        mapped_percent=round(100 * output_records / num_features, 1) if num_features > 0 else None,
    )


def benchmark_dataset(gff_gz: Path, dataset_name: str, output_dir: Path) -> list[BenchmarkResult]:
    num_features = count_gff_features(gff_gz)
    size_mb = get_file_size_mb(gff_gz)
    print(f"  Features: {num_features:,} | Size: {size_mb:.1f} MB")

    gff_file = prepare_input(gff_gz)
    ds_dir = output_dir / gff_gz.stem.replace('.gff3', '').replace('.gff', '')
    ds_dir.mkdir(parents=True, exist_ok=True)
    results = []

    # 1. FastCrossMap
    print(f"  [1/3] FastCrossMap on {dataset_name}")
    out = ds_dir / "fastcrossmap_output.gff3"
    cmd = ["./target/release/fast-crossmap", "gff",
           str(CHAIN_FILE), str(gff_file), str(out)]
    results.append(run_tool_benchmark(
        "FastCrossMap", cmd, num_features, size_mb, dataset_name, str(gff_gz),
        outputs=[out], warm_file=gff_file))

    # 2. CrossMap
    print(f"  [2/3] CrossMap on {dataset_name}")
    out = ds_dir / "crossmap_output.gff3"
    cmd = ["CrossMap", "gff", str(CHAIN_FILE), str(gff_file), str(out)]
    results.append(run_tool_benchmark(
        "CrossMap", cmd, num_features, size_mb, dataset_name, str(gff_gz),
        outputs=[out], warm_file=gff_file))

    # 3. liftOver (supports GFF via -gff flag)
    print(f"  [3/3] liftOver on {dataset_name}")
    out = ds_dir / "liftover_output.gff3"
    unmap = ds_dir / "liftover_output.gff3.unmap"
    cmd = ["liftOver", "-gff", str(gff_file), str(CHAIN_FILE), str(out), str(unmap)]
    results.append(run_tool_benchmark(
        "liftOver", cmd, num_features, size_mb, dataset_name, str(gff_gz),
        outputs=[out], warm_file=gff_file))

    cleanup_outputs(ds_dir)
    return results


def main():
    print("=" * 60)
    print("GFF Format Benchmark (Multi-Dataset)")
    print("=" * 60)
    print("Note: FastRemap does not support GFF")
    print("Note: RefSeq/Ensembl seqids are normalised to chr* before timing")
    print("Note: the input is read into page cache before every timed run\n")

    if not CHAIN_FILE.exists():
        print(f"Error: Chain file not found: {CHAIN_FILE}")
        print("Please run first: bash paper/01_download_data.sh")
        return

    datasets = []
    for ds in GFF_DATASETS:
        f = ds["file"]
        if f.exists():
            datasets.append((f, ds["name"], ds["source"]))
        else:
            print(f"Warning: {ds['name']} not found at {f}, skipping")

    if not datasets:
        print("Error: No GFF datasets found. Run bash paper/01_download_data.sh")
        return

    output_dir = RESULTS_DIR / "gff_benchmark"
    output_dir.mkdir(parents=True, exist_ok=True)

    all_results = []
    for gff_file, name, source in datasets:
        print(f"\n{'─'*60}")
        print(f"Dataset: {name} ({source})")
        print(f"File: {gff_file}")
        print(f"{'─'*60}")
        ds_results = benchmark_dataset(gff_file, name, output_dir)
        all_results.extend(ds_results)

    output_json = RESULTS_DIR / "benchmark_gff.json"
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
    print(f"{'Dataset':<22} {'Tool':<15} {'mean(s)':<10} {'min(s)':<9} {'med(s)':<9} {'Features/s':<14} {'mapped%':<9} {'Mem(MB)':<10} {'OK'}")
    print("-" * 102)
    for r in all_results:
        ok = "✓" if r.success else "✗"
        pct = f"{r.mapped_percent:.1f}" if r.mapped_percent is not None else "-"
        print(f"{r.dataset_name:<22} {r.tool:<15} {r.execution_time_sec:<10.2f} {(r.min_time_sec or 0):<9.2f} {(r.median_time_sec or 0):<9.2f} "
              f"{r.throughput_features_per_sec:<14,.0f} {pct:<9} {r.peak_memory_mb:<10.1f} {ok}")

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

    print("\nNext step: python paper/04_plot_performance.py")


if __name__ == "__main__":
    main()
