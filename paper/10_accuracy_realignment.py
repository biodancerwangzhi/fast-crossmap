#!/usr/bin/env python3
"""
10_accuracy_realignment.py - Compare FastCrossMap liftover vs BWA-MEM re-alignment accuracy

Uses wgsim simulated reads (ground truth in read names) to compare:
1. Re-alignment: FASTQ → BWA-MEM(hg38) → sorted BAM
2. Liftover: FASTQ → BWA-MEM(hg19) → BAM → FastCrossMap → hg38 BAM

Metrics:
- Chromosome concordance rate
- Position offset distribution (median, 95th percentile, max)
- Strand concordance rate
- Overall concordance rate (same chrom + position within threshold + same strand)

Usage: python3 paper/10_accuracy_realignment.py [--scale 1M] [--threads 4]
"""

import argparse
import json
import os
import statistics
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path

DATA_DIR = Path("paper/data")
SIM_DIR = DATA_DIR / "simulated"
RESULTS_DIR = Path("paper/results")
CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
HG19_FA = DATA_DIR / "hg19.fa"
HG38_FA = DATA_DIR / "hg38.fa"
FCM_BIN = Path("./target/release/fast-crossmap")

POSITION_THRESHOLDS = [0, 5, 10, 50, 100, 500]


@dataclass
class AccuracyResult:
    scale: str
    num_reads: int
    total_compared: int
    both_mapped: int
    chrom_concordance: float
    strand_concordance: float
    median_offset: float
    p95_offset: float
    p99_offset: float
    max_offset: int
    concordance_at_threshold: dict
    offset_distribution: dict


# Same reasoning as 09_benchmark_realignment.py: one BWA-MEM + sort run over a
# 10M-read file takes hours on a busy machine, so a 2 h cap fails it wrongly.
DEFAULT_TIMEOUT = 43200  # 12 h


def run_cmd(cmd, desc="", timeout=DEFAULT_TIMEOUT):
    print(f"  {desc}..." if desc else f"  Running: {cmd[:80]}...")
    # Explicit bash -o pipefail: a shell pipeline reports only its LAST command's
    # status, so `bwa mem ... | samtools sort ...` exits 0 even when bwa is
    # missing -- leaving a valid but empty BAM behind, which then looks like a
    # successful run.
    result = subprocess.run(["/bin/bash", "-o", "pipefail", "-c", cmd],
                            capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[:300]}")
        sys.exit(1)
    return result


def _is_empty_bam(path: Path) -> bool:
    """True when `path` exists but holds no alignment records.

    A BAM that exists with 0 reads is worse than a missing one: the `exists()`
    checks below would skip the re-alignment entirely and every downstream
    concordance number would be computed on nothing.
    """
    if not path.exists():
        return True
    try:
        r = subprocess.run(["samtools", "view", "-c", str(path)],
                           capture_output=True, text=True, timeout=1800)
        return int(r.stdout.strip()) == 0
    except (ValueError, FileNotFoundError, subprocess.SubprocessError):
        return False  # cannot tell; let the caller proceed


def prepare_bams(scale: str, num_reads: int, threads: int):
    """Prepare hg19 BAM, re-alignment BAM, and liftover BAM."""
    fq1 = SIM_DIR / f"sim_{scale}_R1.fq"
    fq2 = SIM_DIR / f"sim_{scale}_R2.fq"

    if not fq1.exists() or not fq2.exists():
        print(f"  Simulated reads not found for {scale}. Run 09_benchmark_realignment.py first.")
        sys.exit(1)

    # 1. Align to hg19
    hg19_bam = SIM_DIR / f"sim_{scale}_hg19_accuracy.bam"
    if _is_empty_bam(hg19_bam):
        run_cmd(
            f"bwa mem -t {threads} {HG19_FA} {fq1} {fq2} | samtools sort -@ {threads} -o {hg19_bam}",
            f"Aligning to hg19 ({scale})"
        )
        run_cmd(f"samtools index {hg19_bam}", "Indexing hg19 BAM")
    else:
        print(f"  hg19 BAM exists: {hg19_bam}")

    # 2. Re-align to hg38 (ground truth comparison)
    realign_bam = SIM_DIR / f"sim_{scale}_realign_hg38.bam"
    if _is_empty_bam(realign_bam):
        run_cmd(
            f"bwa mem -t {threads} {HG38_FA} {fq1} {fq2} | samtools sort -@ {threads} -o {realign_bam}",
            f"Re-aligning to hg38 ({scale})"
        )
        run_cmd(f"samtools index {realign_bam}", "Indexing re-alignment BAM")
    else:
        print(f"  Re-alignment BAM exists: {realign_bam}")

    # 3. Liftover hg19 → hg38
    liftover_bam = SIM_DIR / f"sim_{scale}_liftover_hg38.bam"
    if _is_empty_bam(liftover_bam):
        run_cmd(
            f"{FCM_BIN} bam {CHAIN_FILE} {hg19_bam} {liftover_bam}",
            f"FastCrossMap liftover ({scale})"
        )
    else:
        print(f"  Liftover BAM exists: {liftover_bam}")

    return realign_bam, liftover_bam


def parse_bam_to_dict(bam_path: Path) -> dict:
    """Parse BAM into dict: read_name -> (chrom, pos, strand, mapq)."""
    result = subprocess.run(
        f"samtools view {bam_path}",
        shell=True, capture_output=True, text=True, timeout=300
    )
    reads = {}
    for line in result.stdout.split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        if len(fields) < 11:
            continue
        name = fields[0]
        flag = int(fields[1])
        chrom = fields[2]
        pos = int(fields[3])
        mapq = int(fields[4])

        if flag & 4:  # unmapped
            continue
        if flag & 256 or flag & 2048:  # secondary or supplementary
            continue

        strand = '-' if flag & 16 else '+'
        is_read1 = bool(flag & 64)
        key = f"{name}/{'1' if is_read1 else '2'}"
        reads[key] = (chrom, pos, strand, mapq)

    return reads


def compare_bams(realign_reads: dict, liftover_reads: dict) -> dict:
    """Compare liftover vs re-alignment read by read."""
    common_keys = set(realign_reads.keys()) & set(liftover_reads.keys())

    chrom_match = 0
    strand_match = 0
    offsets = []
    concordance = {t: 0 for t in POSITION_THRESHOLDS}
    offset_buckets = Counter()

    for key in common_keys:
        r_chrom, r_pos, r_strand, r_mapq = realign_reads[key]
        l_chrom, l_pos, l_strand, l_mapq = liftover_reads[key]

        same_chrom = (r_chrom == l_chrom)
        same_strand = (r_strand == l_strand)
        offset = abs(r_pos - l_pos) if same_chrom else float('inf')

        if same_chrom:
            chrom_match += 1
            offsets.append(offset)
            for t in POSITION_THRESHOLDS:
                if offset <= t:
                    concordance[t] += 1

            if offset == 0:
                offset_buckets["0"] += 1
            elif offset <= 5:
                offset_buckets["1-5"] += 1
            elif offset <= 10:
                offset_buckets["6-10"] += 1
            elif offset <= 50:
                offset_buckets["11-50"] += 1
            elif offset <= 100:
                offset_buckets["51-100"] += 1
            elif offset <= 500:
                offset_buckets["101-500"] += 1
            else:
                offset_buckets[">500"] += 1

        if same_strand:
            strand_match += 1

    n = len(common_keys)
    return {
        "total_compared": n,
        "both_mapped": n,
        "chrom_match": chrom_match,
        "strand_match": strand_match,
        "offsets": sorted(offsets),
        "concordance": concordance,
        "offset_buckets": dict(offset_buckets),
        "realign_only": len(realign_reads) - n,
        "liftover_only": len(liftover_reads) - n,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scale", default="1M")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    scale = args.scale
    scales_map = {"1M": 1_000_000, "5M": 5_000_000, "10M": 10_000_000,
                  "20M": 20_000_000, "50M": 50_000_000}
    num_reads = scales_map.get(scale, 1_000_000)

    print("=" * 60)
    print(" Liftover vs Re-alignment Accuracy Comparison")
    print("=" * 60)
    print(f"  Scale: {scale} ({num_reads:,} reads)")
    print(f"  Threads: {args.threads}")
    print()

    # Step 1: Prepare BAMs
    print("[1/3] Preparing BAM files...")
    realign_bam, liftover_bam = prepare_bams(scale, num_reads, args.threads)
    print()

    # Step 2: Parse and compare
    print("[2/3] Parsing BAM files...")
    print("  Parsing re-alignment BAM...")
    realign_reads = parse_bam_to_dict(realign_bam)
    print(f"    {len(realign_reads):,} mapped reads")

    print("  Parsing liftover BAM...")
    liftover_reads = parse_bam_to_dict(liftover_bam)
    print(f"    {len(liftover_reads):,} mapped reads")
    print()

    print("[3/3] Comparing...")
    stats = compare_bams(realign_reads, liftover_reads)
    offsets = stats["offsets"]

    n = stats["total_compared"]
    chrom_rate = stats["chrom_match"] / n * 100 if n else 0
    strand_rate = stats["strand_match"] / n * 100 if n else 0
    med_offset = statistics.median(offsets) if offsets else 0
    p95 = sorted(offsets)[int(len(offsets) * 0.95)] if offsets else 0
    p99 = sorted(offsets)[int(len(offsets) * 0.99)] if offsets else 0
    max_offset = max(offsets) if offsets else 0

    conc_rates = {}
    for t in POSITION_THRESHOLDS:
        rate = stats["concordance"][t] / n * 100 if n else 0
        conc_rates[f"<={t}bp"] = round(rate, 3)

    print()
    print("=" * 60)
    print(" Results")
    print("=" * 60)
    print(f"  Re-alignment mapped reads:  {len(realign_reads):,}")
    print(f"  Liftover mapped reads:      {len(liftover_reads):,}")
    print(f"  Common reads compared:      {n:,}")
    print(f"  Re-alignment only:          {stats['realign_only']:,}")
    print(f"  Liftover only:              {stats['liftover_only']:,}")
    print()
    print(f"  Chromosome concordance:     {chrom_rate:.2f}%")
    print(f"  Strand concordance:         {strand_rate:.2f}%")
    print()
    print(f"  Position offset (same-chrom reads):")
    print(f"    Median:   {med_offset:.0f} bp")
    print(f"    95th %%:   {p95} bp")
    print(f"    99th %%:   {p99} bp")
    print(f"    Max:      {max_offset} bp")
    print()
    print(f"  Concordance at thresholds:")
    for label, rate in conc_rates.items():
        print(f"    {label:>10s}: {rate:.3f}%")
    print()
    print(f"  Offset distribution:")
    for bucket in ["0", "1-5", "6-10", "11-50", "51-100", "101-500", ">500"]:
        cnt = stats["offset_buckets"].get(bucket, 0)
        pct = cnt / len(offsets) * 100 if offsets else 0
        print(f"    {bucket:>8s} bp: {cnt:>8,} ({pct:.2f}%)")

    # Save results
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    result = AccuracyResult(
        scale=scale,
        num_reads=num_reads,
        total_compared=n,
        both_mapped=n,
        chrom_concordance=round(chrom_rate, 3),
        strand_concordance=round(strand_rate, 3),
        median_offset=med_offset,
        p95_offset=p95,
        p99_offset=p99,
        max_offset=max_offset,
        concordance_at_threshold=conc_rates,
        offset_distribution=stats["offset_buckets"],
    )
    out_file = RESULTS_DIR / f"accuracy_realignment_{scale}.json"
    with open(out_file, 'w') as f:
        json.dump({
            "benchmark": "liftover_vs_realignment_accuracy",
            "date": datetime.now().isoformat(),
            "result": asdict(result),
        }, f, indent=2)
    print(f"  Results saved to {out_file}")


if __name__ == "__main__":
    main()
