#!/usr/bin/env python3
"""
14_bam_coverage.py - Measure the read count and coverage of every BAM/SAM input

Reviewer Minor 4 asks that BAM inputs not be described by file size alone: the
read count (and, ideally, read length and sequencing depth) is what a reader
needs to judge how heavy a liftover run actually was. This script measures those
numbers rather than transcribing them. Its output, results/bam_coverage.json, is
released with the other result files; it is the record of the measurement and is
not read by 19_generate_supp_table.py, so the statistics below live in the JSON
only.

For each alignment input this measures, with samtools:

  total_reads    -- every record (flagstat's "in total")
  mapped_reads   -- records with the unmapped flag (0x4) clear
  refs_covered   -- the @SQ reference sequences that carry at least one mapped read
  span_bp        -- the summed length (LN) of those references

and then computes:

  avg_read_length  -- from samtools stats (flagstat's "average length:",
                      which ignores clipping)
  depth_x          -- mapped_reads * avg_read_length / span_bp

The depth denominator is deliberately the *covered* span rather than a fixed
3.1 Gb genome: three of the four 1000G inputs are single-chromosome slices
(all reads on one chromosome), and dividing a chr11 slice by the whole genome
would understate its depth by a factor of ~23. Dividing by what is actually
covered is the honest number in both cases, and the span is reported alongside
so the reader can see which was used.

SAM inputs are not scanned separately. Every SAM in the benchmark is derived
from a BAM with `samtools view -h` (see DATA_SOURCES.md), so the alignments --
and therefore the read count, read length and covered span -- are identical by
construction. The script copies the sibling BAM's numbers and marks the row's
`source` as "derived from <bam>" so the provenance is visible in the output
rather than implicit.

Input:   paper/data/*.bam   (paths taken from the benchmark JSONs)
Output:  paper/results/bam_coverage.json
Usage:   python paper/14_bam_coverage.py
"""

import json
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

RESULTS_DIR = Path("paper/results")
BENCHMARK_FILES = [
    RESULTS_DIR / "benchmark_bam.json",
    RESULTS_DIR / "benchmark_sam.json",
]

# dataset names in benchmark_sam.json -> the BAM dataset they were derived from
SAM_FROM_BAM = {
    "K562 CTCF ChIP-seq": "K562 CTCF ChIP-seq",
    "1000G NA12878 chr11": "1000G NA12878 chr11",
    "1000G NA12878 chr20": "1000G NA12878 chr20",
    "1000G HG00096 chr11": "1000G HG00096 chr11",
}


def samtools() -> str:
    exe = shutil.which("samtools")
    if not exe:
        raise SystemExit(
            "samtools not found on PATH. The benchmark environment pins 1.24 "
            "(see zenodo/ENVIRONMENT.md); activate that environment first."
        )
    return exe


def header_lengths(sam: str, bam: Path) -> dict[str, int]:
    """@SQ SN -> LN, i.e. the reference set the file declares."""
    out = subprocess.run([sam, "view", "-H", str(bam)],
                         capture_output=True, text=True, check=True).stdout
    lengths = {}
    for line in out.splitlines():
        if not line.startswith("@SQ"):
            continue
        sn = ln = None
        for field in line.split("\t")[1:]:
            if field.startswith("SN:"):
                sn = field[3:]
            elif field.startswith("LN:"):
                ln = int(field[3:])
        if sn is not None and ln is not None:
            lengths[sn] = ln
    return lengths


def mapped_per_ref(sam: str, bam: Path) -> dict[str, int]:
    """reference -> number of mapped reads on it (unmapped flag 0x4 excluded)."""
    proc = subprocess.Popen([sam, "view", "-F", "4", str(bam)],
                            stdout=subprocess.PIPE, text=True)
    counts: dict[str, int] = {}
    for line in proc.stdout:
        if line.startswith("@"):
            continue
        ref = line.split("\t", 3)[2]
        counts[ref] = counts.get(ref, 0) + 1
    proc.stdout.close()
    if proc.wait() != 0:
        raise SystemExit(f"samtools view failed on {bam}")
    return counts


def avg_read_length(sam: str, bam: Path) -> float | None:
    """flagstat's 'average length' line (samtools stats), in bases."""
    out = subprocess.run([sam, "stats", str(bam)],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.startswith("SN\taverage length:"):
            return float(line.split("\t")[2])
    return None


def total_reads_in_flagstat(sam: str, bam: Path) -> int | None:
    out = subprocess.run([sam, "flagstat", str(bam)],
                         capture_output=True, text=True).stdout
    m = re.search(r"^(\d+) \+ \d+ in total", out, re.MULTILINE)
    return int(m.group(1)) if m else None


def measure(sam: str, bam: Path) -> dict:
    lengths = header_lengths(sam, bam)
    counts = mapped_per_ref(sam, bam)

    mapped_reads = sum(counts.values())
    covered = [c for c in counts if c in lengths]
    span_bp = sum(lengths[c] for c in covered)
    lost = sorted(set(counts) - set(lengths))

    alen = avg_read_length(sam, bam)
    total = total_reads_in_flagstat(sam, bam)

    depth = None
    if alen and span_bp:
        depth = mapped_reads * alen / span_bp

    return {
        "file": str(bam),
        "total_reads": total,
        "mapped_reads": mapped_reads,
        "avg_read_length": alen,
        "refs_declared": len(lengths),
        "refs_covered": len(covered),
        "covered_refs": sorted(covered),
        "span_bp": span_bp,
        "span_mb": round(span_bp / 1e6, 2),
        "depth_x": round(depth, 3) if depth is not None else None,
        "refs_in_reads_but_not_header": lost,
    }


def main():
    sam = samtools()
    print(f"using {sam}")

    # dataset name -> measured coverage, keyed by format so SAM can borrow
    by_fmt: dict[str, dict[str, dict]] = {"BAM": {}, "SAM": {}}
    seen_files: set[Path] = set()

    for bench_file in BENCHMARK_FILES:
        if not bench_file.exists():
            print(f"skip (missing): {bench_file}")
            continue
        fmt = "BAM" if "bam" in bench_file.name else "SAM"
        with open(bench_file) as fh:
            doc = json.load(fh)

        for ds in doc.get("datasets", []):
            name = ds.get("name", "")
            path = Path(ds.get("file", ""))
            if not path.exists():
                print(f"  {fmt} {name}: file missing ({path}), skipped")
                continue

            if fmt == "SAM":
                src = SAM_FROM_BAM.get(name)
                if not src or src not in by_fmt["BAM"]:
                    print(f"  SAM {name}: no sibling BAM measurement, skipped")
                    continue
                row = dict(by_fmt["BAM"][src])
                row["source"] = f"derived from {Path(row['file']).name}"
                by_fmt["SAM"][name] = row
                print(f"  SAM {name}: copied from BAM '{src}' "
                      f"({row['mapped_reads']:,} mapped reads)")
                continue

            print(f"  BAM {name}: scanning {path.name} ...", flush=True)
            row = measure(sam, path)
            row["source"] = "measured"
            by_fmt["BAM"][name] = row
            seen_files.add(path)
            print(f"      total {row['total_reads']:,}  "
                  f"mapped {row['mapped_reads']:,}  "
                  f"span {row['span_mb']:,} Mb over "
                  f"{row['refs_covered']} ref(s)  "
                  f"depth {row['depth_x']}x")

    out = {
        "measured": datetime.now().isoformat(),
        "tool": "samtools flagstat / stats / view (pinned in ENVIRONMENT.md)",
        "depth_formula": "mapped_reads * avg_read_length / span_bp",
        "span_note": ("span_bp is the summed length of the reference sequences "
                      "that actually carry reads -- a single chromosome for the "
                      "1000G slices, the whole genome for the ChIP-seq input"),
        "bam": by_fmt["BAM"],
        "sam": by_fmt["SAM"],
    }
    path = RESULTS_DIR / "bam_coverage.json"
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
