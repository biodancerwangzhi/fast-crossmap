#!/usr/bin/env python3
"""
figS3_threads.py - Figure S3: FastCrossMap multi-threading scalability

CrossMap has no threading, so the paper's cross-tool baseline is -t 1. This
figure is therefore not a comparison against another tool -- it answers a
different question: "what does FastCrossMap itself buy the user who turns on
threads?" It backs the multi-threading row of Table 1.

Eight formats, each run at -t 1/2/4/8/16, five repeats per cell (n=5). The
median wall time is used throughout; the per-cell spread (max-min as a
percentage of the median) is under 16% everywhere in the merged table, so the
curves are not chasing noise.

Two panels, side by side, both against thread count on a log2 axis:
  (a) speed-up relative to -t 1 -- a horizontal dashed line at 1.0 is the
      break-even mark; a format that is never slower than -t 1 stays at or
      above it for every thread count, which is the claim being made.
  (b) parallel efficiency -- speed-up divided by thread count, i.e. the
      fraction of each added core that turns into wall-clock. It decays for
      every format: this is a memory-bandwidth-bound liftover, not a compute
      kernel, so efficiency is the honest framing and the reason the paper
      does not quote a "scales to 16 cores" number.

Thread count is spaced 1/2/4/8/16, so the x-axis is log2 to keep the spacing
even; a linear axis would crowd everything below 4.

DATA PROVENANCE. This script reads the merged summary produced by
paper/11_benchmark_threads.sh (paper/results/threads-merged-20261011/); the
merge itself is documented in that directory's env.txt, including the fact that
the MAF row was re-run with a tuned batch size. That detail is deliberately kept
out of the figure and its caption: multi-threading is a secondary property of
this tool, and a batch-size footnote invites more explanation than the figure is
worth. It stays in env.txt for anyone reproducing the numbers.

check_provenance() below asserts the shape of what was loaded (all eight
formats, every cell at or above 1.0x) so that a future re-run that changes the
input is loud instead of silently absorbed into the curve.

Input:   paper/results/threads-merged-20261011/summary.tsv
Output:  paper/figures/figS3_threads.pdf
         paper/figures/figS3_threads.png

Usage: python paper/plot/figS3_threads.py
"""

import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

RESULTS_DIR = Path("paper/results")
FIGURES_DIR = Path("paper/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY = RESULTS_DIR / "threads-merged-20261011" / "summary.tsv"

# (tsv key, display label, colour, z-order) -- all eight formats the tool handles
FORMATS = [
    ("bed",    "BED",     "#1f77b4"),
    ("vcf",    "VCF",     "#ff7f0e"),
    ("gff",    "GFF",     "#2ca02c"),
    ("wig",    "Wiggle",  "#d62728"),
    ("bigwig", "BigWig",  "#9467bd"),
    ("maf",    "MAF",     "#8c564b"),
    ("bam",    "BAM",     "#e377c2"),
    ("sam",    "SAM",     "#7f7f7f"),
]

THREADS = [1, 2, 4, 8, 16]


def load_summary():
    """summary.tsv -> {format: {threads: row_dict}}"""
    data = {}
    with open(SUMMARY, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fmt = row["format"].strip()
            if fmt not in {k for k, _, _ in FORMATS}:
                continue
            data.setdefault(fmt, {})[int(row["threads"])] = {
                "median_s": float(row["median_s"]),
                "min_s": float(row["min_s"]),
                "max_s": float(row["max_s"]),
                "speedup": float(row["speedup_vs_t1"]),
                "eff": float(row["eff"]),
                "n": int(row["n"]),
            }
    return data


def check_provenance(data):
    """
    Assert the shape this figure depends on: every cell at or above 1.00x vs
    -t 1 (the tuned-batch MAF run is monotone by construction), and all eight
    formats present. A regression here means the input is not the run the
    figure assumes -- fail loudly rather than draw a curve that dips.
    """
    missing = [k for k, _, _ in FORMATS if k not in data]
    if missing:
        raise SystemExit(
            f"summary.tsv is missing formats {missing}; expected all eight. "
            "The figure's claim ('no format is slower at any thread count') "
            "cannot be made from a partial table."
        )

    regressions = []
    for key, label, _ in FORMATS:
        for t, r in sorted(data.get(key, {}).items()):
            if r["speedup"] < 1.0:
                regressions.append((label, t, r["speedup"]))
    if regressions:
        raise SystemExit(
            "regression(s) below 1.00x vs -t 1 in summary.tsv:\n  "
            + "\n  ".join(f"{lab} -t {t}: {s:.2f}x" for lab, t, s in regressions)
            + "\n\nThis is not the tuned-batch merged run the figure documents "
            "(see env.txt: MAF was re-run at FCM_BATCH=65536). Do not draw the "
            "curve; re-check the input."
        )

    print("provenance OK: all eight formats present, no cell below 1.00x.")
    print("(MAF here is the FCM_BATCH=65536 run; the default-config -t 2")
    print(" regression is documented separately and is not in this table.)")
    return data


def panel_speedup(ax, data):
    ax.axhline(1.0, color="#888888", ls="--", lw=1.2, zorder=1)
    # right-aligned just above the break-even line: every curve is above this
    # band at high thread counts, so the label sits in empty space
    ax.text(16, 1.015, "break-even (-t 1)", color="#666666",
            fontsize=8, va="bottom", ha="right")

    for key, label, colour in FORMATS:
        if key not in data:
            continue
        ts = [t for t in THREADS if t in data[key]]
        ys = [data[key][t]["speedup"] for t in ts]
        ax.plot(ts, ys, "-o", color=colour, label=label,
                lw=1.8, ms=5, mec="white", mew=0.8)

    ax.set_xscale("log", base=2)
    ax.set_xticks(THREADS)
    ax.set_xticklabels([str(t) for t in THREADS])
    ax.set_xlabel("threads (-t)")
    ax.set_ylabel("speed-up vs. -t 1  (median wall time)")
    ax.set_title("(a) what turning on threads buys", fontsize=10, loc="left")
    ax.grid(alpha=0.25, which="both", lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8, ncol=2, framealpha=0.9, loc="upper left")


def panel_efficiency(ax, data):
    for key, label, colour in FORMATS:
        if key not in data:
            continue
        ts = [t for t in THREADS if t in data[key]]
        ys = [data[key][t]["eff"] for t in ts]
        ax.plot(ts, ys, "-o", color=colour, lw=1.8, ms=5,
                mec="white", mew=0.8)

    ax.set_xscale("log", base=2)
    ax.set_xticks(THREADS)
    ax.set_xticklabels([str(t) for t in THREADS])
    ax.set_xlabel("threads (-t)")
    ax.set_ylabel("parallel efficiency  (speed-up / threads)")
    ax.set_title("(b) how much of each added core lands in wall time",
                 fontsize=10, loc="left")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.25, which="both", lw=0.5)
    ax.set_axisbelow(True)


def main():
    if not SUMMARY.exists():
        raise SystemExit(f"missing input: {SUMMARY}\n"
                         "run paper/11_benchmark_threads.sh first")

    data = load_summary()
    check_provenance(data)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.6))
    panel_speedup(ax1, data)
    panel_efficiency(ax2, data)

    fig.suptitle("Figure S3. FastCrossMap multi-threading scalability "
                 "(8 formats, median of n=5)", fontsize=11, y=1.00)
    fig.tight_layout()

    out_pdf = FIGURES_DIR / "figS3_threads.pdf"
    out_png = FIGURES_DIR / "figS3_threads.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"wrote {out_pdf}")
    print(f"wrote {out_png}")

    for key, label, _ in FORMATS:
        if key in data:
            best_t = max(THREADS, key=lambda t: data[key].get(t, {}).get("speedup", 0))
            print(f"  {label:<8} peak speed-up {data[key][best_t]['speedup']:.2f}x "
                  f"at -t {best_t}")


if __name__ == "__main__":
    main()
