#!/usr/bin/env python3
"""
fig1_performance.py - Figure 1: FastCrossMap vs CrossMap, all eight formats

Two panels, one bar per supported format, both read straight from the eight
released `benchmark_*.json` files:

  (a) speed-up   = CrossMap time / FastCrossMap time       (both -t 1)
  (b) memory cut = CrossMap peak RSS / FastCrossMap peak RSS (both -t 1)

Why bars of a *ratio* rather than boxes of absolute time and memory
------------------------------------------------------------------
Each format was benchmarked on four or five datasets, so there is a
distribution behind every format and the first instinct is a box plot. The
distributions do not support one: across the eight formats the per-dataset
ratios span at most 1.8x (VCF) and are 1.0-1.2x for most, so a box would
collapse to a hairline while implying uncertainty the data does not have.
What matters here is the *summary* the reader takes away -- "how much faster,
how much lighter, in each format" -- so each bar is the median of the
per-dataset PAIRED ratios and a thin whisker carries the observed min and max.
Pairing first (CrossMap_ds / FastCrossMap_ds on the same dataset) is the
stricter reading and the one Table S1 uses: it asks "on a given dataset, how
much better?", so one heavy dataset in a format cannot move that format's bar.

Why log
-------
The memory cut runs from ~2.0x (MAF) to ~63x (BAM). On a linear axis the four
light formats would be four slivers at the baseline while BAM fills the panel.
A log axis also makes the ratio readable directly: equal factors are equal
distances, and the 1.0x line -- the break-even mark, "no better than CrossMap"
-- is a fixed reference the eye can measure every bar against.

Both panels use the -t 1 baseline, the paper's uniform cross-tool baseline.
CrossMap has no threading, and the BAM/SAM JSONs also carry
`FastCrossMap (4t)` / `(8t)` rows; those are dropped here, as they are in
Figure S4 and Table S1.

BAM vs SAM, and what CrossMap's BAM number contains
---------------------------------------------------
BAM shows a ~63x memory cut while SAM shows ~3x on the same alignments, which
invites the question of whether that is a measurement artifact. It is not an
artifact, and it is not a difference in what was lifted over. Two things:

  * CrossMap converts a BAM *input* by shelling out to `pysam.sort` + `pysam.index`
    on its output (cmmodule/mapbam.py), because liftover moves reads to new
    coordinates and the output is then no longer coordinate-sorted. That sort is
    a real child process with its own ~768 MB heap floor, and it is counted in
    CrossMap's peak RSS because it is part of running CrossMap on BAM. CrossMap's
    SAM-input path does not sort, so SAM carries no such cost.
  * FastCrossMap's `bam` output preserves input order and is not indexed, so it
    does no equivalent work.

The consequence for this figure is that the BAM bars compare "CrossMap,
including the sort+index it must perform on a BAM input" against "FastCrossMap,
as it is". That is the honest reading of what a user runs, and it is stated in
the caption; SAM is the format on which both tools do the same amount of work,
and there the cut is ~3x.

Claim checked before drawing
----------------------------
The figure asserts that FastCrossMap is both faster and lighter than CrossMap
in every one of the eight formats. `check_claim()` verifies exactly that from
the data and refuses to draw if any format's median ratio is <= 1.0, so the
chart cannot be published asserting something the JSONs stopped showing.

Input:   paper/results/benchmark_{bed,bam,vcf,gff,wig,bigwig,maf,sam}.json
Output:  paper/figures/fig1_performance.pdf
         paper/figures/fig1_performance.png
         paper/results/fig1_summary.json   (exactly what was plotted)

Usage: python paper/plot/fig1_performance.py
"""

import json
import statistics as stats
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np

RESULTS_DIR = Path("paper/results")
FIGURES_DIR = Path("paper/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# (display label, benchmark JSON) -- all eight formats the tool handles.
# Panel order is set at draw time (descending by speed-up), see main().
FORMATS = [
    ("BED",     "benchmark_bed.json"),
    ("BAM",     "benchmark_bam.json"),
    ("VCF",     "benchmark_vcf.json"),
    ("GFF",     "benchmark_gff.json"),
    ("Wiggle",  "benchmark_wig.json"),
    ("BigWig",  "benchmark_bigwig.json"),
    ("MAF",     "benchmark_maf.json"),
    ("SAM",     "benchmark_sam.json"),
]

REFERENCE_TOOL = "CrossMap"          # the tool every bar is a ratio against
FCM = "FastCrossMap"                 # the unthreaded row, -t 1
BAR_COLOR = "#4d4d4d"                # greyscale, matching the submitted Figure 1;
                                     # the bars are a single series, so colour
                                     # carried no information here
WHISKER_COLOR = "#555555"


def load_series():
    """
    {format: {tool: {dataset: {"time": s, "mem": MB}}}} for -t 1 rows only.

    `execution_time_sec` and `peak_memory_mb` are the same two fields Table S1
    prints, so the figure and the table cannot disagree about a number.
    Threaded FastCrossMap rows (`FastCrossMap (4t)` / `(8t)`) are dropped: this
    figure's baseline is -t 1 throughout.
    """
    out = {}
    for label, fname in FORMATS:
        path = RESULTS_DIR / fname
        if not path.exists():
            raise SystemExit(f"missing input: {path}")

        with open(path) as fh:
            doc = json.load(fh)

        per_tool: dict[str, dict[str, dict]] = {}
        for r in doc.get("results", []):
            tool = r.get("tool")
            if not tool or not r.get("success"):
                continue
            if tool.startswith("FastCrossMap ("):
                continue                      # -t 4 / -t 8: wrong baseline
            t = r.get("execution_time_sec")
            m = r.get("peak_memory_mb")
            if not t or t <= 0 or not m or m <= 0:
                continue
            per_tool.setdefault(tool, {})[r["dataset_name"]] = {"time": t, "mem": m}

        if FCM not in per_tool or REFERENCE_TOOL not in per_tool:
            raise SystemExit(
                f"{path.name}: needs both {FCM} and {REFERENCE_TOOL} rows, "
                f"found {sorted(per_tool)}"
            )
        out[label] = per_tool
    return out


def pair_stats(series, field: str) -> dict:
    """
    Per-dataset paired ratios, per format, over the datasets BOTH tools mapped.

    {format: {"median", "min", "max", "n", "cross", "fcm"}} where "cross" and
    "fcm" are the two tools' per-dataset medians (absolute units), kept so the
    summary JSON records the underlying numbers and not only the ratio.
    """
    out = {}
    for label, _ in FORMATS:
        ref = series[label][REFERENCE_TOOL]
        fcm = series[label][FCM]
        shared = sorted(set(ref) & set(fcm))
        if not shared:
            raise SystemExit(f"{label}: no dataset was mapped by both tools")
        ratios = [ref[ds][field] / fcm[ds][field] for ds in shared]
        out[label] = {
            "median": stats.median(ratios),
            "min": min(ratios),
            "max": max(ratios),
            "n": len(ratios),
            "datasets": shared,
            "cross_median": stats.median(ref[ds][field] for ds in shared),
            "fcm_median": stats.median(fcm[ds][field] for ds in shared),
        }
    return out


def check_claim(speed, memory) -> list[str]:
    """
    The figure's claim: FastCrossMap is faster AND lighter than CrossMap in
    every format. Return the formats where it is not, so the caller can refuse
    to draw rather than publish a figure whose caption is false.
    """
    violations = []
    for label, _ in FORMATS:
        for name, table in (("speed-up", speed), ("memory cut", memory)):
            med = table[label]["median"]
            if med <= 1.0:
                violations.append(f"{label}: {name} is {med:.2f}x (<= 1.0x)")
    return violations


def draw_panel(ax, label, stats_map, order, ylabel, title, note=None):
    """One panel: horizontal-free grouped bars over `order`, log y, 1.0x ref."""
    pos = np.arange(len(order))
    medians = [stats_map[f]["median"] for f in order]
    lows = [stats_map[f]["median"] - stats_map[f]["min"] for f in order]
    highs = [stats_map[f]["max"] - stats_map[f]["median"] for f in order]

    ax.bar(pos, medians, width=0.62, color=BAR_COLOR, alpha=0.85,
           edgecolor="#333333", linewidth=0.8, zorder=3)
    ax.errorbar(pos, medians, yerr=[lows, highs], fmt="none",
                ecolor=WHISKER_COLOR, elinewidth=1.0, capsize=3.5, zorder=4)

    for p, (f, med) in enumerate(zip(order, medians)):
        # above the upper whisker cap, so a tight whisker (BAM: max/median 1.04)
        # cannot put the cap through the label. Two decimals, matching the
        # caption and the manuscript text (13.09x, not 13.1x).
        ax.text(p, stats_map[f]["max"] * 1.07, f"{med:,.2f}x",
                ha="center", va="bottom", fontsize=11, color="#222222", zorder=5)

    ax.set_yscale("log")
    hi = max(stats_map[f]["max"] for f in order)
    # Floor below 1.0 on both panels, so the break-even line is always inside
    # the axes (a log axis cannot reach 0, and bars are drawn from the floor).
    # The same fixed floor in both panels also lets the two be read alike.
    ax.set_ylim(0.85, hi * 2.4)              # room for the whisker caps + labels
    ax.yaxis.set_major_locator(
        ticker.LogLocator(base=10.0, subs=(1.0, 2.0, 5.0), numticks=12))
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:g}x"))

    # break-even: "no better than CrossMap"; the note sits in the headroom at
    # the top of the panel, the only spot no bar reaches
    ax.axhline(1.0, color="#888888", lw=1.0, ls="--", zorder=1)
    ax.text(0.015, 0.97, "dashed line = 1.0x, on par with CrossMap",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=10, color="#666666", zorder=2)

    ax.set_xticks(pos)
    ax.set_xticklabels([f"{f}\nn={stats_map[f]['n']}" for f in order],
                       fontsize=11)
    ax.set_ylabel(ylabel, fontsize=12.5)
    ax.set_title(title, fontsize=14, loc="left", pad=6)
    if note:
        ax.text(0.0, 1.0, note, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=10, color="#666666")
    ax.grid(alpha=0.25, axis="y", which="major", lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def main():
    series = load_series()
    speed = pair_stats(series, "time")
    memory = pair_stats(series, "mem")

    violations = check_claim(speed, memory)
    if violations:
        raise SystemExit(
            "the figure's claim (FastCrossMap is faster and lighter than "
            "CrossMap in every format) does not hold:\n  "
            + "\n  ".join(violations)
            + "\n\nDo not draw the chart; the caption would be false."
        )
    print("claim holds: FastCrossMap is faster and lighter than CrossMap in all "
          f"{len(FORMATS)} formats")

    # one order for both panels, descending by speed-up, so a format sits in the
    # same column in (a) and (b) and the eye can read across
    order = sorted((f for f, _ in FORMATS), key=lambda f: -speed[f]["median"])

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11.5, 5.0))

    draw_panel(ax_a, "a", speed, order,
               ylabel="Speedup factor over CrossMap (wall-clock)",
               title="(a) Speed-up over CrossMap")
    draw_panel(ax_b, "b", memory, order,
               ylabel="Peak memory reduction factor (RSS)",
               title="(b) Peak-memory reduction vs CrossMap")

    fig.text(0.5, 0.012,
             "both panels: -t 1, all datasets of each format; bar = median of "
             "the per-dataset paired ratios, whisker = observed min-max",
             ha="center", fontsize=11, color="#333333")

    fig.tight_layout(rect=(0, 0.05, 1, 0.95))

    out_pdf = FIGURES_DIR / "fig1_performance.pdf"
    out_png = FIGURES_DIR / "fig1_performance.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"wrote {out_pdf}")
    print(f"wrote {out_png}")

    payload = {
        "baseline": "-t 1",
        "reference": REFERENCE_TOOL,
        "ratio_definition": "median of per-dataset paired ratios "
                            "(crossmap / fastcrossmap) on the same dataset",
        "whisker": "observed min-max of the per-dataset paired ratios",
        "panel_order": order,
        "speedup": {f: speed[f] for f in order},
        "memory_cut": {f: memory[f] for f in order},
    }
    out_json = RESULTS_DIR / "fig1_summary.json"
    with open(out_json, "w") as fh:
        json.dump(payload, fh, indent=2)
        fh.write("\n")
    print(f"wrote {out_json}")

    print("\npanel (a) speed-up, panel (b) memory cut")
    print(f"  {'format':8s} {'n':>2s} {'speedup':>9s} {'memory':>9s}")
    for f in order:
        print(f"  {f:8s} {speed[f]['n']:2d} "
              f"{speed[f]['median']:8.2f}x {memory[f]['median']:8.2f}x")


if __name__ == "__main__":
    main()
