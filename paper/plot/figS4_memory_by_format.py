#!/usr/bin/env python3
"""
figS4_memory_by_format.py - Figure S4: peak memory by format, one subplot each

Eight small multiples, one per supported format. Inside a subplot the y-axis
lists the tools benchmarked on that format and the x-axis is peak RSS on a log
scale, so each tool gets a full-width box instead of being squeezed into a
single group of an eight-row chart.

Every box is built from that format's datasets -- 4 for BAM/SAM, 5 for the
rest, 38 in total. No sampling, and no rows dropped except the threaded
FastCrossMap runs (see below).

Why small multiples rather than one chart
-----------------------------------------
Peak memory spans four orders of magnitude, from ~12 MB (BED) to ~5,150 MB
(CrossMap on BigWig). On one chart that range squeezes every light format into
an indistinguishable cluster at the left edge -- which hides the result rather
than showing it.

Each subplot therefore sets its own x limits, padded around that format's own
data, so a format whose tools sit between 11 and 48 MB is drawn at the same
size as one whose tools sit between 665 and 6,555 MB. This is deliberate:
the figure's job is to show the spread and the gap *within* a format, and eight
independent scales show more of each format than two shared ones would.

The cost is that a box further right in one panel is not larger than a box
further right in another panel, so the panels must be read one at a time. Every
panel prints its own MB ticks, and the value labels are the medians, so no
comparison has to be made by eye. Comparisons *across* formats belong in
Table S1, which carries every per-dataset number.

Why log: within a panel the boxes span roughly one order of magnitude, and on a
linear axis they would bunch at one end; a log axis also makes the *ratio*
readable -- equal factors are equal distances anywhere in any panel.

Reading the boxes
-----------------
Whiskers run to the data extremes (whis=(0, 100)), not to 1.5 x IQR: with four
or five points per box there are no outliers to separate and the reader is
better served by the actual min and max.

The width of a box is itself a result. For the light formats the boxes collapse
almost to a vertical line -- BED across its five datasets is 11.67-11.75 MB --
while for the track formats they are wide (Wiggle 665-1,266 MB). Same tool,
same code path; the difference is what the format makes the tool hold at once.

Only the single-threaded FastCrossMap is drawn. The benchmark JSONs also hold
`FastCrossMap (4t)` / `(8t)` rows for BAM and SAM; those are a threading result,
not a memory result, and including them would put three boxes under one tool
name. -t 1 is also the paper's uniform cross-tool baseline.

The "Nx less" annotation in each subplot is the median of the PER-DATASET
paired ratios (CrossMap_mem / FastCrossMap_mem on the same dataset), with the
observed range, not the ratio of two medians. Pairing first is the stricter
reading: it asks "on a given dataset, how much less?" and it cannot be moved by
one format's datasets being heavier than another's. CrossMap is the reference
because it is the tool this work compares against throughout.

Input:   paper/results/benchmark_{bed,vcf,gff,wig,bigwig,maf,bam,sam}.json
Output:  paper/figures/figS4_memory_by_format.pdf
         paper/figures/figS4_memory_by_format.png
         paper/results/memory_by_format.json

Usage: python paper/plot/figS4_memory_stability.py
"""

import json
import statistics as stats
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.patches import Patch
import numpy as np

RESULTS_DIR = Path("paper/results")
FIGURES_DIR = Path("paper/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# (format label, benchmark JSON) in reading order
FORMATS = [
    ("BED",     "benchmark_bed.json"),
    ("VCF",     "benchmark_vcf.json"),
    ("GFF",     "benchmark_gff.json"),
    ("MAF",     "benchmark_maf.json"),
    ("BAM",     "benchmark_bam.json"),
    ("SAM",     "benchmark_sam.json"),
    ("Wiggle",  "benchmark_wig.json"),
    ("BigWig",  "benchmark_bigwig.json"),
]

# per-subplot x padding, in log space: limits are (min/factor, max*factor).
# Each panel is scaled to its own data, so the boxes fill the panel rather than
# crowding one edge.
X_PAD = 1.8

TOOL_COLORS = {
    "FastCrossMap": "#1f77b4",
    "CrossMap": "#ff7f0e",
    "liftOver": "#2ca02c",
    "FastRemap": "#d62728",
}
# bottom-to-top inside each subplot; FastCrossMap first so it is always in the
# same place relative to the others
TOOL_ORDER = ["FastCrossMap", "CrossMap", "liftOver", "FastRemap"]

REFERENCE_TOOL = "CrossMap"   # denominator of the fold annotation

NCOL = 2


def load_series():
    """
    {format: {tool: {dataset: peak_mb}}} over that format's datasets.

    Threaded FastCrossMap rows are dropped: this figure is -t 1 only.
    """
    out = {}
    for label, fname in FORMATS:
        path = RESULTS_DIR / fname
        if not path.exists():
            raise SystemExit(f"missing input: {path}")

        with open(path) as fh:
            doc = json.load(fh)

        per_tool: dict[str, dict[str, float]] = {}
        for r in doc.get("results", []):
            tool = r.get("tool")
            peak = r.get("peak_memory_mb")
            if not tool or not peak or not r.get("success"):
                continue
            if tool.startswith("FastCrossMap ("):
                continue          # -t 4 / -t 8: wrong figure
            per_tool.setdefault(tool, {})[r["dataset_name"]] = peak

        out[label] = per_tool
    return out


def fold_stats(series) -> dict:
    """
    Per-dataset paired ratios against REFERENCE_TOOL, per format.

    {format: {"median", "min", "max", "n"}}. Datasets where either tool is
    missing are skipped, so every ratio is truly paired.
    """
    out = {}
    for label, _ in FORMATS:
        tool_map = series[label]
        ref = tool_map.get(REFERENCE_TOOL, {})
        fcm = tool_map.get("FastCrossMap", {})
        shared = sorted(set(ref) & set(fcm))
        ratios = [ref[ds] / fcm[ds] for ds in shared if fcm[ds] > 0]
        if ratios:
            out[label] = {
                "median": stats.median(ratios),
                "min": min(ratios),
                "max": max(ratios),
                "n": len(ratios),
            }
    return out


def check_claim(series) -> list[str]:
    """
    Verify the figure's claim: FastCrossMap is the least-memory tool in every
    format, compared on medians. Returns the formats where it is not -- the
    caller fails on a non-empty list rather than publishing a figure that
    asserts something the data does not show.
    """
    violations = []
    for label, _ in FORMATS:
        tool_map = series[label]
        if "FastCrossMap" not in tool_map:
            violations.append(f"{label}: no FastCrossMap row")
            continue
        fcm_med = stats.median(tool_map["FastCrossMap"].values())
        for tool, ds_map in tool_map.items():
            if tool == "FastCrossMap":
                continue
            med = stats.median(ds_map.values())
            if med < fcm_med:
                violations.append(
                    f"{label}: {tool} is lower ({med:.1f} MB vs {fcm_med:.1f} MB)"
                )
    return violations


def x_limits(tool_map) -> tuple[float, float]:
    """
    Per-panel log-x limits: this format's own range, padded by X_PAD on both
    sides so the whisker caps and the median labels stay inside the axes.
    """
    vals = [v for ds_map in tool_map.values() for v in ds_map.values()]
    lo, hi = min(vals), max(vals)
    return lo / X_PAD, hi * X_PAD


def draw_subplot(ax, label, tool_map, fold):
    """One format: horizontal box per tool, log x."""
    present = [t for t in TOOL_ORDER if t in tool_map]
    positions = list(range(len(present)))
    lo, hi = x_limits(tool_map)

    for pos, tool in zip(positions, present):
        vals = list(tool_map[tool].values())
        bp = ax.boxplot(
            [vals], positions=[pos], vert=False,
            widths=0.62, whis=(0, 100),
            patch_artist=True, manage_ticks=False, zorder=3,
            medianprops=dict(color="#111111", lw=1.5),
            whiskerprops=dict(color="#555555", lw=1.0),
            capprops=dict(color="#555555", lw=1.0),
            flierprops=dict(marker="o", ms=3, mec="#555555"),
        )
        for box in bp["boxes"]:
            box.set(facecolor=TOOL_COLORS[tool], alpha=0.85,
                    edgecolor="#333333", lw=0.8)

        # median printed at the box: with one scale per panel the eye cannot
        # carry a value across panels, so the number travels with the box
        med = stats.median(vals)
        ax.text(med, pos + 0.40, f"{med:,.1f}", ha="center", va="bottom",
                fontsize=7, color="#333333", zorder=4)

    ax.set_xscale("log")
    ax.set_xlim(lo, hi)
    # each panel spans roughly one decade, where matplotlib's default log ticks
    # print "2 x 10^1"-style labels that are unreadable at this size. Put the
    # ticks at 1/2/3/5 of each decade and print them as plain MB numbers.
    ax.xaxis.set_major_locator(
        ticker.LogLocator(base=10.0, subs=(1.0, 2.0, 3.0, 5.0), numticks=12))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.set_yticks(positions)
    ax.set_yticklabels(present, fontsize=8)
    ax.set_ylim(-0.7, len(present) - 0.3)
    ax.grid(alpha=0.25, axis="x", which="major", lw=0.5, zorder=0)
    ax.set_axisbelow(True)

    title = label
    if fold:
        title += f"  —  {fold['median']:,.1f}x less than {REFERENCE_TOOL}"
    ax.set_title(title, fontsize=9.5, loc="left", pad=4)


def main():
    series = load_series()

    violations = check_claim(series)
    if violations:
        raise SystemExit(
            "the figure's claim (FastCrossMap has the lowest peak memory in "
            "every format) does not hold:\n  " + "\n  ".join(violations)
            + "\n\nDo not draw the chart; the caption would be false."
        )
    print("claim holds: FastCrossMap has the lowest median peak memory in all "
          f"{len(FORMATS)} format groups")

    folds = fold_stats(series)

    nrow = (len(FORMATS) + NCOL - 1) // NCOL
    fig, axes = plt.subplots(nrow, NCOL, figsize=(10.5, 2.05 * nrow))
    axes = np.atleast_1d(axes).ravel()

    any_tool: set[str] = set()
    for ax, (label, _) in zip(axes, FORMATS):
        tool_map = series[label]
        draw_subplot(ax, label, tool_map, folds.get(label))
        any_tool.update(tool_map)

    for ax in axes[len(FORMATS):]:
        ax.axis("off")

    # one figure-level x label: every panel shows the same quantity and unit,
    # each on its own range (the numeric ticks carry the range per panel)
    fig.text(0.5, 0.045,
             "peak RSS (MB), log scale, one range per panel   —   box = median "
             "and IQR over the format's datasets; whiskers = full range",
             ha="center", fontsize=9)

    handles = [Patch(facecolor=TOOL_COLORS[t], label=t)
               for t in TOOL_ORDER if t in any_tool]
    fig.legend(handles=handles, fontsize=8.5, loc="lower center",
               bbox_to_anchor=(0.5, -0.055), ncol=len(handles),
               framealpha=0.9, title="peak memory, -t 1", title_fontsize=8)

    fig.suptitle("Figure S4. Peak memory by format: FastCrossMap against "
                 "the other tools", fontsize=11, x=0.02, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0.075, 1, 0.97))

    out_pdf = FIGURES_DIR / "figS4_memory_by_format.pdf"
    out_png = FIGURES_DIR / "figS4_memory_by_format.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"wrote {out_pdf}")
    print(f"wrote {out_png}")

    # machine-readable copy of exactly what was plotted
    payload = {
        "unit": "MB",
        "statistic": "box = median + IQR, whisker = full range over the "
                     "format's datasets",
        "baseline": "-t 1",
        "x_axis": "log, one range per panel, padded by x%.1f on each side "
                  "around that format's own min and max" % X_PAD,
        "reference_for_fold": REFERENCE_TOOL,
        "fold_definition": "median of per-dataset paired ratios "
                           "(reference_mb / fastcrossmap_mb)",
        "formats": {
            label: {
                "datasets_used": sorted({
                    ds for ds_map in series[label].values() for ds in ds_map
                }),
                "tools": {
                    tool: {
                        "per_dataset_mb": {k: round(v, 2)
                                           for k, v in sorted(ds_map.items())},
                        "median_mb": round(stats.median(ds_map.values()), 2),
                        "min_mb": round(min(ds_map.values()), 2),
                        "max_mb": round(max(ds_map.values()), 2),
                        "n_datasets": len(ds_map),
                    }
                    for tool, ds_map in series[label].items()
                },
                "fold_vs_reference": (
                    {k: round(v, 2) for k, v in folds[label].items()}
                    if label in folds else None
                ),
            }
            for label, _ in FORMATS
        },
    }
    path = RESULTS_DIR / "memory_by_format.json"
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)
    print(f"wrote {path}")

    print(f"\n{'format':<9} {'tool':<14} {'n':>3} {'median':>9} {'min':>9} {'max':>9}")
    for label, _ in FORMATS:
        for tool in TOOL_ORDER:
            ds_map = series[label].get(tool)
            if not ds_map:
                continue
            v = list(ds_map.values())
            print(f"{label:<9} {tool:<14} {len(v):>3} {stats.median(v):>9.1f} "
                  f"{min(v):>9.1f} {max(v):>9.1f}")
        f = folds.get(label)
        if f:
            print(f"{'':<9} {'-> vs CrossMap':<14} {f['n']:>3} "
                  f"{f['median']:>8.1f}x ({f['min']:.1f}-{f['max']:.1f})")


if __name__ == "__main__":
    main()
