#!/usr/bin/env python3
"""
figS2_realignment.py - Figure S2: BAM liftover vs. re-alignment

Answers the reviewer's question of why a user would lift a BAM over rather than
re-align the reads against the new reference. The same reads go down both paths:

  liftover     : fast-crossmap bam <chain> <bam_hg19> <bam_hg38>
  re-alignment : samtools bam2fq | bwa mem -t N <hg38.fa> | samtools sort | samtools markdup

Both are O(N) in the read count, so the two curves never cross: liftover is
faster and lighter at every measured scale, by a roughly constant factor. The
figure therefore quantifies the gap rather than locating a crossover point.

Two panels, side by side, both against read count:
  (a) runtime      --- log-log; the two curves sit ~1.9 decades apart vertically
  (b) peak memory  --- log; liftover is flat near 16 MB, re-alignment flat near
      5.6 GB, ~2.6 decades apart. Memory is on a log axis for the same reason
      runtime is: on a linear axis the liftover trace would be a flat line on
      the x-axis, indistinguishable from zero beside 5.6 GB, and the two panels
      would no longer share one visual grammar.

Each panel draws two ratio arrows per measured scale: a bold solid one across
the -t 1 boundaries (labelled "109x (t=1)", the paper's headline number) and a
light dotted one across the -t 8 boundaries ("34x (t=8)"), answering "and if the
re-alignment is allowed every core?". The near-constant factors read straight
off the figure: 83-109x at -t 1 on time, 31-34x at -t 8, 300-373x on memory.

Threading is shown as dashed extensions of the same two colours: -t 4 and -t 8
for both paths. They answer "does the gap survive if the re-alignment is
allowed to use threads?" -- it narrows from ~83-109x to ~31-34x but never
closes, and the ratio arrows stay on the -t 1 pair, since -t 1 is the paper's
uniform cross-tool baseline (CrossMap has no threading).

Input:   paper/results/benchmark_realignment.json
Output:  paper/figures/figS2_realignment.pdf
         paper/figures/figS2_realignment.png
         paper/results/realignment_fit.json

Usage: python paper/plot/figS2_realignment.py
"""

import json
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter

RESULTS_DIR = Path("paper/results")
FIGURES_DIR = Path("paper/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

COLORS = {
    "liftover": "#1f77b4",
    "realignment": "#ff7f0e",
}
LABELS = {
    "liftover": "FastCrossMap liftover",
    "realignment": "BWA-MEM re-alignment",
}
METHODS = ["liftover", "realignment"]

# the paper's cross-tool baseline is -t 1 (CrossMap has no threading). The
# threaded runs are drawn as a ribbon: the -t 1 boundary solid, the -t 8
# boundary dotted, and everything between them shaded, so the band itself reads
# as "what threading buys" instead of as four near-parallel curves. -t 4 is
# inside the band and is reported in the fit JSON, not drawn (it adds a marker
# without adding a story).
FIT_THREADS = [1, 4, 8]
BAND_LO, BAND_HI = 1, 8
BAND_ALPHA = 0.16

# axis frame shared by both panels; the right margin keeps the arrow ratio
# labels clear of the axes edge
X_MIN, X_MAX = 6.5e5, 2.6e7
XTICKS = [1e6, 5e6, 1e7]


def load_points(path: Path):
    """Return {threads: {method: [(reads, time_sec, peak_mb), ...]}}."""
    with open(path) as f:
        data = json.load(f)
    out: dict = {}
    for r in data.get("results", []):
        if not r.get("success"):
            continue
        n = r.get("num_reads")
        t = r.get("mean_time_sec")
        if not n or t is None or t <= 0:
            continue
        mem = r.get("peak_memory_mb")
        out.setdefault(r["threads"], {}).setdefault(r["method"], []).append(
            (float(n), float(t), float(mem) if mem is not None else float("nan")))
    for th in out:
        for m in out[th]:
            out[th][m].sort()
    return out


def fit_line(ns, ts):
    """Least squares t = a + k * n (free intercept: both paths carry a fixed
    start-up cost -- chain/index load for liftover, the bwa seed/index pipeline
    for re-alignment -- that bends the line at the small end, so a through-origin
    fit would misstate it). Returns (k seconds/read, a seconds)."""
    x = np.asarray(ns, dtype=float)
    y = np.asarray(ts, dtype=float)
    k, a = np.polyfit(x, y, 1)
    return float(k), float(a)


def render_axes(ax, ylim, ylabel, title):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Number of reads", fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=10.5, fontweight='bold')
    ax.grid(True, which='major', linestyle='--', alpha=0.3)
    ax.grid(True, which='minor', linestyle=':', alpha=0.15)
    ax.set_xticks(XTICKS)
    # plain "1M / 5M / 10M" reads better than 10^6 here: the three scales are
    # concrete dataset sizes, not abstract decades
    ax.set_xticklabels(["1M", "5M", "10M"], fontsize=8.5)
    ax.xaxis.set_minor_formatter(NullFormatter())


def plot_comparison(ax, points, value_of, fmt, ylim, ylabel, title, label_pos,
                    t1_place):
    """One panel: a ribbon per method spanning the -t 1 .. -t 8 range, plus the
    ratio arrows between the two -t 1 curves and the two -t 8 curves.

    The -t 1 boundary is drawn solid (the paper's cross-tool baseline) and the
    -t 8 boundary dotted; everything between them is shaded. Value labels go on
    the -t 1 boundary only -- labelling both edges of both ribbons would bury
    the figure. The legend is figure-level (see main).

    label_pos maps method -> (dx, dy, ha), the per-panel placement that keeps the
    -t 1 value labels off the ratio arrows and off the ribbon fills.
    t1_place is the per-scale (offset, ha) pair that steers each -t 1 ratio label
    clear of the neighbouring column's text; the -t 8 ratio label is stacked
    directly under it, so it inherits the same anchor and needs no placement of
    its own."""
    ns = [p[0] for p in points[BAND_LO]["liftover"]]
    vals = {}

    for method in METHODS:
        color = COLORS[method]
        lo = points[BAND_LO][method]
        hi = points[BAND_HI][method]
        xl = np.asarray([p[0] for p in lo], dtype=float)
        yl = np.asarray([value_of(p) for p in lo], dtype=float)
        xh = np.asarray([p[0] for p in hi], dtype=float)
        yh = np.asarray([value_of(p) for p in hi], dtype=float)
        vals[method] = (list(yl), list(yh))

        # the ribbon first, so both boundary lines sit on top of the fill
        ax.fill_between(xl, yl, yh, color=color, alpha=BAND_ALPHA,
                        linewidth=0, zorder=2)

        ax.plot(xl, yl, '-', color=color, linewidth=1.9, zorder=5,
                label=f"{LABELS[method]}  (t=1)")
        ax.plot(xl, yl, 'o', color=color, markersize=7,
                markeredgecolor='white', markeredgewidth=1.0, zorder=6)
        ax.plot(xh, yh, ':', color=color, linewidth=1.6, zorder=5,
                label=f"{LABELS[method]}  (t=8)")
        ax.plot(xh, yh, '^', color=color, markersize=6, markerfacecolor='white',
                markeredgewidth=1.2, zorder=6)

        dx, dy, ha = label_pos[method]
        for x, y in zip(xl, yl):
            ax.annotate(fmt(y), xy=(x, y), xytext=(dx, dy),
                        textcoords='offset points', ha=ha,
                        va='bottom' if dy > 0 else 'top',
                        color=color, fontsize=7.5, fontweight='bold', zorder=7)

    # Ratio arrows, one pair per measured scale: the solid arrow spans the two
    # -t 1 boundaries (the paper's headline number), and a light dotted arrow
    # spans the two -t 8 boundaries -- it answers "and if the re-alignment is
    # allowed every core?". The two labels are stacked at the -t 1 anchor, the
    # -t 8 value directly under the -t 1 value, so they read as one block instead
    # of the -t 8 text floating loose in the gap between the ribbons. It keeps
    # its lighter grey and smaller size, so the hierarchy survives the stacking.
    lo_t1, lo_t8 = vals["liftover"]      # (t=1 values, t=8 values)
    hi_t1, hi_t8 = vals["realignment"]
    for (x, a, b, c, e), (f1, ha1) in zip(
            zip(ns, lo_t1, hi_t1, lo_t8, hi_t8), t1_place):
        ax.annotate("", xy=(x, b), xytext=(x, a),
                    arrowprops=dict(arrowstyle='<->', color='#555555',
                                    lw=1.0, shrinkA=4, shrinkB=4), zorder=5)
        ax.text(x * f1, (a * b) ** 0.5, f"{b / a:,.0f}× (t=1)",
                color='#333333', fontsize=8.5, fontweight='bold',
                va='center', ha=ha1, zorder=7)

        ax.annotate("", xy=(x, e), xytext=(x, c),
                    arrowprops=dict(arrowstyle='<->', color='#999999',
                                    lw=0.9, linestyle=':', shrinkA=3, shrinkB=3),
                    zorder=5)
        ax.annotate(f"{e / c:,.0f}× (t=8)", xy=(x * f1, (a * b) ** 0.5),
                    xytext=(0, -11), textcoords='offset points',
                    ha=ha1, va='top', color='#777777', fontsize=7.5, zorder=7)

    render_axes(ax, ylim, ylabel, title)
    return vals


def main():
    print("=" * 64)
    print("Figure S2: BAM liftover vs. re-alignment")
    print("=" * 64)

    path = RESULTS_DIR / "benchmark_realignment.json"
    if not path.exists():
        print(f"  Missing {path} — nothing to plot")
        return

    points = load_points(path)
    print(f"\n[realignment] {path.name}")
    for th in sorted(points):
        for m in METHODS:
            if m in points[th]:
                print(f"  t={th:<2} {m:<12} " + "  ".join(
                    f"{n:,.0f}reads->{t:.1f}s/{mem:,.0f}MB"
                    for n, t, mem in points[th][m]))

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12.0, 5.1))

    # panel (a): runtime, seconds
    times = plot_comparison(
        ax_a, points, value_of=lambda p: p[1],
        fmt=lambda s: f"{s:.1f} s" if s < 100 else f"{s:,.0f} s",
        ylim=(4, 4e4), ylabel="Runtime (seconds)",
        title="(a) Runtime — log scale",
        # liftover sits under two threaded traces here, so its -t 1 labels go
        # up-left (clear of the arrow and of the dashed traces); the re-alignment
        # curve has open space above it. The stacked ratio block goes to the
        # right of each arrow, in the gap between the two ribbons.
        label_pos={"liftover": (-9, 8, 'right'), "realignment": (0, 9, 'center')},
        t1_place=[(1.08, 'left')] * 3)

    # panel (b): peak memory, MB
    mems = plot_comparison(
        ax_b, points, value_of=lambda p: p[2],
        fmt=lambda m: f"{m:,.0f} MB" if m < 1000 else f"{m / 1024:.1f} GB",
        ylim=(6, 3e4), ylabel="Peak memory, RSS (MB)",
        title="(b) Peak memory — log scale",
        # here liftover is at the floor with open space below it and the traces
        # above, so its labels go down; re-alignment's go up.
        label_pos={"liftover": (0, -13, 'center'), "realignment": (0, 9, 'center')},
        t1_place=[(1.08, 'left')] * 3)

    # one legend for the whole figure, below both panels: the solid/dotted pair
    # is read once per colour, and the ribbon fill is named in the panel text
    # above, so four entries cover two ribbons.
    handles, labels = ax_a.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, fontsize=8.5,
               frameon=False, bbox_to_anchor=(0.5, 0.005))

    # name the band once, on the panel where it is widest
    ax_a.annotate("shaded band = the whole 1–8-thread range",
                  xy=(0.98, 0.02), xycoords='axes fraction', ha='right',
                  va='bottom', fontsize=7.5, color='#555555', style='italic')

    # fits + ratios, for Supplementary Note S1 and the response text. The fitted
    # model is reported for every thread count; the ratio arrows, like the
    # paper's headline numbers, are the -t 1 pair.
    ns = [p[0] for p in points[1]["liftover"]]
    fits = {}
    for th in FIT_THREADS:
        if th not in points:
            continue
        fits[th] = {}
        for m in METHODS:
            if m not in points[th]:
                continue
            k, a = fit_line(ns, [p[1] for p in points[th][m]])
            fits[th][m] = {
                "n": len(ns),
                "slope_seconds_per_read": k,
                "microseconds_per_read": k * 1e6,
                "intercept_seconds": a,
            }
    time_ratio = [hi / lo for lo, hi in
                  zip(times["liftover"][0], times["realignment"][0])]
    mem_ratio = [hi / lo for lo, hi in
                 zip(mems["liftover"][0], mems["realignment"][0])]
    time_ratio8 = [hi / lo for lo, hi in
                   zip(times["liftover"][1], times["realignment"][1])]

    print("\n  -t 1 ratios")
    for n, tr, mr in zip(ns, time_ratio, mem_ratio):
        print(f"    {n:>12,.0f} reads   time {tr:6.1f}x   memory {mr:6.1f}x")
    print(f"    time   {min(time_ratio):.0f}-{max(time_ratio):.0f}x")
    print(f"    memory {min(mem_ratio):.0f}-{max(mem_ratio):.0f}x")
    print("  -t 8 time ratios  " +
          "  ".join(f"{n:,.0f}->{v:.1f}x" for n, v in zip(ns, time_ratio8)))
    for th in sorted(fits):
        for m in METHODS:
            if m in fits[th]:
                f = fits[th][m]
                print(f"  t={th} {m:<12} {f['microseconds_per_read']:8.2f} us/read  "
                      f"+ {f['intercept_seconds']:7.2f} s  (t = a + k*N)")
    print("\n  threaded time ratios")
    for th in FIT_THREADS:
        if th == 1 or th not in points:
            continue
        r = [hi / lo for lo, hi in
             zip([p[1] for p in points[th]["liftover"]],
                 [p[1] for p in points[th]["realignment"]])]
        print(f"    t={th}  " + "  ".join(f"{n:,.0f}->{v:.1f}x" for n, v in zip(ns, r)))

    fig.suptitle("Figure S2: BAM liftover vs. BWA-MEM re-alignment",
                 fontsize=13, fontweight='bold', y=1.00)
    # reserve the bottom tenth for the figure legend (tight_layout ignores a
    # figure-level legend, so the band has to be carved out explicitly)
    plt.tight_layout(rect=(0, 0.09, 1, 1))
    out_pdf = FIGURES_DIR / "figS2_realignment.pdf"
    out_png = FIGURES_DIR / "figS2_realignment.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  Saved {out_pdf}")

    fit_json = RESULTS_DIR / "realignment_fit.json"
    with open(fit_json, 'w') as f:
        json.dump({
            "generated": datetime.now().isoformat(timespec="seconds"),
            "method": "free-intercept least squares t = a + k*N per (method, threads); "
                      "time in seconds, N in reads",
            "source": str(path.as_posix()),
            "fits": fits,
            "baseline_threads": 1,
            "ratios": {
                "reads": ns,
                "time_ratio": time_ratio,
                "time_ratio_t8": time_ratio8,
                "memory_ratio": mem_ratio,
                "time_ratio_range": [min(time_ratio), max(time_ratio)],
                "memory_ratio_range": [min(mem_ratio), max(mem_ratio)],
            },
        }, f, indent=2)
    print(f"  Saved {fit_json}")
    print("\nDone!")


if __name__ == "__main__":
    main()
