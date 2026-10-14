#!/usr/bin/env python3
"""
figS1_scaling.py - Figure S1: runtime scaling with dataset size

Both tools are O(N) linear scans over the input, so runtime is fitted through
the origin (t = k * N). A free-intercept OLS is also computed and printed as a
diagnostic, but the through-origin fit is what the figure draws: it is the model
that matches the mechanism, and it avoids a small-scale timing wobble inflating
the intercept and faking curvature on a log-log axis.

Panels:
  (a) VCF  - 1000 Genomes Phase 3, five chromosomes (1.10M - 7.08M variants)
  (b) MAF  - synthetic scale series (50K - 1M mutations)

Per-dataset time is the median of the repeats (the small datasets carry an
occasional I/O stall; the median is steadier than the mean), matching the
caliber used by the main tables.

The y-axis uses the unit that keeps the numbers legible for that panel — hours
for VCF (whose extrapolation to 600M variants reaches ~5-18 h), seconds for MAF (which tops out in
the tens of seconds) — and each extrapolated line is labelled with its value.

Input:   paper/results/benchmark_vcf.json
         paper/results/benchmark_maf.json
Output:  paper/figures/figS1_scaling.pdf
         paper/figures/figS1_scaling.png
         paper/results/scaling_fit.json   (coefficients, for Supplementary Note S1)

Usage: python paper/plot/figS1_scaling.py
"""

import json
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

RESULTS_DIR = Path("paper/results")
FIGURES_DIR = Path("paper/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

COLORS = {
    "FastCrossMap": "#1f77b4",
    "CrossMap": "#ff7f0e",
}
TOOLS = ["FastCrossMap", "CrossMap"]

PANELS = [
    {
        "key": "VCF",
        "path": RESULTS_DIR / "benchmark_vcf.json",
        "unit": "input_variants",
        "xlabel": "Number of variants",
        "title": "(a) VCF — GWAS-scale meta-analysis (1000 Genomes Phase 3)",
        "x_min": 8e5,
        "x_max": 1.5e9,                  # room right of the 600M line end for its value label
        "extrapolate_to": 6e8,
        "time_divisor": 3600.0,          # seconds -> hours
        "ylabel": "Runtime (hours)",
        "extrap_label": "{:.1f} h",
        # reference marks carry only the scale: the previous "(N cohorts)"
        # annotations had no citable source and contradicted the cohort count
        # used in 20_large_scale_estimation.py, so they are dropped rather than
        # given a number the manuscript cannot support.
        "marks": [(1e6, "1M"),
                  (1e7, "10M"),
                  (6e8, "600M")],
    },
    {
        "key": "MAF",
        "path": RESULTS_DIR / "benchmark_maf.json",
        "unit": "input_mutations",
        "xlabel": "Number of mutations",
        "title": "(b) MAF — mutation annotation (TCGA-scale)",
        "x_min": 3e4,
        "x_max": 2.5e7,                  # room right of the 10M line end for its value label
        "extrapolate_to": 1e7,
        "time_divisor": 1.0,             # seconds stay seconds
        "ylabel": "Runtime (seconds)",
        "extrap_label": "{:.0f} s",
        "marks": [(5e5, "TCGA pancancer\n(~0.5M)"),
                  (1e7, "10M")],
    },
]


def load_points(path: Path, unit: str):
    """Return {tool: [(n, median_time_sec), ...]} sorted by n."""
    with open(path) as f:
        data = json.load(f)
    out = {t: [] for t in TOOLS}
    for r in data.get("results", []):
        if not r.get("success"):
            continue
        n = r.get(unit)
        if not n:
            continue
        t = r.get("median_time_sec") or r.get("execution_time_sec")
        if t is None or t <= 0:
            continue
        out.setdefault(r["tool"], []).append((float(n), float(t)))
    for k in out:
        out[k].sort()
    return {k: v for k, v in out.items() if v}


def fit_through_origin(ns, ts):
    """t = k * n, least squares. Returns (k, r2)."""
    k = sum(n * t for n, t in zip(ns, ts)) / sum(n * n for n in ns)
    ss_res = sum((t - k * n) ** 2 for n, t in zip(ns, ts))
    mean_t = sum(ts) / len(ts)
    ss_tot = sum((t - mean_t) ** 2 for t in ts)
    return k, 1 - ss_res / ss_tot


def fit_free_intercept(ns, ts):
    """t = a + k * n, least squares. Returns (k, a, r2). Diagnostic only."""
    x = np.asarray(ns, dtype=float)
    y = np.asarray(ts, dtype=float)
    k, a = np.polyfit(x, y, 1)
    resid = y - (a + k * x)
    ss_res = float((resid ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return float(k), float(a), 1 - ss_res / ss_tot


def plot_panel(ax, panel, points):
    fits = {}
    divisor = panel["time_divisor"]
    top = 0.0
    bottom = None
    for tool in TOOLS:
        pts = points.get(tool)
        if not pts:
            continue
        ns = [p[0] for p in pts]
        ts = [p[1] for p in pts]
        ys = [t / divisor for t in ts]
        k, r2 = fit_through_origin(ns, ts)
        k_ols, a_ols, r2_ols = fit_free_intercept(ns, ts)
        fits[tool] = {
            "n": len(pts),
            "span": max(ns) / min(ns),
            "seconds_per_unit": k,
            "microseconds_per_unit": k * 1e6,
            "r2_through_origin": r2,
            "ols_slope_seconds_per_unit": k_ols,
            "ols_intercept_seconds": a_ols,
            "r2_free_intercept": r2_ols,
            "extrapolated_seconds": k * panel["extrapolate_to"],
        }
        color = COLORS[tool]
        kd = k / divisor  # slope in the panel's display unit

        # measured points
        ax.plot(ns, ys, 'o', color=color, markersize=6,
                markeredgecolor='white', markeredgewidth=0.8, zorder=5)

        # fitted line: solid over the measured range, dashed where extrapolated.
        # The dashed segment stops at extrapolate_to, leaving an empty right
        # margin for the value label.
        n_lo, n_hi = min(ns), max(ns)
        n_end = panel["extrapolate_to"]
        ax.plot([n_lo, n_hi], [kd * n_lo, kd * n_hi], '-', color=color, linewidth=1.8,
                label=f"{tool}  ({k * 1e6:.1f} µs/unit, R²={r2:.4f})", zorder=4)
        ax.plot([n_hi, n_end], [kd * n_hi, kd * n_end], '--',
                color=color, linewidth=1.5, alpha=0.75, zorder=3)

        # label the extrapolated endpoint with its value in the display unit,
        # in the empty margin to the right of where the line stops
        y_end = kd * n_end
        ax.annotate(panel["extrap_label"].format(y_end),
                    xy=(n_end, y_end), xytext=(5, 0),
                    textcoords='offset points', ha='left', va='center',
                    color=color, fontsize=7.5, fontweight='bold', zorder=6)

        top = max(top, y_end)
        bottom = min(ys) if bottom is None else min(bottom, min(ys))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(panel["x_min"], panel["x_max"])
    ax.set_ylim(bottom * 0.4, top * 2.2)

    # shade the extrapolated region (limits are final by now)
    measured_max = max(max(p[0] for p in v) for v in points.values())
    ax.axvspan(measured_max, panel["extrapolate_to"], color='gray', alpha=0.08, zorder=1)
    ax.text(measured_max * 1.15, 0.80, "extrapolated", transform=ax.get_xaxis_transform(),
            fontsize=7, color='gray', va='bottom', style='italic')

    # scenario markers along the top band. The leftmost and rightmost labels are
    # anchored inward so they cannot run off the axes (1M sits ~3% from the left
    # edge on the log axis, so a centred label would lose its leading characters).
    first_mark = min(n for n, _ in panel["marks"]) if panel["marks"] else None
    last_mark = max(n for n, _ in panel["marks"]) if panel["marks"] else None
    for n, label in panel["marks"]:
        ax.axvline(n, color='gray', linestyle=':', linewidth=0.8, alpha=0.5, zorder=2)
        ha = 'left' if n == first_mark else 'right' if n == last_mark else 'center'
        ax.text(n, 0.97, label, transform=ax.get_xaxis_transform(),
                fontsize=6.5, color='dimgray', ha=ha, va='top')

    ax.set_xlabel(panel["xlabel"], fontsize=10)
    ax.set_ylabel(panel["ylabel"], fontsize=10)
    ax.set_title(panel["title"], fontsize=10.5, fontweight='bold')
    ax.grid(True, which='major', linestyle='--', alpha=0.3)
    ax.grid(True, which='minor', linestyle=':', alpha=0.15)
    ax.legend(fontsize=7.5, loc='lower right')

    ticks = [t for t in [1e4, 1e5, 1e6, 1e7, 1e8, 1e9] if panel["x_min"] <= t <= panel["x_max"]]
    ax.set_xticks(ticks)
    # 10^n (mathtext), matching the y-axis formatter, so both axes read alike
    ax.set_xticklabels([f"$10^{{{int(round(np.log10(t)))}}}$" for t in ticks], fontsize=8)

    return fits


def main():
    print("=" * 64)
    print("Figure S1: runtime scaling with dataset size")
    print("=" * 64)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5))
    all_fits = {}
    sources = {}

    for ax, panel in zip(axes, PANELS):
        if not panel["path"].exists():
            print(f"  Missing {panel['path']} — skipping panel {panel['key']}")
            ax.set_visible(False)
            continue
        points = load_points(panel["path"], panel["unit"])
        if not points:
            print(f"  No usable points in {panel['path']} — skipping {panel['key']}")
            ax.set_visible(False)
            continue

        print(f"\n[{panel['key']}] {panel['path'].name}")
        for tool in TOOLS:
            pts = points.get(tool)
            if pts:
                print(f"  {tool:<13} " +
                      "  ".join(f"{n:,.0f}->{t:.1f}s" for n, t in pts))

        fits = plot_panel(ax, panel, points)
        for tool, f in fits.items():
            print(f"  {tool:<13} {f['microseconds_per_unit']:6.2f} µs/unit  "
                  f"R²(origin)={f['r2_through_origin']:.5f}  "
                  f"OLS: {f['ols_slope_seconds_per_unit'] * 1e6:.2f} µs/unit "
                  f"+ {f['ols_intercept_seconds']:.2f}s  R²={f['r2_free_intercept']:.5f}")
        if "FastCrossMap" in fits and "CrossMap" in fits:
            ratio = (fits["CrossMap"]["seconds_per_unit"] /
                     fits["FastCrossMap"]["seconds_per_unit"])
            at = panel["extrapolate_to"]
            print(f"  → {at:,.0f} units: FCM {fits['FastCrossMap']['extrapolated_seconds']:.0f}s "
                  f"vs CM {fits['CrossMap']['extrapolated_seconds']:.0f}s  ({ratio:.1f}×)")

        all_fits[panel["key"]] = {
            "unit": panel["unit"],
            "extrapolate_to": panel["extrapolate_to"],
            "y_unit": panel["ylabel"],
            "time_divisor_seconds": panel["time_divisor"],
            "models": fits,
        }
        sources[panel["key"]] = str(panel["path"].as_posix())

    fig.suptitle("Figure S1: Runtime scales linearly with dataset size",
                 fontsize=13, fontweight='bold', y=1.00)
    plt.tight_layout()
    out_pdf = FIGURES_DIR / "figS1_scaling.pdf"
    out_png = FIGURES_DIR / "figS1_scaling.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"\n  Saved {out_pdf}")

    fit_json = RESULTS_DIR / "scaling_fit.json"
    with open(fit_json, 'w') as f:
        json.dump({
            "generated": datetime.now().isoformat(timespec="seconds"),
            "method": "least squares through the origin (t = k * N); "
                      "median of repeats per dataset; free-intercept OLS reported as a diagnostic",
            "sources": sources,
            "panels": all_fits,
        }, f, indent=2)
    print(f"  Saved {fit_json}")
    print("\nDone!")


if __name__ == "__main__":
    main()
