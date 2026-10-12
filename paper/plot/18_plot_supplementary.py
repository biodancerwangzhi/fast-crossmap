#!/usr/bin/env python3
"""
18_plot_supplementary.py - Generate Supplementary Figures

Figure S1: Performance comparison for all formats (VCF/GFF/Wiggle/BigWig/MAF/SAM)
Figure S2: Speedup consistency across datasets (all formats box plot)

Usage: python paper/18_plot_supplementary.py
Output: paper/figures/fig_s1_all_formats.pdf
        paper/figures/fig_s2_speedup_consistency.pdf
"""

import json
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
    "liftOver": "#2ca02c",
    "FastRemap": "#d62728",
}

BENCHMARK_FILES = {
    "BED": RESULTS_DIR / "benchmark_bed.json",
    "BAM": RESULTS_DIR / "benchmark_bam.json",
    "VCF": RESULTS_DIR / "benchmark_vcf.json",
    "GFF": RESULTS_DIR / "benchmark_gff.json",
    "Wiggle": RESULTS_DIR / "benchmark_wig.json",
    "BigWig": RESULTS_DIR / "benchmark_bigwig.json",
    "MAF": RESULTS_DIR / "benchmark_maf.json",
    "SAM": RESULTS_DIR / "benchmark_sam.json",
}


def load_all_benchmarks() -> dict:
    data = {}
    for fmt, path in BENCHMARK_FILES.items():
        if path.exists():
            with open(path) as f:
                data[fmt] = json.load(f)
    return data


def get_speedups(data: dict, fmt: str) -> list[float]:
    results = data.get("results", [])
    datasets = set(r["dataset_name"] for r in results)
    speedups = []
    for ds in datasets:
        fc = next((r for r in results if r["dataset_name"] == ds
                   and r["tool"] == "FastCrossMap" and r.get("success")), None)
        cm = next((r for r in results if r["dataset_name"] == ds
                   and r["tool"] == "CrossMap" and r.get("success")), None)
        if fc and cm and fc["execution_time_sec"] > 0 and cm["execution_time_sec"] > 0:
            speedups.append(cm["execution_time_sec"] / fc["execution_time_sec"])
    return speedups


def plot_fig_s1(all_data: dict):
    """Performance comparison grouped bar chart for non-BED/BAM formats."""
    extra_formats = ["VCF", "GFF", "Wiggle", "BigWig", "MAF", "SAM"]
    available = [fmt for fmt in extra_formats if fmt in all_data]

    if not available:
        print("  No extra format data available for Figure S1")
        return

    n_formats = len(available)
    cols = min(3, n_formats)
    rows = (n_formats + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4.5 * rows), squeeze=False)

    for idx, fmt in enumerate(available):
        r, c = divmod(idx, cols)
        ax = axes[r][c]
        data = all_data[fmt]
        results = data.get("results", [])
        datasets = []
        seen = set()
        for res in results:
            dn = res["dataset_name"]
            if dn not in seen:
                datasets.append(dn)
                seen.add(dn)

        tools_in_fmt = []
        seen_t = set()
        for res in results:
            t = res["tool"]
            if t not in seen_t:
                tools_in_fmt.append(t)
                seen_t.add(t)

        x = np.arange(len(datasets))
        width = 0.8 / max(len(tools_in_fmt), 1)

        for ti, tool in enumerate(tools_in_fmt):
            times = []
            for ds in datasets:
                match = next((res for res in results if res["dataset_name"] == ds
                              and res["tool"] == tool and res.get("success")), None)
                times.append(match["execution_time_sec"] if match else 0)
            offset = (ti - len(tools_in_fmt) / 2 + 0.5) * width
            bars = ax.bar(x + offset, times, width * 0.9,
                          color=COLORS.get(tool, "#999"), label=tool, alpha=0.8)

        ax.set_title(f"{fmt}", fontsize=11, fontweight='bold')
        ax.set_ylabel("Time (s)", fontsize=9)
        short_names = [d[:18] + "…" if len(d) > 18 else d for d in datasets]
        ax.set_xticks(x)
        ax.set_xticklabels(short_names, fontsize=7, rotation=25, ha='right')
        ax.legend(fontsize=7, loc='upper right')

    # Hide unused subplots
    for idx in range(len(available), rows * cols):
        r, c = divmod(idx, cols)
        axes[r][c].set_visible(False)

    fig.suptitle("Figure S1: Performance Comparison Across All Formats",
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    out_pdf = FIGURES_DIR / "fig_s1_all_formats.pdf"
    out_png = FIGURES_DIR / "fig_s1_all_formats.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out_pdf}")


def plot_fig_s2(all_data: dict):
    """Speedup consistency box plot across all formats."""
    format_order = ["BED", "BAM", "VCF", "GFF", "Wiggle", "BigWig", "MAF", "SAM"]
    available = [fmt for fmt in format_order if fmt in all_data]

    if not available:
        print("  No data for Figure S2")
        return

    speedup_data = []
    labels = []
    for fmt in available:
        sp = get_speedups(all_data[fmt], fmt)
        if sp:
            speedup_data.append(sp)
            labels.append(f"{fmt}\n(n={len(sp)})")

    if not speedup_data:
        print("  No speedup data for Figure S2")
        return

    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 1.2), 5))
    bp = ax.boxplot(speedup_data, patch_artist=True, widths=0.5)

    colors_cycle = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(colors_cycle[i % len(colors_cycle)])
        patch.set_alpha(0.7)

    for i, sp in enumerate(speedup_data):
        for val in sp:
            ax.plot(i + 1, val, 'ko', markersize=5, alpha=0.5)

    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5, label='No speedup')
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Speedup (FastCrossMap vs CrossMap)", fontsize=11)
    ax.set_title("Figure S2: Speedup Consistency Across Formats and Datasets",
                 fontsize=12, fontweight='bold')

    medians = [np.median(sp) for sp in speedup_data]
    for i, med in enumerate(medians):
        ax.text(i + 1, med + 0.3, f"{med:.1f}x", ha='center', fontsize=8, fontweight='bold')

    plt.tight_layout()
    out_pdf = FIGURES_DIR / "fig_s2_speedup_consistency.pdf"
    out_png = FIGURES_DIR / "fig_s2_speedup_consistency.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {out_pdf}")


def main():
    print("=" * 60)
    print("Generating Supplementary Figures")
    print("=" * 60)

    all_data = load_all_benchmarks()
    if not all_data:
        print("Error: No benchmark results found. Run benchmarks first.")
        return

    print(f"\nLoaded data for: {', '.join(all_data.keys())}")

    print("\nFigure S1: All-format performance comparison")
    plot_fig_s1(all_data)

    print("\nFigure S2: Speedup consistency")
    plot_fig_s2(all_data)

    print("\nDone!")


if __name__ == "__main__":
    main()
