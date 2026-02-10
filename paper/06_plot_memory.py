#!/usr/bin/env python3
"""
06_plot_memory.py - Generate memory efficiency comparison figure (Figure 2)

Figure 2 layout (1x2):
  Left: Memory usage curves - showing memory changes over time for three tools
  Right: Peak memory comparison - bar chart showing peak memory comparison

Design rationale:
- Left plot shows dynamic memory usage changes, highlighting FastCrossMap's streaming advantage
- Right plot shows direct peak memory comparison for quantitative comparison
- Uses color scheme consistent with Figure 1

Usage: python paper/06_plot_memory.py
Output: paper/figures/fig2_memory.pdf
"""

import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# 配置
# =============================================================================
RESULTS_DIR = Path("paper/results")
FIGURES_DIR = Path("paper/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Tool colors (consistent with Figure 1)
COLORS = {
    "FastCrossMap": "#1f77b4",  # Blue
    "CrossMap": "#ff7f0e",       # Orange
    "FastRemap": "#d62728"       # Red
}

# Tool order
TOOL_ORDER = ["FastCrossMap", "CrossMap", "FastRemap"]


def load_memory_data():
    """Load memory profiling data"""
    memory_file = RESULTS_DIR / "memory_profile.json"
    if not memory_file.exists():
        return None
    
    with open(memory_file) as f:
        return json.load(f)


def plot_memory_curves(data, ax):
    """
    绘制内存使用曲线
    
    参数:
        data: 内存采样数据
        ax: matplotlib axes
    """
    if not data:
        ax.text(0.5, 0.5, 'No memory data', ha='center', va='center', 
                transform=ax.transAxes, fontsize=12)
        ax.set_title('Memory Usage Over Time', fontsize=11, fontweight='bold')
        return
    
    results = {r["tool"]: r for r in data["results"]}
    
    # Plot each tool's memory curve
    for tool in TOOL_ORDER:
        if tool in results and results[tool]["success"]:
            r = results[tool]
            sample_times = r["sample_times"]
            memory_samples = r["memory_samples"]
            
            if sample_times and memory_samples:
                ax.plot(sample_times, memory_samples, 
                       label=tool, 
                       color=COLORS[tool],
                       linewidth=2,
                       alpha=0.8)
    
    ax.set_xlabel('Time (seconds)', fontsize=10)
    ax.set_ylabel('Memory Usage (MB)', fontsize=10)
    ax.set_title('Memory Usage Over Time', fontsize=11, fontweight='bold')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3, linestyle='--')


def plot_peak_memory_comparison(data, ax):
    """
    绘制峰值内存对比条形图
    
    参数:
        data: 内存采样数据
        ax: matplotlib axes
    """
    if not data:
        ax.text(0.5, 0.5, 'No memory data', ha='center', va='center', 
                transform=ax.transAxes, fontsize=12)
        ax.set_title('Peak Memory Comparison', fontsize=11, fontweight='bold')
        return
    
    results = {r["tool"]: r for r in data["results"]}
    
    # Prepare data
    tools = []
    peak_memories = []
    
    for tool in TOOL_ORDER:
        if tool in results and results[tool]["success"]:
            tools.append(tool)
            peak_memories.append(results[tool]["peak_memory_mb"])
    
    if not tools:
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', 
                transform=ax.transAxes)
        return
    
    colors = [COLORS[t] for t in tools]
    
    # Plot bar chart
    bars = ax.bar(range(len(tools)), peak_memories, color=colors, alpha=0.7, edgecolor='black')
    
    # Add value labels above bars
    for i, (bar, mem) in enumerate(zip(bars, peak_memories)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{mem:.1f} MB',
                ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax.set_ylabel('Peak Memory (MB)', fontsize=10)
    ax.set_title('Peak Memory Comparison', fontsize=11, fontweight='bold')
    ax.set_xticks(range(len(tools)))
    
    # Build x-axis labels, highlight FastCrossMap in red
    ax.set_xticklabels([])  # Clear first
    
    for i, tool in enumerate(tools):
        if tool == "FastCrossMap":
            ax.text(i, -0.12, tool, ha='center', va='top', 
                    transform=ax.get_xaxis_transform(), fontsize=9, 
                    color='red', fontweight='bold')
        else:
            ax.text(i, -0.12, tool, ha='center', va='top', 
                    transform=ax.get_xaxis_transform(), fontsize=9, 
                    color='black')
    
    # Add grid lines
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')


def main():
    print("=" * 60)
    print("Generating Memory Efficiency Comparison Figure (Figure 2)")
    print("=" * 60)
    
    # Load data
    memory_data = load_memory_data()
    
    if not memory_data:
        print("Error: No memory profiling data found")
        print("Please run first: python paper/05_memory_profile.py")
        return
    
    print(f"Input file: {memory_data['input_file']}")
    print(f"File size: {memory_data['input_size_mb']:.2f} MB")
    print(f"Sampling interval: {memory_data['sample_interval_sec']} seconds")
    
    # Create 1x2 figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Figure 2: Memory Efficiency Comparison', 
                 fontsize=14, fontweight='bold', y=1.00)
    
    # Left: Memory usage curves
    plot_memory_curves(memory_data, axes[0])
    
    # Right: Peak memory comparison
    plot_peak_memory_comparison(memory_data, axes[1])
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # Save combined figure
    output_pdf = FIGURES_DIR / "fig2_memory.pdf"
    output_png = FIGURES_DIR / "fig2_memory.png"
    
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(output_png, dpi=300, bbox_inches='tight')
    
    print(f"\nCombined figure saved to:")
    print(f"  {output_pdf}")
    print(f"  {output_png}")
    
    # Save individual subplots
    print(f"\nSaving individual subplots...")
    
    # Left: Memory usage curves
    fig_left, ax_left = plt.subplots(figsize=(7, 5))
    plot_memory_curves(memory_data, ax_left)
    plt.tight_layout()
    fig_left.savefig(FIGURES_DIR / "fig2a_memory_curves.pdf", dpi=300, bbox_inches='tight')
    fig_left.savefig(FIGURES_DIR / "fig2a_memory_curves.png", dpi=300, bbox_inches='tight')
    plt.close(fig_left)
    print(f"  {FIGURES_DIR / 'fig2a_memory_curves.pdf'}")
    
    # Right: Peak memory comparison
    fig_right, ax_right = plt.subplots(figsize=(6, 5))
    plot_peak_memory_comparison(memory_data, ax_right)
    plt.tight_layout()
    fig_right.savefig(FIGURES_DIR / "fig2b_peak_memory.pdf", dpi=300, bbox_inches='tight')
    fig_right.savefig(FIGURES_DIR / "fig2b_peak_memory.png", dpi=300, bbox_inches='tight')
    plt.close(fig_right)
    print(f"  {FIGURES_DIR / 'fig2b_peak_memory.pdf'}")
    
    # Print memory summary
    print("\n" + "=" * 60)
    print("Memory Efficiency Summary")
    print("=" * 60)
    
    results = {r["tool"]: r for r in memory_data["results"]}
    
    for tool in TOOL_ORDER:
        if tool in results and results[tool]["success"]:
            r = results[tool]
            print(f"{tool}:")
            print(f"  Execution time: {r['execution_time_sec']:.2f}s")
            print(f"  Peak memory: {r['peak_memory_mb']:.2f} MB")
            print(f"  Samples: {len(r['memory_samples'])}")
    
    # Calculate memory savings
    if "FastCrossMap" in results and "CrossMap" in results:
        fc = results["FastCrossMap"]
        cm = results["CrossMap"]
        if fc["success"] and cm["success"]:
            mem_ratio = cm["peak_memory_mb"] / fc["peak_memory_mb"]
            mem_saved = cm["peak_memory_mb"] - fc["peak_memory_mb"]
            print(f"\nMemory efficiency comparison:")
            print(f"  FastCrossMap vs CrossMap:")
            print(f"    Memory savings: {mem_saved:.2f} MB ({(1 - 1/mem_ratio)*100:.1f}%)")
            print(f"    CrossMap uses {mem_ratio:.2f}x memory")
    
    print("\n" + "=" * 60)
    print("Figure 2 Design Notes:")
    print("=" * 60)
    print("Left: Memory usage curves - showing memory changes over time for three tools")
    print("      Highlights FastCrossMap's streaming advantage (stable memory)")
    print("Right: Peak memory comparison - bar chart for direct quantitative comparison")
    print("      Facilitates quantitative memory efficiency comparison")
    print("\nNext step: python paper/07_accuracy_analysis.py")


if __name__ == "__main__":
    main()
