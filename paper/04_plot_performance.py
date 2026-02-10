#!/usr/bin/env python3
"""
04_plot_performance.py - Generate performance comparison figure (Figure 1)

Figure 1 layout (2x2):
  (a) BED single-thread performance comparison - fair comparison across all tools
  (b) BED multi-thread scalability - FastCrossMap 1/2/4/8 threads
  (c) BAM single-thread performance comparison - fair comparison across all tools
  (d) BAM multi-thread scalability - FastCrossMap 1/2/4/8 threads

Design rationale:
- Left column (a,c): Single-thread fair comparison, showing algorithm efficiency
- Right column (b,d): Multi-thread scalability, showing FastCrossMap's parallel capability
- Only show execution time, avoid redundancy with throughput

Usage: python paper/04_plot_performance.py
Output: paper/figures/fig1_performance.pdf
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

# Tool colors
COLORS = {
    "FastCrossMap": "#1f77b4",  # Blue
    "CrossMap": "#ff7f0e",       # Orange
    "liftOver": "#2ca02c",       # Green
    "FastRemap": "#d62728"       # Red
}

# Multi-thread color gradient (blue series)
THREAD_COLORS = {
    1: "#1f77b4",   # Dark blue
    2: "#4a9fd4",   # Medium blue
    4: "#7ec8e3",   # Light blue
    8: "#aee0f0",   # Lighter blue
    16: "#d4eef7"   # Lightest blue
}

# Tool order
TOOL_ORDER = ["FastCrossMap", "CrossMap", "liftOver", "FastRemap"]
BAM_TOOL_ORDER = ["FastCrossMap", "CrossMap", "FastRemap"]


def load_benchmark_data():
    """Load benchmark data"""
    bed_data = None
    bam_data = None
    bed_mt_data = None  # Multi-thread data
    bam_mt_data = None
    
    bed_file = RESULTS_DIR / "benchmark_bed.json"
    if bed_file.exists():
        with open(bed_file) as f:
            bed_data = json.load(f)
    
    bam_file = RESULTS_DIR / "benchmark_bam.json"
    if bam_file.exists():
        with open(bam_file) as f:
            bam_data = json.load(f)
    
    # Multi-thread data (if available)
    bed_mt_file = RESULTS_DIR / "benchmark_bed_multithread.json"
    if bed_mt_file.exists():
        with open(bed_mt_file) as f:
            bed_mt_data = json.load(f)
    
    bam_mt_file = RESULTS_DIR / "benchmark_bam_multithread.json"
    if bam_mt_file.exists():
        with open(bam_mt_file) as f:
            bam_mt_data = json.load(f)
    
    return bed_data, bam_data, bed_mt_data, bam_mt_data


def plot_single_thread_comparison(data, ax, format_type="BED", tool_order=None):
    """
    绘制单线程性能对比 (箱线图或条形图)
    
    参数:
        data: 基准测试数据
        ax: matplotlib axes
        format_type: "BED" 或 "BAM"
        tool_order: 工具顺序列表
    """
    if not data:
        ax.text(0.5, 0.5, f'No {format_type} data', ha='center', va='center', 
                transform=ax.transAxes, fontsize=12)
        ax.set_title(f'{format_type} Single-Thread', 
                     fontsize=11, fontweight='bold')
        return
    
    if tool_order is None:
        tool_order = TOOL_ORDER
    
    results = {r["tool"]: r for r in data["results"]}
    
    # Prepare data
    tools = []
    all_times_data = []
    
    for tool in tool_order:
        if tool in results and results[tool]["success"]:
            tools.append(tool)
            # Get all run times (for box plot)
            if "all_times" in results[tool] and results[tool]["all_times"]:
                all_times_data.append(results[tool]["all_times"])
            else:
                all_times_data.append([results[tool]["execution_time_sec"]])
    
    if not tools:
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=ax.transAxes)
        return
    
    colors = [COLORS[t] for t in tools]
    
    # Use box plot to show distribution of multiple runs
    bp = ax.boxplot(all_times_data, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # No longer adding median labels (to avoid obscuring the plot)
    
    ax.set_ylabel('Execution Time (seconds)', fontsize=10)
    ax.set_title(f'{format_type} Single-Thread Comparison', 
                 fontsize=11, fontweight='bold')
    
    # Build x-axis labels, highlight FastCrossMap in red, FastRemap's 4T* in red
    labels = []
    label_colors = []
    for tool in tools:
        if tool == "FastCrossMap":
            labels.append(f"{tool}\n(1T)")
            label_colors.append('red')
        elif tool == "FastRemap":
            labels.append(f"{tool}\n(4T*)")
            label_colors.append('black')  # 工具名黑色，但 4T* 需要特殊处理
        else:
            labels.append(f"{tool}\n(1T)")
            label_colors.append('black')
    
    # Set x-axis labels
    ax.set_xticks(range(1, len(tools) + 1))
    ax.set_xticklabels([])  # Clear first, then manually add colored labels
    
    # Manually add colored x-axis labels
    for i, (tool, label) in enumerate(zip(tools, labels)):
        if tool == "FastCrossMap":
            # FastCrossMap entirely in red
            ax.text(i + 1, -0.12, label, ha='center', va='top', 
                    transform=ax.get_xaxis_transform(), fontsize=9, 
                    color='red', fontweight='bold')
        elif tool == "FastRemap":
            # FastRemap: tool name in black, 4T* in red
            ax.text(i + 1, -0.06, tool, ha='center', va='top', 
                    transform=ax.get_xaxis_transform(), fontsize=9, color='black')
            ax.text(i + 1, -0.14, "(4T*)", ha='center', va='top', 
                    transform=ax.get_xaxis_transform(), fontsize=9, 
                    color='red', fontweight='bold')
        else:
            ax.text(i + 1, -0.12, label, ha='center', va='top', 
                    transform=ax.get_xaxis_transform(), fontsize=9, color='black')


def plot_multithread_scaling(mt_data, ax, format_type="BED"):
    """
    绘制多线程扩展性图 (使用箱线图)
    
    参数:
        mt_data: 多线程基准测试数据
        ax: matplotlib axes
        format_type: "BED" 或 "BAM"
    """
    panel_idx = 1 if format_type == "BED" else 3
    
    if not mt_data:
        # 如果没有多线程数据，显示占位信息
        ax.text(0.5, 0.5, f'Run benchmark with\n-t 1,2,4,8,16 to generate\nmulti-thread data', 
                ha='center', va='center', transform=ax.transAxes, fontsize=10)
        ax.set_title(f'FastCrossMap {format_type} Scaling', 
                     fontsize=11, fontweight='bold')
        return
    
    # Parse multi-thread data
    thread_data = []  # [(threads, all_times), ...]
    
    for result in mt_data.get("results", []):
        if result.get("success", False):
            threads = result["threads"]
            # Get all run times (for box plot)
            if "all_times" in result and result["all_times"]:
                all_times = result["all_times"]
            else:
                all_times = [result["execution_time_sec"]]
            thread_data.append((threads, all_times))
    
    if not thread_data:
        ax.text(0.5, 0.5, 'No valid multi-thread data', ha='center', va='center', 
                transform=ax.transAxes)
        return
    
    # Sort by thread count
    thread_data.sort(key=lambda x: x[0])
    threads = [d[0] for d in thread_data]
    all_times_data = [d[1] for d in thread_data]
    
    # Use box plot
    colors = [THREAD_COLORS.get(t, "#1f77b4") for t in threads]
    bp = ax.boxplot(all_times_data, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Calculate speedup and display below X-axis
    baseline_median = np.median(all_times_data[0])  # 1 线程作为基准
    
    # Set x-axis labels
    ax.set_xticks(range(1, len(threads) + 1))
    ax.set_xticklabels([])  # Clear first
    
    # Manually add x-axis labels and speedup
    for i, (th, times) in enumerate(zip(threads, all_times_data)):
        median = np.median(times)
        speedup = baseline_median / median
        
        # Thread count label
        ax.text(i + 1, -0.06, f'{th}T', ha='center', va='top', 
                transform=ax.get_xaxis_transform(), fontsize=9, color='black')
        
        # Speedup label (below thread count) - all in black
        if th == 1:
            speedup_text = "(baseline)"
        else:
            speedup_text = f"({speedup:.2f}x)"
        
        ax.text(i + 1, -0.14, speedup_text, ha='center', va='top', 
                transform=ax.get_xaxis_transform(), fontsize=8, 
                color='black', fontweight='bold')
    
    ax.set_ylabel('Execution Time (seconds)', fontsize=10)
    ax.set_title(f'FastCrossMap {format_type} Scaling', 
                 fontsize=11, fontweight='bold')


def create_mock_multithread_data(single_data, format_type="BED"):
    """
    如果没有多线程数据，基于单线程数据创建模拟数据用于演示
    实际使用时应该运行真实的多线程基准测试
    """
    if not single_data:
        return None
    
    results = {r["tool"]: r for r in single_data["results"]}
    if "FastCrossMap" not in results or not results["FastCrossMap"]["success"]:
        return None
    
    base_time = results["FastCrossMap"]["execution_time_sec"]
    
    # 模拟数据 - 假设接近线性扩展
    # 实际应该运行真实测试
    mock_data = {
        "format": format_type,
        "note": "MOCK DATA - Run actual multi-thread benchmark for real results",
        "results": [
            {"threads": 1, "execution_time_sec": base_time, "success": True},
            {"threads": 2, "execution_time_sec": base_time * 0.55, "success": True},
            {"threads": 4, "execution_time_sec": base_time * 0.30, "success": True},
            {"threads": 8, "execution_time_sec": base_time * 0.18, "success": True},
        ]
    }
    return mock_data


def main():
    print("=" * 60)
    print("Generating Performance Comparison Figure (Figure 1)")
    print("=" * 60)
    
    # Load data
    bed_data, bam_data, bed_mt_data, bam_mt_data = load_benchmark_data()
    
    if not bed_data and not bam_data:
        print("Error: No benchmark data found")
        print("Please run first:")
        print("  python paper/02_benchmark_bed.py")
        print("  python paper/03_benchmark_bam.py")
        return
    
    # If no multi-thread data, create mock data for demonstration
    use_mock = False
    if not bed_mt_data and bed_data:
        print("Warning: No BED multi-thread data found, using mock data")
        bed_mt_data = create_mock_multithread_data(bed_data, "BED")
        use_mock = True
    
    if not bam_mt_data and bam_data:
        print("Warning: No BAM multi-thread data found, using mock data")
        bam_mt_data = create_mock_multithread_data(bam_data, "BAM")
        use_mock = True
    
    # Create 2x2 figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Figure 1: FastCrossMap Performance Benchmark', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    # (a) BED 单线程对比
    if bed_data:
        print(f"BED data: {bed_data['input_records']:,} records")
        plot_single_thread_comparison(bed_data, axes[0, 0], "BED", TOOL_ORDER)
    else:
        axes[0, 0].text(0.5, 0.5, 'No BED data', ha='center', va='center', 
                        transform=axes[0, 0].transAxes)
        axes[0, 0].set_title('BED Single-Thread Comparison', fontsize=11, fontweight='bold')
    
    # (b) BED 多线程扩展性
    plot_multithread_scaling(bed_mt_data, axes[0, 1], "BED")
    
    # (c) BAM 单线程对比
    if bam_data:
        print(f"BAM data: {bam_data['input_size_mb']:.2f} MB")
        plot_single_thread_comparison(bam_data, axes[1, 0], "BAM", BAM_TOOL_ORDER)
    else:
        axes[1, 0].text(0.5, 0.5, 'No BAM data', ha='center', va='center', 
                        transform=axes[1, 0].transAxes)
        axes[1, 0].set_title('BAM Single-Thread Comparison', fontsize=11, fontweight='bold')
    
    # (d) BAM 多线程扩展性
    plot_multithread_scaling(bam_mt_data, axes[1, 1], "BAM")
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=COLORS[t], edgecolor='black', label=t) for t in TOOL_ORDER]
    fig.legend(handles=legend_elements, loc='upper center', ncol=4, 
               bbox_to_anchor=(0.5, 0.94), fontsize=10)
    
    # Add annotation
    if use_mock:
        fig.text(0.5, 0.01, '* Multi-thread data is simulated. Run actual benchmarks for real results.', 
                 ha='center', fontsize=9, style='italic', color='gray')
    
    # fig.text(0.99, 0.01, '* FastRemap uses internal 4 threads (not user-controllable)\n'
    #          '** BED multi-threading shows limited speedup due to I/O-bound nature', 
    #          ha='right', fontsize=8, style='italic', color='gray')
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.92])  # Leave more space at bottom for x-axis labels
    
    # Save combined figure
    output_pdf = FIGURES_DIR / "fig1_performance.pdf"
    output_png = FIGURES_DIR / "fig1_performance.png"
    
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(output_png, dpi=300, bbox_inches='tight')
    
    print(f"\nCombined figure saved to:")
    print(f"  {output_pdf}")
    print(f"  {output_png}")
    
    # Save individual subplots
    print(f"\nSaving individual subplots...")
    
    # (a) BED 单线程对比
    fig_a, ax_a = plt.subplots(figsize=(6, 5))
    if bed_data:
        plot_single_thread_comparison(bed_data, ax_a, "BED", TOOL_ORDER)
    plt.tight_layout()
    fig_a.savefig(FIGURES_DIR / "fig1a_bed_single_thread.pdf", dpi=300, bbox_inches='tight')
    fig_a.savefig(FIGURES_DIR / "fig1a_bed_single_thread.png", dpi=300, bbox_inches='tight')
    plt.close(fig_a)
    print(f"  {FIGURES_DIR / 'fig1a_bed_single_thread.pdf'}")
    
    # (b) BED 多线程扩展性
    fig_b, ax_b = plt.subplots(figsize=(6, 5))
    plot_multithread_scaling(bed_mt_data, ax_b, "BED")
    plt.tight_layout()
    fig_b.savefig(FIGURES_DIR / "fig1b_bed_multithread.pdf", dpi=300, bbox_inches='tight')
    fig_b.savefig(FIGURES_DIR / "fig1b_bed_multithread.png", dpi=300, bbox_inches='tight')
    plt.close(fig_b)
    print(f"  {FIGURES_DIR / 'fig1b_bed_multithread.pdf'}")
    
    # (c) BAM 单线程对比
    fig_c, ax_c = plt.subplots(figsize=(6, 5))
    if bam_data:
        plot_single_thread_comparison(bam_data, ax_c, "BAM", BAM_TOOL_ORDER)
    plt.tight_layout()
    fig_c.savefig(FIGURES_DIR / "fig1c_bam_single_thread.pdf", dpi=300, bbox_inches='tight')
    fig_c.savefig(FIGURES_DIR / "fig1c_bam_single_thread.png", dpi=300, bbox_inches='tight')
    plt.close(fig_c)
    print(f"  {FIGURES_DIR / 'fig1c_bam_single_thread.pdf'}")
    
    # (d) BAM 多线程扩展性
    fig_d, ax_d = plt.subplots(figsize=(6, 5))
    plot_multithread_scaling(bam_mt_data, ax_d, "BAM")
    plt.tight_layout()
    fig_d.savefig(FIGURES_DIR / "fig1d_bam_multithread.pdf", dpi=300, bbox_inches='tight')
    fig_d.savefig(FIGURES_DIR / "fig1d_bam_multithread.png", dpi=300, bbox_inches='tight')
    plt.close(fig_d)
    print(f"  {FIGURES_DIR / 'fig1d_bam_multithread.pdf'}")
    
    # Print performance summary
    print("\n" + "=" * 60)
    print("Performance Summary")
    print("=" * 60)
    
    if bed_data:
        results = {r["tool"]: r for r in bed_data["results"]}
        print("\nBED format (single-thread comparison):")
        for tool in TOOL_ORDER:
            if tool in results and results[tool]["success"]:
                t = results[tool]["execution_time_sec"]
                threads = "4T*" if tool == "FastRemap" else "1T"
                print(f"  {tool} ({threads}): {t:.2f}s")
        
        if "FastCrossMap" in results and "CrossMap" in results:
            fc = results["FastCrossMap"]
            cm = results["CrossMap"]
            if fc["success"] and cm["success"]:
                speedup = cm["execution_time_sec"] / fc["execution_time_sec"]
                print(f"  → FastCrossMap vs CrossMap: {speedup:.1f}x speedup")
    
    if bam_data:
        results = {r["tool"]: r for r in bam_data["results"]}
        print("\nBAM format (single-thread comparison):")
        for tool in BAM_TOOL_ORDER:
            if tool in results and results[tool]["success"]:
                t = results[tool]["execution_time_sec"]
                threads = "4T*" if tool == "FastRemap" else "1T"
                print(f"  {tool} ({threads}): {t:.2f}s")
        
        if "FastCrossMap" in results and "CrossMap" in results:
            fc = results["FastCrossMap"]
            cm = results["CrossMap"]
            if fc["success"] and cm["success"]:
                speedup = cm["execution_time_sec"] / fc["execution_time_sec"]
                print(f"  → FastCrossMap vs CrossMap: {speedup:.1f}x speedup")
    
    print("\n" + "=" * 60)
    print("Figure 1 Design Notes:")
    print("=" * 60)
    print("(a) BED single-thread comparison: Fair comparison across all tools (FastRemap uses internal 4 threads)")
    print("(b) BED multi-thread scaling: FastCrossMap 1/2/4/8/16 thread performance")
    print("    Note: Limited BED multi-thread speedup is expected (I/O-bound)")
    print("(c) BAM single-thread comparison: Fair comparison across all tools")
    print("(d) BAM multi-thread scaling: FastCrossMap 1/2/4/8/16 thread performance")
    print("\nNext step: python paper/05_memory_profile.py")


if __name__ == "__main__":
    main()
