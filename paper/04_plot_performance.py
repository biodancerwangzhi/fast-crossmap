#!/usr/bin/env python3
"""
04_plot_performance.py - 生成性能对比图 (Figure 1)

Figure 1 设计 (2x2 布局):
  (a) BED 单线程性能对比 - 所有工具公平比较
  (b) BED 多线程扩展性 - FastCrossMap 1/2/4/8 线程
  (c) BAM 单线程性能对比 - 所有工具公平比较  
  (d) BAM 多线程扩展性 - FastCrossMap 1/2/4/8 线程

设计理念:
- 左列 (a,c): 单线程公平比较，展示算法效率
- 右列 (b,d): 多线程扩展性，展示 FastCrossMap 的并行能力
- 只展示执行时间，避免与吞吐量冗余

用法: python paper/04_plot_performance.py
输出: paper/figures/fig1_performance.pdf
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

# 工具颜色
COLORS = {
    "FastCrossMap": "#1f77b4",  # 蓝色
    "CrossMap": "#ff7f0e",       # 橙色
    "liftOver": "#2ca02c",       # 绿色
    "FastRemap": "#d62728"       # 红色
}

# 多线程颜色渐变 (蓝色系)
THREAD_COLORS = {
    1: "#1f77b4",   # 深蓝
    2: "#4a9fd4",   # 中蓝
    4: "#7ec8e3",   # 浅蓝
    8: "#aee0f0",   # 浅蓝
    16: "#d4eef7"   # 最浅蓝
}

# 工具顺序
TOOL_ORDER = ["FastCrossMap", "CrossMap", "liftOver", "FastRemap"]
BAM_TOOL_ORDER = ["FastCrossMap", "CrossMap", "FastRemap"]


def load_benchmark_data():
    """加载基准测试数据"""
    bed_data = None
    bam_data = None
    bed_mt_data = None  # 多线程数据
    bam_mt_data = None
    
    bed_file = RESULTS_DIR / "benchmark_bed.json"
    if bed_file.exists():
        with open(bed_file) as f:
            bed_data = json.load(f)
    
    bam_file = RESULTS_DIR / "benchmark_bam.json"
    if bam_file.exists():
        with open(bam_file) as f:
            bam_data = json.load(f)
    
    # 多线程数据 (如果存在)
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
    
    # 准备数据
    tools = []
    all_times_data = []
    
    for tool in tool_order:
        if tool in results and results[tool]["success"]:
            tools.append(tool)
            # 获取所有运行时间 (用于箱线图)
            if "all_times" in results[tool] and results[tool]["all_times"]:
                all_times_data.append(results[tool]["all_times"])
            else:
                all_times_data.append([results[tool]["execution_time_sec"]])
    
    if not tools:
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', transform=ax.transAxes)
        return
    
    colors = [COLORS[t] for t in tools]
    
    # 使用箱线图展示多次运行的分布
    bp = ax.boxplot(all_times_data, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # 不再添加中位数标签 (避免遮挡图片)
    
    ax.set_ylabel('Execution Time (seconds)', fontsize=10)
    ax.set_title(f'{format_type} Single-Thread Comparison', 
                 fontsize=11, fontweight='bold')
    
    # 构建 x 轴标签，FastCrossMap 用红色突出，FastRemap 的 4T* 用红色
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
    
    # 设置 x 轴标签
    ax.set_xticks(range(1, len(tools) + 1))
    ax.set_xticklabels([])  # 先清空，然后手动添加带颜色的标签
    
    # 手动添加带颜色的 x 轴标签
    for i, (tool, label) in enumerate(zip(tools, labels)):
        if tool == "FastCrossMap":
            # FastCrossMap 整体红色
            ax.text(i + 1, -0.12, label, ha='center', va='top', 
                    transform=ax.get_xaxis_transform(), fontsize=9, 
                    color='red', fontweight='bold')
        elif tool == "FastRemap":
            # FastRemap: 工具名黑色，4T* 红色
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
    
    # 解析多线程数据
    thread_data = []  # [(threads, all_times), ...]
    
    for result in mt_data.get("results", []):
        if result.get("success", False):
            threads = result["threads"]
            # 获取所有运行时间 (用于箱线图)
            if "all_times" in result and result["all_times"]:
                all_times = result["all_times"]
            else:
                all_times = [result["execution_time_sec"]]
            thread_data.append((threads, all_times))
    
    if not thread_data:
        ax.text(0.5, 0.5, 'No valid multi-thread data', ha='center', va='center', 
                transform=ax.transAxes)
        return
    
    # 按线程数排序
    thread_data.sort(key=lambda x: x[0])
    threads = [d[0] for d in thread_data]
    all_times_data = [d[1] for d in thread_data]
    
    # 使用箱线图
    colors = [THREAD_COLORS.get(t, "#1f77b4") for t in threads]
    bp = ax.boxplot(all_times_data, patch_artist=True)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # 计算加速比并在 X 轴下方显示
    baseline_median = np.median(all_times_data[0])  # 1 线程作为基准
    
    # 设置 x 轴标签
    ax.set_xticks(range(1, len(threads) + 1))
    ax.set_xticklabels([])  # 先清空
    
    # 手动添加 x 轴标签和加速比
    for i, (th, times) in enumerate(zip(threads, all_times_data)):
        median = np.median(times)
        speedup = baseline_median / median
        
        # 线程数标签
        ax.text(i + 1, -0.06, f'{th}T', ha='center', va='top', 
                transform=ax.get_xaxis_transform(), fontsize=9, color='black')
        
        # 加速比标签 (在线程数下方) - 全部使用黑色
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
    print("生成性能对比图 (Figure 1)")
    print("=" * 60)
    
    # 加载数据
    bed_data, bam_data, bed_mt_data, bam_mt_data = load_benchmark_data()
    
    if not bed_data and not bam_data:
        print("错误: 没有找到基准测试数据")
        print("请先运行:")
        print("  python paper/02_benchmark_bed.py")
        print("  python paper/03_benchmark_bam.py")
        return
    
    # 如果没有多线程数据，创建模拟数据用于演示
    use_mock = False
    if not bed_mt_data and bed_data:
        print("警告: 没有找到 BED 多线程数据，使用模拟数据")
        bed_mt_data = create_mock_multithread_data(bed_data, "BED")
        use_mock = True
    
    if not bam_mt_data and bam_data:
        print("警告: 没有找到 BAM 多线程数据，使用模拟数据")
        bam_mt_data = create_mock_multithread_data(bam_data, "BAM")
        use_mock = True
    
    # 创建 2x2 图表
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Figure 1: FastCrossMap Performance Benchmark', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    # (a) BED 单线程对比
    if bed_data:
        print(f"BED 数据: {bed_data['input_records']:,} 记录")
        plot_single_thread_comparison(bed_data, axes[0, 0], "BED", TOOL_ORDER)
    else:
        axes[0, 0].text(0.5, 0.5, 'No BED data', ha='center', va='center', 
                        transform=axes[0, 0].transAxes)
        axes[0, 0].set_title('BED Single-Thread Comparison', fontsize=11, fontweight='bold')
    
    # (b) BED 多线程扩展性
    plot_multithread_scaling(bed_mt_data, axes[0, 1], "BED")
    
    # (c) BAM 单线程对比
    if bam_data:
        print(f"BAM 数据: {bam_data['input_size_mb']:.2f} MB")
        plot_single_thread_comparison(bam_data, axes[1, 0], "BAM", BAM_TOOL_ORDER)
    else:
        axes[1, 0].text(0.5, 0.5, 'No BAM data', ha='center', va='center', 
                        transform=axes[1, 0].transAxes)
        axes[1, 0].set_title('BAM Single-Thread Comparison', fontsize=11, fontweight='bold')
    
    # (d) BAM 多线程扩展性
    plot_multithread_scaling(bam_mt_data, axes[1, 1], "BAM")
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=COLORS[t], edgecolor='black', label=t) for t in TOOL_ORDER]
    fig.legend(handles=legend_elements, loc='upper center', ncol=4, 
               bbox_to_anchor=(0.5, 0.94), fontsize=10)
    
    # 添加注释
    if use_mock:
        fig.text(0.5, 0.01, '* Multi-thread data is simulated. Run actual benchmarks for real results.', 
                 ha='center', fontsize=9, style='italic', color='gray')
    
    # fig.text(0.99, 0.01, '* FastRemap uses internal 4 threads (not user-controllable)\n'
    #          '** BED multi-threading shows limited speedup due to I/O-bound nature', 
    #          ha='right', fontsize=8, style='italic', color='gray')
    
    plt.tight_layout(rect=[0, 0.05, 1, 0.92])  # 底部留更多空间给 x 轴标签
    
    # 保存组合图表
    output_pdf = FIGURES_DIR / "fig1_performance.pdf"
    output_png = FIGURES_DIR / "fig1_performance.png"
    
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(output_png, dpi=300, bbox_inches='tight')
    
    print(f"\n组合图表已保存到:")
    print(f"  {output_pdf}")
    print(f"  {output_png}")
    
    # 保存单独的子图
    print(f"\n保存单独子图...")
    
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
    
    # 打印性能摘要
    print("\n" + "=" * 60)
    print("性能摘要")
    print("=" * 60)
    
    if bed_data:
        results = {r["tool"]: r for r in bed_data["results"]}
        print("\nBED 格式 (单线程对比):")
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
                print(f"  → FastCrossMap vs CrossMap: {speedup:.1f}x 加速")
    
    if bam_data:
        results = {r["tool"]: r for r in bam_data["results"]}
        print("\nBAM 格式 (单线程对比):")
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
                print(f"  → FastCrossMap vs CrossMap: {speedup:.1f}x 加速")
    
    print("\n" + "=" * 60)
    print("Figure 1 设计说明:")
    print("=" * 60)
    print("(a) BED 单线程对比: 所有工具公平比较 (FastRemap 内部使用 4 线程)")
    print("(b) BED 多线程扩展: FastCrossMap 1/2/4/8/16 线程性能")
    print("    注: BED 多线程加速有限是预期行为 (I/O 密集型)")
    print("(c) BAM 单线程对比: 所有工具公平比较")
    print("(d) BAM 多线程扩展: FastCrossMap 1/2/4/8/16 线程性能")
    print("\n下一步: python paper/05_memory_profile.py")


if __name__ == "__main__":
    main()
