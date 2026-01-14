#!/usr/bin/env python3
"""
06_plot_memory.py - 生成内存效率对比图 (Figure 2)

Figure 2 设计 (1x2 布局):
  左图: 内存使用曲线 - 展示三个工具随时间的内存变化
  右图: 峰值内存对比 - 条形图展示峰值内存对比

设计理念:
- 左图展示内存使用的动态变化，体现 FastCrossMap 的流式处理优势
- 右图展示峰值内存的直接对比，便于量化比较
- 使用与 Figure 1 一致的配色方案

用法: python paper/06_plot_memory.py
输出: paper/figures/fig2_memory.pdf
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

# 工具颜色 (与 Figure 1 保持一致)
COLORS = {
    "FastCrossMap": "#1f77b4",  # 蓝色
    "CrossMap": "#ff7f0e",       # 橙色
    "FastRemap": "#d62728"       # 红色
}

# 工具顺序
TOOL_ORDER = ["FastCrossMap", "CrossMap", "FastRemap"]


def load_memory_data():
    """加载内存采样数据"""
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
    
    # 绘制每个工具的内存曲线
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
    
    # 准备数据
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
    
    # 绘制条形图
    bars = ax.bar(range(len(tools)), peak_memories, color=colors, alpha=0.7, edgecolor='black')
    
    # 在条形图上方添加数值标签
    for i, (bar, mem) in enumerate(zip(bars, peak_memories)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{mem:.1f} MB',
                ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax.set_ylabel('Peak Memory (MB)', fontsize=10)
    ax.set_title('Peak Memory Comparison', fontsize=11, fontweight='bold')
    ax.set_xticks(range(len(tools)))
    
    # 构建 x 轴标签，FastCrossMap 用红色突出
    ax.set_xticklabels([])  # 先清空
    
    for i, tool in enumerate(tools):
        if tool == "FastCrossMap":
            ax.text(i, -0.12, tool, ha='center', va='top', 
                    transform=ax.get_xaxis_transform(), fontsize=9, 
                    color='red', fontweight='bold')
        else:
            ax.text(i, -0.12, tool, ha='center', va='top', 
                    transform=ax.get_xaxis_transform(), fontsize=9, 
                    color='black')
    
    # 添加网格线
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')


def main():
    print("=" * 60)
    print("生成内存效率对比图 (Figure 2)")
    print("=" * 60)
    
    # 加载数据
    memory_data = load_memory_data()
    
    if not memory_data:
        print("错误: 没有找到内存采样数据")
        print("请先运行: python paper/05_memory_profile.py")
        return
    
    print(f"输入文件: {memory_data['input_file']}")
    print(f"文件大小: {memory_data['input_size_mb']:.2f} MB")
    print(f"采样间隔: {memory_data['sample_interval_sec']} 秒")
    
    # 创建 1x2 图表
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Figure 2: Memory Efficiency Comparison', 
                 fontsize=14, fontweight='bold', y=1.00)
    
    # 左图: 内存使用曲线
    plot_memory_curves(memory_data, axes[0])
    
    # 右图: 峰值内存对比
    plot_peak_memory_comparison(memory_data, axes[1])
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # 保存组合图表
    output_pdf = FIGURES_DIR / "fig2_memory.pdf"
    output_png = FIGURES_DIR / "fig2_memory.png"
    
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(output_png, dpi=300, bbox_inches='tight')
    
    print(f"\n组合图表已保存到:")
    print(f"  {output_pdf}")
    print(f"  {output_png}")
    
    # 保存单独的子图
    print(f"\n保存单独子图...")
    
    # 左图: 内存使用曲线
    fig_left, ax_left = plt.subplots(figsize=(7, 5))
    plot_memory_curves(memory_data, ax_left)
    plt.tight_layout()
    fig_left.savefig(FIGURES_DIR / "fig2a_memory_curves.pdf", dpi=300, bbox_inches='tight')
    fig_left.savefig(FIGURES_DIR / "fig2a_memory_curves.png", dpi=300, bbox_inches='tight')
    plt.close(fig_left)
    print(f"  {FIGURES_DIR / 'fig2a_memory_curves.pdf'}")
    
    # 右图: 峰值内存对比
    fig_right, ax_right = plt.subplots(figsize=(6, 5))
    plot_peak_memory_comparison(memory_data, ax_right)
    plt.tight_layout()
    fig_right.savefig(FIGURES_DIR / "fig2b_peak_memory.pdf", dpi=300, bbox_inches='tight')
    fig_right.savefig(FIGURES_DIR / "fig2b_peak_memory.png", dpi=300, bbox_inches='tight')
    plt.close(fig_right)
    print(f"  {FIGURES_DIR / 'fig2b_peak_memory.pdf'}")
    
    # 打印内存摘要
    print("\n" + "=" * 60)
    print("内存效率摘要")
    print("=" * 60)
    
    results = {r["tool"]: r for r in memory_data["results"]}
    
    for tool in TOOL_ORDER:
        if tool in results and results[tool]["success"]:
            r = results[tool]
            print(f"{tool}:")
            print(f"  执行时间: {r['execution_time_sec']:.2f}s")
            print(f"  峰值内存: {r['peak_memory_mb']:.2f} MB")
            print(f"  采样数: {len(r['memory_samples'])}")
    
    # 计算内存节省
    if "FastCrossMap" in results and "CrossMap" in results:
        fc = results["FastCrossMap"]
        cm = results["CrossMap"]
        if fc["success"] and cm["success"]:
            mem_ratio = cm["peak_memory_mb"] / fc["peak_memory_mb"]
            mem_saved = cm["peak_memory_mb"] - fc["peak_memory_mb"]
            print(f"\n内存效率对比:")
            print(f"  FastCrossMap vs CrossMap:")
            print(f"    内存节省: {mem_saved:.2f} MB ({(1 - 1/mem_ratio)*100:.1f}%)")
            print(f"    CrossMap 使用 {mem_ratio:.2f}x 内存")
    
    print("\n" + "=" * 60)
    print("Figure 2 设计说明:")
    print("=" * 60)
    print("左图: 内存使用曲线 - 展示三个工具随时间的内存变化")
    print("      体现 FastCrossMap 的流式处理优势 (内存稳定)")
    print("右图: 峰值内存对比 - 条形图展示峰值内存的直接对比")
    print("      便于量化比较内存效率")
    print("\n下一步: python paper/07_accuracy_analysis.py")


if __name__ == "__main__":
    main()
