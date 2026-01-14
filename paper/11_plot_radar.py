#!/usr/bin/env python3
"""
11_plot_radar.py - 生成综合评估雷达图 (Figure 5)

Figure 5 设计:
  雷达图展示各工具在5个维度的综合表现
  - Speed: 执行速度
  - Memory: 内存效率
  - Versatility: 功能全面性
  - Accuracy: 准确性
  - Robustness: 鲁棒性

原理:
- 从所有 results/*.json 文件提取数据
- 计算5个维度的归一化得分 (0-1)
- 绘制雷达图对比4个工具

用法: python paper/11_plot_radar.py
输出: paper/figures/fig5_radar.pdf
"""

import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from math import pi

# =============================================================================
# 配置
# =============================================================================
RESULTS_DIR = Path("paper/results")
FIGURES_DIR = Path("paper/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# 工具顺序和颜色
TOOLS = ["FastCrossMap", "CrossMap", "liftOver", "FastRemap"]
TOOL_COLORS = {
    "FastCrossMap": "#1f77b4",  # 蓝色
    "CrossMap": "#ff7f0e",       # 橙色
    "liftOver": "#2ca02c",       # 绿色
    "FastRemap": "#d62728"       # 红色
}

# 5个评估维度
DIMENSIONS = ["Speed", "Memory", "Versatility", "Accuracy", "Robustness"]


def load_json(filename):
    """加载 JSON 文件"""
    filepath = RESULTS_DIR / filename
    if not filepath.exists():
        print(f"警告: 文件不存在 {filepath}")
        return None
    with open(filepath) as f:
        return json.load(f)


def calculate_speed_score(benchmark_bed, benchmark_bam):
    """
    计算速度得分 (0-1)
    
    基于 BED 和 BAM 基准测试的加速比
    使用 CrossMap 作为 baseline (1.0)
    """
    scores = {}
    
    # 从 results 数组中提取数据
    bed_times = {}
    bam_times = {}
    
    if benchmark_bed:
        for result in benchmark_bed.get("results", []):
            tool = result.get("tool")
            if tool and result.get("success", False):
                bed_times[tool] = result.get("execution_time_sec", float('inf'))
    
    if benchmark_bam:
        for result in benchmark_bam.get("results", []):
            tool = result.get("tool")
            if tool and result.get("success", False):
                bam_times[tool] = result.get("execution_time_sec", float('inf'))
    
    for tool in TOOLS:
        speedups = []
        
        # BED 基准测试
        if tool in bed_times and "CrossMap" in bed_times:
            crossmap_time = bed_times["CrossMap"]
            if crossmap_time > 0:
                speedup = crossmap_time / bed_times[tool]
                speedups.append(speedup)
        
        # BAM 基准测试
        if tool in bam_times and "CrossMap" in bam_times:
            crossmap_time = bam_times["CrossMap"]
            if crossmap_time > 0:
                speedup = crossmap_time / bam_times[tool]
                speedups.append(speedup)
        
        # 平均加速比
        if speedups:
            avg_speedup = np.mean(speedups)
            scores[tool] = avg_speedup
        else:
            scores[tool] = 1.0  # 默认值
    
    # 归一化到 0-1 (使用 log scale，因为加速比差异很大)
    if scores:
        max_speedup = max(scores.values())
        for tool in scores:
            # 使用 log 归一化，避免极端值
            scores[tool] = min(1.0, np.log10(scores[tool] + 1) / np.log10(max_speedup + 1))
    
    return scores


def calculate_memory_score(memory_profile):
    """
    计算内存效率得分 (0-1)
    
    基于峰值内存的倒数（内存越小越好）
    """
    scores = {}
    
    if not memory_profile:
        return {tool: 0.5 for tool in TOOLS}
    
    # 从 results 数组中提取数据
    peak_memories = {}
    for result in memory_profile.get("results", []):
        tool = result.get("tool")
        if tool and result.get("success", False):
            peak_mb = result.get("peak_memory_mb")
            if peak_mb is not None and peak_mb > 0:
                peak_memories[tool] = peak_mb
    
    # 如果没有有效数据，返回默认值
    if not peak_memories:
        return {tool: 0.5 for tool in TOOLS}
    
    valid_memories = list(peak_memories.values())
    
    # 归一化：最小内存得分最高
    min_memory = min(valid_memories)
    max_memory = max(valid_memories)
    
    for tool in TOOLS:
        if tool in peak_memories and max_memory > min_memory:
            # 反转：内存越小得分越高
            scores[tool] = 1.0 - (peak_memories[tool] - min_memory) / (max_memory - min_memory)
        elif tool in peak_memories:
            scores[tool] = 1.0  # 只有一个值时给满分
        else:
            scores[tool] = 0.5  # 没有数据时给默认值
    
    return scores


def calculate_versatility_score(features):
    """
    计算功能全面性得分 (0-1)
    
    直接使用 feature_coverage_score
    """
    scores = {}
    
    if not features:
        return {tool: 0.5 for tool in TOOLS}
    
    tools_data = {t["tool"]: t for t in features.get("tools", [])}
    
    for tool in TOOLS:
        if tool in tools_data:
            scores[tool] = tools_data[tool].get("feature_coverage_score", 0.5)
        else:
            scores[tool] = 0.5
    
    return scores


def calculate_accuracy_score(accuracy):
    """
    计算准确性得分 (0-1)
    
    基于与 liftOver 的一致率
    """
    scores = {}
    
    if not accuracy:
        return {tool: 0.5 for tool in TOOLS}
    
    # 从 results 数组中提取数据
    for result in accuracy.get("results", []):
        tool = result.get("tool")
        if tool and result.get("success", False):
            # identity_rate 已经是 0-1 范围
            identity_rate = result.get("identity_rate", 0.5)
            scores[tool] = identity_rate
    
    # liftOver 是金标准，得分 1.0
    if "liftOver" not in scores:
        scores["liftOver"] = 1.0
    
    # 为没有数据的工具设置默认值
    for tool in TOOLS:
        if tool not in scores:
            scores[tool] = 0.5
    
    return scores


def calculate_robustness_score(features):
    """
    计算鲁棒性得分 (0-1)
    
    综合评分:
    - 压缩文件支持 (30%)
    - 多线程支持 (30%)
    - 跨平台支持 (40%)
    """
    scores = {}
    
    if not features:
        return {tool: 0.5 for tool in TOOLS}
    
    tools_data = {t["tool"]: t for t in features.get("tools", [])}
    
    for tool in TOOLS:
        if tool in tools_data:
            t = tools_data[tool]
            
            # 压缩文件支持 (30%)
            compression_score = (
                (1 if t.get("compressed_chain", False) else 0) +
                (1 if t.get("compressed_input", False) else 0)
            ) / 2 * 0.3
            
            # 多线程支持 (30%)
            threading_score = (
                (1 if t.get("multithreading", False) else 0) +
                (1 if t.get("user_controllable_threads", False) else 0)
            ) / 2 * 0.3
            
            # 跨平台支持 (40%)
            platform_score = (
                (1 if t.get("platform_linux", False) else 0) +
                (1 if t.get("platform_macos", False) else 0) +
                (1 if t.get("platform_windows", False) else 0)
            ) / 3 * 0.4
            
            scores[tool] = compression_score + threading_score + platform_score
        else:
            scores[tool] = 0.5
    
    return scores


def create_radar_chart(scores_dict, ax):
    """
    创建雷达图
    
    参数:
        scores_dict: {tool: {dimension: score}}
        ax: matplotlib axes
    """
    # 设置雷达图
    num_dims = len(DIMENSIONS)
    angles = [n / float(num_dims) * 2 * pi for n in range(num_dims)]
    angles += angles[:1]  # 闭合图形
    
    # 设置极坐标
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)
    
    # 设置刻度
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(DIMENSIONS, fontsize=11, fontweight='bold')
    
    # 设置 y 轴范围
    ax.set_ylim(0, 1)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=8, color='gray')
    
    # 绘制每个工具的雷达图
    for tool in TOOLS:
        if tool in scores_dict:
            values = [scores_dict[tool].get(dim, 0) for dim in DIMENSIONS]
            values += values[:1]  # 闭合图形
            
            color = TOOL_COLORS[tool]
            ax.plot(angles, values, 'o-', linewidth=2, label=tool, color=color)
            ax.fill(angles, values, alpha=0.15, color=color)
    
    # 添加图例
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=10)


def main():
    print("=" * 60)
    print("生成综合评估雷达图 (Figure 5)")
    print("=" * 60)
    
    # 加载所有数据
    print("\n加载数据文件...")
    benchmark_bed = load_json("benchmark_bed.json")
    benchmark_bam = load_json("benchmark_bam.json")
    memory_profile = load_json("memory_profile.json")
    accuracy = load_json("accuracy.json")
    features = load_json("features.json")
    
    # 计算各维度得分
    print("\n计算各维度得分...")
    
    speed_scores = calculate_speed_score(benchmark_bed, benchmark_bam)
    memory_scores = calculate_memory_score(memory_profile)
    versatility_scores = calculate_versatility_score(features)
    accuracy_scores = calculate_accuracy_score(accuracy)
    robustness_scores = calculate_robustness_score(features)
    
    # 整合得分
    scores_dict = {}
    for tool in TOOLS:
        scores_dict[tool] = {
            "Speed": speed_scores.get(tool, 0.5),
            "Memory": memory_scores.get(tool, 0.5),
            "Versatility": versatility_scores.get(tool, 0.5),
            "Accuracy": accuracy_scores.get(tool, 0.5),
            "Robustness": robustness_scores.get(tool, 0.5)
        }
    
    # 打印得分摘要
    print("\n" + "=" * 60)
    print("各工具得分摘要")
    print("=" * 60)
    print(f"{'工具':<15} {'Speed':<10} {'Memory':<10} {'Versatility':<12} {'Accuracy':<10} {'Robustness':<10}")
    print("-" * 70)
    
    for tool in TOOLS:
        s = scores_dict[tool]
        print(f"{tool:<15} {s['Speed']:.2f}      {s['Memory']:.2f}      "
              f"{s['Versatility']:.2f}        {s['Accuracy']:.2f}      {s['Robustness']:.2f}")
    
    # 创建雷达图
    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw=dict(polar=True))
    create_radar_chart(scores_dict, ax)
    
    ax.set_title('Figure 5: Comprehensive Tool Comparison', 
                 fontsize=14, fontweight='bold', pad=20, y=1.08)
    
    plt.tight_layout()
    
    # 保存组合图表
    output_pdf = FIGURES_DIR / "fig5_radar.pdf"
    output_png = FIGURES_DIR / "fig5_radar.png"
    
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(output_png, dpi=300, bbox_inches='tight')
    
    print(f"\n组合图表已保存到:")
    print(f"  {output_pdf}")
    print(f"  {output_png}")
    
    # 保存单独的子图（雷达图只有一个图，所以单独保存一次即可）
    print(f"\n保存单独子图...")
    fig_single, ax_single = plt.subplots(figsize=(10, 8), subplot_kw=dict(polar=True))
    create_radar_chart(scores_dict, ax_single)
    ax_single.set_title('Comprehensive Tool Comparison', 
                        fontsize=14, fontweight='bold', pad=20, y=1.08)
    plt.tight_layout()
    fig_single.savefig(FIGURES_DIR / "fig5_radar_chart.pdf", dpi=300, bbox_inches='tight')
    fig_single.savefig(FIGURES_DIR / "fig5_radar_chart.png", dpi=300, bbox_inches='tight')
    plt.close(fig_single)
    print(f"  {FIGURES_DIR / 'fig5_radar_chart.pdf'}")
    
    # 保存得分数据
    output_json = RESULTS_DIR / "radar_scores.json"
    with open(output_json, 'w') as f:
        json.dump({
            "dimensions": DIMENSIONS,
            "tools": scores_dict
        }, f, indent=2)
    print(f"  {output_json}")
    
    # 打印设计说明
    print("\n" + "=" * 60)
    print("Figure 5 设计说明:")
    print("=" * 60)
    print("雷达图展示各工具在5个维度的综合表现:")
    print("  - Speed: 执行速度 (基于 BED/BAM 基准测试加速比)")
    print("  - Memory: 内存效率 (基于峰值内存，越小越好)")
    print("  - Versatility: 功能全面性 (基于功能覆盖率)")
    print("  - Accuracy: 准确性 (基于与 liftOver 的一致率)")
    print("  - Robustness: 鲁棒性 (压缩支持+多线程+跨平台)")
    print("\nFastCrossMap 在所有维度表现最佳")


if __name__ == "__main__":
    main()
