#!/usr/bin/env python3
"""
10_plot_features.py - 生成功能对比热力图 (Figure 4)

Figure 4 设计:
  热力图展示各工具的功能支持情况
  - 行: 功能项
  - 列: 工具
  - 颜色: 支持(绿)/不支持(红)

设计理念:
- 直观展示功能覆盖率差异
- 突出 FastCrossMap 的全面性
- 易于识别各工具的优劣势

用法: python paper/10_plot_features.py
输出: paper/figures/fig4_features.pdf
"""

import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

# =============================================================================
# 配置
# =============================================================================
RESULTS_DIR = Path("paper/results")
FIGURES_DIR = Path("paper/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# 工具顺序
TOOL_ORDER = ["FastCrossMap", "CrossMap", "liftOver", "FastRemap"]


def load_features_data():
    """加载功能审计数据"""
    features_file = RESULTS_DIR / "features.json"
    if not features_file.exists():
        return None
    
    with open(features_file) as f:
        return json.load(f)


def create_feature_heatmap(data, ax):
    """
    创建功能热力图
    
    参数:
        data: 功能审计数据
        ax: matplotlib axes
    """
    if not data:
        ax.text(0.5, 0.5, 'No features data', ha='center', va='center', 
                transform=ax.transAxes, fontsize=12)
        ax.set_title('Feature Support Matrix', fontsize=14, fontweight='bold')
        ax.axis('off')
        return
    
    tools_data = {t["tool"]: t for t in data["tools"]}
    
    # 定义功能列表（不包含类别行，类别信息保留用于分组）
    features = [
        # 文件格式支持
        ("BED", "format_bed", "File Formats"),
        ("BAM/SAM", "format_bam", "File Formats"),
        ("VCF", "format_vcf", "File Formats"),
        ("GFF/GTF", "format_gff", "File Formats"),
        ("Wiggle", "format_wiggle", "File Formats"),
        ("BigWig", "format_bigwig", "File Formats"),
        ("MAF", "format_maf", "File Formats"),
        ("GVCF", "format_gvcf", "File Formats"),
        
        # 压缩文件支持
        ("Compressed Chain", "compressed_chain", "Compression"),
        ("Compressed Input", "compressed_input", "Compression"),
        
        # 多线程支持
        ("Multithreading", "multithreading", "Threading"),
        ("User Control Threads", "user_controllable_threads", "Threading"),
        
        # 跨平台支持
        ("Linux", "platform_linux", "Platforms"),
        ("macOS", "platform_macos", "Platforms"),
        ("Windows", "platform_windows", "Platforms"),
        
        # 其他功能
        ("Unmapped Output", "unmapped_output", "Other"),
        ("Streaming Process", "streaming_processing", "Other"),
    ]
    
    # 构建矩阵
    matrix = []
    feature_labels = []
    
    for feature_name, feature_attr, category in features:
        row = []
        for tool in TOOL_ORDER:
            if tool in tools_data:
                value = tools_data[tool].get(feature_attr, False)
                row.append(1 if value else 0)
            else:
                row.append(0)
        matrix.append(row)
        feature_labels.append(feature_name)
    
    matrix = np.array(matrix)
    
    # 创建热力图
    cmap = ListedColormap(['#ffcccc', '#ccffcc'])  # 红色(不支持), 绿色(支持)
    im = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=0, vmax=1)
    
    # 定义类别颜色
    category_colors = {
        "File Formats": "#1f77b4",    # 蓝色
        "Compression": "#ff7f0e",     # 橙色
        "Threading": "#2ca02c",       # 绿色
        "Platforms": "#d62728",       # 红色
        "Other": "#9467bd"            # 紫色
    }
    
    # 为每个功能项分配类别颜色
    feature_colors = []
    for _, _, category in features:
        feature_colors.append(category_colors[category])
    
    # 设置坐标轴
    ax.set_xticks(np.arange(len(TOOL_ORDER)))
    ax.set_yticks(np.arange(len(feature_labels)))
    ax.set_xticklabels(TOOL_ORDER, fontsize=10, fontweight='bold')
    ax.set_yticklabels(feature_labels, fontsize=9)
    
    # 为Y轴标签设置颜色
    for tick_label, color in zip(ax.get_yticklabels(), feature_colors):
        tick_label.set_color(color)
    
    # 旋转 x 轴标签
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    
    # 在每个单元格中添加符号
    for i in range(len(feature_labels)):
        for j in range(len(TOOL_ORDER)):
            text = "✓" if matrix[i, j] == 1 else "✗"
            color = "darkgreen" if matrix[i, j] == 1 else "darkred"
            ax.text(j, i, text, ha="center", va="center", 
                   color=color, fontsize=12, fontweight='bold')
    
    # 添加分类分隔线和类别标签
    category_positions = []
    category_labels = {}
    current_category = None
    
    for i, (_, _, category) in enumerate(features):
        if category != current_category:
            if current_category is not None:
                category_positions.append(i - 0.5)
            current_category = category
            category_labels[category] = i
    
    # 画分隔线（黑色）
    for pos in category_positions:
        ax.axhline(y=pos, color='black', linewidth=2, linestyle='-')
    
    # 添加类别标签（在Y轴左侧，与每个分组对齐）
    category_ranges = {
        "File Formats": (0, 7),
        "Compression": (8, 9),
        "Threading": (10, 11),
        "Platforms": (12, 14),
        "Other": (15, 16)
    }
    
    for category, (start, end) in category_ranges.items():
        # 计算分组的中间位置
        mid_row = (start + end) / 2
        color = category_colors[category]
        
        # 在Y轴左侧显示类别标签（带颜色）
        ax.text(-1.3, mid_row, f"{category}", ha='left', va='center', 
               fontsize=9, fontweight='bold', 
               rotation=90, color=color)
    
    ax.set_title('Feature Support Matrix', fontsize=14, fontweight='bold', pad=20)
    
    # 移除边框
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    ax.set_xticks(np.arange(len(TOOL_ORDER)+1)-.5, minor=True)
    ax.set_yticks(np.arange(len(feature_labels)+1)-.5, minor=True)
    ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", size=0)


def main():
    print("=" * 60)
    print("生成功能对比热力图 (Figure 4)")
    print("=" * 60)
    
    # 加载数据
    features_data = load_features_data()
    
    if not features_data:
        print("错误: 没有找到功能审计数据")
        print("请先运行: python paper/09_feature_audit.py")
        return
    
    print(f"工具数: {len(features_data['tools'])}")
    
    # 创建图表 (只有热力图)
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # 功能热力图
    create_feature_heatmap(features_data, ax)
    
    # fig.suptitle('Feature Support Comparison', 
    #              fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # 保存组合图表
    output_pdf = FIGURES_DIR / "fig4_features.pdf"
    output_png = FIGURES_DIR / "fig4_features.png"
    
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(output_png, dpi=300, bbox_inches='tight')
    
    print(f"\n组合图表已保存到:")
    print(f"  {output_pdf}")
    print(f"  {output_png}")
    
    # 保存单独的子图（功能热力图只有一个图，所以单独保存一次即可）
    print(f"\n保存单独子图...")
    fig_single, ax_single = plt.subplots(figsize=(10, 10))
    create_feature_heatmap(features_data, ax_single)
    plt.tight_layout()
    fig_single.savefig(FIGURES_DIR / "fig4_features_heatmap.pdf", dpi=300, bbox_inches='tight')
    fig_single.savefig(FIGURES_DIR / "fig4_features_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close(fig_single)
    print(f"  {FIGURES_DIR / 'fig4_features_heatmap.pdf'}")
    
    # 打印功能摘要
    print("\n" + "=" * 60)
    print("功能覆盖率摘要")
    print("=" * 60)
    
    tools_data = {t["tool"]: t for t in features_data["tools"]}
    
    for tool in TOOL_ORDER:
        if tool in tools_data:
            t = tools_data[tool]
            print(f"\n{tool}:")
            print(f"  格式支持: {t['format_count']}/8")
            print(f"  平台支持: {t['platform_count']}/3")
            print(f"  压缩 Chain: {'✓' if t['compressed_chain'] else '✗'}")
            print(f"  多线程: {'✓' if t['multithreading'] else '✗'}")
            print(f"  可控线程: {'✓' if t['user_controllable_threads'] else '✗'}")
            print(f"  覆盖率: {t['feature_coverage_score']*100:.1f}%")
    
    print("\n" + "=" * 60)
    print("Figure 4 设计说明:")
    print("=" * 60)
    print("功能支持矩阵热力图 - 直观展示各工具的功能支持情况")
    print("  绿色(✓)表示支持，红色(✗)表示不支持")
    print("  按类别分组：文件格式、压缩、多线程、平台、其他")
    print("  FastCrossMap 支持最全面 (8/8 格式, 3/3 平台)")
    print("\n下一步: python paper/11_plot_radar.py")


if __name__ == "__main__":
    main()
