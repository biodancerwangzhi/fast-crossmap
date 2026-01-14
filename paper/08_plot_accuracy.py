#!/usr/bin/env python3
"""
08_plot_accuracy.py - 生成准确性对比图 (Figure 3)

Figure 3 设计:
  表格形式展示准确性指标对比
  - 映射率 (Mapping Rate)
  - 一致率 (Identity Rate with liftOver)
  - 坐标偏差统计 (Coordinate Drift)
  - 染色体偏差 (Chromosome Mismatch)

设计理念:
- 使用表格清晰展示数值指标
- 突出 FastCrossMap 与 liftOver 的高一致性
- 量化各工具的准确性差异

用法: python paper/08_plot_accuracy.py
输出: paper/figures/fig3_accuracy.pdf
"""

import json
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.table import Table

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


def load_accuracy_data():
    """加载准确性分析数据"""
    accuracy_file = RESULTS_DIR / "accuracy.json"
    if not accuracy_file.exists():
        return None
    
    with open(accuracy_file) as f:
        return json.load(f)


def create_accuracy_table(data, ax):
    """
    创建准确性对比表格
    
    参数:
        data: 准确性分析数据
        ax: matplotlib axes
    """
    if not data:
        ax.text(0.5, 0.5, 'No accuracy data', ha='center', va='center', 
                transform=ax.transAxes, fontsize=12)
        ax.set_title('Accuracy Comparison', fontsize=14, fontweight='bold')
        ax.axis('off')
        return
    
    results = {r["tool"]: r for r in data["results"]}
    
    # 准备表格数据
    headers = ["Tool", "Mapping\nRate (%)", "Identity\nRate (%)", 
               "Partial\nMatch", "Coord.\nMismatch"]
    
    table_data = []
    for tool in TOOL_ORDER:
        if tool in results and results[tool]["success"]:
            r = results[tool]
            
            row = [
                tool,
                f"{r['mapping_rate']*100:.2f}",
                f"{r['identity_rate']*100:.2f}",
                str(r.get('partial_match', 0)),
                str(r.get('coordinate_mismatch', 0))
            ]
            table_data.append(row)
    
    if not table_data:
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', 
                transform=ax.transAxes)
        ax.axis('off')
        return
    
    # 创建表格
    ax.axis('off')
    
    # 计算表格位置和大小
    n_rows = len(table_data) + 1  # +1 for header
    n_cols = len(headers)
    
    table = ax.table(cellText=table_data, colLabels=headers,
                    cellLoc='center', loc='center',
                    bbox=[0.05, 0.2, 0.9, 0.6])
    
    # 设置表格样式
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2.5)
    
    # 设置表头样式
    for i in range(n_cols):
        cell = table[(0, i)]
        cell.set_facecolor('#4CAF50')
        cell.set_text_props(weight='bold', color='white')
    
    # 设置工具名称列的颜色
    for i, tool in enumerate(TOOL_ORDER):
        if i < len(table_data):
            cell = table[(i+1, 0)]
            cell.set_facecolor(COLORS.get(tool, 'white'))
            cell.set_text_props(weight='bold', color='white')
            
            # 如果是 FastCrossMap，整行背景稍微高亮
            if tool == "FastCrossMap":
                for j in range(1, n_cols):
                    table[(i+1, j)].set_facecolor('#E3F2FD')
    
    # 添加标题
    ax.text(0.5, 0.95, 'Accuracy Comparison (Gold Standard: liftOver)', 
            ha='center', va='top', transform=ax.transAxes,
            fontsize=14, fontweight='bold')
    
    # 添加说明
    note_text = (
        f"Total input records: {data['total_input_records']:,}\n"
        f"liftOver mapped: {data['gold_standard_mapped']:,} "
        f"({data['gold_standard_mapped']/data['total_input_records']*100:.2f}%)\n\n"
        "Identity Rate: Percentage of records with identical coordinates to liftOver\n"
        "Partial Match: Records with different number of output regions (split)\n"
        "Coord. Mismatch: Records with same count but different coordinates"
    )
    
    # ax.text(0.5, 0.05, note_text, ha='center', va='bottom', 
    #         transform=ax.transAxes, fontsize=8, style='italic',
    #         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))


def create_comparison_chart(data, ax):
    """
    创建准确性对比图（条形图）
    
    参数:
        data: 准确性分析数据
        ax: matplotlib axes
    """
    if not data:
        ax.text(0.5, 0.5, 'No accuracy data', ha='center', va='center', 
                transform=ax.transAxes, fontsize=12)
        ax.set_title('Accuracy Metrics Comparison', fontsize=11, fontweight='bold')
        return
    
    results = {r["tool"]: r for r in data["results"]}
    
    # 准备数据
    tools = []
    identity_rates = []
    mapping_rates = []
    
    for tool in TOOL_ORDER:
        if tool in results and results[tool]["success"]:
            tools.append(tool)
            identity_rates.append(results[tool]["identity_rate"] * 100)
            mapping_rates.append(results[tool]["mapping_rate"] * 100)
    
    if not tools:
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center', 
                transform=ax.transAxes, fontsize=10)
        ax.set_title('Accuracy Metrics Comparison', fontsize=11, fontweight='bold')
        return
    
    # 绘制分组条形图
    x = range(len(tools))
    width = 0.35
    
    bars1 = ax.bar([i - width/2 for i in x], mapping_rates, width, 
                   label='Mapping Rate', color='#4CAF50', alpha=0.7)
    bars2 = ax.bar([i + width/2 for i in x], identity_rates, width, 
                   label='Identity Rate', color='#2196F3', alpha=0.7)
    
    # 在条形图上方添加数值标签
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.2f}%',
                    ha='center', va='bottom', fontsize=8)
    
    ax.set_ylabel('Percentage (%)', fontsize=10)
    ax.set_title('Accuracy Metrics Comparison', fontsize=11, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(tools, fontsize=9)
    ax.set_ylim([95, 101])  # 聚焦在高准确率区域
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')


def main():
    print("=" * 60)
    print("生成准确性对比图 (Figure 3)")
    print("=" * 60)
    
    # 加载数据
    accuracy_data = load_accuracy_data()
    
    if not accuracy_data:
        print("错误: 没有找到准确性分析数据")
        print("请先运行: python paper/07_accuracy_analysis.py")
        return
    
    print(f"输入文件: {accuracy_data['input_file']}")
    print(f"总记录数: {accuracy_data['total_input_records']:,}")
    print(f"金标准: {accuracy_data['gold_standard']}")
    
    # 创建图表 (2行1列布局)
    fig = plt.figure(figsize=(12, 10))
    
    # 上半部分：准确性表格
    ax1 = plt.subplot2grid((2, 1), (0, 0))
    create_accuracy_table(accuracy_data, ax1)
    
    # 下半部分：准确性对比图
    ax2 = plt.subplot2grid((2, 1), (1, 0))
    create_comparison_chart(accuracy_data, ax2)
    
    fig.suptitle('Figure 3: Accuracy Analysis', 
                 fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # 保存组合图表
    output_pdf = FIGURES_DIR / "fig3_accuracy.pdf"
    output_png = FIGURES_DIR / "fig3_accuracy.png"
    
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(output_png, dpi=300, bbox_inches='tight')
    
    print(f"\n组合图表已保存到:")
    print(f"  {output_pdf}")
    print(f"  {output_png}")
    
    # 保存单独的子图
    print(f"\n保存单独子图...")
    
    # 上图: 准确性表格
    fig_table, ax_table = plt.subplots(figsize=(12, 6))
    create_accuracy_table(accuracy_data, ax_table)
    plt.tight_layout()
    fig_table.savefig(FIGURES_DIR / "fig3a_accuracy_table.pdf", dpi=300, bbox_inches='tight')
    fig_table.savefig(FIGURES_DIR / "fig3a_accuracy_table.png", dpi=300, bbox_inches='tight')
    plt.close(fig_table)
    print(f"  {FIGURES_DIR / 'fig3a_accuracy_table.pdf'}")
    
    # 下图: 准确性对比图
    fig_chart, ax_chart = plt.subplots(figsize=(8, 5))
    create_comparison_chart(accuracy_data, ax_chart)
    plt.tight_layout()
    fig_chart.savefig(FIGURES_DIR / "fig3b_accuracy_chart.pdf", dpi=300, bbox_inches='tight')
    fig_chart.savefig(FIGURES_DIR / "fig3b_accuracy_chart.png", dpi=300, bbox_inches='tight')
    plt.close(fig_chart)
    print(f"  {FIGURES_DIR / 'fig3b_accuracy_chart.pdf'}")
    
    # 打印准确性摘要
    print("\n" + "=" * 60)
    print("准确性摘要")
    print("=" * 60)
    
    results = {r["tool"]: r for r in accuracy_data["results"]}
    
    for tool in TOOL_ORDER:
        if tool in results and results[tool]["success"]:
            r = results[tool]
            print(f"\n{tool}:")
            print(f"  映射率: {r['mapping_rate']*100:.2f}%")
            print(f"  一致率: {r['identity_rate']*100:.2f}%")
            print(f"  部分匹配: {r.get('partial_match', 0)}")
            print(f"  坐标偏差: {r.get('coordinate_mismatch', 0)}")
    
    print("\n" + "=" * 60)
    print("Figure 3 设计说明:")
    print("=" * 60)
    print("上图: 准确性对比表格 - 展示各工具与 liftOver 的一致性")
    print("下图: 准确性指标对比 - 映射率和一致率的条形图对比")
    print("      聚焦在 95-100% 区域，突出高准确率")
    print("\n下一步: python paper/09_feature_audit.py")


if __name__ == "__main__":
    main()
