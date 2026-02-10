#!/usr/bin/env python3
"""
10_plot_features.py - Generate feature comparison heatmap (Figure 4)

Figure 4 design:
  Heatmap showing feature support for each tool
  - Rows: Features
  - Columns: Tools
  - Colors: Supported(green)/Not supported(red)

Design rationale:
- Visually show feature coverage differences
- Highlight FastCrossMap's comprehensiveness
- Easy to identify strengths and weaknesses of each tool

Usage: python paper/10_plot_features.py
Output: paper/figures/fig4_features.pdf
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
    """Load feature audit data"""
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
    
    # Define feature list (categories preserved for grouping)
    features = [
        # File format support
        ("BED", "format_bed", "File Formats"),
        ("BAM/SAM", "format_bam", "File Formats"),
        ("VCF", "format_vcf", "File Formats"),
        ("GFF/GTF", "format_gff", "File Formats"),
        ("Wiggle", "format_wiggle", "File Formats"),
        ("BigWig", "format_bigwig", "File Formats"),
        ("MAF", "format_maf", "File Formats"),
        ("GVCF", "format_gvcf", "File Formats"),
        
        # Compressed file support
        ("Compressed Chain", "compressed_chain", "Compression"),
        ("Compressed Input", "compressed_input", "Compression"),
        
        # Multi-threading support
        ("Multithreading", "multithreading", "Threading"),
        ("User Control Threads", "user_controllable_threads", "Threading"),
        
        # Cross-platform support
        ("Linux", "platform_linux", "Platforms"),
        ("macOS", "platform_macos", "Platforms"),
        ("Windows", "platform_windows", "Platforms"),
        
        # Other features
        ("Unmapped Output", "unmapped_output", "Other"),
        ("Streaming Process", "streaming_processing", "Other"),
    ]
    
    # Build matrix
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
    
    # Create heatmap
    cmap = ListedColormap(['#ffcccc', '#ccffcc'])  # Red(not supported), Green(supported)
    im = ax.imshow(matrix, cmap=cmap, aspect='auto', vmin=0, vmax=1)
    
    # Define category colors
    category_colors = {
        "File Formats": "#1f77b4",    # 蓝色
        "Compression": "#ff7f0e",     # 橙色
        "Threading": "#2ca02c",       # 绿色
        "Platforms": "#d62728",       # 红色
        "Other": "#9467bd"            # 紫色
    }
    
    # Assign category colors to each feature
    feature_colors = []
    for _, _, category in features:
        feature_colors.append(category_colors[category])
    
    # Set axes
    ax.set_xticks(np.arange(len(TOOL_ORDER)))
    ax.set_yticks(np.arange(len(feature_labels)))
    ax.set_xticklabels(TOOL_ORDER, fontsize=10, fontweight='bold')
    ax.set_yticklabels(feature_labels, fontsize=9)
    
    # Set Y-axis label colors
    for tick_label, color in zip(ax.get_yticklabels(), feature_colors):
        tick_label.set_color(color)
    
    # Rotate x-axis labels
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    
    # Add symbol in each cell
    for i in range(len(feature_labels)):
        for j in range(len(TOOL_ORDER)):
            text = "✓" if matrix[i, j] == 1 else "✗"
            color = "darkgreen" if matrix[i, j] == 1 else "darkred"
            ax.text(j, i, text, ha="center", va="center", 
                   color=color, fontsize=12, fontweight='bold')
    
    # Add category separator lines and labels
    category_positions = []
    category_labels = {}
    current_category = None
    
    for i, (_, _, category) in enumerate(features):
        if category != current_category:
            if current_category is not None:
                category_positions.append(i - 0.5)
            current_category = category
            category_labels[category] = i
    
    # Draw separator lines (black)
    for pos in category_positions:
        ax.axhline(y=pos, color='black', linewidth=2, linestyle='-')
    
    # Add category labels (on left side of Y-axis, aligned with each group)
    category_ranges = {
        "File Formats": (0, 7),
        "Compression": (8, 9),
        "Threading": (10, 11),
        "Platforms": (12, 14),
        "Other": (15, 16)
    }
    
    for category, (start, end) in category_ranges.items():
        # Calculate group center position
        mid_row = (start + end) / 2
        color = category_colors[category]
        
        # Display category label on left side of Y-axis (with color)
        ax.text(-1.3, mid_row, f"{category}", ha='left', va='center', 
               fontsize=9, fontweight='bold', 
               rotation=90, color=color)
    
    ax.set_title('Feature Support Matrix', fontsize=14, fontweight='bold', pad=20)
    
    # Remove borders
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    ax.set_xticks(np.arange(len(TOOL_ORDER)+1)-.5, minor=True)
    ax.set_yticks(np.arange(len(feature_labels)+1)-.5, minor=True)
    ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", size=0)


def main():
    print("=" * 60)
    print("Generating Feature Comparison Heatmap (Figure 4)")
    print("=" * 60)
    
    # Load data
    features_data = load_features_data()
    
    if not features_data:
        print("Error: No feature audit data found")
        print("Please run first: python paper/09_feature_audit.py")
        return
    
    print(f"Tools: {len(features_data['tools'])}")
    
    # Create figure (heatmap only)
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Feature heatmap
    create_feature_heatmap(features_data, ax)
    
    # fig.suptitle('Feature Support Comparison', 
    #              fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    # Save combined figure
    output_pdf = FIGURES_DIR / "fig4_features.pdf"
    output_png = FIGURES_DIR / "fig4_features.png"
    
    fig.savefig(output_pdf, dpi=300, bbox_inches='tight')
    fig.savefig(output_png, dpi=300, bbox_inches='tight')
    
    print(f"\nCombined figure saved to:")
    print(f"  {output_pdf}")
    print(f"  {output_png}")
    
    # Save individual subplots (only one heatmap, save once)
    print(f"\nSaving individual subplots...")
    fig_single, ax_single = plt.subplots(figsize=(10, 10))
    create_feature_heatmap(features_data, ax_single)
    plt.tight_layout()
    fig_single.savefig(FIGURES_DIR / "fig4_features_heatmap.pdf", dpi=300, bbox_inches='tight')
    fig_single.savefig(FIGURES_DIR / "fig4_features_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close(fig_single)
    print(f"  {FIGURES_DIR / 'fig4_features_heatmap.pdf'}")
    
    # Print feature summary
    print("\n" + "=" * 60)
    print("Feature Coverage Summary")
    print("=" * 60)
    
    tools_data = {t["tool"]: t for t in features_data["tools"]}
    
    for tool in TOOL_ORDER:
        if tool in tools_data:
            t = tools_data[tool]
            print(f"\n{tool}:")
            print(f"  Format support: {t['format_count']}/8")
            print(f"  Platform support: {t['platform_count']}/3")
            print(f"  Compressed Chain: {'✓' if t['compressed_chain'] else '✗'}")
            print(f"  Multi-threading: {'✓' if t['multithreading'] else '✗'}")
            print(f"  Controllable threads: {'✓' if t['user_controllable_threads'] else '✗'}")
            print(f"  Coverage: {t['feature_coverage_score']*100:.1f}%")
    
    print("\n" + "=" * 60)
    print("Figure 4 Design Notes:")
    print("=" * 60)
    print("Feature support matrix heatmap - visual overview of tool capabilities")
    print("  Green(✓) = supported, Red(✗) = not supported")
    print("  Grouped by category: File Formats, Compression, Threading, Platforms, Other")
    print("  FastCrossMap has the most comprehensive support (8/8 formats, 3/3 platforms)")
    print("\nNext step: python paper/11_plot_radar.py")


if __name__ == "__main__":
    main()
