#!/usr/bin/env python3
"""
09_feature_audit.py - 功能审计

检测各工具的功能支持情况，生成功能对比矩阵

检测项目：
1. 文件格式支持 (BED, BAM, VCF, GFF, etc.)
2. 压缩文件支持 (gzip chain, gzip input)
3. 多线程支持
4. 输出未映射记录
5. 命令行易用性

原理：
- 基于工具文档和实际测试结果
- 生成功能矩阵用于热力图可视化
- 计算 Feature Coverage Score

用法: python paper/09_feature_audit.py
输出: paper/results/features.json
"""

import json
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List

# =============================================================================
# 配置
# =============================================================================
RESULTS_DIR = Path("paper/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


@dataclass
class FeatureMatrix:
    """功能矩阵"""
    tool: str
    
    # 文件格式支持 (8 种)
    format_bed: bool
    format_bam: bool
    format_vcf: bool
    format_gff: bool
    format_wiggle: bool
    format_bigwig: bool
    format_maf: bool
    format_gvcf: bool
    
    # 压缩文件支持
    compressed_chain: bool      # 支持 .gz chain 文件
    compressed_input: bool      # 支持 .gz 输入文件
    
    # 多线程支持
    multithreading: bool        # 支持多线程
    user_controllable_threads: bool  # 用户可控制线程数
    
    # 跨平台支持
    platform_linux: bool        # Linux 支持
    platform_macos: bool        # macOS 支持
    platform_windows: bool      # Windows 支持
    
    # 其他功能
    unmapped_output: bool       # 输出未映射记录
    streaming_processing: bool  # 流式处理（低内存）
    cli_simplicity: int         # 命令行简洁度 (1-5, 5最简洁)
    
    # 统计
    format_count: int           # 支持的格式数
    platform_count: int         # 支持的平台数
    feature_coverage_score: float  # 功能覆盖率 (0-1)


def audit_fastcrossmap() -> FeatureMatrix:
    """审计 FastCrossMap 功能"""
    return FeatureMatrix(
        tool="FastCrossMap",
        
        # 文件格式支持 (8/8)
        format_bed=True,
        format_bam=True,
        format_vcf=True,
        format_gff=True,
        format_wiggle=True,
        format_bigwig=True,
        format_maf=True,
        format_gvcf=True,
        
        # 压缩文件支持
        compressed_chain=True,      # 支持 .gz chain
        compressed_input=True,      # 支持 .gz 输入 (通过 flate2)
        
        # 多线程支持
        multithreading=True,        # 支持多线程
        user_controllable_threads=True,  # -t 参数控制
        
        # 跨平台支持 (3/3)
        platform_linux=True,        # Rust 原生支持
        platform_macos=True,        # Rust 原生支持
        platform_windows=True,      # Rust 原生支持
        
        # 其他功能
        unmapped_output=True,       # 自动生成 .unmap 文件
        streaming_processing=True,  # 流式处理，低内存
        cli_simplicity=5,           # 命令简洁: fast-crossmap bed chain.gz input.bed output.bed
        
        # 统计
        format_count=8,
        platform_count=3,
        feature_coverage_score=0.0  # 稍后计算
    )


def audit_crossmap() -> FeatureMatrix:
    """审计 CrossMap 功能"""
    return FeatureMatrix(
        tool="CrossMap",
        
        # 文件格式支持 (8/8)
        format_bed=True,
        format_bam=True,
        format_vcf=True,
        format_gff=True,
        format_wiggle=True,
        format_bigwig=True,
        format_maf=True,
        format_gvcf=True,
        
        # 压缩文件支持
        compressed_chain=True,      # 支持 .gz chain
        compressed_input=True,      # 支持 .gz 输入
        
        # 多线程支持
        multithreading=False,       # 单线程
        user_controllable_threads=False,
        
        # 跨平台支持 (2.5/3)
        platform_linux=True,        # Python 原生支持
        platform_macos=True,        # Python 原生支持
        platform_windows=False,     # Windows 需要 WSL，不算完全支持
        
        # 其他功能
        unmapped_output=True,       # 生成 .unmap 文件
        streaming_processing=False, # Python 实现，内存占用较高
        cli_simplicity=4,           # 命令: CrossMap bed chain.gz input.bed output.bed
        
        # 统计
        format_count=8,
        platform_count=2,
        feature_coverage_score=0.0
    )


def audit_liftover() -> FeatureMatrix:
    """审计 liftOver 功能"""
    return FeatureMatrix(
        tool="liftOver",
        
        # 文件格式支持 (2/8) - 只支持 BED 和 GFF
        format_bed=True,
        format_bam=False,
        format_vcf=False,
        format_gff=True,
        format_wiggle=False,
        format_bigwig=False,
        format_maf=False,
        format_gvcf=False,
        
        # 压缩文件支持
        compressed_chain=True,      # 支持 .gz chain
        compressed_input=False,     # 不支持 .gz 输入
        
        # 多线程支持
        multithreading=False,       # 单线程
        user_controllable_threads=False,
        
        # 跨平台支持 (2/3)
        platform_linux=True,        # C 程序，Linux 原生
        platform_macos=True,        # macOS 原生
        platform_windows=False,     # Windows 不支持
        
        # 其他功能
        unmapped_output=True,       # 需要指定 unmap 文件
        streaming_processing=True,  # C 实现，流式处理
        cli_simplicity=3,           # 命令: liftOver input.bed chain.gz output.bed unmap.bed
        
        # 统计
        format_count=2,
        platform_count=2,
        feature_coverage_score=0.0
    )


def audit_fastremap() -> FeatureMatrix:
    """审计 FastRemap 功能"""
    return FeatureMatrix(
        tool="FastRemap",
        
        # 文件格式支持 (2/8) - 只支持 BED 和 BAM
        format_bed=True,
        format_bam=True,
        format_vcf=False,
        format_gff=False,
        format_wiggle=False,
        format_bigwig=False,
        format_maf=False,
        format_gvcf=False,
        
        # 压缩文件支持
        compressed_chain=False,     # 不支持 .gz chain (致命缺陷)
        compressed_input=False,     # 不支持 .gz 输入
        
        # 多线程支持
        multithreading=True,        # 内部使用 4 线程
        user_controllable_threads=False,  # 用户无法控制线程数
        
        # 跨平台支持 (2/3)
        platform_linux=True,        # C++ SeqAn2
        platform_macos=True,        # macOS 支持
        platform_windows=False,     # Windows 编译困难
        
        # 其他功能
        unmapped_output=True,       # -u 参数指定
        streaming_processing=True,  # C++ 实现
        cli_simplicity=2,           # 命令复杂: FastRemap -f bed -c chain -i input -o output -u unmap
        
        # 统计
        format_count=2,
        platform_count=2,
        feature_coverage_score=0.0
    )


def calculate_feature_score(matrix: FeatureMatrix) -> float:
    """
    计算功能覆盖率分数 (0-1)
    
    权重分配:
    - 文件格式支持: 35% (8 种格式)
    - 压缩文件支持: 15% (2 项)
    - 多线程支持: 15% (2 项)
    - 跨平台支持: 15% (3 个平台)
    - 其他功能: 20% (3 项)
    """
    # 文件格式支持 (8 种，每种 4.375%)
    format_score = sum([
        matrix.format_bed,
        matrix.format_bam,
        matrix.format_vcf,
        matrix.format_gff,
        matrix.format_wiggle,
        matrix.format_bigwig,
        matrix.format_maf,
        matrix.format_gvcf
    ]) / 8 * 0.35
    
    # 压缩文件支持 (2 项，每项 7.5%)
    compression_score = sum([
        matrix.compressed_chain,
        matrix.compressed_input
    ]) / 2 * 0.15
    
    # 多线程支持 (2 项，每项 7.5%)
    threading_score = sum([
        matrix.multithreading,
        matrix.user_controllable_threads
    ]) / 2 * 0.15
    
    # 跨平台支持 (3 个平台，每个 5%)
    platform_score = sum([
        matrix.platform_linux,
        matrix.platform_macos,
        matrix.platform_windows
    ]) / 3 * 0.15
    
    # 其他功能 (3 项)
    other_score = (
        (1 if matrix.unmapped_output else 0) * 0.07 +
        (1 if matrix.streaming_processing else 0) * 0.07 +
        (matrix.cli_simplicity / 5) * 0.06
    )
    
    total_score = format_score + compression_score + threading_score + platform_score + other_score
    return round(total_score, 3)


def main():
    print("=" * 60)
    print("功能审计")
    print("=" * 60)
    
    # 审计各工具
    tools = [
        audit_fastcrossmap(),
        audit_crossmap(),
        audit_liftover(),
        audit_fastremap()
    ]
    
    # 计算功能覆盖率分数
    for tool in tools:
        tool.feature_coverage_score = calculate_feature_score(tool)
    
    # 打印摘要
    print("\n" + "=" * 60)
    print("功能审计摘要")
    print("=" * 60)
    print(f"{'工具':<15} {'格式':<8} {'压缩':<6} {'多线程':<8} {'平台':<8} {'覆盖率':<10}")
    print("-" * 60)
    
    for tool in tools:
        compressed_chain = "✓" if tool.compressed_chain else "✗"
        multithreading = "✓" if tool.multithreading else "✗"
        print(f"{tool.tool:<15} {tool.format_count}/8{'':<4} {compressed_chain:<6} "
              f"{multithreading:<8} {tool.platform_count}/3{'':<4} {tool.feature_coverage_score*100:.1f}%")
    
    # 详细功能矩阵
    print("\n" + "=" * 60)
    print("详细功能矩阵")
    print("=" * 60)
    
    features = [
        ("BED", "format_bed"),
        ("BAM/SAM", "format_bam"),
        ("VCF", "format_vcf"),
        ("GFF/GTF", "format_gff"),
        ("Wiggle", "format_wiggle"),
        ("BigWig", "format_bigwig"),
        ("MAF", "format_maf"),
        ("GVCF", "format_gvcf"),
        ("压缩 Chain", "compressed_chain"),
        ("压缩输入", "compressed_input"),
        ("多线程", "multithreading"),
        ("可控线程数", "user_controllable_threads"),
        ("Linux 支持", "platform_linux"),
        ("macOS 支持", "platform_macos"),
        ("Windows 支持", "platform_windows"),
        ("未映射输出", "unmapped_output"),
        ("流式处理", "streaming_processing"),
    ]
    
    # 打印表头
    print(f"{'功能':<20}", end="")
    for tool in tools:
        print(f"{tool.tool:<15}", end="")
    print()
    print("-" * 80)
    
    # 打印每个功能
    for feature_name, feature_attr in features:
        print(f"{feature_name:<20}", end="")
        for tool in tools:
            value = getattr(tool, feature_attr)
            symbol = "✓" if value else "✗"
            print(f"{symbol:<15}", end="")
        print()
    
    # 保存结果
    output_json = RESULTS_DIR / "features.json"
    with open(output_json, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "tools": [asdict(tool) for tool in tools],
            "feature_categories": {
                "file_formats": ["BED", "BAM", "VCF", "GFF", "Wiggle", "BigWig", "MAF", "GVCF"],
                "compression": ["Compressed Chain", "Compressed Input"],
                "threading": ["Multithreading", "User Controllable Threads"],
                "platforms": ["Linux", "macOS", "Windows"],
                "other": ["Unmapped Output", "Streaming Processing", "CLI Simplicity"]
            }
        }, f, indent=2)
    
    print(f"\n结果已保存到: {output_json}")
    print("\n下一步: python paper/10_plot_features.py")


if __name__ == "__main__":
    main()
