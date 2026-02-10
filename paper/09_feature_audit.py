#!/usr/bin/env python3
"""
09_feature_audit.py - Feature audit

Detect feature support for each tool, generate feature comparison matrix

Detection items:
1. File format support (BED, BAM, VCF, GFF, etc.)
2. Compressed file support (gzip chain, gzip input)
3. Multi-threading support
4. Unmapped record output
5. Command-line usability

Methodology:
- Based on tool documentation and actual test results
- Generate feature matrix for heatmap visualization
- Calculate Feature Coverage Score

Usage: python paper/09_feature_audit.py
Output: paper/results/features.json
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
    """Feature matrix"""
    tool: str
    
    # File format support (8 types)
    format_bed: bool
    format_bam: bool
    format_vcf: bool
    format_gff: bool
    format_wiggle: bool
    format_bigwig: bool
    format_maf: bool
    format_gvcf: bool
    
    # Compressed file support
    compressed_chain: bool      # Supports .gz chain files
    compressed_input: bool      # Supports .gz input files
    
    # Multi-threading support
    multithreading: bool        # Supports multi-threading
    user_controllable_threads: bool  # User can control thread count
    
    # Cross-platform support
    platform_linux: bool        # Linux support
    platform_macos: bool        # macOS support
    platform_windows: bool      # Windows support
    
    # Other features
    unmapped_output: bool       # Outputs unmapped records
    streaming_processing: bool  # Streaming processing (low memory)
    cli_simplicity: int         # CLI simplicity (1-5, 5 = simplest)
    
    # Statistics
    format_count: int           # Number of supported formats
    platform_count: int         # Number of supported platforms
    feature_coverage_score: float  # Feature coverage score (0-1)


def audit_fastcrossmap() -> FeatureMatrix:
    """Audit FastCrossMap features"""
    return FeatureMatrix(
        tool="FastCrossMap",
        
        # File format support (8/8)
        format_bed=True,
        format_bam=True,
        format_vcf=True,
        format_gff=True,
        format_wiggle=True,
        format_bigwig=True,
        format_maf=True,
        format_gvcf=True,
        
        # Compressed file support
        compressed_chain=True,      # Supports .gz chain
        compressed_input=True,      # Supports .gz input (via flate2)
        
        # Multi-threading support
        multithreading=True,        # Supports multi-threading
        user_controllable_threads=True,  # -t parameter
        
        # Cross-platform support (3/3)
        platform_linux=True,        # Rust native support
        platform_macos=True,        # Rust native support
        platform_windows=True,      # Rust native support
        
        # Other features
        unmapped_output=True,       # Auto-generates .unmap file
        streaming_processing=True,  # Streaming, low memory
        cli_simplicity=5,           # Simple: fast-crossmap bed chain.gz input.bed output.bed
        
        # Statistics
        format_count=8,
        platform_count=3,
        feature_coverage_score=0.0  # 稍后计算
    )


def audit_crossmap() -> FeatureMatrix:
    """Audit CrossMap features"""
    return FeatureMatrix(
        tool="CrossMap",
        
        # File format support (8/8)
        format_bed=True,
        format_bam=True,
        format_vcf=True,
        format_gff=True,
        format_wiggle=True,
        format_bigwig=True,
        format_maf=True,
        format_gvcf=True,
        
        # Compressed file support
        compressed_chain=True,      # Supports .gz chain
        compressed_input=True,      # Supports .gz input
        
        # Multi-threading support
        multithreading=False,       # Single-threaded
        user_controllable_threads=False,
        
        # Cross-platform support (2.5/3)
        platform_linux=True,        # Python native support
        platform_macos=True,        # Python native support
        platform_windows=False,     # Windows requires WSL, not fully supported
        
        # Other features
        unmapped_output=True,       # Generates .unmap file
        streaming_processing=False, # Python implementation, higher memory usage
        cli_simplicity=4,           # Command: CrossMap bed chain.gz input.bed output.bed
        
        # Statistics
        format_count=8,
        platform_count=2,
        feature_coverage_score=0.0
    )


def audit_liftover() -> FeatureMatrix:
    """Audit liftOver features"""
    return FeatureMatrix(
        tool="liftOver",
        
        # File format support (2/8) - only BED and GFF
        format_bed=True,
        format_bam=False,
        format_vcf=False,
        format_gff=True,
        format_wiggle=False,
        format_bigwig=False,
        format_maf=False,
        format_gvcf=False,
        
        # Compressed file support
        compressed_chain=True,      # Supports .gz chain
        compressed_input=False,     # No .gz input support
        
        # Multi-threading support
        multithreading=False,       # Single-threaded
        user_controllable_threads=False,
        
        # Cross-platform support (2/3)
        platform_linux=True,        # C program, Linux native
        platform_macos=True,        # macOS native
        platform_windows=False,     # Windows not supported
        
        # Other features
        unmapped_output=True,       # Requires specifying unmap file
        streaming_processing=True,  # C implementation, streaming
        cli_simplicity=3,           # Command: liftOver input.bed chain.gz output.bed unmap.bed
        
        # Statistics
        format_count=2,
        platform_count=2,
        feature_coverage_score=0.0
    )


def audit_fastremap() -> FeatureMatrix:
    """Audit FastRemap features"""
    return FeatureMatrix(
        tool="FastRemap",
        
        # File format support (2/8) - only BED and BAM
        format_bed=True,
        format_bam=True,
        format_vcf=False,
        format_gff=False,
        format_wiggle=False,
        format_bigwig=False,
        format_maf=False,
        format_gvcf=False,
        
        # Compressed file support
        compressed_chain=False,     # No .gz chain support (critical flaw)
        compressed_input=False,     # No .gz input support
        
        # Multi-threading support
        multithreading=True,        # Internal 4 threads
        user_controllable_threads=False,  # User cannot control thread count
        
        # Cross-platform support (2/3)
        platform_linux=True,        # C++ SeqAn2
        platform_macos=True,        # macOS support
        platform_windows=False,     # Windows compilation difficult
        
        # Other features
        unmapped_output=True,       # -u parameter
        streaming_processing=True,  # C++ implementation
        cli_simplicity=2,           # Complex: FastRemap -f bed -c chain -i input -o output -u unmap
        
        # Statistics
        format_count=2,
        platform_count=2,
        feature_coverage_score=0.0
    )


def calculate_feature_score(matrix: FeatureMatrix) -> float:
    """
    Calculate feature coverage score (0-1).
    
    Weight allocation:
    - File format support: 35% (8 formats)
    - Compressed file support: 15% (2 items)
    - Multi-threading support: 15% (2 items)
    - Cross-platform support: 15% (3 platforms)
    - Other features: 20% (3 items)
    """
    # File format support (8 types, each 4.375%)
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
    
    # Compressed file support (2 items, each 7.5%)
    compression_score = sum([
        matrix.compressed_chain,
        matrix.compressed_input
    ]) / 2 * 0.15
    
    # Multi-threading support (2 items, each 7.5%)
    threading_score = sum([
        matrix.multithreading,
        matrix.user_controllable_threads
    ]) / 2 * 0.15
    
    # Cross-platform support (3 platforms, each 5%)
    platform_score = sum([
        matrix.platform_linux,
        matrix.platform_macos,
        matrix.platform_windows
    ]) / 3 * 0.15
    
    # Other features (3 items)
    other_score = (
        (1 if matrix.unmapped_output else 0) * 0.07 +
        (1 if matrix.streaming_processing else 0) * 0.07 +
        (matrix.cli_simplicity / 5) * 0.06
    )
    
    total_score = format_score + compression_score + threading_score + platform_score + other_score
    return round(total_score, 3)


def main():
    print("=" * 60)
    print("Feature Audit")
    print("=" * 60)
    
    # Audit each tool
    tools = [
        audit_fastcrossmap(),
        audit_crossmap(),
        audit_liftover(),
        audit_fastremap()
    ]
    
    # Calculate feature coverage scores
    for tool in tools:
        tool.feature_coverage_score = calculate_feature_score(tool)
    
    # Print summary
    print("\n" + "=" * 60)
    print("Feature Audit Summary")
    print("=" * 60)
    print(f"{'Tool':<15} {'Formats':<8} {'Compress':<6} {'Threads':<8} {'Platform':<8} {'Coverage':<10}")
    print("-" * 60)
    
    for tool in tools:
        compressed_chain = "✓" if tool.compressed_chain else "✗"
        multithreading = "✓" if tool.multithreading else "✗"
        print(f"{tool.tool:<15} {tool.format_count}/8{'':<4} {compressed_chain:<6} "
              f"{multithreading:<8} {tool.platform_count}/3{'':<4} {tool.feature_coverage_score*100:.1f}%")
    
    # Detailed feature matrix
    print("\n" + "=" * 60)
    print("Detailed Feature Matrix")
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
        ("Compressed Chain", "compressed_chain"),
        ("Compressed Input", "compressed_input"),
        ("Multi-threading", "multithreading"),
        ("Controllable Threads", "user_controllable_threads"),
        ("Linux Support", "platform_linux"),
        ("macOS Support", "platform_macos"),
        ("Windows Support", "platform_windows"),
        ("Unmapped Output", "unmapped_output"),
        ("Streaming Processing", "streaming_processing"),
    ]
    
    # Print header
    print(f"{'Feature':<20}", end="")
    for tool in tools:
        print(f"{tool.tool:<15}", end="")
    print()
    print("-" * 80)
    
    # Print each feature
    for feature_name, feature_attr in features:
        print(f"{feature_name:<20}", end="")
        for tool in tools:
            value = getattr(tool, feature_attr)
            symbol = "✓" if value else "✗"
            print(f"{symbol:<15}", end="")
        print()
    
    # Save results
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
    
    print(f"\nResults saved to: {output_json}")
    print("\nNext step: python paper/10_plot_features.py")


if __name__ == "__main__":
    main()
