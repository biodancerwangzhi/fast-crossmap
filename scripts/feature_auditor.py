#!/usr/bin/env python3
"""功能审计脚本 - 对比各工具的功能支持情况"""

import argparse
import json
import subprocess
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List

FORMATS = ["bed", "vcf", "gff", "gtf", "bam", "sam", "wig", "bigwig", "maf", "gvcf"]
INPUT_FEATURES = ["gz_compression", "bz2_compression", "url_input", "stdin_input"]
ADVANCED_FEATURES = ["chain_gz", "multi_threading", "memory_efficient", "batch_processing"]
USABILITY_FEATURES = ["progress_bar", "verbose_output", "error_recovery", "unmapped_output"]

@dataclass
class FeatureSupport:
    supported: bool
    partial: bool = False
    notes: str = ""

@dataclass
class ToolFeatures:
    tool: str
    version: str = ""
    formats: Dict[str, FeatureSupport] = field(default_factory=dict)
    input_features: Dict[str, FeatureSupport] = field(default_factory=dict)
    advanced_features: Dict[str, FeatureSupport] = field(default_factory=dict)
    usability_features: Dict[str, FeatureSupport] = field(default_factory=dict)
    
    def get_format_count(self) -> int:
        return sum(1 for f in self.formats.values() if f.supported)
    
    def get_total_features(self) -> int:
        count = 0
        for features in [self.formats, self.input_features, self.advanced_features, self.usability_features]:
            count += sum(1 for f in features.values() if f.supported)
        return count
    
    def get_coverage_score(self) -> float:
        total = len(FORMATS) + len(INPUT_FEATURES) + len(ADVANCED_FEATURES) + len(USABILITY_FEATURES)
        return round(self.get_total_features() / total * 100, 1) if total > 0 else 0

@dataclass
class FeatureAuditReport:
    timestamp: str
    tools: List[ToolFeatures] = field(default_factory=list)


# 预定义功能矩阵
KNOWN_FEATURES = {
    "fastcrossmap": {
        "formats": {"bed": True, "vcf": True, "gff": True, "gtf": True, "bam": True, 
                   "sam": True, "wig": True, "bigwig": True, "maf": True, "gvcf": True},
        "input": {"gz_compression": True, "bz2_compression": False, "url_input": False, "stdin_input": False},
        "advanced": {"chain_gz": True, "multi_threading": True, "memory_efficient": True, "batch_processing": True},
        "usability": {"progress_bar": False, "verbose_output": True, "error_recovery": True, "unmapped_output": True},
    },
    "crossmap": {
        "formats": {"bed": True, "vcf": True, "gff": True, "gtf": True, "bam": True,
                   "sam": True, "wig": True, "bigwig": True, "maf": True, "gvcf": True},
        "input": {"gz_compression": True, "bz2_compression": True, "url_input": True, "stdin_input": False},
        "advanced": {"chain_gz": True, "multi_threading": False, "memory_efficient": False, "batch_processing": False},
        "usability": {"progress_bar": False, "verbose_output": True, "error_recovery": True, "unmapped_output": True},
    },
    "liftover": {
        "formats": {"bed": True, "vcf": False, "gff": False, "gtf": False, "bam": False,
                   "sam": False, "wig": False, "bigwig": False, "maf": False, "gvcf": False},
        "input": {"gz_compression": True, "bz2_compression": False, "url_input": False, "stdin_input": False},
        "advanced": {"chain_gz": True, "multi_threading": False, "memory_efficient": True, "batch_processing": False},
        "usability": {"progress_bar": False, "verbose_output": False, "error_recovery": False, "unmapped_output": True},
    },
    "fastremap": {
        "formats": {"bed": True, "vcf": False, "gff": False, "gtf": False, "bam": True,
                   "sam": True, "wig": False, "bigwig": False, "maf": False, "gvcf": False},
        "input": {"gz_compression": False, "bz2_compression": False, "url_input": False, "stdin_input": False},
        "advanced": {"chain_gz": False, "multi_threading": False, "memory_efficient": True, "batch_processing": False},
        "usability": {"progress_bar": False, "verbose_output": False, "error_recovery": False, "unmapped_output": True},
    },
}


class FeatureAuditor:
    def __init__(self, output_dir: Path = Path("results")):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def get_tool_version(self, tool: str) -> str:
        try:
            if tool == "fastcrossmap":
                r = subprocess.run(["./target/release/fast-crossmap", "--version"], capture_output=True, text=True, timeout=10)
                return r.stdout.strip() or "unknown"
            elif tool == "crossmap":
                r = subprocess.run(["CrossMap", "--version"], capture_output=True, text=True, timeout=10)
                return r.stdout.strip() or r.stderr.strip() or "unknown"
            elif tool == "liftover":
                return "UCSC liftOver"
            elif tool == "fastremap":
                return "FastRemap 1.0"
        except Exception:
            pass
        return "unknown"
    
    def audit_tool(self, tool: str) -> ToolFeatures:
        features = ToolFeatures(tool=tool, version=self.get_tool_version(tool))
        known = KNOWN_FEATURES.get(tool, {})
        
        for fmt in FORMATS:
            features.formats[fmt] = FeatureSupport(supported=known.get("formats", {}).get(fmt, False))
        for feat in INPUT_FEATURES:
            features.input_features[feat] = FeatureSupport(supported=known.get("input", {}).get(feat, False))
        for feat in ADVANCED_FEATURES:
            features.advanced_features[feat] = FeatureSupport(supported=known.get("advanced", {}).get(feat, False))
        for feat in USABILITY_FEATURES:
            features.usability_features[feat] = FeatureSupport(supported=known.get("usability", {}).get(feat, False))
        return features
    
    def audit_all_tools(self) -> FeatureAuditReport:
        report = FeatureAuditReport(timestamp=datetime.now().isoformat())
        for tool in ["fastcrossmap", "crossmap", "liftover", "fastremap"]:
            report.tools.append(self.audit_tool(tool))
        return report


def generate_markdown(report: FeatureAuditReport) -> str:
    lines = ["# 功能对比表 (Feature Comparison)", "", f"生成时间: {report.timestamp}", ""]
    
    # 版本信息
    lines.extend(["## 工具版本", "", "| 工具 | 版本 | 格式支持 | 功能覆盖率 |", "|------|------|----------|------------|"])
    for t in report.tools:
        lines.append(f"| {t.tool} | {t.version} | {t.get_format_count()}/10 | {t.get_coverage_score()}% |")
    lines.append("")
    
    # 格式支持
    lines.extend(["## 格式支持", "", "| 格式 | FastCrossMap | CrossMap | liftOver | FastRemap |", "|------|--------------|----------|----------|-----------|"])
    for fmt in FORMATS:
        row = [fmt]
        for t in report.tools:
            s = t.formats.get(fmt, FeatureSupport(False))
            row.append("✓" if s.supported else "✗")
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")
    
    # 输入功能
    feat_names = {"gz_compression": "GZ压缩", "bz2_compression": "BZ2压缩", "url_input": "URL输入", "stdin_input": "标准输入"}
    lines.extend(["## 输入功能", "", "| 功能 | FastCrossMap | CrossMap | liftOver | FastRemap |", "|------|--------------|----------|----------|-----------|"])
    for feat in INPUT_FEATURES:
        row = [feat_names.get(feat, feat)]
        for t in report.tools:
            s = t.input_features.get(feat, FeatureSupport(False))
            row.append("✓" if s.supported else "✗")
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")
    
    # 高级功能
    adv_names = {"chain_gz": "Chain GZ支持", "multi_threading": "多线程", "memory_efficient": "内存高效", "batch_processing": "批量处理"}
    lines.extend(["## 高级功能", "", "| 功能 | FastCrossMap | CrossMap | liftOver | FastRemap |", "|------|--------------|----------|----------|-----------|"])
    for feat in ADVANCED_FEATURES:
        row = [adv_names.get(feat, feat)]
        for t in report.tools:
            s = t.advanced_features.get(feat, FeatureSupport(False))
            row.append("✓" if s.supported else "✗")
        lines.append("| " + " | ".join(row) + " |")
    lines.append("")
    
    # 可用性功能
    use_names = {"progress_bar": "进度条", "verbose_output": "详细输出", "error_recovery": "错误恢复", "unmapped_output": "未映射输出"}
    lines.extend(["## 可用性功能", "", "| 功能 | FastCrossMap | CrossMap | liftOver | FastRemap |", "|------|--------------|----------|----------|-----------|"])
    for feat in USABILITY_FEATURES:
        row = [use_names.get(feat, feat)]
        for t in report.tools:
            s = t.usability_features.get(feat, FeatureSupport(False))
            row.append("✓" if s.supported else "✗")
        lines.append("| " + " | ".join(row) + " |")
    
    return "\n".join(lines)


def print_summary(report: FeatureAuditReport):
    print("\n" + "="*70)
    print("FEATURE AUDIT SUMMARY")
    print("="*70)
    print(f"{'Tool':<15} {'Formats':<12} {'Features':<12} {'Coverage':<12}")
    print("-"*70)
    for t in report.tools:
        print(f"{t.tool:<15} {t.get_format_count()}/10{'':<6} {t.get_total_features()}/22{'':<5} {t.get_coverage_score()}%")
    print("="*70)

def main():
    parser = argparse.ArgumentParser(description="Feature audit for coordinate conversion tools")
    parser.add_argument("--output", "-o", type=Path, default=Path("results/feature_comparison.md"))
    parser.add_argument("--json", action="store_true", help="Also output JSON")
    args = parser.parse_args()
    
    auditor = FeatureAuditor()
    report = auditor.audit_all_tools()
    
    # 生成 Markdown
    md_content = generate_markdown(report)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(md_content)
    print(f"Markdown report saved to: {args.output}")
    
    # 可选 JSON 输出
    if args.json:
        json_path = args.output.with_suffix('.json')
        data = {"timestamp": report.timestamp, "tools": []}
        for t in report.tools:
            data["tools"].append({
                "tool": t.tool, "version": t.version,
                "format_count": t.get_format_count(),
                "total_features": t.get_total_features(),
                "coverage_score": t.get_coverage_score()
            })
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"JSON report saved to: {json_path}")
    
    print_summary(report)

if __name__ == "__main__":
    main()
