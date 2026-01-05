#!/usr/bin/env python3
"""
准确性分析脚本

比较 FastCrossMap 与参考工具的输出，计算 Identity Rate 和坐标漂移。

参考标准：
- BED 格式：以 UCSC liftOver 为"真值"
- VCF/BAM 格式：以 CrossMap (Python) 为"行为标准"

Usage:
    python scripts/accuracy_analyzer.py --input test_data/benchmark_test.bed --format bed
    python scripts/accuracy_analyzer.py --input test_data/test.vcf --format vcf
    python scripts/accuracy_analyzer.py --results-dir results/  # 使用已有结果

Features:
    - Identity Rate 计算（映射记录和未映射记录分别计算）
    - 坐标漂移距离统计（Coordinate Drift）
    - 差异报告生成
    - 低于 99% 一致率的警告
"""

import argparse
import csv
import json
import os
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set


@dataclass
class BedRecord:
    """BED 记录"""
    chrom: str
    start: int
    end: int
    name: str = ""
    score: str = ""
    strand: str = ""
    extra: List[str] = field(default_factory=list)
    
    @classmethod
    def parse(cls, line: str) -> Optional['BedRecord']:
        """解析 BED 行"""
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('track') or line.startswith('browser'):
            return None
        
        fields = line.split('\t')
        if len(fields) < 3:
            return None
        
        try:
            return cls(
                chrom=fields[0],
                start=int(fields[1]),
                end=int(fields[2]),
                name=fields[3] if len(fields) > 3 else "",
                score=fields[4] if len(fields) > 4 else "",
                strand=fields[5] if len(fields) > 5 else "",
                extra=fields[6:] if len(fields) > 6 else []
            )
        except (ValueError, IndexError):
            return None
    
    def key(self) -> str:
        """生成用于比较的键（基于 name 或坐标）"""
        if self.name:
            return self.name
        return f"{self.chrom}:{self.start}-{self.end}"
    
    def coords_key(self) -> str:
        """坐标键"""
        return f"{self.chrom}:{self.start}-{self.end}"


@dataclass
class Discrepancy:
    """差异记录"""
    record_key: str
    field: str
    reference_value: str
    test_value: str
    drift_distance: int = 0  # 坐标漂移距离


@dataclass
class AccuracyResult:
    """准确性分析结果"""
    test_tool: str
    reference_tool: str
    format: str
    total_records: int = 0
    matched_records: int = 0
    mismatched_records: int = 0
    
    # 分离的 Identity Rate
    mapped_total: int = 0
    mapped_matched: int = 0
    unmapped_total: int = 0
    unmapped_matched: int = 0
    
    # Identity Rate
    identity_rate: float = 0.0
    mapped_identity_rate: float = 0.0
    unmapped_identity_rate: float = 0.0
    
    # 坐标漂移统计
    drift_distances: List[int] = field(default_factory=list)
    drift_zero: int = 0      # 完全一致
    drift_one: int = 0       # 坐标系差异 (off-by-one)
    drift_small: int = 0     # 小漂移 (2-100)
    drift_large: int = 0     # 大漂移 (>100)
    
    # 差异详情
    discrepancies: List[Discrepancy] = field(default_factory=list)
    
    # 警告
    warnings: List[str] = field(default_factory=list)


@dataclass
class AccuracyReport:
    """准确性分析报告"""
    timestamp: str
    input_file: str
    format: str
    reference_tool: str
    results: List[AccuracyResult] = field(default_factory=list)


class AccuracyAnalyzer:
    """准确性分析器"""
    
    # 参考工具配置
    REFERENCE_TOOLS = {
        "bed": "liftover",      # BED 使用 liftOver 作为参考
        "vcf": "crossmap",      # VCF 使用 CrossMap 作为参考
        "bam": "crossmap",      # BAM 使用 CrossMap 作为参考
        "gff": "crossmap",      # GFF 使用 CrossMap 作为参考
    }
    
    def __init__(
        self,
        chain_file: Path,
        output_dir: Path = Path("results/accuracy"),
        threads: int = 4
    ):
        self.chain_file = chain_file
        self.output_dir = output_dir
        self.threads = threads
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 准备解压的 chain 文件
        self.chain_file_uncompressed = self._prepare_uncompressed_chain()
    
    def _prepare_uncompressed_chain(self) -> Path:
        """准备解压的 chain 文件"""
        if not str(self.chain_file).endswith('.gz'):
            return self.chain_file
        
        uncompressed = Path(str(self.chain_file)[:-3])
        if uncompressed.exists():
            return uncompressed
        
        try:
            import gzip
            import shutil
            with gzip.open(self.chain_file, 'rb') as f_in:
                with open(uncompressed, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            return uncompressed
        except Exception:
            return self.chain_file
    
    def run_tool(self, tool: str, input_file: Path, format: str) -> Tuple[Path, Path]:
        """运行工具并返回输出文件路径"""
        output_file = self.output_dir / f"{tool}_output.{format}"
        unmapped_file = self.output_dir / f"{tool}_output.{format}.unmap"
        
        if tool == "fastcrossmap":
            cmd = f"./target/release/fast-crossmap {format} {self.chain_file} {input_file} {output_file} -t {self.threads}"
        elif tool == "crossmap":
            cmd = f"CrossMap {format} {self.chain_file} {input_file} {output_file}"
        elif tool == "liftover":
            cmd = f"liftOver {input_file} {self.chain_file} {output_file} {unmapped_file}"
        elif tool == "fastremap":
            cmd = f"FastRemap -f {format} -c {self.chain_file_uncompressed} -i {input_file} -u {unmapped_file} -o {output_file}"
        else:
            raise ValueError(f"Unknown tool: {tool}")
        
        print(f"  Running {tool}...")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0 and tool != "liftover":
            print(f"    Warning: {tool} returned {result.returncode}")
        
        return output_file, unmapped_file
    
    def load_bed_records(self, file_path: Path) -> Dict[str, BedRecord]:
        """加载 BED 文件为字典"""
        records = {}
        if not file_path.exists():
            return records
        
        with open(file_path, 'r') as f:
            for line in f:
                rec = BedRecord.parse(line)
                if rec:
                    records[rec.key()] = rec
        
        return records
    
    @staticmethod
    def normalize_chrom(chrom: str) -> str:
        """标准化染色体名称（移除或添加 chr 前缀）"""
        if chrom.startswith('chr'):
            return chrom[3:]  # 移除 chr 前缀
        return chrom
    
    def compare_bed_records(
        self,
        test_records: Dict[str, BedRecord],
        ref_records: Dict[str, BedRecord],
        test_tool: str,
        ref_tool: str
    ) -> AccuracyResult:
        """比较两组 BED 记录"""
        result = AccuracyResult(
            test_tool=test_tool,
            reference_tool=ref_tool,
            format="bed"
        )
        
        # 找到共同的记录
        common_keys = set(test_records.keys()) & set(ref_records.keys())
        result.total_records = len(common_keys)
        
        for key in common_keys:
            test_rec = test_records[key]
            ref_rec = ref_records[key]
            
            # 比较染色体（标准化后比较，忽略 chr 前缀差异）
            test_chrom_norm = self.normalize_chrom(test_rec.chrom)
            ref_chrom_norm = self.normalize_chrom(ref_rec.chrom)
            
            if test_chrom_norm != ref_chrom_norm:
                result.mismatched_records += 1
                result.discrepancies.append(Discrepancy(
                    record_key=key,
                    field="chrom",
                    reference_value=ref_rec.chrom,
                    test_value=test_rec.chrom
                ))
                continue
            
            # 计算坐标漂移
            start_drift = abs(test_rec.start - ref_rec.start)
            end_drift = abs(test_rec.end - ref_rec.end)
            max_drift = max(start_drift, end_drift)
            
            result.drift_distances.append(max_drift)
            
            if max_drift == 0:
                result.drift_zero += 1
                result.matched_records += 1
            elif max_drift == 1:
                result.drift_one += 1
                result.matched_records += 1  # off-by-one 也算匹配
            elif max_drift <= 100:
                result.drift_small += 1
                result.mismatched_records += 1
                result.discrepancies.append(Discrepancy(
                    record_key=key,
                    field="coordinates",
                    reference_value=f"{ref_rec.start}-{ref_rec.end}",
                    test_value=f"{test_rec.start}-{test_rec.end}",
                    drift_distance=max_drift
                ))
            else:
                result.drift_large += 1
                result.mismatched_records += 1
                result.discrepancies.append(Discrepancy(
                    record_key=key,
                    field="coordinates",
                    reference_value=f"{ref_rec.start}-{ref_rec.end}",
                    test_value=f"{test_rec.start}-{test_rec.end}",
                    drift_distance=max_drift
                ))
        
        # 计算 Identity Rate
        if result.total_records > 0:
            result.identity_rate = round(
                result.matched_records / result.total_records * 100, 2
            )
        
        # 添加警告
        if result.identity_rate < 99.0:
            result.warnings.append(
                f"Identity Rate ({result.identity_rate}%) is below 99% threshold"
            )
        
        return result


    def compare_unmapped_records(
        self,
        test_unmapped: Set[str],
        ref_unmapped: Set[str],
        result: AccuracyResult
    ):
        """比较未映射记录（分离 Identity Rate）"""
        # 计算未映射记录的一致性
        common_unmapped = test_unmapped & ref_unmapped
        only_test_unmapped = test_unmapped - ref_unmapped
        only_ref_unmapped = ref_unmapped - test_unmapped
        
        result.unmapped_total = len(test_unmapped | ref_unmapped)
        result.unmapped_matched = len(common_unmapped)
        
        if result.unmapped_total > 0:
            result.unmapped_identity_rate = round(
                result.unmapped_matched / result.unmapped_total * 100, 2
            )
        
        # 记录差异
        for key in only_test_unmapped:
            result.discrepancies.append(Discrepancy(
                record_key=key,
                field="unmapped_status",
                reference_value="mapped",
                test_value="unmapped"
            ))
        
        for key in only_ref_unmapped:
            result.discrepancies.append(Discrepancy(
                record_key=key,
                field="unmapped_status",
                reference_value="unmapped",
                test_value="mapped"
            ))
    
    def analyze_bed_accuracy(self, input_file: Path) -> AccuracyReport:
        """分析 BED 格式的准确性"""
        print(f"\n{'='*60}")
        print(f"BED Accuracy Analysis")
        print(f"{'='*60}")
        print(f"Input: {input_file}")
        print(f"Reference tool: liftOver")
        print(f"{'='*60}\n")
        
        report = AccuracyReport(
            timestamp=datetime.now().isoformat(),
            input_file=str(input_file),
            format="bed",
            reference_tool="liftover"
        )
        
        # 运行所有工具
        tools = ["fastcrossmap", "crossmap", "liftover", "fastremap"]
        tool_outputs = {}
        tool_unmapped = {}
        
        for tool in tools:
            try:
                output_file, unmapped_file = self.run_tool(tool, input_file, "bed")
                tool_outputs[tool] = self.load_bed_records(output_file)
                
                # 加载未映射记录
                unmapped_records = self.load_bed_records(unmapped_file)
                tool_unmapped[tool] = set(unmapped_records.keys())
                
                print(f"    {tool}: {len(tool_outputs[tool])} mapped, {len(tool_unmapped[tool])} unmapped")
            except Exception as e:
                print(f"    {tool}: Error - {e}")
                tool_outputs[tool] = {}
                tool_unmapped[tool] = set()
        
        # 以 liftOver 为参考，比较其他工具
        ref_records = tool_outputs.get("liftover", {})
        ref_unmapped = tool_unmapped.get("liftover", set())
        
        if not ref_records:
            print("\nWarning: liftOver output is empty, cannot perform comparison")
            return report
        
        # 比较每个工具与参考
        for tool in ["fastcrossmap", "crossmap", "fastremap"]:
            if tool not in tool_outputs:
                continue
            
            test_records = tool_outputs[tool]
            test_unmapped = tool_unmapped[tool]
            
            print(f"\nComparing {tool} vs liftOver...")
            
            # 比较映射记录
            result = self.compare_bed_records(
                test_records, ref_records, tool, "liftover"
            )
            
            # 更新映射记录统计
            result.mapped_total = result.total_records
            result.mapped_matched = result.matched_records
            if result.mapped_total > 0:
                result.mapped_identity_rate = round(
                    result.mapped_matched / result.mapped_total * 100, 2
                )
            
            # 比较未映射记录
            self.compare_unmapped_records(test_unmapped, ref_unmapped, result)
            
            report.results.append(result)
            
            print(f"    Mapped Identity Rate: {result.mapped_identity_rate}%")
            print(f"    Unmapped Identity Rate: {result.unmapped_identity_rate}%")
            print(f"    Overall Identity Rate: {result.identity_rate}%")
        
        return report
    
    def generate_discrepancy_report(self, result: AccuracyResult) -> str:
        """生成差异报告"""
        lines = []
        lines.append(f"# 差异报告: {result.test_tool} vs {result.reference_tool}")
        lines.append(f"格式: {result.format}")
        lines.append(f"生成时间: {datetime.now().isoformat()}")
        lines.append("")
        
        # 摘要
        lines.append("## 摘要")
        lines.append(f"- 总记录数: {result.total_records}")
        lines.append(f"- 匹配记录: {result.matched_records}")
        lines.append(f"- 不匹配记录: {result.mismatched_records}")
        lines.append(f"- Identity Rate: {result.identity_rate}%")
        lines.append("")
        
        # 分离 Identity Rate
        lines.append("## 分离 Identity Rate")
        lines.append(f"- 映射记录 Identity Rate: {result.mapped_identity_rate}%")
        lines.append(f"  - 总数: {result.mapped_total}")
        lines.append(f"  - 匹配: {result.mapped_matched}")
        lines.append(f"- 未映射记录 Identity Rate: {result.unmapped_identity_rate}%")
        lines.append(f"  - 总数: {result.unmapped_total}")
        lines.append(f"  - 匹配: {result.unmapped_matched}")
        lines.append("")
        
        # 坐标漂移统计
        lines.append("## 坐标漂移分布")
        lines.append(f"- 完全一致 (Dist=0): {result.drift_zero}")
        lines.append(f"- 坐标系差异 (Dist=1): {result.drift_one}")
        lines.append(f"- 小漂移 (2-100): {result.drift_small}")
        lines.append(f"- 大漂移 (>100): {result.drift_large}")
        lines.append("")
        
        # 警告
        if result.warnings:
            lines.append("## ⚠️ 警告")
            for warning in result.warnings:
                lines.append(f"- {warning}")
            lines.append("")
        
        # 差异详情（限制数量）
        if result.discrepancies:
            lines.append("## 差异详情（前 50 条）")
            lines.append("| 记录 | 字段 | 参考值 | 测试值 | 漂移距离 |")
            lines.append("|------|------|--------|--------|----------|")
            for disc in result.discrepancies[:50]:
                drift_str = str(disc.drift_distance) if disc.drift_distance > 0 else "-"
                lines.append(f"| {disc.record_key[:30]} | {disc.field} | {disc.reference_value} | {disc.test_value} | {drift_str} |")
            
            if len(result.discrepancies) > 50:
                lines.append(f"\n... 还有 {len(result.discrepancies) - 50} 条差异未显示")
        
        return "\n".join(lines)



def export_report(report: AccuracyReport, output_dir: Path):
    """导出准确性分析报告"""
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # 导出 JSON
    json_path = output_dir / f"accuracy_{report.format}_{timestamp}.json"
    
    # 转换为可序列化的字典
    report_dict = {
        "timestamp": report.timestamp,
        "input_file": report.input_file,
        "format": report.format,
        "reference_tool": report.reference_tool,
        "results": []
    }
    
    for result in report.results:
        result_dict = {
            "test_tool": result.test_tool,
            "reference_tool": result.reference_tool,
            "format": result.format,
            "total_records": result.total_records,
            "matched_records": result.matched_records,
            "mismatched_records": result.mismatched_records,
            "mapped_total": result.mapped_total,
            "mapped_matched": result.mapped_matched,
            "unmapped_total": result.unmapped_total,
            "unmapped_matched": result.unmapped_matched,
            "identity_rate": result.identity_rate,
            "mapped_identity_rate": result.mapped_identity_rate,
            "unmapped_identity_rate": result.unmapped_identity_rate,
            "drift_zero": result.drift_zero,
            "drift_one": result.drift_one,
            "drift_small": result.drift_small,
            "drift_large": result.drift_large,
            "warnings": result.warnings,
            "discrepancy_count": len(result.discrepancies)
        }
        report_dict["results"].append(result_dict)
    
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(report_dict, f, indent=2, ensure_ascii=False)
    print(f"\nJSON report saved to: {json_path}")
    
    # 导出 CSV
    csv_path = output_dir / f"accuracy_{report.format}_{timestamp}.csv"
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            "test_tool", "reference_tool", "format",
            "total_records", "matched_records", "mismatched_records",
            "identity_rate", "mapped_identity_rate", "unmapped_identity_rate",
            "drift_zero", "drift_one", "drift_small", "drift_large"
        ])
        for result in report.results:
            writer.writerow([
                result.test_tool, result.reference_tool, result.format,
                result.total_records, result.matched_records, result.mismatched_records,
                result.identity_rate, result.mapped_identity_rate, result.unmapped_identity_rate,
                result.drift_zero, result.drift_one, result.drift_small, result.drift_large
            ])
    print(f"CSV report saved to: {csv_path}")
    
    # 保存最新版本
    json_latest = output_dir / f"accuracy_{report.format}_latest.json"
    with open(json_latest, 'w', encoding='utf-8') as f:
        json.dump(report_dict, f, indent=2, ensure_ascii=False)


def print_summary(report: AccuracyReport):
    """打印准确性分析摘要"""
    print(f"\n{'='*80}")
    print("ACCURACY ANALYSIS SUMMARY")
    print(f"{'='*80}")
    print(f"Reference tool: {report.reference_tool}")
    print(f"Input file: {report.input_file}")
    print(f"{'='*80}")
    
    print(f"\n{'Tool':<15} {'Identity%':<12} {'Mapped%':<12} {'Unmapped%':<12} {'Drift=0':<10} {'Drift=1':<10}")
    print(f"{'-'*80}")
    
    for result in report.results:
        print(f"{result.test_tool:<15} {result.identity_rate:<12.2f} {result.mapped_identity_rate:<12.2f} "
              f"{result.unmapped_identity_rate:<12.2f} {result.drift_zero:<10} {result.drift_one:<10}")
    
    print(f"{'='*80}")
    
    # 检查警告
    for result in report.results:
        if result.warnings:
            print(f"\n⚠️ Warnings for {result.test_tool}:")
            for warning in result.warnings:
                print(f"  - {warning}")
    
    # FastCrossMap 特别检查
    fcm_result = next((r for r in report.results if r.test_tool == "fastcrossmap"), None)
    if fcm_result:
        if fcm_result.identity_rate >= 99.0:
            print(f"\n✓ FastCrossMap achieves {fcm_result.identity_rate}% identity with liftOver (≥99% threshold)")
        else:
            print(f"\n✗ FastCrossMap identity rate ({fcm_result.identity_rate}%) is below 99% threshold")



def main():
    parser = argparse.ArgumentParser(
        description="Accuracy analysis for genome coordinate conversion tools"
    )
    parser.add_argument("--input", "-i", type=Path,
                        help="Input file for testing")
    parser.add_argument("--format", "-f", choices=["bed", "vcf", "gff", "bam"],
                        default="bed", help="Input file format")
    parser.add_argument("--chain", "-c", type=Path,
                        default=Path("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz"),
                        help="Chain file for coordinate conversion")
    parser.add_argument("--output-dir", "-o", type=Path,
                        default=Path("results/accuracy"),
                        help="Output directory for results")
    parser.add_argument("--threads", "-t", type=int, default=4,
                        help="Number of threads for FastCrossMap")
    parser.add_argument("--generate-test-data", action="store_true",
                        help="Generate test data before analysis")
    parser.add_argument("--test-records", type=int, default=10000,
                        help="Number of test records to generate")
    
    args = parser.parse_args()
    
    # 生成测试数据
    if args.generate_test_data:
        # 内联生成测试数据
        import random
        test_data_dir = Path("test_data")
        test_data_dir.mkdir(parents=True, exist_ok=True)
        bed_file = test_data_dir / "benchmark_test.bed"
        
        print(f"Generating {args.test_records} test BED records...")
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        
        with open(bed_file, 'w') as f:
            for i in range(args.test_records):
                chrom = random.choice(chromosomes)
                start = random.randint(1000000, 200000000)
                end = start + random.randint(100, 10000)
                name = f"region_{i}"
                score = random.randint(0, 1000)
                strand = random.choice(['+', '-'])
                f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
        
        print(f"Test data generated: {bed_file}")
        args.input = bed_file
    
    # 验证输入
    if not args.input or not args.input.exists():
        print("Error: Input file not specified or does not exist")
        print("Use --generate-test-data to create test data, or specify --input")
        sys.exit(1)
    
    if not args.chain.exists():
        print(f"Error: Chain file not found: {args.chain}")
        sys.exit(1)
    
    # 确保 FastCrossMap 已编译
    if not Path("./target/release/fast-crossmap").exists():
        print("Building FastCrossMap...")
        subprocess.run(["cargo", "build", "--release"], check=True)
    
    # 创建分析器
    analyzer = AccuracyAnalyzer(
        chain_file=args.chain,
        output_dir=args.output_dir,
        threads=args.threads
    )
    
    # 运行分析
    if args.format == "bed":
        report = analyzer.analyze_bed_accuracy(args.input)
    else:
        print(f"Format {args.format} analysis not yet implemented")
        sys.exit(1)
    
    # 导出报告
    export_report(report, args.output_dir)
    
    # 打印摘要
    print_summary(report)
    
    # 生成详细差异报告
    for result in report.results:
        if result.discrepancies:
            report_content = analyzer.generate_discrepancy_report(result)
            report_path = args.output_dir / f"discrepancy_{result.test_tool}_vs_{result.reference_tool}.md"
            with open(report_path, 'w', encoding='utf-8') as f:
                f.write(report_content)
            print(f"Discrepancy report saved to: {report_path}")


if __name__ == "__main__":
    main()
