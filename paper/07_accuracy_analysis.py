#!/usr/bin/env python3
"""
07_accuracy_analysis.py - 准确性分析

以 UCSC liftOver 为金标准，评估 FastCrossMap、CrossMap、FastRemap 的准确性

注意：坐标转换可能产生一对多映射（区间分割），本脚本通过 name 字段追踪记录

用法: python paper/07_accuracy_analysis.py
输出: paper/results/accuracy.json
"""

import json
import subprocess
from collections import defaultdict
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

# =============================================================================
# 配置
# =============================================================================
DATA_DIR = Path("paper/data")
RESULTS_DIR = Path("paper/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# 输入文件
CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
BED_FILE = DATA_DIR / "encode_dnase_peaks.bed.gz"

# 解压的 chain 文件 (FastRemap 需要)
CHAIN_FILE_UNZIPPED = DATA_DIR / "hg19ToHg38.over.chain"


@dataclass
class BedRecord:
    """BED 记录"""
    chrom: str
    start: int
    end: int
    name: str = "."
    
    @classmethod
    def from_line(cls, line: str):
        """从 BED 行解析"""
        fields = line.strip().split('\t')
        if len(fields) < 3:
            return None
        
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3] if len(fields) > 3 else "."
        
        return cls(chrom, start, end, name)
    
    def __eq__(self, other) -> bool:
        """判断两个记录是否完全一致（忽略 name）"""
        if not isinstance(other, BedRecord):
            return False
        return (self.chrom == other.chrom and 
                self.start == other.start and 
                self.end == other.end)


@dataclass
class AccuracyMetrics:
    """准确性指标"""
    tool: str
    total_input_records: int
    mapped_records: int
    unmapped_records: int
    mapping_rate: float
    
    # 与 liftOver 对比
    identical_records: int          # 完全一致
    partial_match: int              # 部分匹配（区间分割）
    coordinate_mismatch: int        # 坐标不一致
    missing_in_tool: int            # 工具未映射但 liftOver 映射了
    
    identity_rate: float            # 完全一致率
    
    success: bool
    error_message: Optional[str] = None


def create_indexed_bed(input_bed: Path, output_bed: Path) -> int:
    """
    创建带索引的 BED 文件，将行号作为 name 字段
    返回记录数
    """
    import gzip
    
    if str(input_bed).endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    count = 0
    with opener(input_bed, mode) as fin:
        with open(output_bed, 'w') as fout:
            for line in fin:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    # 使用行号作为 name (从 0 开始)
                    chrom, start, end = fields[0], fields[1], fields[2]
                    fout.write(f"{chrom}\t{start}\t{end}\tID_{count}\t0\t.\n")
                    count += 1
    
    return count


def run_tool_and_load_output(tool: str, indexed_bed: Path, chain_file: Path, 
                             output_dir: Path) -> Dict[int, List[BedRecord]]:
    """
    运行工具并加载输出
    返回: {record_id: [BedRecord, ...]}
    
    record_id 是输入记录的索引（从 name 字段解析）
    一条输入记录可能对应多条输出记录（区间分割）
    """
    output_file = output_dir / f"{tool.lower()}_accuracy.bed"
    unmap_file = Path(str(output_file) + ".unmap")
    
    # 根据工具选择命令
    if tool == "FastCrossMap":
        cmd = [
            "./target/release/fast-crossmap",
            "bed",
            str(chain_file),
            str(indexed_bed),
            str(output_file)
        ]
    elif tool == "CrossMap":
        cmd = [
            "CrossMap", "bed",
            str(chain_file),
            str(indexed_bed),
            str(output_file)
        ]
    elif tool == "liftOver":
        cmd = [
            "liftOver",
            str(indexed_bed),
            str(chain_file),
            str(output_file),
            str(unmap_file)
        ]
    elif tool == "FastRemap":
        # FastRemap 需要解压的 chain 文件
        chain_unzipped = CHAIN_FILE_UNZIPPED
        if not chain_unzipped.exists():
            subprocess.run(["gunzip", "-k", str(chain_file)], check=True)
        
        # FastRemap 会自动在输出文件名后加 .bed 后缀
        # 所以我们需要去掉 .bed 后缀
        output_base = str(output_file).replace('.bed', '')
        unmap_base = str(unmap_file).replace('.bed.unmap', '')
        
        cmd = [
            "FastRemap",
            "-f", "bed",
            "-c", str(chain_unzipped),
            "-i", str(indexed_bed),
            "-o", output_base,
            "-u", unmap_base
        ]
    else:
        raise ValueError(f"Unknown tool: {tool}")
    
    # 运行命令
    print(f"  Running {tool}...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"  Warning: {tool} failed: {result.stderr[:200]}")
        return {}
    
    # 加载输出 - 按 name 字段分组
    mapped_records = defaultdict(list)
    
    if output_file.exists():
        with open(output_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    record = BedRecord.from_line(line)
                    if record and record.name.startswith("ID_"):
                        try:
                            record_id = int(record.name.split("_")[1])
                            mapped_records[record_id].append(record)
                        except:
                            pass
    
    return dict(mapped_records)


def compare_with_gold_standard(tool_mapped: Dict[int, List[BedRecord]], 
                               gold_mapped: Dict[int, List[BedRecord]],
                               total_records: int) -> Dict:
    """
    与金标准对比，计算准确性指标
    """
    identical = 0
    partial_match = 0
    coord_mismatch = 0
    missing_in_tool = 0
    
    for record_id in range(total_records):
        gold_records = gold_mapped.get(record_id, [])
        tool_records = tool_mapped.get(record_id, [])
        
        if gold_records and tool_records:
            # 两者都映射了
            # 检查是否完全一致（所有输出记录都相同）
            if len(gold_records) == len(tool_records):
                # 排序后比较
                gold_sorted = sorted(gold_records, key=lambda r: (r.chrom, r.start, r.end))
                tool_sorted = sorted(tool_records, key=lambda r: (r.chrom, r.start, r.end))
                
                if all(g == t for g, t in zip(gold_sorted, tool_sorted)):
                    identical += 1
                else:
                    coord_mismatch += 1
            else:
                # 输出记录数不同，算作部分匹配
                partial_match += 1
        elif gold_records and not tool_records:
            # 金标准映射了，但工具未映射
            missing_in_tool += 1
        # 如果工具映射了但金标准未映射，这种情况很少见，暂不单独统计
    
    return {
        "identical": identical,
        "partial_match": partial_match,
        "coordinate_mismatch": coord_mismatch,
        "missing_in_tool": missing_in_tool
    }


def analyze_accuracy(tool: str, indexed_bed: Path, chain_file: Path, 
                    gold_mapped: Dict[int, List[BedRecord]], 
                    total_records: int,
                    output_dir: Path) -> AccuracyMetrics:
    """分析工具的准确性"""
    print(f"\n[{tool}]")
    
    # 运行工具并加载输出
    tool_mapped = run_tool_and_load_output(tool, indexed_bed, chain_file, output_dir)
    
    if not tool_mapped:
        return AccuracyMetrics(
            tool=tool,
            total_input_records=total_records,
            mapped_records=0,
            unmapped_records=total_records,
            mapping_rate=0.0,
            identical_records=0,
            partial_match=0,
            coordinate_mismatch=0,
            missing_in_tool=0,
            identity_rate=0.0,
            success=False,
            error_message="Failed to run tool or load output"
        )
    
    # 与金标准对比
    comparison = compare_with_gold_standard(
        tool_mapped, gold_mapped, total_records
    )
    
    # 计算指标
    mapped_count = len(tool_mapped)
    unmapped_count = total_records - mapped_count
    mapping_rate = mapped_count / total_records if total_records > 0 else 0.0
    
    identical = comparison["identical"]
    identity_rate = identical / total_records if total_records > 0 else 0.0
    
    print(f"  Mapped: {mapped_count} / {total_records} ({mapping_rate*100:.2f}%)")
    print(f"  Identical with liftOver: {identical} ({identity_rate*100:.2f}%)")
    print(f"  Partial match (split regions): {comparison['partial_match']}")
    print(f"  Coordinate mismatch: {comparison['coordinate_mismatch']}")
    print(f"  Missing in tool: {comparison['missing_in_tool']}")
    
    return AccuracyMetrics(
        tool=tool,
        total_input_records=total_records,
        mapped_records=mapped_count,
        unmapped_records=unmapped_count,
        mapping_rate=round(mapping_rate, 4),
        identical_records=identical,
        partial_match=comparison["partial_match"],
        coordinate_mismatch=comparison["coordinate_mismatch"],
        missing_in_tool=comparison["missing_in_tool"],
        identity_rate=round(identity_rate, 4),
        success=True
    )


def main():
    print("=" * 60)
    print("准确性分析 (以 liftOver 为金标准)")
    print("=" * 60)
    
    # 检查输入文件
    if not BED_FILE.exists():
        print(f"错误: BED 文件不存在: {BED_FILE}")
        print("请先运行: bash paper/01_download_data.sh")
        return
    
    if not CHAIN_FILE.exists():
        print(f"错误: Chain 文件不存在: {CHAIN_FILE}")
        print("请先运行: bash paper/01_download_data.sh")
        return
    
    # 创建输出目录
    output_dir = RESULTS_DIR / "accuracy_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 创建带索引的 BED 文件
    print(f"\n创建带索引的 BED 文件...")
    indexed_bed = output_dir / "indexed_input.bed"
    total_records = create_indexed_bed(BED_FILE, indexed_bed)
    print(f"输入记录数: {total_records:,}")
    
    # 运行 liftOver 作为金标准
    print("\n" + "=" * 60)
    print("运行 liftOver (金标准)")
    print("=" * 60)
    gold_mapped = run_tool_and_load_output(
        "liftOver", indexed_bed, CHAIN_FILE, output_dir
    )
    
    gold_mapped_count = len(gold_mapped)
    gold_unmapped_count = total_records - gold_mapped_count
    
    print(f"  liftOver mapped: {gold_mapped_count}")
    print(f"  liftOver unmapped: {gold_unmapped_count}")
    
    if not gold_mapped:
        print("错误: liftOver 未能生成输出")
        return
    
    # 分析各工具的准确性
    print("\n" + "=" * 60)
    print("分析工具准确性")
    print("=" * 60)
    
    results = []
    
    for tool in ["FastCrossMap", "CrossMap", "FastRemap"]:
        metrics = analyze_accuracy(
            tool, indexed_bed, CHAIN_FILE,
            gold_mapped, total_records,
            output_dir
        )
        results.append(metrics)
    
    # 保存结果
    output_json = RESULTS_DIR / "accuracy.json"
    with open(output_json, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "input_file": str(BED_FILE),
            "total_input_records": total_records,
            "gold_standard": "liftOver",
            "gold_standard_mapped": gold_mapped_count,
            "gold_standard_unmapped": gold_unmapped_count,
            "results": [asdict(r) for r in results]
        }, f, indent=2)
    
    print(f"\n结果已保存到: {output_json}")
    
    # 打印摘要
    print("\n" + "=" * 60)
    print("准确性分析摘要")
    print("=" * 60)
    print(f"{'工具':<15} {'映射率':<10} {'一致率':<10} {'部分匹配':<10} {'坐标偏差':<10}")
    print("-" * 60)
    
    for r in results:
        if r.success:
            print(f"{r.tool:<15} {r.mapping_rate*100:<10.2f}% {r.identity_rate*100:<10.2f}% "
                  f"{r.partial_match:<10} {r.coordinate_mismatch:<10}")
    
    print("\n说明:")
    print("- 一致率: 与 liftOver 输出完全相同的记录比例")
    print("- 部分匹配: 区间被分割，输出记录数不同")
    print("- 坐标偏差: 输出记录数相同但坐标不同")
    print("\n下一步: python paper/08_plot_accuracy.py")


if __name__ == "__main__":
    main()
