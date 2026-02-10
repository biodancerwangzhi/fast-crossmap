#!/usr/bin/env python3
"""
07_accuracy_analysis.py - Accuracy analysis

Evaluate accuracy of FastCrossMap, CrossMap, FastRemap using UCSC liftOver as gold standard

Note: Coordinate conversion may produce one-to-many mappings (interval splitting);
this script tracks records via the name field

Usage: python paper/07_accuracy_analysis.py
Output: paper/results/accuracy.json
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

# Input files
CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
BED_FILE = DATA_DIR / "encode_dnase_peaks.bed.gz"

# Uncompressed chain file (required by FastRemap)
CHAIN_FILE_UNZIPPED = DATA_DIR / "hg19ToHg38.over.chain"


@dataclass
class BedRecord:
    """BED record"""
    chrom: str
    start: int
    end: int
    name: str = "."
    
    @classmethod
    def from_line(cls, line: str):
        """Parse from BED line"""
        fields = line.strip().split('\t')
        if len(fields) < 3:
            return None
        
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3] if len(fields) > 3 else "."
        
        return cls(chrom, start, end, name)
    
    def __eq__(self, other) -> bool:
        """Check if two records are identical (ignoring name)"""
        if not isinstance(other, BedRecord):
            return False
        return (self.chrom == other.chrom and 
                self.start == other.start and 
                self.end == other.end)


@dataclass
class AccuracyMetrics:
    """Accuracy metrics"""
    tool: str
    total_input_records: int
    mapped_records: int
    unmapped_records: int
    mapping_rate: float
    
    # 与 liftOver 对比
    identical_records: int          # Exact match
    partial_match: int              # Partial match (interval splitting)
    coordinate_mismatch: int        # Coordinate mismatch
    missing_in_tool: int            # Tool unmapped but liftOver mapped
    
    identity_rate: float            # Exact match rate
    
    success: bool
    error_message: Optional[str] = None


def create_indexed_bed(input_bed: Path, output_bed: Path) -> int:
    """
    Create indexed BED file with line number as name field.
    Returns record count.
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
                    # Use line number as name (starting from 0)
                    chrom, start, end = fields[0], fields[1], fields[2]
                    fout.write(f"{chrom}\t{start}\t{end}\tID_{count}\t0\t.\n")
                    count += 1
    
    return count


def run_tool_and_load_output(tool: str, indexed_bed: Path, chain_file: Path, 
                             output_dir: Path) -> Dict[int, List[BedRecord]]:
    """
    Run tool and load output.
    Returns: {record_id: [BedRecord, ...]}
    
    record_id is the input record index (parsed from name field)
    One input record may correspond to multiple output records (interval splitting)
    """
    output_file = output_dir / f"{tool.lower()}_accuracy.bed"
    unmap_file = Path(str(output_file) + ".unmap")
    
    # Based on tool, choose command
    if tool == "FastCrossMap":
        cmd = [
            "./fast-crossmap-linux-x64/fast-crossmap",
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
        # FastRemap needs uncompressed chain file
        chain_unzipped = CHAIN_FILE_UNZIPPED
        if not chain_unzipped.exists():
            subprocess.run(["gunzip", "-k", str(chain_file)], check=True)
        
        # FastRemap automatically appends .bed suffix to output filename
        # So we need to remove the .bed suffix
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
    
    # Run command
    print(f"  Running {tool}...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"  Warning: {tool} failed: {result.stderr[:200]}")
        return {}
    
    # Load output - group by name field
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
    Compare with gold standard, calculate accuracy metrics.
    """
    identical = 0
    partial_match = 0
    coord_mismatch = 0
    missing_in_tool = 0
    
    for record_id in range(total_records):
        gold_records = gold_mapped.get(record_id, [])
        tool_records = tool_mapped.get(record_id, [])
        
        if gold_records and tool_records:
            # Both mapped
            # Check if completely identical (all output records match)
            if len(gold_records) == len(tool_records):
                # Sort and compare
                gold_sorted = sorted(gold_records, key=lambda r: (r.chrom, r.start, r.end))
                tool_sorted = sorted(tool_records, key=lambda r: (r.chrom, r.start, r.end))
                
                if all(g == t for g, t in zip(gold_sorted, tool_sorted)):
                    identical += 1
                else:
                    coord_mismatch += 1
            else:
                # Different number of output records, count as partial match
                partial_match += 1
        elif gold_records and not tool_records:
            # Gold standard mapped, but tool did not
            missing_in_tool += 1
        # If tool mapped but gold standard did not, this is rare, not counted separately
    
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
    """Analyze tool accuracy"""
    print(f"\n[{tool}]")
    
    # Run tool and load output
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
    
    # Compare with gold standard
    comparison = compare_with_gold_standard(
        tool_mapped, gold_mapped, total_records
    )
    
    # Calculate metrics
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
    print("Accuracy Analysis (Gold Standard: liftOver)")
    print("=" * 60)
    
    # Check input files
    if not BED_FILE.exists():
        print(f"Error: BED file not found: {BED_FILE}")
        print("Please run first: bash paper/01_download_data.sh")
        return
    
    if not CHAIN_FILE.exists():
        print(f"Error: Chain file not found: {CHAIN_FILE}")
        print("Please run first: bash paper/01_download_data.sh")
        return
    
    # Create output directory
    output_dir = RESULTS_DIR / "accuracy_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create indexed BED file
    print(f"\nCreating indexed BED file...")
    indexed_bed = output_dir / "indexed_input.bed"
    total_records = create_indexed_bed(BED_FILE, indexed_bed)
    print(f"Input records: {total_records:,}")
    
    # Run liftOver as gold standard
    print("\n" + "=" * 60)
    print("Running liftOver (Gold Standard)")
    print("=" * 60)
    gold_mapped = run_tool_and_load_output(
        "liftOver", indexed_bed, CHAIN_FILE, output_dir
    )
    
    gold_mapped_count = len(gold_mapped)
    gold_unmapped_count = total_records - gold_mapped_count
    
    print(f"  liftOver mapped: {gold_mapped_count}")
    print(f"  liftOver unmapped: {gold_unmapped_count}")
    
    if not gold_mapped:
        print("Error: liftOver failed to generate output")
        return
    
    # Analyze accuracy of each tool
    print("\n" + "=" * 60)
    print("Analyzing Tool Accuracy")
    print("=" * 60)
    
    results = []
    
    for tool in ["FastCrossMap", "CrossMap", "FastRemap"]:
        metrics = analyze_accuracy(
            tool, indexed_bed, CHAIN_FILE,
            gold_mapped, total_records,
            output_dir
        )
        results.append(metrics)
    
    # Save results
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
    
    print(f"\nResults saved to: {output_json}")
    
    # Print summary
    print("\n" + "=" * 60)
    print("Accuracy Analysis Summary")
    print("=" * 60)
    print(f"{'Tool':<15} {'Map Rate':<10} {'Identity':<10} {'Partial':<10} {'Coord Diff':<10}")
    print("-" * 60)
    
    for r in results:
        if r.success:
            print(f"{r.tool:<15} {r.mapping_rate*100:<10.2f}% {r.identity_rate*100:<10.2f}% "
                  f"{r.partial_match:<10} {r.coordinate_mismatch:<10}")
    
    print("\nNotes:")
    print("- Identity rate: Percentage of records with identical coordinates to liftOver")
    print("- Partial match: Intervals split into different number of output regions")
    print("- Coord mismatch: Same number of output records but different coordinates")
    print("\nNext step: python paper/08_plot_accuracy.py")


if __name__ == "__main__":
    main()
