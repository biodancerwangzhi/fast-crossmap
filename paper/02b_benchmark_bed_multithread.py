#!/usr/bin/env python3
"""
02b_benchmark_bed_multithread.py - BED 格式多线程扩展性测试

测试 FastCrossMap 在不同线程数下的性能
用于生成 Figure 1(b) 的数据

用法: python paper/02b_benchmark_bed_multithread.py
输出: paper/results/benchmark_bed_multithread.json
"""

import subprocess
import time
import json
from pathlib import Path
from datetime import datetime

# =============================================================================
# 配置
# =============================================================================
DATA_DIR = Path("paper/data")
RESULTS_DIR = Path("paper/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# 测试文件
CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
BED_FILE = DATA_DIR / "encode_dnase_peaks.bed.gz"

# 测试线程数
THREAD_COUNTS = [1, 2, 4, 8, 16]

# 每个配置运行次数
NUM_RUNS = 5


def count_bed_records(bed_file):
    """统计 BED 文件记录数 (支持 .gz 压缩)"""
    import gzip
    
    count = 0
    bed_path = Path(bed_file)
    
    if str(bed_path).endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    with opener(bed_path, mode) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                count += 1
    return count


def run_fastcrossmap(chain_file, input_file, output_file, threads=1):
    """运行 FastCrossMap 并返回执行时间"""
    cmd = [
        "./target/release/fast-crossmap", "bed",
        "-t", str(threads),
        str(chain_file),
        str(input_file),
        str(output_file)
    ]
    
    start = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - start
    
    return {
        "success": result.returncode == 0,
        "time": elapsed,
        "stderr": result.stderr
    }


def main():
    print("=" * 60)
    print("FastCrossMap BED 多线程扩展性测试")
    print("=" * 60)
    
    # 检查文件
    if not CHAIN_FILE.exists():
        print(f"错误: Chain 文件不存在: {CHAIN_FILE}")
        print("请先运行: bash paper/01_download_data.sh")
        return
    
    if not BED_FILE.exists():
        print(f"错误: BED 文件不存在: {BED_FILE}")
        print("请先运行: bash paper/01_download_data.sh")
        return
    
    # 如果 BED 文件是 .gz 格式，先解压
    import gzip
    import shutil
    
    bed_file_to_use = BED_FILE
    if str(BED_FILE).endswith('.gz'):
        bed_file_unzipped = Path(str(BED_FILE)[:-3])  # 去掉 .gz 后缀
        if not bed_file_unzipped.exists():
            print(f"解压 BED 文件: {BED_FILE} -> {bed_file_unzipped}")
            with gzip.open(BED_FILE, 'rb') as f_in:
                with open(bed_file_unzipped, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        bed_file_to_use = bed_file_unzipped
    
    # 统计记录数
    num_records = count_bed_records(BED_FILE)
    print(f"输入文件: {bed_file_to_use}")
    print(f"记录数: {num_records:,}")
    print(f"测试线程数: {THREAD_COUNTS}")
    print(f"每配置运行次数: {NUM_RUNS}")
    print()
    
    results = []
    
    for threads in THREAD_COUNTS:
        print(f"\n测试 {threads} 线程...")
        output_file = RESULTS_DIR / f"fastcrossmap_mt{threads}_output.bed"
        
        times = []
        for run in range(NUM_RUNS):
            result = run_fastcrossmap(CHAIN_FILE, bed_file_to_use, output_file, threads)
            if result["success"]:
                times.append(result["time"])
                print(f"  Run {run+1}: {result['time']:.3f}s")
            else:
                print(f"  Run {run+1}: FAILED - {result['stderr'][:100]}")
        
        if times:
            avg_time = sum(times) / len(times)
            min_time = min(times)
            max_time = max(times)
            throughput = num_records / avg_time
            
            results.append({
                "threads": threads,
                "execution_time_sec": avg_time,
                "min_time_sec": min_time,
                "max_time_sec": max_time,
                "all_times": times,
                "throughput_rec_per_sec": throughput,
                "success": True
            })
            
            print(f"  平均: {avg_time:.3f}s (min: {min_time:.3f}s, max: {max_time:.3f}s)")
            print(f"  吞吐量: {throughput:,.0f} records/sec")
        else:
            results.append({
                "threads": threads,
                "success": False,
                "error": "All runs failed"
            })
    
    # 计算加速比
    if results and results[0]["success"]:
        baseline = results[0]["execution_time_sec"]
        print("\n" + "=" * 60)
        print("扩展性分析")
        print("=" * 60)
        for r in results:
            if r["success"]:
                speedup = baseline / r["execution_time_sec"]
                efficiency = speedup / r["threads"] * 100
                print(f"{r['threads']}T: {r['execution_time_sec']:.3f}s, "
                      f"加速比: {speedup:.2f}x, 效率: {efficiency:.1f}%")
    
    # 保存结果
    output_data = {
        "timestamp": datetime.now().isoformat(),
        "format": "BED",
        "input_file": str(BED_FILE),
        "input_records": num_records,
        "chain_file": str(CHAIN_FILE),
        "num_runs": NUM_RUNS,
        "results": results
    }
    
    output_file = RESULTS_DIR / "benchmark_bed_multithread.json"
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\n结果已保存到: {output_file}")
    print("\n下一步: python paper/04_plot_performance.py")


if __name__ == "__main__":
    main()
