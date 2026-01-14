#!/usr/bin/env python3
"""
03b_benchmark_bam_multithread.py - BAM 格式多线程扩展性测试

测试 FastCrossMap 在不同线程数下的性能
用于生成 Figure 1(d) 的数据

用法: python paper/03b_benchmark_bam_multithread.py
输出: paper/results/benchmark_bam_multithread.json
"""

import subprocess
import time
import json
import os
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
BAM_FILE = DATA_DIR / "encode_chipseq.bam"

# 测试线程数
THREAD_COUNTS = [1, 2, 4, 8, 16]

# 每个配置运行次数
NUM_RUNS = 5


def get_file_size_mb(filepath):
    """获取文件大小 (MB)"""
    return os.path.getsize(filepath) / (1024 * 1024)


def run_fastcrossmap_bam(chain_file, input_file, output_file, threads=1):
    """运行 FastCrossMap BAM 转换并返回执行时间"""
    cmd = [
        "./target/release/fast-crossmap", "bam",
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
    print("FastCrossMap BAM 多线程扩展性测试")
    print("=" * 60)
    
    # 检查文件
    if not CHAIN_FILE.exists():
        print(f"错误: Chain 文件不存在: {CHAIN_FILE}")
        print("请先运行: bash paper/01_download_data.sh")
        return
    
    if not BAM_FILE.exists():
        print(f"错误: BAM 文件不存在: {BAM_FILE}")
        print("请先运行: bash paper/01_download_data.sh")
        return
    
    # 获取文件大小
    file_size_mb = get_file_size_mb(BAM_FILE)
    print(f"输入文件: {BAM_FILE}")
    print(f"文件大小: {file_size_mb:.2f} MB")
    print(f"测试线程数: {THREAD_COUNTS}")
    print(f"每配置运行次数: {NUM_RUNS}")
    print()
    
    results = []
    
    for threads in THREAD_COUNTS:
        print(f"\n测试 {threads} 线程...")
        output_file = RESULTS_DIR / f"fastcrossmap_bam_mt{threads}_output.bam"
        
        times = []
        for run in range(NUM_RUNS):
            result = run_fastcrossmap_bam(CHAIN_FILE, BAM_FILE, output_file, threads)
            if result["success"]:
                times.append(result["time"])
                print(f"  Run {run+1}: {result['time']:.2f}s")
            else:
                print(f"  Run {run+1}: FAILED - {result['stderr'][:100]}")
        
        if times:
            avg_time = sum(times) / len(times)
            min_time = min(times)
            max_time = max(times)
            throughput = file_size_mb / avg_time
            
            results.append({
                "threads": threads,
                "execution_time_sec": avg_time,
                "min_time_sec": min_time,
                "max_time_sec": max_time,
                "all_times": times,
                "throughput_mb_per_sec": throughput,
                "success": True
            })
            
            print(f"  平均: {avg_time:.2f}s (min: {min_time:.2f}s, max: {max_time:.2f}s)")
            print(f"  吞吐量: {throughput:.2f} MB/sec")
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
                print(f"{r['threads']}T: {r['execution_time_sec']:.2f}s, "
                      f"加速比: {speedup:.2f}x, 效率: {efficiency:.1f}%")
    
    # 保存结果
    output_data = {
        "timestamp": datetime.now().isoformat(),
        "format": "BAM",
        "input_file": str(BAM_FILE),
        "input_size_mb": file_size_mb,
        "chain_file": str(CHAIN_FILE),
        "num_runs": NUM_RUNS,
        "results": results
    }
    
    output_file = RESULTS_DIR / "benchmark_bam_multithread.json"
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\n结果已保存到: {output_file}")
    print("\n下一步: python paper/04_plot_performance.py")


if __name__ == "__main__":
    main()
