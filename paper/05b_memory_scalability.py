#!/usr/bin/env python3
"""
05b_memory_scalability.py - 测试 FastCrossMap 内存可扩展性

目的: 证明流式处理架构的优势 - 内存占用与文件大小无关

测试流程:
1. 从原始 BAM 文件生成不同大小的子集
2. 对每个文件运行 FastCrossMap 并采样内存
3. 保存结果到 paper/results/memory_scalability.json

用法: python paper/05b_memory_scalability.py
输出: paper/results/memory_scalability.json
"""

import json
import subprocess
import time
import psutil
import threading
from pathlib import Path
from datetime import datetime

# =============================================================================
# 配置
# =============================================================================
RESULTS_DIR = Path("paper/results")
DATA_DIR = Path("paper/data")
TEMP_DIR = Path("paper/data/temp_scalability")

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)
TEMP_DIR.mkdir(parents=True, exist_ok=True)

# FastCrossMap 可执行文件路径
FASTCROSSMAP_BIN = "target/release/fast-crossmap"

# Chain 文件
CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"

# 原始 BAM 文件（需要足够大，至少 1GB+）
SOURCE_BAM = DATA_DIR / "encode_chipseq.bam"

# 测试的文件大小（MB）
# 从小到大测试，证明内存占用不随文件大小增长
TEST_SIZES_MB = [50, 100, 200, 500, 1000, 2000]

# 内存采样间隔（秒）
SAMPLE_INTERVAL = 0.5


def check_dependencies():
    """检查依赖工具"""
    print("检查依赖工具...")
    
    # 检查 FastCrossMap
    if not Path(FASTCROSSMAP_BIN).exists():
        print(f"错误: FastCrossMap 未找到: {FASTCROSSMAP_BIN}")
        print("请先构建: cargo build --release")
        return False
    
    # 检查 samtools
    try:
        subprocess.run(["samtools", "--version"], 
                      capture_output=True, check=True)
        print("  ✓ samtools")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("错误: samtools 未安装")
        print("请安装: conda install -c bioconda samtools")
        return False
    
    # 检查 Chain 文件
    if not CHAIN_FILE.exists():
        print(f"错误: Chain 文件未找到: {CHAIN_FILE}")
        print("请先运行: python paper/01_download_data.sh")
        return False
    
    print("  ✓ Chain 文件")
    
    # 检查源 BAM 文件
    if not SOURCE_BAM.exists():
        print(f"错误: 源 BAM 文件未找到: {SOURCE_BAM}")
        print("请先运行: python paper/01_download_data.sh")
        return False
    
    print("  ✓ 源 BAM 文件")
    
    return True


def get_file_size_mb(filepath):
    """获取文件大小（MB）"""
    return filepath.stat().st_size / (1024 * 1024)


def create_bam_subset(source_bam, output_bam, target_size_mb):
    """
    从源 BAM 文件创建指定大小的子集
    
    参数:
        source_bam: 源 BAM 文件路径
        output_bam: 输出 BAM 文件路径
        target_size_mb: 目标文件大小（MB）
    
    返回:
        实际文件大小（MB）
    """
    print(f"\n生成 {target_size_mb} MB BAM 子集...")
    
    # 获取源文件大小
    source_size_mb = get_file_size_mb(source_bam)
    print(f"  源文件大小: {source_size_mb:.2f} MB")
    
    if target_size_mb >= source_size_mb:
        print(f"  目标大小 >= 源文件，直接复制")
        subprocess.run(["cp", str(source_bam), str(output_bam)], check=True)
        return source_size_mb
    
    # 计算需要提取的比例
    ratio = target_size_mb / source_size_mb
    print(f"  提取比例: {ratio:.2%}")
    
    # 使用 samtools view 提取子集
    # -s 参数指定采样比例（需要加上随机种子）
    seed = 42  # 固定种子保证可重复性
    subsample_fraction = f"{seed}.{int(ratio * 100)}"
    
    cmd = [
        "samtools", "view",
        "-b",  # 输出 BAM 格式
        "-s", subsample_fraction,  # 采样比例
        "-o", str(output_bam),  # 输出文件
        str(source_bam)
    ]
    
    print(f"  运行: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    # 索引 BAM 文件
    print(f"  索引 BAM 文件...")
    subprocess.run(["samtools", "index", str(output_bam)], check=True)
    
    actual_size_mb = get_file_size_mb(output_bam)
    print(f"  ✓ 生成完成: {actual_size_mb:.2f} MB")
    
    return actual_size_mb


def memory_sampler(process, sample_interval, results):
    """
    内存采样线程
    
    参数:
        process: psutil.Process 对象
        sample_interval: 采样间隔（秒）
        results: 结果字典（用于存储采样数据）
    """
    results["memory_samples"] = []
    results["sample_times"] = []
    start_time = time.time()
    
    try:
        while process.is_running() and process.status() != psutil.STATUS_ZOMBIE:
            try:
                # 获取内存使用（RSS）
                mem_info = process.memory_info()
                memory_mb = mem_info.rss / (1024 * 1024)
                
                elapsed = time.time() - start_time
                results["memory_samples"].append(memory_mb)
                results["sample_times"].append(elapsed)
                
                time.sleep(sample_interval)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
    except Exception as e:
        print(f"    内存采样错误: {e}")


def run_fastcrossmap_with_memory_profiling(input_bam, output_bam, chain_file):
    """
    运行 FastCrossMap 并采样内存
    
    参数:
        input_bam: 输入 BAM 文件
        output_bam: 输出 BAM 文件
        chain_file: Chain 文件
    
    返回:
        {
            "execution_time_sec": float,
            "peak_memory_mb": float,
            "memory_samples": [float],
            "sample_times": [float]
        }
    """
    print(f"  运行 FastCrossMap...")
    
    # 构建命令
    cmd = [
        FASTCROSSMAP_BIN,
        "bam",
        str(chain_file),
        str(input_bam),
        str(output_bam)
    ]
    
    # 启动进程
    start_time = time.time()
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    # 获取 psutil Process 对象
    ps_process = psutil.Process(process.pid)
    
    # 启动内存采样线程
    results = {}
    sampler_thread = threading.Thread(
        target=memory_sampler,
        args=(ps_process, SAMPLE_INTERVAL, results)
    )
    sampler_thread.start()
    
    # 等待进程完成
    stdout, stderr = process.communicate()
    execution_time = time.time() - start_time
    
    # 等待采样线程完成
    sampler_thread.join()
    
    # 检查是否成功
    if process.returncode != 0:
        print(f"    错误: FastCrossMap 失败")
        print(f"    stderr: {stderr.decode()}")
        return None
    
    # 计算峰值内存
    peak_memory_mb = max(results["memory_samples"]) if results["memory_samples"] else 0
    
    print(f"  ✓ 完成: {execution_time:.2f}s, 峰值内存: {peak_memory_mb:.2f} MB")
    
    return {
        "execution_time_sec": execution_time,
        "peak_memory_mb": peak_memory_mb,
        "memory_samples": results["memory_samples"],
        "sample_times": results["sample_times"]
    }


def main():
    print("=" * 60)
    print("FastCrossMap 内存可扩展性测试")
    print("=" * 60)
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 检查依赖
    if not check_dependencies():
        return
    
    # 测试结果
    test_results = []
    
    # 对每个目标大小进行测试
    for target_size_mb in TEST_SIZES_MB:
        print(f"\n{'=' * 60}")
        print(f"测试文件大小: {target_size_mb} MB")
        print(f"{'=' * 60}")
        
        # 生成 BAM 子集
        subset_bam = TEMP_DIR / f"subset_{target_size_mb}mb.bam"
        output_bam = TEMP_DIR / f"output_{target_size_mb}mb.bam"
        
        try:
            actual_size_mb = create_bam_subset(SOURCE_BAM, subset_bam, target_size_mb)
            
            # 运行 FastCrossMap 并采样内存
            result = run_fastcrossmap_with_memory_profiling(
                subset_bam, output_bam, CHAIN_FILE
            )
            
            if result:
                test_results.append({
                    "target_size_mb": target_size_mb,
                    "actual_size_mb": actual_size_mb,
                    "execution_time_sec": result["execution_time_sec"],
                    "peak_memory_mb": result["peak_memory_mb"],
                    "memory_samples": result["memory_samples"],
                    "sample_times": result["sample_times"]
                })
            
            # 清理临时文件（保留子集文件，删除输出文件）
            if output_bam.exists():
                output_bam.unlink()
            if (output_bam.parent / f"{output_bam.name}.unmap").exists():
                (output_bam.parent / f"{output_bam.name}.unmap").unlink()
        
        except Exception as e:
            print(f"  错误: {e}")
            continue
    
    # 保存结果
    output_file = RESULTS_DIR / "memory_scalability.json"
    
    output_data = {
        "tool": "FastCrossMap",
        "test_date": datetime.now().isoformat(),
        "source_bam": str(SOURCE_BAM),
        "chain_file": str(CHAIN_FILE),
        "sample_interval_sec": SAMPLE_INTERVAL,
        "test_results": test_results
    }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\n{'=' * 60}")
    print("测试完成")
    print(f"{'=' * 60}")
    print(f"结果已保存到: {output_file}")
    
    # 打印摘要
    print(f"\n{'=' * 60}")
    print("测试摘要")
    print(f"{'=' * 60}")
    print(f"{'文件大小 (MB)':<15} {'执行时间 (s)':<15} {'峰值内存 (MB)':<15}")
    print("-" * 45)
    
    for result in test_results:
        print(f"{result['actual_size_mb']:<15.2f} "
              f"{result['execution_time_sec']:<15.2f} "
              f"{result['peak_memory_mb']:<15.2f}")
    
    # 计算内存占用的变化
    if len(test_results) >= 2:
        min_mem = min(r["peak_memory_mb"] for r in test_results)
        max_mem = max(r["peak_memory_mb"] for r in test_results)
        mem_variation = max_mem - min_mem
        
        print(f"\n内存占用变化:")
        print(f"  最小: {min_mem:.2f} MB")
        print(f"  最大: {max_mem:.2f} MB")
        print(f"  变化: {mem_variation:.2f} MB ({mem_variation/min_mem*100:.1f}%)")
        
        if mem_variation < 10:
            print(f"  ✓ 内存占用几乎恒定，证明流式处理架构有效")
    
    print(f"\n下一步: python paper/06b_plot_memory_scalability.py")


if __name__ == "__main__":
    main()
