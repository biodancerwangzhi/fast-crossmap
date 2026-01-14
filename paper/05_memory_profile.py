#!/usr/bin/env python3
"""
05_memory_profile.py - 内存效率测试

对 FastCrossMap, CrossMap, FastRemap 进行内存采样
每秒采样一次内存使用 (RSS)，记录完整的内存曲线

用法: python paper/05_memory_profile.py
输出: paper/results/memory_profile.json
"""

import json
import os
import subprocess
import time
import threading
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import List, Optional
import psutil

# =============================================================================
# 配置
# =============================================================================
DATA_DIR = Path("paper/data")
RESULTS_DIR = Path("paper/results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# 输入文件 - 使用 BAM 文件进行内存测试
CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
BAM_FILE = DATA_DIR / "encode_chipseq.bam"

# 解压的 chain 文件 (FastRemap 需要)
CHAIN_FILE_UNZIPPED = DATA_DIR / "hg19ToHg38.over.chain"

# 采样间隔 (秒)
SAMPLE_INTERVAL = 1.0


@dataclass
class MemoryProfile:
    """内存采样结果"""
    tool: str
    format: str
    input_file: str
    execution_time_sec: float
    peak_memory_mb: float
    memory_samples: List[float]  # 每秒采样的内存值 (MB)
    sample_times: List[float]    # 采样时间点 (相对于开始时间的秒数)
    success: bool
    error_message: Optional[str] = None


class MemorySampler:
    """内存采样器 - 在后台线程中定期采样进程内存"""
    
    def __init__(self, pid: int, interval: float = 1.0):
        self.pid = pid
        self.interval = interval
        self.samples = []
        self.sample_times = []
        self.start_time = None
        self.stop_flag = threading.Event()
        self.thread = None
        
    def _sample_loop(self):
        """采样循环"""
        try:
            process = psutil.Process(self.pid)
            while not self.stop_flag.is_set():
                try:
                    # 获取 RSS (Resident Set Size) 内存
                    mem_info = process.memory_info()
                    rss_mb = mem_info.rss / (1024 * 1024)
                    
                    elapsed = time.time() - self.start_time
                    self.samples.append(rss_mb)
                    self.sample_times.append(elapsed)
                    
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    # 进程已结束或无法访问
                    break
                
                # 等待下一次采样
                self.stop_flag.wait(self.interval)
                
        except Exception as e:
            print(f"    采样线程错误: {e}")
    
    def start(self):
        """开始采样"""
        self.start_time = time.time()
        self.stop_flag.clear()
        self.thread = threading.Thread(target=self._sample_loop, daemon=True)
        self.thread.start()
    
    def stop(self):
        """停止采样"""
        self.stop_flag.set()
        if self.thread:
            self.thread.join(timeout=2.0)
    
    def get_results(self):
        """获取采样结果"""
        return self.samples, self.sample_times


def run_with_memory_profile(cmd: List[str], output_file: Path) -> tuple[float, List[float], List[float], bool, str]:
    """
    运行命令并采样内存
    返回: (执行时间秒, 内存采样列表MB, 采样时间点列表, 是否成功, 错误信息)
    """
    try:
        # 启动进程
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # 启动内存采样器
        sampler = MemorySampler(process.pid, interval=SAMPLE_INTERVAL)
        sampler.start()
        
        start_time = time.time()
        
        # 等待进程完成
        try:
            stdout, stderr = process.communicate(timeout=1800)  # 30 分钟超时
            elapsed = time.time() - start_time
            
            # 停止采样
            sampler.stop()
            
            # 获取采样结果
            memory_samples, sample_times = sampler.get_results()
            
            if process.returncode != 0:
                return elapsed, memory_samples, sample_times, False, stderr[:500]
            
            return elapsed, memory_samples, sample_times, True, ""
            
        except subprocess.TimeoutExpired:
            process.kill()
            sampler.stop()
            return 1800, [], [], False, "Timeout after 1800 seconds"
            
    except Exception as e:
        return 0, [], [], False, str(e)


def profile_fastcrossmap(bam_file: Path, chain_file: Path, output_dir: Path) -> MemoryProfile:
    """测试 FastCrossMap 内存使用"""
    print("  Profiling FastCrossMap...")
    output_file = output_dir / "fastcrossmap_memory.bam"
    
    cmd = [
        "./target/release/fast-crossmap",
        "bam",
        str(chain_file),
        str(bam_file),
        str(output_file)
    ]
    
    elapsed, mem_samples, sample_times, success, error_msg = run_with_memory_profile(cmd, output_file)
    
    if not success or not mem_samples:
        return MemoryProfile(
            tool="FastCrossMap",
            format="BAM",
            input_file=str(bam_file),
            execution_time_sec=elapsed,
            peak_memory_mb=0,
            memory_samples=[],
            sample_times=[],
            success=False,
            error_message=error_msg if not success else "No memory samples collected"
        )
    
    peak_memory = max(mem_samples)
    
    return MemoryProfile(
        tool="FastCrossMap",
        format="BAM",
        input_file=str(bam_file),
        execution_time_sec=round(elapsed, 2),
        peak_memory_mb=round(peak_memory, 2),
        memory_samples=[round(m, 2) for m in mem_samples],
        sample_times=[round(t, 2) for t in sample_times],
        success=True
    )


def profile_crossmap(bam_file: Path, chain_file: Path, output_dir: Path) -> MemoryProfile:
    """测试 CrossMap 内存使用"""
    print("  Profiling CrossMap...")
    output_file = output_dir / "crossmap_memory.bam"
    
    cmd = [
        "CrossMap", "bam",
        str(chain_file),
        str(bam_file),
        str(output_file)
    ]
    
    elapsed, mem_samples, sample_times, success, error_msg = run_with_memory_profile(cmd, output_file)
    
    if not success or not mem_samples:
        return MemoryProfile(
            tool="CrossMap",
            format="BAM",
            input_file=str(bam_file),
            execution_time_sec=elapsed,
            peak_memory_mb=0,
            memory_samples=[],
            sample_times=[],
            success=False,
            error_message=error_msg if not success else "No memory samples collected"
        )
    
    peak_memory = max(mem_samples)
    
    return MemoryProfile(
        tool="CrossMap",
        format="BAM",
        input_file=str(bam_file),
        execution_time_sec=round(elapsed, 2),
        peak_memory_mb=round(peak_memory, 2),
        memory_samples=[round(m, 2) for m in mem_samples],
        sample_times=[round(t, 2) for t in sample_times],
        success=True
    )


def profile_fastremap(bam_file: Path, chain_file: Path, output_dir: Path) -> MemoryProfile:
    """测试 FastRemap 内存使用"""
    print("  Profiling FastRemap...")
    output_file = output_dir / "fastremap_memory.bam"
    unmap_file = output_dir / "fastremap_memory.bam.unmap"
    
    # FastRemap 不支持 .gz，需要解压的 chain 文件
    chain_unzipped = CHAIN_FILE_UNZIPPED
    if not chain_unzipped.exists():
        print("    解压 chain 文件供 FastRemap 使用...")
        subprocess.run(["gunzip", "-k", str(chain_file)], check=True)
    
    cmd = [
        "FastRemap",
        "-f", "bam",
        "-c", str(chain_unzipped),
        "-i", str(bam_file),
        "-o", str(output_file),
        "-u", str(unmap_file)
    ]
    
    elapsed, mem_samples, sample_times, success, error_msg = run_with_memory_profile(cmd, output_file)
    
    if not success or not mem_samples:
        return MemoryProfile(
            tool="FastRemap",
            format="BAM",
            input_file=str(bam_file),
            execution_time_sec=elapsed,
            peak_memory_mb=0,
            memory_samples=[],
            sample_times=[],
            success=False,
            error_message=error_msg if not success else "No memory samples collected"
        )
    
    peak_memory = max(mem_samples)
    
    return MemoryProfile(
        tool="FastRemap",
        format="BAM",
        input_file=str(bam_file),
        execution_time_sec=round(elapsed, 2),
        peak_memory_mb=round(peak_memory, 2),
        memory_samples=[round(m, 2) for m in mem_samples],
        sample_times=[round(t, 2) for t in sample_times],
        success=True
    )


def main():
    print("=" * 60)
    print("内存效率测试")
    print("=" * 60)
    
    # 检查输入文件
    if not BAM_FILE.exists():
        print(f"错误: BAM 文件不存在: {BAM_FILE}")
        print("请先运行: bash paper/01_download_data.sh")
        return
    
    if not CHAIN_FILE.exists():
        print(f"错误: Chain 文件不存在: {CHAIN_FILE}")
        print("请先运行: bash paper/01_download_data.sh")
        return
    
    # 检查 psutil 是否安装
    try:
        import psutil
    except ImportError:
        print("错误: 需要安装 psutil 库")
        print("请运行: pip install psutil")
        return
    
    # 获取文件大小
    file_size_mb = BAM_FILE.stat().st_size / (1024 * 1024)
    print(f"\n输入文件: {BAM_FILE}")
    print(f"文件大小: {file_size_mb:.2f} MB")
    print(f"采样间隔: {SAMPLE_INTERVAL} 秒")
    print()
    
    # 创建输出目录
    output_dir = RESULTS_DIR / "memory_profile"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 运行内存采样
    results = []
    
    print("[1/3] FastCrossMap")
    results.append(profile_fastcrossmap(BAM_FILE, CHAIN_FILE, output_dir))
    
    print("[2/3] CrossMap")
    results.append(profile_crossmap(BAM_FILE, CHAIN_FILE, output_dir))
    
    print("[3/3] FastRemap")
    results.append(profile_fastremap(BAM_FILE, CHAIN_FILE, output_dir))
    
    # 保存结果
    output_json = RESULTS_DIR / "memory_profile.json"
    with open(output_json, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "input_file": str(BAM_FILE),
            "input_size_mb": round(file_size_mb, 2),
            "sample_interval_sec": SAMPLE_INTERVAL,
            "results": [asdict(r) for r in results]
        }, f, indent=2)
    
    print(f"\n结果已保存到: {output_json}")
    
    # 打印摘要
    print("\n" + "=" * 60)
    print("内存采样结果摘要")
    print("=" * 60)
    print(f"{'工具':<15} {'时间(s)':<10} {'峰值内存(MB)':<15} {'采样数':<10} {'状态':<10}")
    print("-" * 60)
    
    for r in results:
        status = "✓" if r.success else "✗"
        sample_count = len(r.memory_samples)
        print(f"{r.tool:<15} {r.execution_time_sec:<10.2f} {r.peak_memory_mb:<15.2f} {sample_count:<10} {status:<10}")
    
    # 打印内存效率对比
    print("\n" + "=" * 60)
    print("内存效率对比")
    print("=" * 60)
    
    successful_results = [r for r in results if r.success]
    if len(successful_results) >= 2:
        # 找到 FastCrossMap 和 CrossMap
        fc = next((r for r in successful_results if r.tool == "FastCrossMap"), None)
        cm = next((r for r in successful_results if r.tool == "CrossMap"), None)
        
        if fc and cm:
            mem_ratio = cm.peak_memory_mb / fc.peak_memory_mb
            print(f"FastCrossMap 峰值内存: {fc.peak_memory_mb:.2f} MB")
            print(f"CrossMap 峰值内存: {cm.peak_memory_mb:.2f} MB")
            print(f"内存节省: {mem_ratio:.2f}x (CrossMap 使用 {mem_ratio:.1f} 倍内存)")
    
    print("\n下一步: python paper/06_plot_memory.py")


if __name__ == "__main__":
    main()
