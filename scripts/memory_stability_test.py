#!/usr/bin/env python3
"""
内存稳定性长程测试脚本

针对大型 BAM 文件进行内存采样测试，检测内存泄漏模式，
验证 FastCrossMap 内存曲线是否水平稳定。

Usage:
    python scripts/memory_stability_test.py --input test_data/external/H-CAF.fq.bam
    python scripts/memory_stability_test.py --input test_data/external/H-CAF.fq.bam --tools fastcrossmap crossmap
    python scripts/memory_stability_test.py --input test_data/external/H-CAF.fq.bam --sample-interval 0.5

Features:
    - 每秒（或自定义间隔）采样内存使用
    - 检测内存泄漏模式（线性增长）
    - 生成内存随时间变化的数据
    - 支持多工具对比测试
"""

import argparse
import json
import os
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import psutil
except ImportError:
    print("Error: psutil is required. Install with: pip install psutil")
    sys.exit(1)


# 工具命令模板（BAM 格式）
TOOL_COMMANDS = {
    "fastcrossmap": "./target/release/fast-crossmap bam {chain} {input} {output} -t {threads}",
    "crossmap": "CrossMap bam {chain} {input} {output}",
    # liftOver 不支持 BAM 格式
    # FastRemap 支持 BAM，但需要解压的 chain 文件
    "fastremap": "FastRemap -f bam -c {chain_uncompressed} -i {input} -u {unmapped} -o {output}",
}


@dataclass
class MemorySample:
    """单个内存采样点"""
    timestamp: float  # 相对于开始时间的秒数
    rss_mb: float     # Resident Set Size (MB)
    vms_mb: float     # Virtual Memory Size (MB)


@dataclass
class MemoryProfile:
    """工具的内存使用概况"""
    tool: str
    input_file: str
    input_size_mb: float
    samples: List[MemorySample] = field(default_factory=list)
    peak_rss_mb: float = 0.0
    avg_rss_mb: float = 0.0
    min_rss_mb: float = 0.0
    max_rss_mb: float = 0.0
    execution_time_sec: float = 0.0
    exit_code: int = 0
    memory_leak_detected: bool = False
    leak_rate_mb_per_sec: float = 0.0  # 内存增长率
    is_stable: bool = True  # 内存曲线是否稳定
    error_message: str = ""


@dataclass
class MemoryTestReport:
    """内存测试报告"""
    timestamp: str
    input_file: str
    input_size_mb: float
    sample_interval_sec: float
    profiles: List[MemoryProfile] = field(default_factory=list)


class MemorySampler:
    """内存采样器 - 在后台线程中采样进程内存"""
    
    def __init__(self, pid: int, interval: float = 1.0):
        self.pid = pid
        self.interval = interval
        self.samples: List[MemorySample] = []
        self._stop_event = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._start_time = 0.0
    
    def start(self):
        """开始采样"""
        self._start_time = time.time()
        self._stop_event.clear()
        self._thread = threading.Thread(target=self._sample_loop, daemon=True)
        self._thread.start()
    
    def stop(self) -> List[MemorySample]:
        """停止采样并返回结果"""
        self._stop_event.set()
        if self._thread:
            self._thread.join(timeout=5.0)
        return self.samples
    
    def _sample_loop(self):
        """采样循环"""
        try:
            process = psutil.Process(self.pid)
        except psutil.NoSuchProcess:
            return
        
        while not self._stop_event.is_set():
            try:
                # 获取内存信息
                mem_info = process.memory_info()
                rss_mb = mem_info.rss / (1024 * 1024)
                vms_mb = mem_info.vms / (1024 * 1024)
                
                # 也检查子进程
                try:
                    for child in process.children(recursive=True):
                        try:
                            child_mem = child.memory_info()
                            rss_mb += child_mem.rss / (1024 * 1024)
                            vms_mb += child_mem.vms / (1024 * 1024)
                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            pass
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
                
                # 记录采样
                timestamp = time.time() - self._start_time
                self.samples.append(MemorySample(
                    timestamp=round(timestamp, 3),
                    rss_mb=round(rss_mb, 2),
                    vms_mb=round(vms_mb, 2)
                ))
                
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
            
            # 等待下一次采样
            self._stop_event.wait(self.interval)


def analyze_memory_stability(samples: List[MemorySample]) -> Tuple[bool, float, bool]:
    """
    分析内存稳定性
    
    Returns:
        (is_stable, leak_rate_mb_per_sec, memory_leak_detected)
    """
    if len(samples) < 10:
        return True, 0.0, False
    
    # 跳过前 10% 的样本（启动阶段）
    skip_count = max(1, len(samples) // 10)
    stable_samples = samples[skip_count:]
    
    if len(stable_samples) < 5:
        return True, 0.0, False
    
    # 计算线性回归斜率（内存增长率）
    n = len(stable_samples)
    sum_x = sum(s.timestamp for s in stable_samples)
    sum_y = sum(s.rss_mb for s in stable_samples)
    sum_xy = sum(s.timestamp * s.rss_mb for s in stable_samples)
    sum_x2 = sum(s.timestamp ** 2 for s in stable_samples)
    
    denominator = n * sum_x2 - sum_x ** 2
    if abs(denominator) < 1e-10:
        return True, 0.0, False
    
    slope = (n * sum_xy - sum_x * sum_y) / denominator  # MB/sec
    
    # 计算内存波动范围
    rss_values = [s.rss_mb for s in stable_samples]
    rss_range = max(rss_values) - min(rss_values)
    rss_mean = sum(rss_values) / len(rss_values)
    
    # 判断标准：
    # 1. 斜率 > 1 MB/sec 认为有内存泄漏
    # 2. 波动范围 > 平均值的 50% 认为不稳定
    memory_leak_detected = slope > 1.0
    is_stable = rss_range < rss_mean * 0.5 and not memory_leak_detected
    
    return is_stable, round(slope, 4), memory_leak_detected


class MemoryStabilityTester:
    """内存稳定性测试器"""
    
    def __init__(
        self,
        chain_file: Path,
        output_dir: Path = Path("results/memory_test"),
        sample_interval: float = 1.0,
        threads: int = 4
    ):
        self.chain_file = chain_file
        self.output_dir = output_dir
        self.sample_interval = sample_interval
        self.threads = threads
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 准备解压的 chain 文件（给 FastRemap 用）
        self.chain_file_uncompressed = self._prepare_uncompressed_chain()
    
    def _prepare_uncompressed_chain(self) -> Path:
        """准备解压的 chain 文件"""
        if not str(self.chain_file).endswith('.gz'):
            return self.chain_file
        
        uncompressed = Path(str(self.chain_file)[:-3])
        if uncompressed.exists():
            return uncompressed
        
        print(f"Decompressing chain file for FastRemap...")
        try:
            import gzip
            import shutil
            with gzip.open(self.chain_file, 'rb') as f_in:
                with open(uncompressed, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            return uncompressed
        except Exception as e:
            print(f"Warning: Could not decompress chain file: {e}")
            return self.chain_file
    
    def build_command(self, tool: str, input_file: Path, output_file: Path) -> Optional[str]:
        """构建工具命令"""
        template = TOOL_COMMANDS.get(tool)
        if not template:
            return None
        
        unmapped_file = output_file.with_suffix('.bam.unmap')
        
        return template.format(
            chain=str(self.chain_file),
            chain_uncompressed=str(self.chain_file_uncompressed),
            input=str(input_file),
            output=str(output_file),
            unmapped=str(unmapped_file),
            threads=self.threads
        )
    
    def run_with_memory_sampling(
        self,
        tool: str,
        input_file: Path,
        timeout: int = 7200  # 2 小时超时
    ) -> MemoryProfile:
        """运行工具并采样内存"""
        
        input_size_mb = input_file.stat().st_size / (1024 * 1024)
        
        profile = MemoryProfile(
            tool=tool,
            input_file=str(input_file),
            input_size_mb=round(input_size_mb, 2)
        )
        
        # 构建命令
        output_file = self.output_dir / f"{tool}_output.bam"
        cmd = self.build_command(tool, input_file, output_file)
        
        if not cmd:
            profile.error_message = f"Tool {tool} not supported for BAM format"
            profile.exit_code = -1
            return profile
        
        print(f"\n{'='*60}")
        print(f"Testing: {tool}")
        print(f"Command: {cmd}")
        print(f"Input size: {input_size_mb:.2f} MB")
        print(f"Sample interval: {self.sample_interval}s")
        print(f"{'='*60}")
        
        start_time = time.time()
        
        try:
            # 启动进程
            process = subprocess.Popen(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            
            # 启动内存采样器
            sampler = MemorySampler(process.pid, self.sample_interval)
            sampler.start()
            
            # 等待进程完成
            try:
                stdout, stderr = process.communicate(timeout=timeout)
                profile.exit_code = process.returncode
            except subprocess.TimeoutExpired:
                process.kill()
                stdout, stderr = process.communicate()
                profile.exit_code = -1
                profile.error_message = f"Timeout after {timeout}s"
            
            # 停止采样
            profile.samples = sampler.stop()
            profile.execution_time_sec = round(time.time() - start_time, 2)
            
            # 分析内存数据
            if profile.samples:
                rss_values = [s.rss_mb for s in profile.samples]
                profile.peak_rss_mb = round(max(rss_values), 2)
                profile.avg_rss_mb = round(sum(rss_values) / len(rss_values), 2)
                profile.min_rss_mb = round(min(rss_values), 2)
                profile.max_rss_mb = round(max(rss_values), 2)
                
                # 分析稳定性
                is_stable, leak_rate, leak_detected = analyze_memory_stability(profile.samples)
                profile.is_stable = is_stable
                profile.leak_rate_mb_per_sec = leak_rate
                profile.memory_leak_detected = leak_detected
            
            # 打印结果摘要
            print(f"\nResults for {tool}:")
            print(f"  Execution time: {profile.execution_time_sec}s")
            print(f"  Peak RSS: {profile.peak_rss_mb} MB")
            print(f"  Avg RSS: {profile.avg_rss_mb} MB")
            print(f"  Memory range: {profile.min_rss_mb} - {profile.max_rss_mb} MB")
            print(f"  Samples collected: {len(profile.samples)}")
            print(f"  Memory stable: {'✓ Yes' if profile.is_stable else '✗ No'}")
            print(f"  Memory leak: {'✗ Detected' if profile.memory_leak_detected else '✓ None'}")
            if profile.memory_leak_detected:
                print(f"  Leak rate: {profile.leak_rate_mb_per_sec} MB/sec")
            
        except Exception as e:
            profile.error_message = str(e)
            profile.exit_code = -1
            print(f"Error running {tool}: {e}")
        
        return profile
    
    def run_all_tools(
        self,
        input_file: Path,
        tools: List[str] = None
    ) -> MemoryTestReport:
        """运行所有工具的内存测试"""
        
        if tools is None:
            tools = ["fastcrossmap", "crossmap"]  # 默认只测试支持 BAM 的工具
        
        input_size_mb = input_file.stat().st_size / (1024 * 1024)
        
        report = MemoryTestReport(
            timestamp=datetime.now().isoformat(),
            input_file=str(input_file),
            input_size_mb=round(input_size_mb, 2),
            sample_interval_sec=self.sample_interval
        )
        
        print(f"\n{'#'*60}")
        print(f"Memory Stability Test")
        print(f"{'#'*60}")
        print(f"Input file: {input_file}")
        print(f"Input size: {input_size_mb:.2f} MB")
        print(f"Tools to test: {', '.join(tools)}")
        print(f"Sample interval: {self.sample_interval}s")
        print(f"{'#'*60}")
        
        for tool in tools:
            profile = self.run_with_memory_sampling(tool, input_file)
            report.profiles.append(profile)
        
        return report


def export_report(report: MemoryTestReport, output_dir: Path):
    """导出测试报告"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # 导出 JSON（包含完整采样数据）
    json_path = output_dir / f"memory_test_{timestamp}.json"
    with open(json_path, 'w') as f:
        json.dump(asdict(report), f, indent=2)
    print(f"\nJSON report saved to: {json_path}")
    
    # 导出简化的 CSV（不含采样数据）
    csv_path = output_dir / f"memory_test_{timestamp}.csv"
    with open(csv_path, 'w') as f:
        f.write("tool,input_size_mb,execution_time_sec,peak_rss_mb,avg_rss_mb,")
        f.write("min_rss_mb,max_rss_mb,is_stable,memory_leak_detected,leak_rate_mb_per_sec\n")
        for p in report.profiles:
            f.write(f"{p.tool},{p.input_size_mb},{p.execution_time_sec},{p.peak_rss_mb},")
            f.write(f"{p.avg_rss_mb},{p.min_rss_mb},{p.max_rss_mb},{p.is_stable},")
            f.write(f"{p.memory_leak_detected},{p.leak_rate_mb_per_sec}\n")
    print(f"CSV report saved to: {csv_path}")
    
    # 导出采样数据（用于绘图）
    for profile in report.profiles:
        if profile.samples:
            samples_path = output_dir / f"memory_samples_{profile.tool}_{timestamp}.csv"
            with open(samples_path, 'w') as f:
                f.write("timestamp_sec,rss_mb,vms_mb\n")
                for s in profile.samples:
                    f.write(f"{s.timestamp},{s.rss_mb},{s.vms_mb}\n")
            print(f"Samples saved to: {samples_path}")
    
    # 保存最新版本
    json_latest = output_dir / "memory_test_latest.json"
    with open(json_latest, 'w') as f:
        json.dump(asdict(report), f, indent=2)


def print_summary(report: MemoryTestReport):
    """打印测试摘要"""
    print(f"\n{'='*80}")
    print("MEMORY STABILITY TEST SUMMARY")
    print(f"{'='*80}")
    print(f"{'Tool':<15} {'Time (s)':<12} {'Peak (MB)':<12} {'Avg (MB)':<12} {'Stable':<10} {'Leak':<10}")
    print(f"{'-'*80}")
    
    for p in report.profiles:
        stable_str = "✓ Yes" if p.is_stable else "✗ No"
        leak_str = "✗ Yes" if p.memory_leak_detected else "✓ No"
        print(f"{p.tool:<15} {p.execution_time_sec:<12.2f} {p.peak_rss_mb:<12.2f} {p.avg_rss_mb:<12.2f} {stable_str:<10} {leak_str:<10}")
    
    print(f"{'='*80}")
    
    # 检查 FastCrossMap 是否稳定
    fcm_profile = next((p for p in report.profiles if p.tool == "fastcrossmap"), None)
    if fcm_profile:
        if fcm_profile.is_stable and not fcm_profile.memory_leak_detected:
            print("\n✓ FastCrossMap memory curve is STABLE (no memory leak detected)")
        else:
            print("\n✗ FastCrossMap memory stability issue detected!")
            if fcm_profile.memory_leak_detected:
                print(f"  Memory leak rate: {fcm_profile.leak_rate_mb_per_sec} MB/sec")


def main():
    parser = argparse.ArgumentParser(
        description="Memory stability test for genome coordinate conversion tools"
    )
    parser.add_argument("--input", "-i", type=Path, required=True,
                        help="Input BAM file for testing")
    parser.add_argument("--chain", "-c", type=Path,
                        default=Path("ref/CrossMap/chain_files/human/GRCh37_to_GRCh38.chain.gz"),
                        help="Chain file for coordinate conversion")
    parser.add_argument("--output-dir", "-o", type=Path,
                        default=Path("results/memory_test"),
                        help="Output directory for results")
    parser.add_argument("--tools", "-t", nargs="+",
                        choices=["fastcrossmap", "crossmap", "fastremap"],
                        default=["fastcrossmap", "crossmap"],
                        help="Tools to test")
    parser.add_argument("--sample-interval", "-s", type=float, default=1.0,
                        help="Memory sampling interval in seconds")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads for FastCrossMap")
    
    args = parser.parse_args()
    
    # 验证输入文件
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    if not args.chain.exists():
        print(f"Error: Chain file not found: {args.chain}")
        sys.exit(1)
    
    # 确保 FastCrossMap 已编译
    if not Path("./target/release/fast-crossmap").exists():
        print("Building FastCrossMap...")
        subprocess.run(["cargo", "build", "--release"], check=True)
    
    # 运行测试
    tester = MemoryStabilityTester(
        chain_file=args.chain,
        output_dir=args.output_dir,
        sample_interval=args.sample_interval,
        threads=args.threads
    )
    
    report = tester.run_all_tools(args.input, args.tools)
    
    # 导出报告
    export_report(report, args.output_dir)
    
    # 打印摘要
    print_summary(report)


if __name__ == "__main__":
    main()