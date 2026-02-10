#!/bin/bash
# =============================================================================
# 12_run_all.sh - 一键运行所有论文复现脚本
# =============================================================================
#
# 用法: bash paper/12_run_all.sh [选项]
#
# 选项:
#   --skip-download    跳过数据下载 (如果数据已存在)
#   --skip-benchmark   跳过基准测试 (只生成图表)
#   --help             显示帮助信息
#
# 输出:
#   paper/results/     中间结果 JSON 文件
#   paper/figures/     输出图表 PDF/PNG 文件
#
# =============================================================================

set -e  # 遇到错误立即退出

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 默认选项
SKIP_DOWNLOAD=false
SKIP_BENCHMARK=false

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-download)
            SKIP_DOWNLOAD=true
            shift
            ;;
        --skip-benchmark)
            SKIP_BENCHMARK=true
            shift
            ;;
        --help)
            echo "用法: bash paper/12_run_all.sh [选项]"
            echo ""
            echo "选项:"
            echo "  --skip-download    跳过数据下载 (如果数据已存在)"
            echo "  --skip-benchmark   跳过基准测试 (只生成图表)"
            echo "  --help             显示帮助信息"
            exit 0
            ;;
        *)
            echo -e "${RED}未知选项: $1${NC}"
            exit 1
            ;;
    esac
done

# 打印带颜色的消息
print_step() {
    echo -e "\n${BLUE}============================================================${NC}"
    echo -e "${GREEN}$1${NC}"
    echo -e "${BLUE}============================================================${NC}\n"
}

print_warning() {
    echo -e "${YELLOW}警告: $1${NC}"
}

print_error() {
    echo -e "${RED}错误: $1${NC}"
}

# 检查依赖
check_dependencies() {
    print_step "检查依赖..."
    
    local missing=()
    
    # 检查 Python
    if ! command -v python3 &> /dev/null; then
        missing+=("python3")
    fi
    
    # 检查 Python 包
    python3 -c "import matplotlib" 2>/dev/null || missing+=("matplotlib")
    python3 -c "import numpy" 2>/dev/null || missing+=("numpy")
    python3 -c "import psutil" 2>/dev/null || missing+=("psutil")
    
    # 检查工具
    if ! command -v ./fast-crossmap-linux-x64/fast-crossmap &> /dev/null; then
        if [ ! -f "./fast-crossmap-linux-x64/fast-crossmap" ]; then
            missing+=("fast-crossmap (请先运行 cargo build --release)")
        fi
    fi
    
    if ! command -v CrossMap.py &> /dev/null; then
        print_warning "CrossMap.py 未找到，相关测试将跳过"
    fi
    
    if ! command -v liftOver &> /dev/null; then
        print_warning "liftOver 未找到，相关测试将跳过"
    fi
    
    if ! command -v FastRemap &> /dev/null; then
        print_warning "FastRemap 未找到，相关测试将跳过"
    fi
    
    if [ ${#missing[@]} -gt 0 ]; then
        print_error "缺少以下依赖: ${missing[*]}"
        echo "请安装后重试"
        exit 1
    fi
    
    echo -e "${GREEN}依赖检查通过${NC}"
}

# 创建目录
create_directories() {
    mkdir -p paper/data
    mkdir -p paper/results
    mkdir -p paper/figures
}

# 记录开始时间
START_TIME=$(date +%s)

echo -e "${BLUE}"
echo "╔════════════════════════════════════════════════════════════╗"
echo "║     FastCrossMap 论文复现脚本 - 一键运行                   ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# 检查依赖
check_dependencies

# 创建目录
create_directories

# =============================================================================
# Step 1: 下载数据
# =============================================================================
if [ "$SKIP_DOWNLOAD" = false ]; then
    print_step "Step 1/7: 下载公共数据"
    bash paper/01_download_data.sh
else
    print_warning "跳过数据下载 (--skip-download)"
fi

# =============================================================================
# Step 2-4: 基准测试
# =============================================================================
if [ "$SKIP_BENCHMARK" = false ]; then
    print_step "Step 2/7: BED 格式基准测试"
    python3 paper/02_benchmark_bed.py
    
    print_step "Step 2b/7: BED 多线程扩展测试"
    python3 paper/02b_benchmark_bed_multithread.py
    
    print_step "Step 3/7: BAM 格式基准测试"
    python3 paper/03_benchmark_bam.py
    
    print_step "Step 3b/7: BAM 多线程扩展测试"
    python3 paper/03b_benchmark_bam_multithread.py
    
    print_step "Step 4/7: 内存分析"
    python3 paper/05_memory_profile.py
    
    print_step "Step 5/7: 准确性分析"
    python3 paper/07_accuracy_analysis.py
    
    print_step "Step 6/7: 功能审计"
    python3 paper/09_feature_audit.py
else
    print_warning "跳过基准测试 (--skip-benchmark)"
fi

# =============================================================================
# Step 7: 生成图表
# =============================================================================
print_step "Step 7/7: 生成图表"

echo "生成 Figure 1: 性能对比图..."
python3 paper/04_plot_performance.py

echo "生成 Figure 2: 内存曲线图..."
python3 paper/06_plot_memory.py

echo "生成 Figure 3: 准确性分析图..."
python3 paper/08_plot_accuracy.py

echo "生成 Figure 4: 功能热力图..."
python3 paper/10_plot_features.py

echo "生成 Figure 5: 综合雷达图..."
python3 paper/11_plot_radar.py

# =============================================================================
# 完成
# =============================================================================
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
MINUTES=$((ELAPSED / 60))
SECONDS=$((ELAPSED % 60))

echo -e "\n${BLUE}"
echo "╔════════════════════════════════════════════════════════════╗"
echo "║                      运行完成!                             ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

echo -e "${GREEN}总耗时: ${MINUTES}分${SECONDS}秒${NC}"
echo ""
echo "输出文件:"
echo "  结果数据: paper/results/*.json"
echo "  图表文件: paper/figures/*.pdf"
echo ""
echo "生成的图表:"
ls -la paper/figures/*.pdf 2>/dev/null || echo "  (无 PDF 文件)"
echo ""
echo -e "${GREEN}论文复现完成!${NC}"
