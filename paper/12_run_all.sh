#!/bin/bash
# =============================================================================
# 12_run_all.sh - Run all paper reproduction scripts
# =============================================================================
#
# Usage: bash paper/12_run_all.sh [options]
#
# Options:
#   --skip-download    Skip data download (if data already exists)
#   --skip-benchmark   Skip benchmarks (only generate figures)
#   --help             Show help information
#
# Output:
#   paper/results/     Intermediate result JSON files
#   paper/figures/     Output figure PDF/PNG files
#
# =============================================================================

set -e  # Exit immediately on error

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default options
SKIP_DOWNLOAD=false
SKIP_BENCHMARK=false

# Parse command line arguments
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
            echo "Usage: bash paper/12_run_all.sh [options]"
            echo ""
            echo "Options:"
            echo "  --skip-download    Skip data download (if data already exists)"
            echo "  --skip-benchmark   Skip benchmarks (only generate figures)"
            echo "  --help             Show help information"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            exit 1
            ;;
    esac
done

# Print colored messages
print_step() {
    echo -e "\n${BLUE}============================================================${NC}"
    echo -e "${GREEN}$1${NC}"
    echo -e "${BLUE}============================================================${NC}\n"
}

print_warning() {
    echo -e "${YELLOW}Warning: $1${NC}"
}

print_error() {
    echo -e "${RED}Error: $1${NC}"
}

# Check dependencies
check_dependencies() {
    print_step "Checking dependencies..."
    
    local missing=()
    
    # Check Python
    if ! command -v python3 &> /dev/null; then
        missing+=("python3")
    fi
    
    # Check Python packages
    python3 -c "import matplotlib" 2>/dev/null || missing+=("matplotlib")
    python3 -c "import numpy" 2>/dev/null || missing+=("numpy")
    python3 -c "import psutil" 2>/dev/null || missing+=("psutil")
    
    # Check tools
    if ! command -v ./fast-crossmap-linux-x64/fast-crossmap &> /dev/null; then
        if [ ! -f "./fast-crossmap-linux-x64/fast-crossmap" ]; then
            missing+=("fast-crossmap (please run cargo build --release first)")
        fi
    fi
    
    if ! command -v CrossMap.py &> /dev/null; then
        print_warning "CrossMap.py not found, related tests will be skipped"
    fi
    
    if ! command -v liftOver &> /dev/null; then
        print_warning "liftOver not found, related tests will be skipped"
    fi
    
    if ! command -v FastRemap &> /dev/null; then
        print_warning "FastRemap not found, related tests will be skipped"
    fi
    
    if [ ${#missing[@]} -gt 0 ]; then
        print_error "Missing dependencies: ${missing[*]}"
        echo "Please install and retry"
        exit 1
    fi
    
    echo -e "${GREEN}Dependency check passed${NC}"
}

# Create directories
create_directories() {
    mkdir -p paper/data
    mkdir -p paper/results
    mkdir -p paper/figures
}

# Record start time
START_TIME=$(date +%s)

echo -e "${BLUE}"
echo "╔════════════════════════════════════════════════════════════╗"
echo "║     FastCrossMap Paper Reproduction - Run All              ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# Check dependencies
check_dependencies

# Create directories
create_directories

# =============================================================================
# Step 1: Download data
# =============================================================================
if [ "$SKIP_DOWNLOAD" = false ]; then
    print_step "Step 1/7: Download public data"
    bash paper/01_download_data.sh
else
    print_warning "Skipping data download (--skip-download)"
fi

# =============================================================================
# Step 2-4: Benchmarks
# =============================================================================
if [ "$SKIP_BENCHMARK" = false ]; then
    print_step "Step 2/7: BED format benchmark"
    python3 paper/02_benchmark_bed.py
    
    print_step "Step 2b/7: BED multi-thread scaling test"
    python3 paper/02b_benchmark_bed_multithread.py
    
    print_step "Step 3/7: BAM format benchmark"
    python3 paper/03_benchmark_bam.py
    
    print_step "Step 3b/7: BAM multi-thread scaling test"
    python3 paper/03b_benchmark_bam_multithread.py
    
    print_step "Step 4/7: Memory profiling"
    python3 paper/05_memory_profile.py
    
    print_step "Step 5/7: Accuracy analysis"
    python3 paper/07_accuracy_analysis.py
    
    print_step "Step 6/7: Feature audit"
    python3 paper/09_feature_audit.py
else
    print_warning "Skipping benchmarks (--skip-benchmark)"
fi

# =============================================================================
# Step 7: Generate figures
# =============================================================================
print_step "Step 7/7: Generating figures"

echo "Generating Figure 1: Performance comparison..."
python3 paper/04_plot_performance.py

echo "Generating Figure 2: Memory curves..."
python3 paper/06_plot_memory.py

echo "Generating Figure 3: Accuracy analysis..."
python3 paper/08_plot_accuracy.py

echo "Generating Figure 4: Feature heatmap..."
python3 paper/10_plot_features.py

echo "Generating Figure 5: Comprehensive radar chart..."
python3 paper/11_plot_radar.py

# =============================================================================
# Complete
# =============================================================================
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
MINUTES=$((ELAPSED / 60))
SECONDS=$((ELAPSED % 60))

echo -e "\n${BLUE}"
echo "╔════════════════════════════════════════════════════════════╗"
echo "║                      Run Complete!                         ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

echo -e "${GREEN}Total time: ${MINUTES}m${SECONDS}s${NC}"
echo ""
echo "Output files:"
echo "  Result data: paper/results/*.json"
echo "  Figure files: paper/figures/*.pdf"
echo ""
echo "Generated figures:"
ls -la paper/figures/*.pdf 2>/dev/null || echo "  (No PDF files)"
echo ""
echo -e "${GREEN}Paper reproduction complete!${NC}"
