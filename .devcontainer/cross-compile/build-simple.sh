#!/bin/bash
# =============================================================================
# FastCrossMap 简化构建脚本
# 使用 cross 工具进行交叉编译 (需要 Docker)
# =============================================================================
#
# 用法: ./build-simple.sh [目标]
#
# 目标:
#   linux       Linux x86_64 (静态链接 musl)
#   linux-gnu   Linux x86_64 (动态链接 glibc)
#   windows     Windows x86_64
#   all         构建所有 (默认)
#
# 前提条件:
#   - 安装 Rust: curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
#   - 安装 cross: cargo install cross --git https://github.com/cross-rs/cross
#   - 安装 Docker 并确保当前用户有权限运行
#
# =============================================================================

set -e

# 颜色
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

# 检查 cross 是否安装
check_cross() {
    if ! command -v cross &> /dev/null; then
        echo -e "${YELLOW}cross 未安装，正在安装...${NC}"
        cargo install cross --git https://github.com/cross-rs/cross
    fi
}

# 检查 Docker
check_docker() {
    if ! docker info &> /dev/null; then
        echo -e "${RED}错误: Docker 未运行或无权限${NC}"
        echo "请确保 Docker 已启动，并且当前用户在 docker 组中"
        exit 1
    fi
}

# 获取版本号
get_version() {
    grep '^version' Cargo.toml | head -1 | sed 's/.*"\(.*\)".*/\1/'
}

VERSION=$(get_version)
OUTPUT_DIR="release"
mkdir -p "$OUTPUT_DIR"

echo -e "${BLUE}FastCrossMap v${VERSION} 构建脚本${NC}"
echo ""

# =============================================================================
# 构建函数
# =============================================================================

build_linux_musl() {
    echo -e "${GREEN}构建 Linux x86_64 (musl 静态链接)...${NC}"
    cross build --release --target x86_64-unknown-linux-musl
    
    mkdir -p "$OUTPUT_DIR/fast-crossmap-${VERSION}-linux-x64"
    cp target/x86_64-unknown-linux-musl/release/fast-crossmap "$OUTPUT_DIR/fast-crossmap-${VERSION}-linux-x64/"
    
    cd "$OUTPUT_DIR"
    tar -czvf "fast-crossmap-${VERSION}-linux-x64.tar.gz" "fast-crossmap-${VERSION}-linux-x64"
    rm -rf "fast-crossmap-${VERSION}-linux-x64"
    cd ..
    
    echo -e "${GREEN}✓ Linux x86_64 完成${NC}"
}

build_linux_gnu() {
    echo -e "${GREEN}构建 Linux x86_64 (glibc 动态链接)...${NC}"
    cross build --release --target x86_64-unknown-linux-gnu
    
    mkdir -p "$OUTPUT_DIR/fast-crossmap-${VERSION}-linux-x64-gnu"
    cp target/x86_64-unknown-linux-gnu/release/fast-crossmap "$OUTPUT_DIR/fast-crossmap-${VERSION}-linux-x64-gnu/"
    
    cd "$OUTPUT_DIR"
    tar -czvf "fast-crossmap-${VERSION}-linux-x64-gnu.tar.gz" "fast-crossmap-${VERSION}-linux-x64-gnu"
    rm -rf "fast-crossmap-${VERSION}-linux-x64-gnu"
    cd ..
    
    echo -e "${GREEN}✓ Linux x86_64 (glibc) 完成${NC}"
}

build_windows() {
    echo -e "${GREEN}构建 Windows x86_64...${NC}"
    cross build --release --target x86_64-pc-windows-gnu
    
    mkdir -p "$OUTPUT_DIR/fast-crossmap-${VERSION}-windows-x64"
    cp target/x86_64-pc-windows-gnu/release/fast-crossmap.exe "$OUTPUT_DIR/fast-crossmap-${VERSION}-windows-x64/"
    
    cd "$OUTPUT_DIR"
    zip -r "fast-crossmap-${VERSION}-windows-x64.zip" "fast-crossmap-${VERSION}-windows-x64"
    rm -rf "fast-crossmap-${VERSION}-windows-x64"
    cd ..
    
    echo -e "${GREEN}✓ Windows x86_64 完成${NC}"
}

# =============================================================================
# 主程序
# =============================================================================

TARGET="${1:-all}"

check_cross
check_docker

case "$TARGET" in
    linux)
        build_linux_musl
        ;;
    linux-gnu)
        build_linux_gnu
        ;;
    windows)
        build_windows
        ;;
    all)
        build_linux_musl
        build_windows
        ;;
    *)
        echo -e "${RED}未知目标: $TARGET${NC}"
        echo "用法: ./build-simple.sh [linux|linux-gnu|windows|all]"
        exit 1
        ;;
esac

echo ""
echo -e "${BLUE}============================================================${NC}"
echo -e "${GREEN}构建完成!${NC}"
echo -e "${BLUE}============================================================${NC}"
echo "输出目录: $OUTPUT_DIR/"
ls -la "$OUTPUT_DIR/"
