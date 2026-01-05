#!/bin/bash
# =============================================================================
# FastCrossMap 多平台构建脚本
# =============================================================================
#
# 用法: build-all.sh [目标平台]
#
# 目标平台:
#   linux-x64      Linux x86_64
#   linux-arm64    Linux aarch64
#   windows        Windows x86_64
#   macos-x64      macOS x86_64 (Intel)
#   macos-arm64    macOS aarch64 (Apple Silicon)
#   all            构建所有平台 (默认)
#
# =============================================================================

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# 输出目录
OUTPUT_DIR="${OUTPUT_DIR:-/workspace/release}"
mkdir -p "$OUTPUT_DIR"

# 版本号 (从 Cargo.toml 读取)
VERSION=$(grep '^version' /workspace/Cargo.toml | head -1 | sed 's/.*"\(.*\)".*/\1/')
echo -e "${BLUE}FastCrossMap 版本: ${VERSION}${NC}"

# =============================================================================
# 构建函数
# =============================================================================

build_linux_x64() {
    echo -e "\n${GREEN}构建 Linux x86_64...${NC}"
    cargo build --release --target x86_64-unknown-linux-gnu
    
    local out="$OUTPUT_DIR/fast-crossmap-${VERSION}-linux-x64"
    mkdir -p "$out"
    cp target/x86_64-unknown-linux-gnu/release/fast-crossmap "$out/"
    
    # 创建压缩包
    cd "$OUTPUT_DIR"
    tar -czvf "fast-crossmap-${VERSION}-linux-x64.tar.gz" "fast-crossmap-${VERSION}-linux-x64"
    rm -rf "fast-crossmap-${VERSION}-linux-x64"
    
    echo -e "${GREEN}✓ Linux x86_64 构建完成${NC}"
}

build_linux_arm64() {
    echo -e "\n${GREEN}构建 Linux aarch64...${NC}"
    
    # 检查交叉编译工具链
    if ! command -v aarch64-linux-gnu-gcc &> /dev/null; then
        echo -e "${RED}错误: 缺少 aarch64-linux-gnu-gcc${NC}"
        return 1
    fi
    
    cargo build --release --target aarch64-unknown-linux-gnu
    
    local out="$OUTPUT_DIR/fast-crossmap-${VERSION}-linux-arm64"
    mkdir -p "$out"
    cp target/aarch64-unknown-linux-gnu/release/fast-crossmap "$out/"
    
    cd "$OUTPUT_DIR"
    tar -czvf "fast-crossmap-${VERSION}-linux-arm64.tar.gz" "fast-crossmap-${VERSION}-linux-arm64"
    rm -rf "fast-crossmap-${VERSION}-linux-arm64"
    
    echo -e "${GREEN}✓ Linux aarch64 构建完成${NC}"
}

build_windows() {
    echo -e "\n${GREEN}构建 Windows x86_64...${NC}"
    
    # 检查 MinGW 工具链
    if ! command -v x86_64-w64-mingw32-gcc &> /dev/null; then
        echo -e "${RED}错误: 缺少 mingw-w64${NC}"
        return 1
    fi
    
    cargo build --release --target x86_64-pc-windows-gnu
    
    local out="$OUTPUT_DIR/fast-crossmap-${VERSION}-windows-x64"
    mkdir -p "$out"
    cp target/x86_64-pc-windows-gnu/release/fast-crossmap.exe "$out/"
    
    cd "$OUTPUT_DIR"
    zip -r "fast-crossmap-${VERSION}-windows-x64.zip" "fast-crossmap-${VERSION}-windows-x64"
    rm -rf "fast-crossmap-${VERSION}-windows-x64"
    
    echo -e "${GREEN}✓ Windows x86_64 构建完成${NC}"
}

build_macos_x64() {
    echo -e "\n${GREEN}构建 macOS x86_64 (Intel)...${NC}"
    
    # 检查 osxcross 是否配置
    if [ ! -f /opt/osxcross/target/bin/x86_64-apple-darwin21.4-clang ]; then
        echo -e "${YELLOW}警告: osxcross 未配置，跳过 macOS x86_64 构建${NC}"
        echo -e "${YELLOW}请参考 README 配置 macOS SDK${NC}"
        return 0
    fi
    
    export PATH="/opt/osxcross/target/bin:$PATH"
    export CC=x86_64-apple-darwin21.4-clang
    export CXX=x86_64-apple-darwin21.4-clang++
    
    cargo build --release --target x86_64-apple-darwin
    
    local out="$OUTPUT_DIR/fast-crossmap-${VERSION}-macos-x64"
    mkdir -p "$out"
    cp target/x86_64-apple-darwin/release/fast-crossmap "$out/"
    
    cd "$OUTPUT_DIR"
    tar -czvf "fast-crossmap-${VERSION}-macos-x64.tar.gz" "fast-crossmap-${VERSION}-macos-x64"
    rm -rf "fast-crossmap-${VERSION}-macos-x64"
    
    echo -e "${GREEN}✓ macOS x86_64 构建完成${NC}"
}

build_macos_arm64() {
    echo -e "\n${GREEN}构建 macOS aarch64 (Apple Silicon)...${NC}"
    
    # 检查 osxcross 是否配置
    if [ ! -f /opt/osxcross/target/bin/aarch64-apple-darwin21.4-clang ]; then
        echo -e "${YELLOW}警告: osxcross 未配置，跳过 macOS aarch64 构建${NC}"
        echo -e "${YELLOW}请参考 README 配置 macOS SDK${NC}"
        return 0
    fi
    
    export PATH="/opt/osxcross/target/bin:$PATH"
    export CC=aarch64-apple-darwin21.4-clang
    export CXX=aarch64-apple-darwin21.4-clang++
    
    cargo build --release --target aarch64-apple-darwin
    
    local out="$OUTPUT_DIR/fast-crossmap-${VERSION}-macos-arm64"
    mkdir -p "$out"
    cp target/aarch64-apple-darwin/release/fast-crossmap "$out/"
    
    cd "$OUTPUT_DIR"
    tar -czvf "fast-crossmap-${VERSION}-macos-arm64.tar.gz" "fast-crossmap-${VERSION}-macos-arm64"
    rm -rf "fast-crossmap-${VERSION}-macos-arm64"
    
    echo -e "${GREEN}✓ macOS aarch64 构建完成${NC}"
}

# =============================================================================
# 主程序
# =============================================================================

cd /workspace

TARGET="${1:-all}"

case "$TARGET" in
    linux-x64)
        build_linux_x64
        ;;
    linux-arm64)
        build_linux_arm64
        ;;
    windows)
        build_windows
        ;;
    macos-x64)
        build_macos_x64
        ;;
    macos-arm64)
        build_macos_arm64
        ;;
    all)
        echo -e "${BLUE}构建所有平台...${NC}"
        build_linux_x64
        build_linux_arm64
        build_windows
        build_macos_x64
        build_macos_arm64
        ;;
    *)
        echo -e "${RED}未知目标: $TARGET${NC}"
        echo "用法: build-all.sh [linux-x64|linux-arm64|windows|macos-x64|macos-arm64|all]"
        exit 1
        ;;
esac

echo -e "\n${BLUE}============================================================${NC}"
echo -e "${GREEN}构建完成!${NC}"
echo -e "${BLUE}============================================================${NC}"
echo "输出目录: $OUTPUT_DIR"
ls -la "$OUTPUT_DIR"
