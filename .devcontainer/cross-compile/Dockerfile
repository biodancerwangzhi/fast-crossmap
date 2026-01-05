# =============================================================================
# FastCrossMap 交叉编译容器
# 支持编译 Linux (x86_64, aarch64), Windows (x86_64), macOS (x86_64, aarch64)
# =============================================================================
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# =============================================================================
# 基础工具安装
# =============================================================================
RUN apt-get update && apt-get install -y \
    build-essential \
    pkg-config \
    cmake \
    git \
    curl \
    wget \
    vim \
    ca-certificates \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    # Windows 交叉编译工具链
    mingw-w64 \
    # macOS 交叉编译依赖
    clang \
    llvm \
    libxml2-dev \
    libz-dev \
    # 其他工具
    zip \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# =============================================================================
# 安装 Rust 和交叉编译目标
# =============================================================================
ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain stable

# 添加交叉编译目标
RUN rustup target add \
    x86_64-unknown-linux-gnu \
    aarch64-unknown-linux-gnu \
    x86_64-pc-windows-gnu \
    x86_64-apple-darwin \
    aarch64-apple-darwin

# =============================================================================
# 安装 osxcross (macOS 交叉编译工具链)
# =============================================================================
# 注意: macOS SDK 需要从 Apple 官方获取，这里提供安装框架
# 用户需要自行下载 MacOSX SDK 并放入容器

# 安装 osxcross 依赖
RUN apt-get update && apt-get install -y \
    libssl-dev \
    lzma-dev \
    libxml2-dev \
    xz-utils \
    cpio \
    libbz2-dev \
    && rm -rf /var/lib/apt/lists/*

# 克隆 osxcross
RUN git clone https://github.com/tpoechtrager/osxcross /opt/osxcross

# 创建 SDK 目录 (用户需要手动添加 SDK)
RUN mkdir -p /opt/osxcross/tarballs

# =============================================================================
# 配置 Cargo 交叉编译
# =============================================================================
RUN mkdir -p /root/.cargo && \
    echo '[target.x86_64-pc-windows-gnu]' >> /root/.cargo/config.toml && \
    echo 'linker = "x86_64-w64-mingw32-gcc"' >> /root/.cargo/config.toml && \
    echo '' >> /root/.cargo/config.toml && \
    echo '[target.x86_64-apple-darwin]' >> /root/.cargo/config.toml && \
    echo 'linker = "x86_64-apple-darwin21.4-clang"' >> /root/.cargo/config.toml && \
    echo '' >> /root/.cargo/config.toml && \
    echo '[target.aarch64-apple-darwin]' >> /root/.cargo/config.toml && \
    echo 'linker = "aarch64-apple-darwin21.4-clang"' >> /root/.cargo/config.toml && \
    echo '' >> /root/.cargo/config.toml && \
    echo '[target.aarch64-unknown-linux-gnu]' >> /root/.cargo/config.toml && \
    echo 'linker = "aarch64-linux-gnu-gcc"' >> /root/.cargo/config.toml

# 安装 aarch64 Linux 交叉编译工具链
RUN apt-get update && apt-get install -y \
    gcc-aarch64-linux-gnu \
    g++-aarch64-linux-gnu \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# =============================================================================
# 构建脚本入口
# =============================================================================
COPY build-all.sh /usr/local/bin/build-all.sh
RUN chmod +x /usr/local/bin/build-all.sh

CMD ["/bin/bash"]
