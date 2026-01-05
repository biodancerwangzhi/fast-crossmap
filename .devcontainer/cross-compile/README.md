# FastCrossMap 交叉编译指南

本目录包含用于构建多平台二进制文件的 Docker 配置。

## 支持的目标平台

| 平台 | 架构 | 目标三元组 | 状态 |
|------|------|-----------|------|
| Linux | x86_64 | x86_64-unknown-linux-gnu | ✅ 开箱即用 |
| Linux | aarch64 | aarch64-unknown-linux-gnu | ✅ 开箱即用 |
| Windows | x86_64 | x86_64-pc-windows-gnu | ✅ 开箱即用 |
| macOS | x86_64 | x86_64-apple-darwin | ⚠️ 需要 SDK |
| macOS | aarch64 | aarch64-apple-darwin | ⚠️ 需要 SDK |

## 快速开始

### 1. 构建 Docker 镜像

```bash
cd .devcontainer/cross-compile
docker build -t fastcrossmap-cross .
```

### 2. 运行容器并构建

```bash
# 挂载项目目录并构建所有平台
docker run --rm -v $(pwd):/workspace fastcrossmap-cross build-all.sh all

# 只构建 Linux 和 Windows (不需要 macOS SDK)
docker run --rm -v $(pwd):/workspace fastcrossmap-cross bash -c "
    build-all.sh linux-x64 && \
    build-all.sh linux-arm64 && \
    build-all.sh windows
"
```

### 3. 输出文件

构建完成后，二进制文件位于 `release/` 目录：

```
release/
├── fast-crossmap-x.x.x-linux-x64.tar.gz
├── fast-crossmap-x.x.x-linux-arm64.tar.gz
├── fast-crossmap-x.x.x-windows-x64.zip
├── fast-crossmap-x.x.x-macos-x64.tar.gz    (需要 SDK)
└── fast-crossmap-x.x.x-macos-arm64.tar.gz  (需要 SDK)
```

## macOS 交叉编译配置

由于 Apple 许可证限制，macOS SDK 不能直接包含在 Docker 镜像中。需要手动配置：

### 方法 1: 使用 osxcross (推荐)

1. 从 Apple Developer 下载 Xcode 或 Command Line Tools
2. 提取 MacOSX SDK (例如 MacOSX12.3.sdk)
3. 打包 SDK:
   ```bash
   tar -cJf MacOSX12.3.sdk.tar.xz MacOSX12.3.sdk
   ```
4. 将 SDK 复制到容器中:
   ```bash
   docker cp MacOSX12.3.sdk.tar.xz <container_id>:/opt/osxcross/tarballs/
   ```
5. 在容器中构建 osxcross:
   ```bash
   cd /opt/osxcross
   UNATTENDED=1 ./build.sh
   ```

### 方法 2: 使用 GitHub Actions (推荐用于 CI/CD)

参考 `.github/workflows/release.yml` 使用 macOS runner 进行原生编译。

## 单独构建某个平台

```bash
# Linux x86_64
docker run --rm -v $(pwd):/workspace fastcrossmap-cross build-all.sh linux-x64

# Linux ARM64
docker run --rm -v $(pwd):/workspace fastcrossmap-cross build-all.sh linux-arm64

# Windows
docker run --rm -v $(pwd):/workspace fastcrossmap-cross build-all.sh windows

# macOS Intel (需要 SDK)
docker run --rm -v $(pwd):/workspace fastcrossmap-cross build-all.sh macos-x64

# macOS Apple Silicon (需要 SDK)
docker run --rm -v $(pwd):/workspace fastcrossmap-cross build-all.sh macos-arm64
```

## 交互式开发

```bash
# 进入容器进行调试
docker run -it --rm -v $(pwd):/workspace fastcrossmap-cross bash

# 在容器内手动构建
cargo build --release --target x86_64-pc-windows-gnu
```

## 常见问题

### Q: Windows 构建失败，提示找不到 mingw

确保 Dockerfile 中安装了 `mingw-w64`:
```dockerfile
RUN apt-get install -y mingw-w64
```

### Q: Linux ARM64 构建失败

确保安装了交叉编译工具链:
```dockerfile
RUN apt-get install -y gcc-aarch64-linux-gnu g++-aarch64-linux-gnu
```

### Q: macOS 构建提示 osxcross 未配置

需要手动配置 macOS SDK，参考上面的 "macOS 交叉编译配置" 部分。

### Q: 如何验证构建的二进制文件？

```bash
# 检查文件类型
file release/fast-crossmap-*-linux-x64/fast-crossmap
# 输出: ELF 64-bit LSB executable, x86-64

file release/fast-crossmap-*-windows-x64/fast-crossmap.exe
# 输出: PE32+ executable (console) x86-64

file release/fast-crossmap-*-macos-x64/fast-crossmap
# 输出: Mach-O 64-bit x86_64 executable
```

## 依赖说明

| 目标平台 | 依赖 |
|----------|------|
| Linux x64 | gcc, glibc |
| Linux ARM64 | gcc-aarch64-linux-gnu |
| Windows | mingw-w64 |
| macOS | osxcross + MacOSX SDK |
