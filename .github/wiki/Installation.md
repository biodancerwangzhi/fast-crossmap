# Installation

## System Requirements

- **OS**: Linux, macOS, or Windows
- **Memory**: 512 MB minimum (2 GB+ recommended for VCF/GVCF/MAF with reference genome)
- **Note**: Windows builds do not support BAM/SAM/CRAM format

---

## Method 1: Pre-built Binaries (Recommended)

Download from [Releases](https://github.com/biodancerwangzhi/fast-crossmap/releases).

### Linux
```bash
# Download and extract
wget https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fast-crossmap-linux-x64.tar.gz
tar -xzf fast-crossmap-linux-x64.tar.gz
chmod +x fast-crossmap

# Test
./fast-crossmap --version
```

### macOS (Apple Silicon)
```bash
wget https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fast-crossmap-macos-arm64.tar.gz
tar -xzf fast-crossmap-macos-arm64.tar.gz
chmod +x fast-crossmap
./fast-crossmap --version
```

### Windows
1. Download `fast-crossmap-windows-x64.zip` from [Releases](https://github.com/biodancerwangzhi/fast-crossmap/releases)
2. Extract the ZIP file
3. Run `fast-crossmap.exe --version`

### Optional: Add to PATH

For frequent use, add FastCrossMap to your system PATH:

```bash
# Option 1: Move to system directory (Linux/macOS)
sudo mv fast-crossmap /usr/local/bin/

# Option 2: Add current directory to PATH (temporary)
export PATH="$PWD:$PATH"

# Option 3: Add to ~/.bashrc or ~/.zshrc (permanent)
echo 'export PATH="/path/to/fast-crossmap:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

---

## Method 2: Build from Source

```bash
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env

# Build
git clone https://github.com/biodancerwangzhi/fast-crossmap.git
cd fast-crossmap
cargo build --release
./target/release/fast-crossmap --version
```

---

## Download Chain Files

Chain files are required for coordinate conversion. Download from UCSC:

### Human Genome
```bash
# hg19 → hg38
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

# hg38 → hg19
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
```

### Mouse Genome
```bash
# mm10 → mm39
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/liftOver/mm10ToMm39.over.chain.gz

# mm39 → mm10
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm39/liftOver/mm39ToMm10.over.chain.gz
```

### Reference Genomes (for VCF/GVCF/MAF)
```bash
# hg38 reference (required for VCF/GVCF/MAF conversion)
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

Browse all chain files: https://hgdownload.cse.ucsc.edu/downloads.html

---

## Next Steps

- [Quick Start](QuickStart) - Get started in 5 minutes
- [Advanced Usage](Advanced) - Detailed format guides
