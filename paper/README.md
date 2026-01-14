# FastCrossMap Benchmark Reproducibility

## Quick Test (1 minute)

**Step 1: Clone repository**
```bash
git clone https://github.com/biodancerwangzhi/fast-crossmap.git
cd fast-crossmap
```

**Step 2: Download FastCrossMap binary**

| Platform | Download |
|----------|----------|
| Linux x64 | [fast-crossmap-linux-x64.tar.gz](https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fast-crossmap-linux-x64.tar.gz) |
| Linux ARM64 | [fast-crossmap-linux-arm64.tar.gz](https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fast-crossmap-linux-arm64.tar.gz) |
| macOS Intel | [fast-crossmap-macos-x64.tar.gz](https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fast-crossmap-macos-x64.tar.gz) |
| macOS Apple Silicon | [fast-crossmap-macos-arm64.tar.gz](https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fast-crossmap-macos-arm64.tar.gz) |
| Windows x64 | [fast-crossmap-windows-x64.zip](https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fast-crossmap-windows-x64.zip) ⚠️ No BAM support |

```bash
# Linux example:
wget https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fast-crossmap-linux-x64.tar.gz
tar -xzf fast-crossmap-linux-x64.tar.gz
```

**Step 3: Run test**
```bash
# BED conversion test (all platforms)
./fast-crossmap-linux-x64/fast-crossmap bed \
    paper/sample_data/hg19ToHg38.over.chain.gz \
    paper/sample_data/sample.bed \
    output.bed

# Check output
head output.bed
wc -l output.bed          
wc -l output.bed.unmap    # Unmapped records

# SAM conversion test (Linux/macOS only)
./fast-crossmap-linux-x64/fast-crossmap bam \
    paper/sample_data/hg19ToHg38.over.chain.gz \
    paper/sample_data/sample.sam \
    output.sam
```

**Expected output:**
- `output.bed`: successfully converted BED records
- `output.bed.unmap`: unmapped records
- `output.sam`: Converted SAM file (Linux/macOS only)

---

## Sample Data (included in repository)

| File | Size | Description |
|------|------|-------------|
| `paper/sample_data/sample.bed` | 500 KB | 10,000 BED records (hg19) |
| `paper/sample_data/sample.sam` | 150 KB | 1,000 SAM reads (hg19) |
| `paper/sample_data/hg19ToHg38.over.chain.gz` | 1 MB | UCSC chain file |

---

## Full Benchmark (~30 min)

Requires Linux and conda environment.

**Step 1: Install comparison tools**
```bash
conda install -c bioconda crossmap ucsc-liftover fastremap-bio
pip install matplotlib numpy pandas seaborn psutil
```

**Step 2: Download ENCODE data (~2GB)**
```bash
bash paper/01_download_data.sh
```

**Step 3: Run benchmarks**
```bash
# All at once
bash paper/12_run_all.sh

# Or step by step
python paper/02_benchmark_bed.py        # BED benchmark
python paper/03_benchmark_bam.py        # BAM benchmark
python paper/05_memory_profile.py       # Memory profiling
python paper/07_accuracy_analysis.py    # Accuracy validation
```

---

## Expected Results

### Performance (BED, 296,898 records)

| Tool | Threads | Time (s) | Speedup |
|------|---------|----------|---------|
| FastCrossMap | 1 | 0.35 | 11x |
| FastCrossMap | 4 | 0.12 | 32x |
| CrossMap | 1 | 3.81 | 1x |

### Memory (BAM, 1.3 GB)

| Tool | Peak Memory |
|------|-------------|
| FastCrossMap | 18 MB |
| CrossMap | 1,100 MB |

### Accuracy

FastCrossMap produces bit-exact identical output to CrossMap (99.8% identical to liftOver, 0.2% partial mappings handled identically by both tools).

---

## Data Sources

| Data | Source | Accession |
|------|--------|-----------|
| BED | ENCODE | [ENCFF001WBV](https://www.encodeproject.org/files/ENCFF001WBV/) |
| BAM | ENCODE | [ENCFF000PED](https://www.encodeproject.org/files/ENCFF000PED/) |
| Chain | UCSC | [hg19ToHg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/) |
