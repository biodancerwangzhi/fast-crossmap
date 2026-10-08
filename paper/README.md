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
| Linux x86_64 (glibc ≥ 2.34) | [fcm-0.5.0-linux-x86_64.tar.gz](https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fcm-0.5.0-linux-x86_64.tar.gz) |
| Linux ARM64 | [fcm-0.5.0-linux-arm64.tar.gz](https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fcm-0.5.0-linux-arm64.tar.gz) |
| macOS Intel | [fcm-0.5.0-macos-x86_64.tar.gz](https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fcm-0.5.0-macos-x86_64.tar.gz) |
| macOS Apple Silicon | [fcm-0.5.0-macos-arm64.tar.gz](https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fcm-0.5.0-macos-arm64.tar.gz) |
| Windows x64 | [fcm-0.5.0-windows-x64.zip](https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fcm-0.5.0-windows-x64.zip) ⚠️ No BAM support |

```bash
# Linux example:
wget https://github.com/biodancerwangzhi/fast-crossmap/releases/latest/download/fcm-0.5.0-linux-x86_64.tar.gz
tar -xzf fcm-0.5.0-linux-x86_64.tar.gz
sha256sum -c fcm-0.5.0-linux-x86_64/SHA256SUMS   # optional
```

**Step 3: Run test**
```bash
# BED conversion test (all platforms)
./fcm-0.5.0-linux-x86_64/fast-crossmap bed \
    paper/sample_data/hg19ToHg38.over.chain.gz \
    paper/sample_data/sample.bed \
    output.bed

# Check output
head output.bed
wc -l output.bed          
wc -l output.bed.unmap    # Unmapped records

# SAM conversion test (Linux/macOS only)
./fcm-0.5.0-linux-x86_64/fast-crossmap bam \
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

## Multi-threading Scalability (~1 hour)

Measures FastCrossMap's own `-t` scaling per format. Prints `timings.tsv`
(raw runs) and `summary.tsv` (median / spread / speedup / parallel efficiency).
This does not change the paper's cross-tool comparison, which is `-t 1`
throughout (CrossMap has no multithreading).

```bash
FCM=/path/to/fast-crossmap bash paper/20_benchmark_threads.sh
FCM=/path/to/fast-crossmap bash paper/21_check_determinism.sh   # output must be byte-identical across -t
```

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

Single-thread vs multi-thread curves for BED/BAM, used for Figure 1(b)/(d):
```bash
python paper/02b_benchmark_bed_multithread.py   # -> paper/results/benchmark_bed_multithread.json
python paper/03b_benchmark_bam_multithread.py   # -> paper/results/benchmark_bam_multithread.json
```

All benchmark scripts expect the binary at `./target/release/fast-crossmap`
(build it with `cargo build --release`, or drop a pre-built binary there).
`02b` / `03b` / `07` also honour `FCM_BIN=/path/to/fast-crossmap`, and the two
shell scripts (`20`, `21`) honour `FCM=`.

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
