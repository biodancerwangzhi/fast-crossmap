# FastCrossMap Benchmark Reproducibility

## Quick Test (1 minute)

**Step 1: Clone repository**
```bash
git clone https://github.com/biodancerwangzhi/fast-crossmap.git
cd fast-crossmap
```

**Step 2: Download FastCrossMap (Linux)**

The scripts in this repository are Linux-only (`bash`, `-t` pinning, `samtools` /
`bwa` pipelines), so only the Linux build is needed here. macOS and Windows builds
are published on the [Releases page](https://github.com/biodancerwangzhi/fast-crossmap/releases/latest).

The release asset name embeds the version (`fcm-<version>-linux-x86_64.tar.gz`),
so a URL that pins `0.5.0` goes stale the moment a new release is cut -- under
`releases/latest/` the asset no longer exists and the download 404s. Ask the
GitHub API which asset the newest release carries instead, so the command never
has to be edited again:

```bash
# x86_64, glibc >= 2.34 -- always the newest release, no version to maintain
curl -sSL https://api.github.com/repos/biodancerwangzhi/fast-crossmap/releases/latest \
    | grep -o 'https://github.com/[^"]*fcm-[^"]*-linux-x86_64\.tar\.gz' \
    | head -1 | xargs -I{} curl -fLO {}

# with the GitHub CLI, same idea and no version either:
#   gh release download --repo biodancerwangzhi/fast-crossmap \
#       --pattern 'fcm-*-linux-x86_64.tar.gz' --dir . --clobber

tar -xzf fcm-*-linux-x86_64.tar.gz
sha256sum -c fcm-*-linux-x86_64/SHA256SUMS   # optional
```

Other Linux builds, same version-agnostic pattern:

| Build | `--pattern` for `gh release download` |
|-------|----------------------------------------|
| x86_64, glibc >= 2.34 | `fcm-*-linux-x86_64.tar.gz` |
| x86_64, statically linked (musl, no glibc) | `fcm-*-linux-x86_64-static.tar.gz` |
| ARM64 | `fcm-*-linux-arm64.tar.gz` |

**Step 3: Run test**

The tarball unpacks to `fcm-<version>-linux-x86_64/`; `FCM` below picks it up
whatever the version happens to be.

```bash
# Point FCM at the binary just unpacked (matches fcm-0.5.0-linux-x86_64, fcm-0.6.0-...)
FCM=$(ls -d fcm-*-linux-x86_64)/fast-crossmap

# BED conversion test
"$FCM" bed \
    paper/sample_data/hg19ToHg38.over.chain.gz \
    paper/sample_data/input/sample.bed \
    output.bed

# SAM conversion test (via the `bam` subcommand)
"$FCM" bam \
    paper/sample_data/hg19ToHg38.over.chain.gz \
    paper/sample_data/input/sample.sam \
    output.sam
```

Each run reports `100` records, all mapped. `paper/sample_data/output/` holds the
committed FastCrossMap output for those same inputs, so the test also checks the
binary:

```bash
diff output.bed paper/sample_data/output/sample_hg38.bed   # both print nothing
diff output.sam paper/sample_data/output/sample_hg38.sam   # on a correct build
```

---

## Sample Data (included in repository)

| Path | Description |
|------|-------------|
| `paper/sample_data/input/` | one small real-data slice per supported format (bed/bam/vcf/gff3/wig/bigWig/maf/sam) |
| `paper/sample_data/output/` | the FastCrossMap output committed for every `input/` file -- the reference the quick test diffs against |
| `paper/sample_data/hg19ToHg38.over.chain.gz` | UCSC chain file used by every conversion above |

Every `input/` file holds 100 records. See `paper/sample_data/README.md` for the
per-format input/output pairing.

---

## Script Inventory

Every `.py` / `.sh` in this directory, what it reads, and what it writes.
All paths are relative to the repository root.

### Data download

| Script | Reads | Writes |
|--------|-------|--------|
| `01_download_data.sh` | public URLs (UCSC / ENCODE / GENCODE / 1000G) | `paper/data/*` |

### Benchmarks - write `paper/results/*.json`

These invoke FastCrossMap plus the comparison tools (`CrossMap`, `liftOver`,
`FastRemap`, `samtools`, and for re-alignment `bwa`), so they need the tool
environment described under *Full Benchmark* below.

| Script | Command form | Writes |
|--------|--------------|--------|
| `02_benchmark_bed.py` | `bed` | `benchmark_bed.json` |
| `02b_benchmark_bed_multithread.py` | `bed -t N` | `benchmark_bed_multithread.json` |
| `02c_benchmark_vcf.py` | `vcf` | `benchmark_vcf.json` |
| `02d_benchmark_gff.py` | `gff` | `benchmark_gff.json` |
| `02e_benchmark_wig.py` | `wig` | `benchmark_wig.json` |
| `02f_benchmark_bigwig.py` | `bigwig` | `benchmark_bigwig.json` |
| `02g_benchmark_maf.py` | `maf` | `benchmark_maf.json` |
| `02h_benchmark_sam.py` | `bam` (SAM input) | `benchmark_sam.json` |
| `03_benchmark_bam.py` | `bam -t N` | `benchmark_bam.json` |
| `03b_benchmark_bam_multithread.py` | `bam -t N` | `benchmark_bam_multithread.json` |
| `05_memory_profile.py` | `bam` (RSS sampling) | `memory_profile.json` |
| `05b_memory_scalability.py` | `bam` (RSS sampling) | `memory_scalability.json` |
| `07_accuracy_analysis.py` | `bed` vs liftOver | `accuracy.json` |
| `09_benchmark_realignment.py` | `bam`, `bwa` | `benchmark_realignment.json` |
| `10_accuracy_realignment.py --scale 1M\|5M\|10M` | `bam`, `bwa` | `accuracy_realignment_<scale>.json` |

### Figures and tables - read the JSON above, write `paper/figures/*` and `paper/results/*`

These need only `matplotlib` + `numpy`; no comparison tools.

| Script | Reads | Writes |
|--------|-------|--------|
| `04_plot_performance.py` | `benchmark_bed.json`, `benchmark_bam.json`, `benchmark_{bed,bam}_multithread.json` | `fig1_performance.*`, `fig1a`-`fig1d` |
| `06_plot_memory.py` | `memory_profile.json` | `fig2_memory.*`, `fig2a`, `fig2b` |
| `06b_plot_memory_scalability.py` | `memory_scalability.json` | `figS2_memory_scalability.*`, `figS2a`-`figS2c` |
| `08_plot_accuracy.py` | `accuracy.json` | `fig3_accuracy.*`, `fig3a`, `fig3b` |
| `09_feature_audit.py` | (feature matrix is entered here) | `features.json` |
| `10_plot_features.py` | `features.json` | `fig4_features.*`, `fig4_features_heatmap` |
| `11_plot_radar.py` | `benchmark_bed.json`, `benchmark_bam.json`, `memory_profile.json`, `features.json`, `accuracy.json` | `fig5_radar.*`, `fig5_radar_chart`, `radar_scores.json` |
| `13_large_scale_estimation.py` | `benchmark_bed.json`, `benchmark_bam.json`, `benchmark_vcf.json` | `large_scale_estimation.json` |
| `18_plot_supplementary.py` | all eight `benchmark_*.json` | `fig_s1_all_formats.*`, `fig_s2_speedup_consistency.*` |
| `19_generate_supp_table.py` | all eight `benchmark_*.json` | `supplementary_tables.json` |

### Threading and determinism tools - need the FastCrossMap binary

| Script | Purpose | Writes |
|--------|---------|--------|
| `20_benchmark_threads.sh` | `-t N` scaling per format | `paper/results/threads-<timestamp>/{timings,summary}.tsv` |
| `21_check_determinism.sh` | output must be byte-identical across `-t` | `paper/results/determinism-<timestamp>/` |
| `22_run_benchmarks_pinned.sh` | run `02*`/`03` with one core pinned each | per-script logs |

### Drivers and environment

| Script | Purpose | Notes |
|--------|---------|-------|
| `12_run_all.sh` | convenience driver | Runs only `01`, `02`, `02b`, `03`, `03b`, `04`, `05`, `06`, `07`, `08`, `09_feature_audit`, `10_plot_features`, `11`. It does **not** cover `02c`-`02h`, `05b`, `06b`, `09_benchmark_realignment`, `10_accuracy_realignment`, `13`, `18`, `19` |
| `run_all_benchmarks.sh` | legacy 8-format driver | Superseded by `22_run_benchmarks_pinned.sh` (which pins CPUs) |
| `setup_benchmark_env.sh` | fresh-server bootstrap (conda env, tools, optional data) | - |

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
# Convenience driver (see the caveat on coverage in Script Inventory above)
bash paper/12_run_all.sh

# Or the per-format scripts directly, one per format:
python paper/02_benchmark_bed.py        # -> benchmark_bed.json
python paper/02c_benchmark_vcf.py       # -> benchmark_vcf.json
python paper/02d_benchmark_gff.py       # -> benchmark_gff.json
python paper/02e_benchmark_wig.py       # -> benchmark_wig.json
python paper/02f_benchmark_bigwig.py    # -> benchmark_bigwig.json
python paper/02g_benchmark_maf.py       # -> benchmark_maf.json
python paper/02h_benchmark_sam.py       # -> benchmark_sam.json
python paper/03_benchmark_bam.py        # -> benchmark_bam.json

# Memory and accuracy
python paper/05_memory_profile.py       # -> memory_profile.json
python paper/05b_memory_scalability.py  # -> memory_scalability.json
python paper/07_accuracy_analysis.py    # -> accuracy.json

# Re-alignment accuracy (needs bwa + the hg19/hg38 reference FASTA)
python paper/09_benchmark_realignment.py
python paper/10_accuracy_realignment.py --scale 10M --threads 4
```

Single-thread vs multi-thread curves for BED/BAM, used for Figure 1(b)/(d):
```bash
python paper/02b_benchmark_bed_multithread.py   # -> paper/results/benchmark_bed_multithread.json
python paper/03b_benchmark_bam_multithread.py   # -> paper/results/benchmark_bam_multithread.json
```

**Step 4: Generate figures and tables** (local, JSON only)
```bash
python paper/04_plot_performance.py
python paper/06_plot_memory.py
python paper/06b_plot_memory_scalability.py
python paper/08_plot_accuracy.py
python paper/10_plot_features.py
python paper/11_plot_radar.py
python paper/13_large_scale_estimation.py
python paper/18_plot_supplementary.py
python paper/19_generate_supp_table.py
```

All benchmark scripts expect the binary at `./target/release/fast-crossmap`
(build it with `cargo build --release`, or drop a pre-built binary there).
`02b` / `03b` / `05b` / `07` also honour `FCM_BIN=/path/to/fast-crossmap`, and the
shell scripts `20` / `21` honour `FCM=`.

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

## Expected Results

### Performance (BED, 296,898 records - K562 DNase-seq)

| Tool | Threads | Mean time (s) | Peak memory (MB) |
|------|---------|---------------|------------------|
| FastCrossMap | 1 | 0.27 | 11.7 |
| liftOver | 1 | 1.93 | 12.1 |
| FastRemap | multi (fixed) | 1.02 | 20.9 |
| CrossMap | 1 | 3.05 | 48.1 |

### Memory (BAM, 1.30 GB - encode_chipseq.bam)

| Tool | Threads | Mean time (s) | Peak memory (MB) |
|------|---------|---------------|------------------|
| FastCrossMap | 1 | 76.7 | 14.5 |
| FastRemap | multi (fixed) | 52.9 | 41.1 |
| CrossMap | 1 | 818.3 | 1022.0 |

### Accuracy

FastCrossMap produces bit-exact identical output to CrossMap (99.8% of records
identical to liftOver, the remaining 0.2% being partial/edge mappings handled
identically by both tools).

> These numbers are illustrative and were produced by the scripts above on the
> current benchmark data; re-running on another machine will differ in the
> absolute timings while preserving the order.

---

## Data Sources

| Data | Source | Accession |
|------|--------|-----------|
| BED | ENCODE | [ENCFF001WBV](https://www.encodeproject.org/files/ENCFF001WBV/) |
| BAM | ENCODE | [ENCFF000PED](https://www.encodeproject.org/files/ENCFF000PED/) |
| Chain | UCSC | [hg19ToHg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/) |
