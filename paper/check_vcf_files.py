#!/usr/bin/env python3
"""Check if all files required by 02c_benchmark_vcf.py exist."""
from pathlib import Path

DATA_DIR = Path("paper/data")
CHAIN_FILE = DATA_DIR / "hg19ToHg38.over.chain.gz"
REF_GENOME = DATA_DIR / "hg38.fa"

VCF_FILES = [
    ("chain file", CHAIN_FILE),
    ("ref genome", REF_GENOME),
    ("ref genome index", DATA_DIR / "hg38.fa.fai"),
    ("1KG chr22", DATA_DIR / "vcf_1kg_chr22.vcf.gz"),
    ("1KG chr22 (legacy)", DATA_DIR / "1000g_chr22.vcf.gz"),
    ("1KG chr1", DATA_DIR / "vcf_1kg_chr1.vcf.gz"),
    ("1KG chr11", DATA_DIR / "vcf_1kg_chr11.vcf.gz"),
    ("1KG chr20", DATA_DIR / "vcf_1kg_chr20.vcf.gz"),
    ("1KG chr2", DATA_DIR / "vcf_1kg_chr2.vcf.gz"),
    ("FastCrossMap binary", Path("./target/release/fast-crossmap")),
]

print("=" * 60)
print("VCF Benchmark File Check")
print("=" * 60)

missing = []
for label, path in VCF_FILES:
    if path.exists():
        size = path.stat().st_size
        if size > 1024 * 1024:
            size_str = f"{size / (1024*1024):.1f} MB"
        elif size > 1024:
            size_str = f"{size / 1024:.1f} KB"
        else:
            size_str = f"{size} B"
        print(f"  OK   {label:25s}  {size_str:>10s}  {path}")
    else:
        print(f"  MISS {label:25s}  {'':>10s}  {path}")
        missing.append((label, path))

print()
if missing:
    print(f"MISSING {len(missing)} file(s):")
    for label, path in missing:
        print(f"  - {label}: {path}")
else:
    print("All files present.")

# Also check CrossMap availability
import shutil
cm = shutil.which("CrossMap")
print(f"\nCrossMap in PATH: {'OK (' + cm + ')' if cm else 'NOT FOUND'}")

import subprocess
try:
    r = subprocess.run(["conda", "run", "-n", "fcm", "CrossMap", "--version"],
                       capture_output=True, text=True, timeout=15)
    print(f"CrossMap in conda fcm: {'OK - ' + r.stdout.strip() if r.returncode == 0 else 'FAILED - ' + r.stderr.strip()[:100]}")
except Exception as e:
    print(f"CrossMap in conda fcm: ERROR - {e}")
