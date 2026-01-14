#!/usr/bin/env python3
"""
generate_sample_data.py - Generate small sample datasets for quick benchmark

This script creates small test files that can be included in the GitHub repository
for reviewers to quickly verify FastCrossMap functionality.

Output directory: paper/sample_data/
"""

import os
import random
import subprocess
import urllib.request
from pathlib import Path

# Configuration
SAMPLE_DIR = Path("paper/sample_data")
SAMPLE_DIR.mkdir(parents=True, exist_ok=True)

# hg19 chromosome sizes (for generating valid coordinates)
HG19_CHROMS = {
    "chr1": 249250621, "chr2": 243199373, "chr3": 198022430,
    "chr4": 191154276, "chr5": 180915260, "chr6": 171115067,
    "chr7": 159138663, "chr8": 146364022, "chr9": 141213431,
    "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
    "chr16": 90354753, "chr17": 81195210, "chr18": 78077248,
    "chr19": 59128983, "chr20": 63025520, "chr21": 48129895,
    "chr22": 51304566, "chrX": 155270560, "chrY": 59373566,
}


def download_chain_file():
    """Download hg19ToHg38 chain file from UCSC"""
    chain_file = SAMPLE_DIR / "hg19ToHg38.over.chain.gz"
    
    if chain_file.exists():
        print(f"  ✓ Chain file already exists: {chain_file}")
        return chain_file
    
    print("  Downloading chain file from UCSC...")
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
    
    try:
        urllib.request.urlretrieve(url, chain_file)
        print(f"  ✓ Downloaded: {chain_file}")
    except Exception as e:
        print(f"  ✗ Failed to download chain file: {e}")
        return None
    
    return chain_file


def generate_sample_bed(num_records=10000):
    """Generate synthetic BED file with random genomic coordinates"""
    bed_file = SAMPLE_DIR / "sample.bed"
    
    if bed_file.exists():
        print(f"  ✓ BED file already exists: {bed_file}")
        return bed_file
    
    print(f"  Generating {num_records} BED records...")
    random.seed(42)  # Reproducibility
    
    chroms = list(HG19_CHROMS.keys())
    
    with open(bed_file, 'w') as f:
        for i in range(num_records):
            chrom = random.choice(chroms)
            chrom_size = HG19_CHROMS[chrom]
            
            # Random start position (leave room for region)
            start = random.randint(0, chrom_size - 10000)
            # Random region length (100-5000 bp)
            length = random.randint(100, 5000)
            end = start + length
            
            # BED format: chrom, start, end, name, score, strand
            name = f"region_{i+1:06d}"
            score = random.randint(0, 1000)
            strand = random.choice(['+', '-'])
            
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
    
    print(f"  ✓ Generated: {bed_file}")
    return bed_file


def generate_sample_sam():
    """Generate a minimal SAM file for testing (no external dependencies)"""
    sam_file = SAMPLE_DIR / "sample.sam"
    
    if sam_file.exists():
        print(f"  ✓ SAM file already exists: {sam_file}")
        return sam_file
    
    print("  Generating sample SAM file...")
    random.seed(42)
    
    # SAM header
    header = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:249250621
@SQ\tSN:chr2\tLN:243199373
@SQ\tSN:chr3\tLN:198022430
@SQ\tSN:chr10\tLN:135534747
@SQ\tSN:chr22\tLN:51304566
@PG\tID:sample_generator\tPN:generate_sample_data.py\tVN:1.0
"""
    
    chroms = ['chr1', 'chr2', 'chr3', 'chr10', 'chr22']
    chrom_sizes = {c: HG19_CHROMS[c] for c in chroms}
    
    # Generate random sequence
    def random_seq(length):
        return ''.join(random.choices('ACGT', k=length))
    
    with open(sam_file, 'w') as f:
        f.write(header)
        
        # Generate 1000 reads
        for i in range(1000):
            chrom = random.choice(chroms)
            pos = random.randint(1, chrom_sizes[chrom] - 200)
            read_len = 100
            
            qname = f"read_{i+1:06d}"
            flag = random.choice([0, 16])  # Forward or reverse
            mapq = random.randint(20, 60)
            cigar = f"{read_len}M"
            seq = random_seq(read_len)
            qual = 'I' * read_len  # Phred 40
            
            f.write(f"{qname}\t{flag}\t{chrom}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\n")
    
    print(f"  ✓ Generated: {sam_file}")
    return sam_file


def main():
    """Main function to generate all sample data"""
    print("=" * 60)
    print("FastCrossMap Sample Data Generator")
    print("=" * 60)
    print()
    
    print("[1/3] Downloading chain file...")
    chain_file = download_chain_file()
    
    print()
    print("[2/3] Generating sample BED file...")
    bed_file = generate_sample_bed(num_records=10000)
    
    print()
    print("[3/3] Generating sample SAM file...")
    sam_file = generate_sample_sam()
    
    print()
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    
    # List generated files
    print("\nGenerated files:")
    for f in SAMPLE_DIR.iterdir():
        if f.is_file():
            size = f.stat().st_size
            if size > 1024 * 1024:
                size_str = f"{size / 1024 / 1024:.1f} MB"
            elif size > 1024:
                size_str = f"{size / 1024:.1f} KB"
            else:
                size_str = f"{size} B"
            print(f"  {f.name}: {size_str}")
    
    print()
    print("Next steps:")
    print("  1. Build FastCrossMap: cargo build --release")
    print("  2. Run quick benchmark: python paper/benchmark_quick.py")
    print()


if __name__ == "__main__":
    main()