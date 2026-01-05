#!/usr/bin/env python3
"""
Generate test data for FastCrossMap benchmarks.

Usage:
    python scripts/generate_test_data.py --format bed --count 100000 --output test_100k.bed
    python scripts/generate_test_data.py --format vcf --count 50000 --output test_50k.vcf
"""

import argparse
import random
from pathlib import Path

# Human chromosome sizes (GRCh37)
CHROM_SIZES = {
    "chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276,
    "chr5": 180915260, "chr6": 171115067, "chr7": 159138663, "chr8": 146364022,
    "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753,
    "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
    "chr21": 48129895, "chr22": 51304566, "chrX": 155270560, "chrY": 59373566,
}

CHROMS = list(CHROM_SIZES.keys())
BASES = "ACGT"


def generate_bed(count: int, output: Path):
    """Generate random BED records."""
    with open(output, 'w') as f:
        for i in range(count):
            chrom = random.choice(CHROMS[:22])  # Autosomes only
            max_pos = CHROM_SIZES[chrom] - 10000
            start = random.randint(10000, max_pos)
            size = random.randint(50, 5000)
            end = start + size
            strand = random.choice(['+', '-'])
            name = f"region_{i}"
            score = random.randint(0, 1000)
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
    print(f"Generated {count} BED records to {output}")


def generate_vcf(count: int, output: Path):
    """Generate random VCF records."""
    with open(output, 'w') as f:
        # Header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        for i in range(count):
            chrom = random.choice(CHROMS[:22])
            max_pos = CHROM_SIZES[chrom] - 1000
            pos = random.randint(10000, max_pos)
            ref = random.choice(BASES)
            alt = random.choice([b for b in BASES if b != ref])
            qual = random.randint(20, 100)
            dp = random.randint(10, 200)
            f.write(f"{chrom}\t{pos}\tvar_{i}\t{ref}\t{alt}\t{qual}\tPASS\tDP={dp}\n")
    print(f"Generated {count} VCF records to {output}")


def generate_gff(count: int, output: Path):
    """Generate random GFF records."""
    features = ["gene", "exon", "CDS", "mRNA", "transcript"]
    sources = ["ensembl", "refseq", "havana"]
    
    with open(output, 'w') as f:
        f.write("##gff-version 3\n")
        
        for i in range(count):
            chrom = random.choice(CHROMS[:22])
            max_pos = CHROM_SIZES[chrom] - 10000
            start = random.randint(10000, max_pos)
            size = random.randint(100, 10000)
            end = start + size
            strand = random.choice(['+', '-'])
            feature = random.choice(features)
            source = random.choice(sources)
            score = "."
            phase = random.choice([".", "0", "1", "2"])
            attrs = f"ID=feature_{i};Name=test_{i}"
            f.write(f"{chrom}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attrs}\n")
    print(f"Generated {count} GFF records to {output}")


def generate_region(count: int, output: Path):
    """Generate random large region records for region mapping."""
    with open(output, 'w') as f:
        for i in range(count):
            chrom = random.choice(CHROMS[:22])
            max_pos = CHROM_SIZES[chrom] - 2000000
            start = random.randint(10000, max(10001, max_pos))
            size = random.randint(100000, 1000000)  # Large regions
            end = min(start + size, CHROM_SIZES[chrom])
            strand = random.choice(['+', '-'])
            name = f"region_{i}"
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
    print(f"Generated {count} region records to {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--format', '-f', required=True, 
                        choices=['bed', 'vcf', 'gff', 'region'],
                        help='Output format')
    parser.add_argument('--count', '-c', type=int, default=10000,
                        help='Number of records to generate')
    parser.add_argument('--output', '-o', type=Path, required=True,
                        help='Output file path')
    parser.add_argument('--seed', '-s', type=int, default=42,
                        help='Random seed for reproducibility')
    
    args = parser.parse_args()
    random.seed(args.seed)
    
    generators = {
        'bed': generate_bed,
        'vcf': generate_vcf,
        'gff': generate_gff,
        'region': generate_region,
    }
    
    generators[args.format](args.count, args.output)


if __name__ == '__main__':
    main()
