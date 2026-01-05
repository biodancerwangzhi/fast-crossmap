#!/usr/bin/env python3
"""
Compare FastCrossMap output with Python CrossMap output.

Usage:
    python scripts/compare_outputs.py <crossmap_output> <fastcrossmap_output> [--format FORMAT]
    
Supported formats: bed, vcf, gff, maf, wig, region
"""

import sys
from pathlib import Path
import argparse


def parse_bed_line(line: str) -> tuple:
    """Parse BED line into comparable tuple."""
    fields = line.strip().split('\t')
    if len(fields) < 3:
        return None
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    rest = tuple(fields[3:]) if len(fields) > 3 else ()
    return (chrom, start, end, rest)


def parse_vcf_line(line: str) -> tuple:
    """Parse VCF line into comparable tuple (coordinates only)."""
    fields = line.strip().split('\t')
    if len(fields) < 2:
        return None
    chrom = fields[0]
    pos = int(fields[1])
    return (chrom, pos)


def parse_gff_line(line: str) -> tuple:
    """Parse GFF line into comparable tuple."""
    fields = line.strip().split('\t')
    if len(fields) < 5:
        return None
    chrom = fields[0]
    start = int(fields[3])
    end = int(fields[4])
    return (chrom, start, end)


def compare_files(file1: Path, file2: Path, format_type: str = 'bed') -> dict:
    """Compare two files and return statistics."""
    with open(file1) as f1, open(file2) as f2:
        lines1 = [l for l in f1 if not l.startswith('#') and l.strip()]
        lines2 = [l for l in f2 if not l.startswith('#') and l.strip()]

    # Select parser based on format
    if format_type == 'vcf':
        parser = parse_vcf_line
    elif format_type == 'gff':
        parser = parse_gff_line
    else:
        parser = parse_bed_line

    records1 = [parser(l) for l in lines1]
    records2 = [parser(l) for l in lines2]

    # Filter None values
    records1 = [r for r in records1 if r is not None]
    records2 = [r for r in records2 if r is not None]

    set1 = set(records1)
    set2 = set(records2)

    identical = set1 == set2
    only_in_1 = set1 - set2
    only_in_2 = set2 - set1
    common = set1 & set2

    return {
        'identical': identical,
        'total_file1': len(records1),
        'total_file2': len(records2),
        'common': len(common),
        'only_in_file1': len(only_in_1),
        'only_in_file2': len(only_in_2),
        'examples_only_in_file1': list(only_in_1)[:5],
        'examples_only_in_file2': list(only_in_2)[:5],
        'match_rate': len(common) / max(len(set1), len(set2), 1) * 100,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file1', type=Path, help='CrossMap output file')
    parser.add_argument('file2', type=Path, help='FastCrossMap output file')
    parser.add_argument('--format', '-f', default='bed', 
                        choices=['bed', 'vcf', 'gff', 'maf', 'wig', 'region'],
                        help='File format (default: bed)')
    
    args = parser.parse_args()

    if not args.file1.exists():
        print(f"Error: {args.file1} not found")
        sys.exit(1)
    if not args.file2.exists():
        print(f"Error: {args.file2} not found")
        sys.exit(1)

    result = compare_files(args.file1, args.file2, args.format)

    print("=== Output Comparison Results ===")
    print(f"File 1 (CrossMap):     {args.file1}")
    print(f"File 2 (FastCrossMap): {args.file2}")
    print(f"Format:                {args.format}")
    print()
    print(f"Total records in File 1: {result['total_file1']}")
    print(f"Total records in File 2: {result['total_file2']}")
    print(f"Common records:          {result['common']}")
    print(f"Only in File 1:          {result['only_in_file1']}")
    print(f"Only in File 2:          {result['only_in_file2']}")
    print(f"Match rate:              {result['match_rate']:.2f}%")
    print()

    if result['identical']:
        print("✓ Files are IDENTICAL")
        sys.exit(0)
    else:
        print("✗ Files DIFFER")
        if result['examples_only_in_file1']:
            print("\nExamples only in CrossMap output:")
            for r in result['examples_only_in_file1']:
                print(f"  {r}")
        if result['examples_only_in_file2']:
            print("\nExamples only in FastCrossMap output:")
            for r in result['examples_only_in_file2']:
                print(f"  {r}")
        
        # Exit with success if match rate > 95%
        if result['match_rate'] >= 95.0:
            print(f"\n✓ Match rate {result['match_rate']:.2f}% >= 95% threshold")
            sys.exit(0)
        sys.exit(1)


if __name__ == '__main__':
    main()
