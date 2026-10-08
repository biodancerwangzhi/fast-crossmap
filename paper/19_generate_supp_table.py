#!/usr/bin/env python3
"""
19_generate_supp_table.py - Generate Supplementary Tables from benchmark results

Reads all benchmark JSON files and generates:
  - Table S1: Complete benchmark results (all tools × formats × datasets)
  - Table S2: Dataset details (source, accession, size, record count)

Usage: python paper/19_generate_supp_table.py
Output: paper/results/supplementary_tables.json
        paper/results/table_s1.tsv
        paper/results/table_s2.tsv
"""

import json
from datetime import datetime
from pathlib import Path

RESULTS_DIR = Path("paper/results")

BENCHMARK_FILES = {
    "BED": RESULTS_DIR / "benchmark_bed.json",
    "BAM": RESULTS_DIR / "benchmark_bam.json",
    "VCF": RESULTS_DIR / "benchmark_vcf.json",
    "GFF": RESULTS_DIR / "benchmark_gff.json",
    "Wiggle": RESULTS_DIR / "benchmark_wig.json",
    "BigWig": RESULTS_DIR / "benchmark_bigwig.json",
    "MAF": RESULTS_DIR / "benchmark_maf.json",
    "SAM": RESULTS_DIR / "benchmark_sam.json",
}

RECORD_KEY_MAP = {
    "BED": ("input_records", "records"),
    "BAM": ("total_reads", "reads"),
    "VCF": ("input_variants", "variants"),
    "GFF": ("input_features", "features"),
    "Wiggle": ("input_datapoints", "datapoints"),
    "BigWig": (None, None),
    "MAF": ("input_mutations", "mutations"),
    "SAM": ("total_reads", "reads"),
}

THROUGHPUT_KEY_MAP = {
    "BED": "throughput_rec_per_sec",
    "BAM": "throughput_mb_per_sec",
    "VCF": "throughput_variants_per_sec",
    "GFF": "throughput_features_per_sec",
    "Wiggle": "throughput_datapoints_per_sec",
    "BigWig": "throughput_mb_per_sec",
    "MAF": "throughput_mutations_per_sec",
    "SAM": "throughput_mb_per_sec",
}


def load_benchmark(json_file: Path) -> dict | None:
    if not json_file.exists():
        return None
    with open(json_file) as f:
        return json.load(f)


def compute_std(times: list[float]) -> float:
    if len(times) < 2:
        return 0.0
    mean = sum(times) / len(times)
    variance = sum((t - mean) ** 2 for t in times) / (len(times) - 1)
    return variance ** 0.5


def generate_table_s1() -> list[dict]:
    rows = []
    for fmt, json_file in BENCHMARK_FILES.items():
        data = load_benchmark(json_file)
        if not data:
            continue
        record_key, record_unit = RECORD_KEY_MAP.get(fmt, (None, None))
        throughput_key = THROUGHPUT_KEY_MAP.get(fmt, "")

        for r in data.get("results", []):
            if not r.get("success"):
                continue
            times = r.get("all_times", [])
            std = compute_std(times)
            record_count = r.get(record_key, 0) if record_key else None
            throughput = r.get(throughput_key, 0)

            row = {
                "Format": fmt,
                "Dataset": r.get("dataset_name", ""),
                "Tool": r.get("tool", ""),
                "Input Size (MB)": r.get("input_size_mb", 0),
                "Records": record_count,
                "Record Unit": record_unit,
                "Time Mean (s)": r.get("execution_time_sec", 0),
                "Time SD (s)": round(std, 3),
                "Throughput": round(throughput, 1),
                "Throughput Unit": throughput_key.replace("throughput_", "").replace("_", "/") if throughput_key else "",
                "Peak Memory (MB)": r.get("peak_memory_mb", 0),
                "Num Runs": len(times),
            }
            rows.append(row)

    # Compute speedup relative to CrossMap for each format+dataset
    for row in rows:
        if row["Tool"] == "CrossMap":
            row["Speedup vs CrossMap"] = "1.0x (ref)"
        else:
            cm = next((r for r in rows if r["Format"] == row["Format"]
                       and r["Dataset"] == row["Dataset"]
                       and r["Tool"] == "CrossMap"), None)
            if cm and cm["Time Mean (s)"] > 0 and row["Time Mean (s)"] > 0:
                speedup = cm["Time Mean (s)"] / row["Time Mean (s)"]
                row["Speedup vs CrossMap"] = f"{speedup:.1f}x"
            else:
                row["Speedup vs CrossMap"] = "N/A"

    return rows


def generate_table_s2() -> list[dict]:
    rows = []
    for fmt, json_file in BENCHMARK_FILES.items():
        data = load_benchmark(json_file)
        if not data:
            continue
        record_key, record_unit = RECORD_KEY_MAP.get(fmt, (None, None))

        for ds in data.get("datasets", []):
            ds_name = ds.get("name", "")
            first_result = next(
                (r for r in data.get("results", []) if r.get("dataset_name") == ds_name),
                {}
            )
            record_count = first_result.get(record_key, 0) if record_key else None
            stats = ds.get("stats", {})

            row = {
                "Format": fmt,
                "Dataset": ds_name,
                "Source": ds.get("source", ""),
                "File": Path(ds.get("file", "")).name,
                "Size (MB)": first_result.get("input_size_mb", 0),
                "Records": record_count,
                "Record Unit": record_unit,
                "Genome Build": "hg19/GRCh37",
            }
            if stats:
                if "total_reads" in stats:
                    row["Total Reads"] = stats["total_reads"]
                if "mapped_reads" in stats:
                    row["Mapped Reads"] = stats["mapped_reads"]
            rows.append(row)

    return rows


def write_tsv(rows: list[dict], output_file: Path):
    if not rows:
        return
    headers = list(rows[0].keys())
    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for row in rows:
            vals = []
            for h in headers:
                v = row.get(h, "")
                vals.append(str(v) if v is not None else "")
            f.write('\t'.join(vals) + '\n')


def main():
    print("=" * 60)
    print("Generating Supplementary Tables")
    print("=" * 60)

    # Table S1
    print("\nTable S1: Complete benchmark results")
    s1_rows = generate_table_s1()
    if s1_rows:
        write_tsv(s1_rows, RESULTS_DIR / "table_s1.tsv")
        print(f"  {len(s1_rows)} rows -> paper/results/table_s1.tsv")

        formats_found = sorted(set(r["Format"] for r in s1_rows))
        tools_found = sorted(set(r["Tool"] for r in s1_rows))
        print(f"  Formats: {', '.join(formats_found)}")
        print(f"  Tools: {', '.join(tools_found)}")
    else:
        print("  No benchmark results found. Run benchmarks first.")

    # Table S2
    print("\nTable S2: Dataset details")
    s2_rows = generate_table_s2()
    if s2_rows:
        write_tsv(s2_rows, RESULTS_DIR / "table_s2.tsv")
        print(f"  {len(s2_rows)} rows -> paper/results/table_s2.tsv")
    else:
        print("  No dataset info found.")

    # Combined JSON
    output_json = RESULTS_DIR / "supplementary_tables.json"
    with open(output_json, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "table_s1": s1_rows,
            "table_s2": s2_rows,
        }, f, indent=2)
    print(f"\nAll tables saved to: {output_json}")

    # Print Table S1 summary
    if s1_rows:
        print(f"\n{'='*90}")
        print("Table S1 Preview (FastCrossMap speedup vs CrossMap)")
        print(f"{'='*90}")
        print(f"{'Format':<10} {'Dataset':<25} {'Tool':<15} {'Time(s)':<10} {'Speedup'}")
        print("-" * 70)
        for r in s1_rows:
            print(f"{r['Format']:<10} {r['Dataset']:<25} {r['Tool']:<15} "
                  f"{r['Time Mean (s)']:<10} {r['Speedup vs CrossMap']}")


if __name__ == "__main__":
    main()
