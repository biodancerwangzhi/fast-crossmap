#!/usr/bin/env python3
"""
19_generate_supp_table.py - Generate Supplementary Tables from benchmark results

Reads all benchmark JSON files and generates:
  - Table S1: Complete benchmark results (all tools x formats x datasets),
              single-threaded only, each row carrying its dataset's provenance
              (source, file, genome build)

This was two tables. Table S2 was one row per dataset; everything it held that
Table S1 did not is now three columns on Table S1 -- Source, File and Genome
Build -- so a dataset's provenance sits on the row that reports its timings
instead of in a second table the reader has to join by name. The rest of Table
S2 (size, record count) was already in Table S1.

The measured read statistics from paper/14_bam_coverage.py are NOT carried by
either table any more. That script remains the record of those measurements
(results/bam_coverage.json), but nothing here reads it.

BigWig carries no record count in its own benchmark JSON, so its datapoint count
is derived from the Wiggle benchmark's row for the same dataset (the .wig files
are the bigWigToWig conversion of these .bigWig files); a dataset with no wig
counterpart is a hard error. Table S1's Throughput column stays a byte rate for
BigWig, the one rate that benchmark actually measured.

Usage: python paper/19_generate_supp_table.py
Output: paper/results/supplementary_tables.json
        paper/results/table_s1.tsv

benchmark_wig.json is required: without it the BigWig rows have no record count
and the script stops rather than write blank cells.
"""

import json
from datetime import datetime
from pathlib import Path

RESULTS_DIR = Path("paper/results")

# stable column order for the TSV: previously the header came from the first
# row's keys, so a column present only on later rows (e.g. the BAM read counts)
# vanished from the output entirely.
S1_COLUMNS = [
    "Format", "Dataset", "Tool", "Source", "File", "Genome Build",
    "Input Size (MB)", "Records", "Record Unit",
    "Time Mean (s)", "Time SD (s)", "Throughput", "Throughput Unit",
    "Peak Memory (MB)", "Num Runs", "Speedup vs CrossMap",
]

# Numeric columns rendered to a fixed two decimals in the released TSV, so the
# supplementary table and the manuscript text quote numbers to the same
# precision (the text asserts e.g. 0.27 s where the old TSV printed 0.274).
# Kept as floats in the JSON dump, which is a machine-readable artifact.
TWO_DP_COLUMNS = {
    "Input Size (MB)", "Time Mean (s)", "Time SD (s)",
    "Throughput", "Peak Memory (MB)",
}

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
    "BigWig": (None, None),          # see load_bigwig_datapoints() below
    "MAF": ("input_mutations", "mutations"),
    "SAM": ("total_reads", "reads"),
}

THROUGHPUT_KEY_MAP = {
    "BED": "throughput_rec_per_sec",
    "BAM": "throughput_reads_per_sec",
    "VCF": "throughput_variants_per_sec",
    "GFF": "throughput_features_per_sec",
    "Wiggle": "throughput_datapoints_per_sec",
    "BigWig": "throughput_mb_per_sec",
    "MAF": "throughput_mutations_per_sec",
    "SAM": "throughput_reads_per_sec",
}

# The column header is units, not a restatement of the key name. Deriving it
# from the key (strip "throughput_", replace "_" with "/") produced "rec/per/sec"
# and "mb/per/sec" -- a slash for every underscore, including the one before
# "per" -- so it is spelled out per format instead.
#
# Seven of the eight are the format's own record rate, which is what the
# Records column counts. BAM and SAM count reads too, so they report reads/s
# rather than the byte rate their benchmark happens to also record (both
# throughput_mb_per_sec and throughput_reads_per_sec are in the JSON): a byte
# rate in that column would be the one row out of eight not measuring a record
# rate, and would restate Input Size / Time.
#
# BigWig is the exception and stays a byte rate. Its Records column is now
# populated (see load_bigwig_datapoints), but the bigWig benchmark timed only
# bytes per second -- it never computed a per-datapoint rate, and deriving one
# would be Records / Time Mean, the same kind of pure unit conversion that got
# the GB/min and Mreads/min columns dropped. So BigWig reports the one rate
# that was actually measured.
THROUGHPUT_UNIT_MAP = {
    "BED": "records/s",
    "BAM": "reads/s",
    "VCF": "variants/s",
    "GFF": "features/s",
    "Wiggle": "datapoints/s",
    "BigWig": "MB/s",
    "MAF": "mutations/s",
    "SAM": "reads/s",
}

BIGWIG_FORMAT = "BigWig"
BIGWIG_RECORD_UNIT = "datapoints"

# The bigWig benchmark records no record count at all -- only bytes and seconds
# -- so BigWig is the one format whose Records cell the JSON cannot fill. The
# five .bigWig inputs are exactly the sources the five .wig inputs were
# converted from (bigWigToWig), and the Wiggle benchmark counts a datapoint for
# every data line, i.e. every interval, so the wig row for the same dataset
# carries the number. The two agree by construction; pyBigWig interval counts
# confirm it for two of the five (GM12878 H3K27ac 45,812,448 and H3K4me3
# 83,889,296, both exact).
#
# The dataset name is the join key and is spelled identically in both files. A
# dataset with no wig counterpart is a hard error rather than a blank cell: a
# blank is the exact defect this fills, so it must not be able to return.


def load_bigwig_datapoints() -> dict[str, int]:
    """{bigWig dataset name: datapoint count}, from the Wiggle benchmark."""
    data = load_benchmark(RESULTS_DIR / "benchmark_wig.json")
    if not data:
        raise SystemExit(
            "benchmark_wig.json is missing, and it is where the BigWig "
            "datapoint counts come from (the bigWig benchmark does not record "
            "one). Refusing to write a table with an empty BigWig Records "
            "column -- re-run the Wiggle benchmark first."
        )
    out: dict[str, int] = {}
    for r in data.get("results", []):
        name, n = r.get("dataset_name"), r.get("input_datapoints")
        if r.get("success") and name and n:
            out.setdefault(name, n)
    return out


def bigwig_record_count(dataset: str, lookup: dict[str, int]):
    """
    The datapoint count for a BigWig dataset, or a hard error naming it.

    Raises rather than returning a blank: a BigWig row whose Records cell is
    empty is precisely the problem this lookup exists to fix, so failing loudly
    is the only outcome that cannot silently reintroduce it.
    """
    n = lookup.get(dataset)
    if n is None:
        raise SystemExit(
            f"no Wiggle row for the BigWig dataset {dataset!r}, so its datapoint "
            "count cannot be derived; the two benchmarks' dataset names must "
            f"match. Wiggle datasets available: {sorted(lookup)}"
        )
    return n

# Minor 4 asked for the read count alongside the file size, and for BAM/SAM the
# read count is now readable two ways off the table: Records / Record Unit
# ("reads") and the reads/s throughput. An earlier draft also carried GB/min and
# Mreads/min columns; both were dropped because neither says anything the row
# does not already say -- Mreads/min is reads/s x 60/1e6, and GB/min is
# Input Size / Time Mean x 60/1024. They were pure unit conversions of columns
# sitting on the same row.
#
# The measured read statistics (mapped reads, read length, covered span, depth)
# that paper/14_bam_coverage.py produced are deliberately not among the columns:
# see that script's docstring for the numbers, which are released as
# results/bam_coverage.json but are not carried by this table.

# every benchmark records the size of the *archived* input -- the .gz for the
# text formats (02d_benchmark_gff.py measures the .gz before decompressing),
# the file itself for BigWig/MAF/BAM/SAM. The BED benchmark predates that field
# and stores nothing, so its 20 rows would print 0 MB; fall back to measuring
# the input file named in the results.
MIB = 1024 * 1024


def resolve_input_size_mb(recorded, file_hint: str) -> float:
    """
    The benchmark's recorded input size, or the on-disk size of the input when
    the format's benchmark never stored one (BED). Returns 0.0 only when the
    file cannot be found either, so a missing number is visible rather than
    silently plausible.
    """
    if recorded:
        return recorded
    if file_hint:
        path = Path(file_hint)
        if path.exists():
            return round(path.stat().st_size / MIB, 2)
    return 0.0


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


def dataset_metadata(data: dict) -> dict[str, dict]:
    """
    {dataset name: {column name: value}} for the table's three provenance
    columns, taken from a benchmark's own `datasets` list.

    These used to be a second table keyed by the same dataset names; joining
    two released TSVs by hand is exactly the step a reader should not have to
    do, so they travel on the row instead. The keys are the S1 column names so
    the result can be merged straight into the row. Every benchmark measures
    against hg19/GRCh37 as the source assembly.
    """
    out = {}
    for ds in data.get("datasets", []):
        name = ds.get("name", "")
        if not name:
            continue
        out[name] = {
            "Source": ds.get("source", ""),
            "File": Path(ds.get("file", "")).name,
            "Genome Build": "hg19/GRCh37",
        }
    return out


def generate_table_s1(bigwig_dp: dict[str, int]) -> list[dict]:
    rows = []
    for fmt, json_file in BENCHMARK_FILES.items():
        data = load_benchmark(json_file)
        if not data:
            continue
        record_key, record_unit = RECORD_KEY_MAP.get(fmt, (None, None))
        throughput_key = THROUGHPUT_KEY_MAP.get(fmt, "")
        meta = dataset_metadata(data)

        for r in data.get("results", []):
            if not r.get("success"):
                continue
            # the threaded FastCrossMap rows (BAM/SAM only) are a threading
            # result, not a single-thread comparison: their "speedup vs
            # CrossMap" would divide -t 4/-t 8 wall time by CrossMap's -t 1,
            # which is not the baseline used anywhere else in this work.
            # Figure S4 drops them too, for the same reason.
            if r.get("tool", "").startswith("FastCrossMap ("):
                continue
            times = r.get("all_times", [])
            std = compute_std(times)
            dataset = r.get("dataset_name", "")
            # a result whose dataset is not in the file's own `datasets` list
            # would print blank provenance: exactly the defect the join is
            # meant to prevent, so it stops rather than blanks the row
            if dataset not in meta:
                raise SystemExit(
                    f"{json_file.name}: result names dataset {dataset!r}, which "
                    f"is not in that file's `datasets` list "
                    f"({sorted(meta)}). The provenance columns cannot be "
                    "filled, and a blank Source/File cell is worse than no row."
                )
            if fmt == BIGWIG_FORMAT:
                # the bigWig benchmark has no record count; the wig row for the
                # same dataset supplies it (load_bigwig_datapoints)
                record_count = bigwig_record_count(dataset, bigwig_dp)
                record_unit = BIGWIG_RECORD_UNIT
            else:
                record_count = r.get(record_key, 0) if record_key else None
            throughput = r.get(throughput_key, 0)

            row = {
                "Format": fmt,
                "Dataset": dataset,
                "Tool": r.get("tool", ""),
                "Input Size (MB)": resolve_input_size_mb(
                    r.get("input_size_mb"), r.get("input_file", "")
                ),
                "Records": record_count,
                "Record Unit": record_unit,
                "Time Mean (s)": r.get("execution_time_sec", 0),
                "Time SD (s)": round(std, 3),
                "Throughput": round(throughput, 1),
                "Throughput Unit": THROUGHPUT_UNIT_MAP.get(fmt, ""),
                "Peak Memory (MB)": r.get("peak_memory_mb", 0),
                "Num Runs": len(times),
            }
            row.update(meta[dataset])
            rows.append(row)

    # Compute speedup relative to CrossMap for each format+dataset
    for row in rows:
        if row["Tool"] == "CrossMap":
            row["Speedup vs CrossMap"] = "1.00x (ref)"
        else:
            cm = next((r for r in rows if r["Format"] == row["Format"]
                       and r["Dataset"] == row["Dataset"]
                       and r["Tool"] == "CrossMap"), None)
            if cm and cm["Time Mean (s)"] > 0 and row["Time Mean (s)"] > 0:
                speedup = cm["Time Mean (s)"] / row["Time Mean (s)"]
                row["Speedup vs CrossMap"] = f"{speedup:.2f}x"
            else:
                row["Speedup vs CrossMap"] = "N/A"

    # the provenance columns are merged in with dict.update, which appends them,
    # so the in-memory order would differ from S1_COLUMNS (the TSV is written in
    # S1_COLUMNS regardless, but the JSON dump follows insertion order and would
    # list them last). Reorder so both released files read the same way.
    return [{k: r[k] for k in S1_COLUMNS} for r in rows]


def write_tsv(rows: list[dict], output_file: Path, columns: list[str] | None = None):
    if not rows:
        return
    if columns is None:
        columns = list(rows[0].keys())
    # keep the declared order first, then append anything unexpected
    headers = columns + [k for k in rows[0] if k not in columns]
    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for row in rows:
            vals = []
            for h in headers:
                v = row.get(h, "")
                if h in TWO_DP_COLUMNS and isinstance(v, (int, float)):
                    vals.append(f"{v:.2f}")
                else:
                    vals.append(str(v) if v is not None else "")
            f.write('\t'.join(vals) + '\n')


def main():
    print("=" * 60)
    print("Generating Supplementary Tables")
    print("=" * 60)

    # resolved up front so a missing wig counterpart fails before any output is
    # written, rather than half way through the first table
    bigwig_dp = load_bigwig_datapoints()

    # Table S1
    print("\nTable S1: Complete benchmark results")
    s1_rows = generate_table_s1(bigwig_dp)
    if s1_rows:
        write_tsv(s1_rows, RESULTS_DIR / "table_s1.tsv", S1_COLUMNS)
        print(f"  {len(s1_rows)} rows -> paper/results/table_s1.tsv")

        formats_found = sorted(set(r["Format"] for r in s1_rows))
        tools_found = sorted(set(r["Tool"] for r in s1_rows))
        print(f"  Formats: {', '.join(formats_found)}")
        print(f"  Tools: {', '.join(tools_found)}")
        with_records = [r for r in s1_rows if r.get("Records")]
        print(f"  {len(with_records)}/{len(s1_rows)} rows carry a record count")
        with_prov = [r for r in s1_rows if r.get("Source") and r.get("File")]
        print(f"  {len(with_prov)}/{len(s1_rows)} rows carry dataset provenance "
              "(source + file); Table S2 has been folded into these columns")
    else:
        print("  No benchmark results found. Run benchmarks first.")

    # Combined JSON
    output_json = RESULTS_DIR / "supplementary_tables.json"
    with open(output_json, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "table_s1": s1_rows,
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
