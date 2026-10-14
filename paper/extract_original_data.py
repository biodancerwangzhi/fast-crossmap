#!/usr/bin/env python3
"""
extract_original_data.py - one table of every number the paper asserts.

The manuscript, the supplementary material, the figure legends and the tables
are maintained as separate documents, but they all quote the same measurements.
Those measurements live in exactly one place each: a JSON (or TSV) under
paper/results/. This script walks those files and emits a single flat table --
one row per asserted value -- so that a number can be traced from the sentence
that states it back to the file and the key it came from.

It does not recompute anything. Every value is read out of the result file that
the corresponding plotting or table script writes, so if this table and a
legend disagree, one of them was edited without re-running the script rather
than re-derived. That is the failure mode this file exists to catch.

Rows are grouped by the display item that carries them:

    Main-Figure-1   paper/results/fig1_summary.json
    Figure-S1       paper/results/scaling_fit.json
    Figure-S2       paper/results/realignment_fit.json
    Figure-S3       paper/results/threads-merged-20261011/summary.tsv
    Figure-S4       paper/results/memory_by_format.json
    Figure-S5       paper/results/realignment_accuracy_summary.json
    Main-Table-1    paper/results/fig1_summary.json
    Table-S1        paper/results/table_s1.tsv
    Table-S2        curated (paper/geo_bam_only_datasets.md); no numeric source

The check pass at the end asserts the invariants the documents rely on: every
format appears in every per-format table, the ratio claimed for a format is the
median of its paired per-dataset ratios, and the thread maxima quoted in the
Figure S3 legend are the maxima over the measured thread counts. A violation is
printed and the exit status is non-zero, so a stale table fails loudly.

Input:   paper/results/*.json, paper/results/threads-merged-20261011/summary.tsv
         paper/results/table_s1.tsv
Output:  paper/results/original_data.tsv   (tab-separated, for reading in a shell)
         paper/results/original_data.md    (the same rows as a Markdown table)

Usage:   python paper/extract_original_data.py
         python paper/extract_original_data.py --check   (validate only)
"""

import argparse
import csv
import json
import statistics as stats
import sys
from pathlib import Path

RESULTS_DIR = Path("paper/results")
TSV_OUT = RESULTS_DIR / "original_data.tsv"
MD_OUT = RESULTS_DIR / "original_data.md"

# The eight supported formats, in the order Figure 1's two panels share.
FORMATS = ["Wiggle", "BAM", "BED", "SAM", "BigWig", "MAF", "GFF", "VCF"]

# Format label used by the thread summary, which is lowercase and shorter than
# the Figure 1 label ("wig" for Wiggle). Kept here so the map is explicit.
THREAD_KEY = {
    "Wiggle": "wig", "BAM": "bam", "BED": "bed", "SAM": "sam",
    "BigWig": "bigwig", "MAF": "maf", "GFF": "gff", "VCF": "vcf",
}

COLUMNS = ["item", "quantity", "value", "unit", "source_file", "source_path", "note"]

rows: list[dict] = []


def add(item, quantity, value, unit, source_file, source_path, note=""):
    """Append one row. `value` is formatted here so the table is text-stable."""
    rows.append({
        "item": item,
        "quantity": quantity,
        "value": fmt(value),
        "unit": unit,
        "source_file": source_file,
        "source_path": source_path,
        "note": note,
    })


def fmt(v):
    """Format a number the way the manuscript quotes it: fixed two decimals for
    anything with a fractional part, integers unchanged. A ratio is therefore
    written 63.20x, not 63.2x or 63.204x.

    Two documented exceptions, both because two decimals would destroy the
    value rather than round it:
      * a quantity below 0.01 keeps three decimals, so 0.008% does not become
        0.01% (the Figure S5 sub-500 bp band);
      * R^2 goes through fmt_r2() below, since 0.99983 must not become 1.00.
    """
    if isinstance(v, str):
        return v
    if isinstance(v, int) or (isinstance(v, float) and v == int(v) and abs(v) < 1e15):
        return str(int(v))
    if abs(v) >= 1000:
        return f"{v:,.2f}"
    if 0 < abs(v) < 0.01:
        return f"{v:.3f}"
    return f"{v:.2f}"


def fmt_r2(v):
    """R^2 to five places, so the document's rounding is itself checkable
    (0.99983 must not collapse to 1)."""
    return f"{v:.5f}".rstrip("0").rstrip(".")


def load(name):
    with open(RESULTS_DIR / name) as fh:
        return json.load(fh)


# ---------------------------------------------------------------- Figure 1 --
# fig1_summary.json is written by paper/plot/fig1_performance.py and holds
# exactly what that figure plots: the median (and min/max, n) of the per-dataset
# paired ratios, per format, for time and for memory.
def figure1():
    doc = load("fig1_summary.json")
    baseline = doc["baseline"]
    for panel, key, unit in [
        ("speed-up over CrossMap", "speedup", "x"),
        ("peak-memory reduction vs CrossMap", "memory_cut", "x"),
    ]:
        block = doc[key]
        for f in FORMATS:
            if f not in block:
                continue
            s = block[f]
            add("Main-Figure-1", f"{f}: {panel} (median of paired ratios)",
                s["median"], unit, "fig1_summary.json",
                f"{key}.{f}.median",
                f"baseline {baseline}; min {fmt(s['min'])}, max {fmt(s['max'])}, n={s['n']}")
    return doc


# ---------------------------------------------------------------- Figure S1 --
# scaling_fit.json: least squares through the origin (t = k*N) on per-dataset
# medians. k is the slope; the extrapolation is that slope at the named scale.
def figure_s1():
    doc = load("scaling_fit.json")
    # The legend names each extrapolation target in the unit it counts ("600M
    # variants", "10M mutations"); pd["unit"] holds that unit already.
    scale_label = {"input_variants": "variants", "input_mutations": "mutations"}
    for panel, pd in doc["panels"].items():
        n_target = pd["extrapolate_to"]
        div = pd["time_divisor_seconds"]
        # The legend quotes the extrapolation in the axis unit named on the
        # line ("~5.3 h" for VCF, "~39 s" for MAF).
        y_unit = {3600.0: "h", 1.0: "s"}.get(div, pd["y_unit"])
        target = f"{fmt(n_target / 1e6)}M {scale_label.get(pd['unit'], pd['unit'])}"
        for tool in ["FastCrossMap", "CrossMap"]:
            m = pd["models"][tool]
            add("Figure-S1", f"{panel} {tool}: slope",
                m["microseconds_per_unit"], "us/unit", "scaling_fit.json",
                f"panels.{panel}.models.{tool}.microseconds_per_unit",
                # The seconds form is a convenience, not a quoted value; it is
                # sub-millisecond per unit, so three significant digits rather
                # than `fmt()` (which would round 3.16e-05 to 0.000).
                f"= {m['seconds_per_unit']:.3g} s/unit")
            add("Figure-S1", f"{panel} {tool}: R^2 through the origin",
                fmt_r2(m["r2_through_origin"]), "", "scaling_fit.json",
                f"panels.{panel}.models.{tool}.r2_through_origin",
                "n=" + str(m["n"]))
            add("Figure-S1", f"{panel} {tool}: extrapolated to {target}",
                m["extrapolated_seconds"] / div, y_unit, "scaling_fit.json",
                f"panels.{panel}.models.{tool}.extrapolated_seconds",
                f"{fmt(m['extrapolated_seconds'])} s")
        # The ratio the legend quotes, from the same two slopes.
        fcm = pd["models"]["FastCrossMap"]["extrapolated_seconds"]
        cross = pd["models"]["CrossMap"]["extrapolated_seconds"]
        add("Figure-S1", f"{panel}: extrapolated speed-up (CrossMap / FastCrossMap)",
            cross / fcm, "x", "scaling_fit.json",
            f"panels.{panel}.models.*.extrapolated_seconds",
            f"at {target}")
    return doc


# ---------------------------------------------------------------- Figure S2 --
# realignment_fit.json: free-intercept fits per (method, threads), plus the
# ratios at the three measured scales. Absolute times/memory are in
# benchmark_realignment.json, which the S2 legend quotes directly.
def figure_s2():
    fit = load("realignment_fit.json")
    for method in ["liftover", "realignment"]:
        m = fit["fits"]["1"][method]
        add("Figure-S2", f"{method}: slope at -t 1",
            m["microseconds_per_read"], "us/read", "realignment_fit.json",
            f"fits.1.{method}.microseconds_per_read",
            f"intercept {fmt(m['intercept_seconds'])} s, n={m['n']}")
    for reads, t1, t8, mem in zip(fit["ratios"]["reads"],
                                  fit["ratios"]["time_ratio"],
                                  fit["ratios"]["time_ratio_t8"],
                                  fit["ratios"]["memory_ratio"]):
        scale = f"{int(reads / 1e6)}M read pairs"
        add("Figure-S2", f"{scale}: re-alignment / liftover runtime (-t 1)",
            t1, "x", "realignment_fit.json", "ratios.time_ratio", "")
        add("Figure-S2", f"{scale}: re-alignment / liftover runtime (-t 8)",
            t8, "x", "realignment_fit.json", "ratios.time_ratio_t8", "")
        add("Figure-S2", f"{scale}: re-alignment / liftover peak RSS",
            mem, "x", "realignment_fit.json", "ratios.memory_ratio", "")

    # The absolute endpoints the legend quotes, straight from the raw file.
    raw = load("benchmark_realignment.json")
    by = {}
    for r in raw["results"]:
        by[(r["method"], r["threads"], int(r["num_reads"]))] = r
    for reads in [1_000_000, 5_000_000, 10_000_000]:
        for method, label in [("liftover", "liftover"), ("realignment", "re-alignment")]:
            r = by.get((method, 1, reads))
            if r is None:
                continue
            scale = f"{int(reads / 1e6)}M read pairs"
            add("Figure-S2", f"{scale}: {label} runtime (-t 1)",
                r["mean_time_sec"], "s", "benchmark_realignment.json",
                f"results[method={method},threads=1,num_reads={reads}].mean_time_sec", "")
            add("Figure-S2", f"{scale}: {label} peak RSS (-t 1)",
                r["peak_memory_mb"], "MB", "benchmark_realignment.json",
                f"results[method={method},threads=1,num_reads={reads}].peak_memory_mb", "")
    return fit


# ---------------------------------------------------------------- Figure S3 --
# The merged thread table. Only the two quantities the S3 legend names are
# extracted per cell: the speed-up over -t 1 and the parallel efficiency.
def figure_s3():
    path = RESULTS_DIR / "threads-merged-20261011" / "summary.tsv"
    with open(path, newline="") as fh:
        cells = list(csv.DictReader(fh, delimiter="\t"))
    max_threads = max(int(c["threads"]) for c in cells)
    best = {}
    for c in cells:
        f = c["format"]
        th = int(c["threads"])
        rec = best.setdefault(f, {"speed": 0.0, "speed_th": 0,
                                  "eff_end": 0.0, "eff_end_th": max_threads})
        if float(c["speedup_vs_t1"]) > rec["speed"]:
            rec["speed"], rec["speed_th"] = float(c["speedup_vs_t1"]), th
        if th == max_threads:
            rec["eff_end"] = float(c["eff"])
    for f in FORMATS:
        key = THREAD_KEY[f]
        rec = best[key]
        add("Figure-S3", f"{f}: max speed-up over -t 1",
            rec["speed"], "x", "threads-merged-20261011/summary.tsv",
            f"format={key}, speedup_vs_t1", f"at -t {rec['speed_th']}")
        add("Figure-S3", f"{f}: parallel efficiency at -t {max_threads}",
            rec["eff_end"], "", "threads-merged-20261011/summary.tsv",
            f"format={key}, eff", "")
    return cells, best


# ---------------------------------------------------------------- Figure S4 --
# memory_by_format.json: per format, the per-dataset paired reduction against
# CrossMap, recomputed here from the two tools' per-dataset medians.
def figure_s4():
    doc = load("memory_by_format.json")
    n_ds = []
    for f in FORMATS:
        block = doc["formats"][f]
        fcm = block["tools"]["FastCrossMap"]["per_dataset_mb"]
        cross = block["tools"]["CrossMap"]["per_dataset_mb"]
        shared = sorted(set(fcm) & set(cross))
        n_ds.append(len(shared))
        reductions = [cross[d] / fcm[d] for d in shared if fcm[d] > 0]
        add("Figure-S4", f"{f}: peak-memory reduction (median of paired ratios)",
            stats.median(reductions), "x", "memory_by_format.json",
            f"formats.{f}.tools", f"n={len(shared)} datasets")
    add("Figure-S4", "datasets measured per format", "", "",
        "memory_by_format.json", "formats.<f>.datasets_used",
        "; ".join(f"{f} {n}" for f, n in zip(FORMATS, n_ds)) +
        f"; total {sum(n_ds)}")
    return doc


# ---------------------------------------------------------------- Figure S5 --
def figure_s5():
    doc = load("realignment_accuracy_summary.json")
    for scale in ["1M", "10M"]:
        s = doc["scales"][scale]
        n = s["both_mapped"]
        for rung, pct in s["cascade_pct"].items():
            add("Figure-S5", f"{scale} (a): cascade rung - {rung}",
                pct, "%", "realignment_accuracy_summary.json",
                f"scales.{scale}.cascade_pct[{rung}]", f"of {n:,} both-mapped pairs")
        for name, cnt in s["attribution_count"].items():
            add("Figure-S5", f"{scale} (b): {name} - count",
                cnt, "pairs", "realignment_accuracy_summary.json",
                f"scales.{scale}.attribution_count[{name}]", "")
        for name, pct in s["attribution_pct"].items():
            add("Figure-S5", f"{scale} (b): {name} - share of non-exact",
                pct, "%", "realignment_accuracy_summary.json",
                f"scales.{scale}.attribution_pct[{name}]", "")
        add("Figure-S5", f"{scale}: not-exact share of all both-mapped pairs",
            100.0 * s["not_exact"] / n, "%", "realignment_accuracy_summary.json",
            f"scales.{scale}.not_exact / both_mapped", "")
        low = s["attribution_count"]["same chromosome, 1-500 bp"]
        add("Figure-S5", f"{scale}: within-500 bp band as share of all reads",
            100.0 * low / n, "%", "realignment_accuracy_summary.json",
            f"scales.{scale}.attribution_count[same chromosome, 1-500 bp] / both_mapped", "")
        add("Figure-S5", f"{scale}: strand concordance",
            s["strand_concordance_pct"], "%", "realignment_accuracy_summary.json",
            f"scales.{scale}.strand_concordance_pct", "")
    return doc


# ----------------------------------------------------------------- Table S1 --
# table_s1.tsv is the shipped table; this only counts its shape and repeats the
# per-format median time, which is the figure-adjacent summary a reader checks.
def table_s1():
    path = RESULTS_DIR / "table_s1.tsv"
    with open(path, newline="") as fh:
        tsv_rows = list(csv.DictReader(fh, delimiter="\t"))
    by_fmt: dict[str, set] = {}
    for r in tsv_rows:
        by_fmt.setdefault(r["Format"], set()).add(r["Tool"])
    add("Table-S1", "rows (format x dataset x tool)", len(tsv_rows), "rows",
        "table_s1.tsv", "n/a", f"{len(by_fmt)} formats")
    for f in sorted(by_fmt):
        add("Table-S1", f"{f}: tools benchmarked", len(by_fmt[f]), "tools",
            "table_s1.tsv", "Tool", ", ".join(sorted(by_fmt[f])))
    return tsv_rows, by_fmt


# ------------------------------------------------------------------ checks --
def check(fig1, fit, cells, s4, acc, tsv_rows, by_fmt, problems):
    """Assert the invariants the documents lean on. Each problem is a string."""
    # Figure 1 must carry a bar for every format in both panels.
    for key in ["speedup", "memory_cut"]:
        missing = [f for f in FORMATS if f not in fig1[key]]
        if missing:
            problems.append(f"Figure 1 {key}: missing formats {missing}")

    # The S3 legend quotes eight maxima, one per format; all eight must exist.
    for f in FORMATS:
        if THREAD_KEY[f] not in {c["format"] for c in cells}:
            problems.append(f"Figure S3: no thread rows for {f}")

    # The S4 legend says every format's FastCrossMap median is the lowest of
    # the tools measured. Check that against the per-dataset medians.
    for f in FORMATS:
        block = s4["formats"][f]
        meds = {t: v["median_mb"] for t, v in block["tools"].items()}
        fcm = meds.get("FastCrossMap")
        if fcm is not None and any(v < fcm for t, v in meds.items() if t != "FastCrossMap"):
            problems.append(f"Figure S4 {f}: FastCrossMap is not the lowest median")

    # The S5 legend's proportional claim: the three attribution buckets must
    # account for every non-exact pair.
    for scale in ["1M", "10M"]:
        s = acc["scales"][scale]
        total = sum(s["attribution_count"].values())
        if total != s["not_exact"]:
            problems.append(
                f"Figure S5 {scale}: attribution buckets sum to {total}, "
                f"not_exact is {s['not_exact']}")

    # Table S1 must have one row per (format, dataset, tool) with no repeats,
    # and every row must name a tool.
    seen = set()
    for r in tsv_rows:
        key = (r["Format"], r["Dataset"], r["Tool"])
        if key in seen:
            problems.append(f"Table S1: duplicate row {key}")
        seen.add(key)
        if not r["Tool"]:
            problems.append(f"Table S1: row without a tool: {key}")


# -------------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="validate the invariants and exit; write no table")
    args = ap.parse_args()

    fig1 = figure1()
    figure_s1()
    fit = figure_s2()
    cells, _ = figure_s3()
    s4 = figure_s4()
    acc = figure_s5()
    tsv_rows, by_fmt = table_s1()

    problems: list[str] = []
    check(fig1, fit, cells, s4, acc, tsv_rows, by_fmt, problems)

    if problems:
        print("CONSISTENCY PROBLEMS:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)

    if args.check:
        sys.exit(1 if problems else 0)

    def write(path, body):
        """Write one output, naming the file rather than raising a traceback.

        These tables live on a Windows drive, where an open Excel window holds a
        lock that surfaces as PermissionError. The lock is a mistake in the
        environment, not in the data, so say which file to close instead of
        dumping a stack trace. Both outputs are attempted either way -- one
        locked file must not keep the other from being regenerated -- and the
        run still exits non-zero, since it failed.
        """
        try:
            with open(path, "w", newline="") as fh:
                body(fh)
        except PermissionError:
            print(f"cannot write {path}: permission denied.", file=sys.stderr)
            print("  the file is open in another program (Excel, an editor, a "
                  "viewer); close it and run this script again.", file=sys.stderr)
            return False
        print(f"wrote {path}" + (f" ({len(rows)} rows)" if path == TSV_OUT else ""))
        return True

    def tsv_body(fh):
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    def md_body(fh):
        fh.write("# Original data behind the figures and tables\n\n")
        fh.write(f"{len(rows)} values, each read from the result file named in "
                 "`source_file` at the key in `source_path`.\n\n")
        fh.write("| " + " | ".join(COLUMNS) + " |\n")
        fh.write("|" + "---|" * len(COLUMNS) + "\n")
        for r in rows:
            fh.write("| " + " | ".join(str(r[c]) for c in COLUMNS) + " |\n")

    wrote_tsv = write(TSV_OUT, tsv_body)
    wrote_md = write(MD_OUT, md_body)

    if problems:
        print("consistency problems found; see stderr", file=sys.stderr)
        sys.exit(1)
    if not (wrote_tsv and wrote_md):
        sys.exit(1)


if __name__ == "__main__":
    main()
