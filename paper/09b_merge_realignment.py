#!/usr/bin/env python3
"""
09b_merge_realignment.py - Merge separate `09_benchmark_realignment.py` runs
into one results JSON.

Why this exists: 09 writes its output with `output_file.write_text(...)`, so a
run covers only the scales it was asked for and *replaces* the previous file
rather than adding to it. That is fine for a single full run, but the published
`benchmark_realignment.json` (10M reads, t=1/4/8, 2026-10-07) is hours of
server time, and adding the missing 5M scale with `--scales 5M` alone would
overwrite it. The safe sequence is therefore:

    cp paper/results/benchmark_realignment.json \\
       paper/results/benchmark_realignment.json.bak
    python3 paper/09_benchmark_realignment.py --scales 5M --threads 1,4,8
    python3 paper/09b_merge_realignment.py

which reads the 10M backup and the fresh 5M file and writes the union back to
`benchmark_realignment.json`.

Usage:
    python3 paper/09b_merge_realignment.py
    python3 paper/09b_merge_realignment.py --base OLD.json --add NEW.json -o OUT.json
    python3 paper/09b_merge_realignment.py --prefer add     # new wins on conflict
    python3 paper/09b_merge_realignment.py --allow-failures # merge failed (zero) rows too

Exit status is non-zero when a conflict or a failed row is found and the merge
was refused, so a `set -e` pipeline stops instead of silently publishing.
"""

import argparse
import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

RESULTS_DIR = Path("paper/results")
DEFAULT_BASE = RESULTS_DIR / "benchmark_realignment.json.bak"
DEFAULT_ADD = RESULTS_DIR / "benchmark_realignment.json"
DEFAULT_OUT = RESULTS_DIR / "benchmark_realignment.json"

# Liftover is the derived number (it is what the speedup is computed against),
# so list the reference method first in each (scale, threads) group.
METHOD_ORDER = {"realignment": 0, "liftover": 1}


def load(path: Path) -> dict:
    if not path.exists():
        sys.exit(f"ERROR: {path} not found")
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as e:
        sys.exit(f"ERROR: {path} is not valid JSON: {e}")


def entry_key(e: dict) -> tuple:
    """Identity of a measurement: one method, one scale, one thread count."""
    return (e.get("method"), e.get("scale"), e.get("threads"))


def describe(e: dict) -> str:
    return (f"{e.get('scale')} t={e.get('threads')} {e.get('method')}: "
            f"mean={e.get('mean_time_sec')}s sd={e.get('std_time_sec')} "
            f"mem={e.get('peak_memory_mb')}MB n={len(e.get('all_times') or [])}")


def check_failures(name: str, doc: dict) -> list:
    """Rows recorded as zeros by a run whose configurations all failed.

    `09` prints `!! Some configurations FAILED and were recorded as zeros.` and
    writes `success: false` with a zero mean. Those are not measurements, and
    letting one through the merge is how a zero ends up in a figure.
    """
    bad = [e for e in doc.get("results", []) if not e.get("success", True)]
    if bad:
        print(f"WARNING: {name} has {len(bad)} failed row(s) (zeros):")
        for e in bad:
            print(f"    {describe(e)}  error={str(e.get('error'))[:70]}")
    return bad


def sort_key(e: dict) -> tuple:
    # num_reads is written by 09, so scales order as 1M/5M/10M without parsing
    # the label; fall back to the label for a hand-edited file.
    reads = e.get("num_reads") or 0
    scale = e.get("scale") or ""
    return (reads, scale, e.get("threads") or 0,
            METHOD_ORDER.get(e.get("method"), 9))


def merge(base: dict, add: dict, prefer: str) -> tuple[list, list]:
    """Union of results keyed by (method, scale, threads). Returns (merged, conflicts).

    A scale that only one side measured is simply carried over; a key both sides
    measured is a conflict unless the two rows are byte-identical, and is
    reported rather than silently picked, because the two runs are hours apart
    and a disagreement there is worth a look.
    """
    merged: dict[tuple, dict] = {}
    conflicts: list = []

    for doc in (base, add):
        for e in doc.get("results", []):
            k = entry_key(e)
            if k not in merged:
                merged[k] = e
                continue
            if json.dumps(merged[k], sort_keys=True) != json.dumps(e, sort_keys=True):
                conflicts.append((k, merged[k], e))
            if prefer == "add":       # same key in both: the later run wins if asked
                merged[k] = e

    return sorted(merged.values(), key=sort_key), conflicts


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Merge two benchmark_realignment.json runs (e.g. 10M + 5M)")
    ap.add_argument("--base", type=Path, default=DEFAULT_BASE,
                    help=f"existing/authoritative JSON (default: {DEFAULT_BASE})")
    ap.add_argument("--add", type=Path, default=DEFAULT_ADD,
                    help=f"freshly produced JSON (default: {DEFAULT_ADD})")
    ap.add_argument("-o", "--out", type=Path, default=DEFAULT_OUT,
                    help=f"where to write the union (default: {DEFAULT_OUT})")
    ap.add_argument("--prefer", choices=("base", "add"), default="base",
                    help="which side wins when the same (method, scale, threads) "
                         "appears in both (default: base)")
    ap.add_argument("--allow-failures", action="store_true",
                    help="merge even if a side has `success: false` (zero) rows")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the union without writing anything")
    args = ap.parse_args()

    base, add = load(args.base), load(args.add)
    print(f"base: {args.base}  ({len(base.get('results', []))} rows, "
          f"scales={base.get('scales')}, threads={base.get('thread_counts')}, "
          f"date={base.get('date')})")
    print(f"add : {args.add}  ({len(add.get('results', []))} rows, "
          f"scales={add.get('scales')}, threads={add.get('thread_counts')}, "
          f"date={add.get('date')})")

    bad = check_failures("base", base) + check_failures("add", add)
    if bad and not args.allow_failures:
        print("\nRefusing to merge: a failed row would be published as a measurement.")
        print("Re-run that configuration, or pass --allow-failures if the zeros are "
              "known and you are merging deliberately.")
        return 2

    merged, conflicts = merge(base, add, args.prefer)

    if conflicts:
        print(f"\n{len(conflicts)} conflicting row(s) (same method/scale/threads, "
              f"different values); --prefer {args.prefer} decided:")
        for k, b, a in conflicts:
            print(f"  {k}:")
            print(f"    base: {describe(b)}")
            print(f"    add : {describe(a)}")

    print(f"\nmerged: {len(merged)} rows")
    print(f"{'scale':<7}{'t':<4}{'method':<14}{'mean(s)':>11}{'sd':>8}{'mem(MB)':>10}  all_times")
    for e in merged:
        print(f"{str(e.get('scale')):<7}{e.get('threads'):<4}{e.get('method'):<14}"
              f"{e.get('mean_time_sec', 0):>11}{e.get('std_time_sec', 0):>8}"
              f"{e.get('peak_memory_mb', 0):>10}  {e.get('all_times')}")

    # Speedup per (scale, threads), the number the paper actually quotes.
    print(f"\n{'scale':<7}{'t':<4}{'realign(s)':>12}{'liftover(s)':>13}{'speedup':>10}")
    by_key = {(e.get("scale"), e.get("threads"), e.get("method")): e for e in merged}
    scales = sorted({e.get("scale") for e in merged},
                    key=lambda s: next((e.get("num_reads", 0)
                                        for e in merged if e.get("scale") == s), 0))
    for s in scales:
        for t in sorted({e.get("threads") for e in merged if e.get("scale") == s}):
            r = by_key.get((s, t, "realignment"))
            l = by_key.get((s, t, "liftover"))
            # A failed configuration (no bwa, disk full, timeout) is written as
            # a zero row, and 09 skips the liftover side when re-alignment
            # failed. Print the gap explicitly -- a missing line is easy to
            # read as "this scale was never run".
            why = None
            if r is None or l is None:
                why = f"missing {'liftover' if l is None else 'realignment'} row"
            elif not (r.get("success", True) and l.get("success", True)):
                why = "failed row (zeros)"
            elif not l.get("mean_time_sec"):
                why = "liftover mean is zero"
            if why:
                print(f"{str(s):<7}{t:<4}{'N/A':>12}{'N/A':>13}{'N/A':>10}   <- {why}")
                continue
            sp = r["mean_time_sec"] / l["mean_time_sec"]
            print(f"{str(s):<7}{t:<4}{r['mean_time_sec']:>12.1f}"
                  f"{l['mean_time_sec']:>13.2f}{sp:>9.1f}x")

    if args.dry_run:
        print("\n--dry-run: nothing written")
        return 0

    if args.out.exists() and not args.out.with_suffix(args.out.suffix + ".premerge").exists():
        shutil.copy2(args.out, args.out.with_suffix(args.out.suffix + ".premerge"))
        print(f"\n(kept the pre-merge file as {args.out.with_suffix(args.out.suffix + '.premerge')})")

    out_doc = {
        "benchmark": base.get("benchmark", "realignment_vs_liftover"),
        "date": datetime.now().isoformat(),
        "merged_from": [str(args.base), str(args.add)],
        "num_runs": max(base.get("num_runs", 0), add.get("num_runs", 0)),
        "scales": scales,
        "thread_counts": sorted({e.get("threads") for e in merged}),
        "results": merged,
    }
    args.out.write_text(json.dumps(out_doc, indent=2))
    print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
