#!/usr/bin/env bash
# =============================================================================
# 22_run_benchmarks_pinned.sh
#
# Run the single-format timing benchmarks with each benchmark pinned to its own
# physical core, so the tools inside them never end up sharing a core.
#
# Why this exists
# ---------------
# All eight of these scripts time single-threaded tools (FastCrossMap is invoked
# with `-t 1`; CrossMap and liftOver have no threading at all). Each running
# script therefore needs one core -- but a script also does other work between
# timed runs (counting features, decompressing, warming the page cache), so a
# bare parallel launch oversubscribes the machine.
#
# That is not hypothetical. The 2026-10-07 19:52 GFF run was launched as eight
# concurrent scripts on a 24-core / 125 GB box and its timings were unusable:
#
#   * FastCrossMap took 43.4-59.9 s on GENCODE full, where the same machine had
#     measured 5.20 s by hand on the same file with the page cache warm;
#   * within one five-run group, FastCrossMap spread [36.5 -> 81.2] s and
#     CrossMap [58.3 -> 99.5] s;
#   * FastCrossMap came out *slower* on the 1.11 GB basic file (58.4 s) than on
#     the 1.87 GB full one (55.1 s), and CrossMap was ordered the opposite way.
#
# `uptime` reported a sustained load of ~8.9 -- the eight scripts, one core
# apiece -- and the scripts' own `all_disk_read_kb` was 0 for every run, so the
# page cache was not the cause and the machine was not short of memory. Pinning
# removes the contention *without* serialising the run: the wall-clock cost is
# the same as the unpinned parallel launch, because each tool still gets a whole
# core to itself instead of sharing one with seven neighbours.
#
# Usage
# -----
#   bash paper/22_run_benchmarks_pinned.sh                 # all eight, pinned
#   bash paper/22_run_benchmarks_pinned.sh --dry-run       # print the plan only
#   bash paper/22_run_benchmarks_pinned.sh 02d_benchmark_gff 03_benchmark_bam
#   bash paper/22_run_benchmarks_pinned.sh --with-multithread
#   bash paper/22_run_benchmarks_pinned.sh --jobs 4        # at most 4 at a time
#
# --jobs caps how many run concurrently (using that many physical cores); the
# rest wait their turn. Only needed on a machine whose cores are wanted for
# something else -- the default runs one per physical core.
#
# The two multithread-scalability scripts (02b, 03b) are deliberately NOT pinned:
# they measure how a tool scales across N cores, so confining one to a single
# core would defeat the measurement. Pass --with-multithread to run them after
# the pinned batch, one at a time, once the machine has gone quiet.
#
# Run from anywhere; the script cd's to the repository root itself, because the
# benchmark scripts resolve their inputs as `paper/data/...` relative paths.
# =============================================================================

set -euo pipefail

# --- the eight single-format timing scripts ---------------------------------
DEFAULT_SCRIPTS=(
    02_benchmark_bed
    03_benchmark_bam
    02c_benchmark_vcf
    02d_benchmark_gff
    02e_benchmark_wig
    02f_benchmark_bigwig
    02g_benchmark_maf
    02h_benchmark_sam
)

# Multithread scalability: run unpinned, after the pinned batch.
MULTITHREAD_SCRIPTS=(
    02b_benchmark_bed_multithread
    03b_benchmark_bam_multithread
)

WITH_MULTITHREAD=0
DRY_RUN=0
JOBS=0                 # 0 = one per physical core
SCRIPTS=()
LOGDIR=""

usage() { sed -n '2,50p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

# --- arguments ---------------------------------------------------------------
while [ $# -gt 0 ]; do
    case "$1" in
        --with-multithread) WITH_MULTITHREAD=1 ;;
        --dry-run|--plan)   DRY_RUN=1 ;;
        --jobs)             JOBS="${2:?--jobs needs a number}"; shift ;;
        --logdir)           LOGDIR="${2:?--logdir needs a directory}"; shift ;;
        -h|--help)          usage 0 ;;
        -*)                 echo "unknown option: $1" >&2; usage 1 ;;
        *)                  SCRIPTS+=("$1") ;;
    esac
    shift
done

case "$JOBS" in
    ''|*[!0-9]*) echo "ERROR: --jobs must be a positive integer" >&2; exit 1 ;;
esac
[ "$JOBS" -ge 1 ] || JOBS=0

# --- location ----------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"
[ -n "$LOGDIR" ] || LOGDIR="$REPO_ROOT/paper/results/logs"

if [ "${#SCRIPTS[@]}" -eq 0 ]; then
    SCRIPTS=("${DEFAULT_SCRIPTS[@]}")
fi

# Accept either "02d_benchmark_gff" or "paper/02d_benchmark_gff.py".
norm() {
    local s="${1##*/}"
    s="${s%.py}"
    printf '%s' "$s"
}

# --- prerequisites -----------------------------------------------------------
fail=0
for tool in python3 taskset; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: '$tool' not found." >&2
        [ "$tool" = taskset ] && echo "       taskset ships with util-linux: apt-get install util-linux" >&2
        fail=1
    fi
done
for s in "${SCRIPTS[@]}"; do
    if [ ! -f "paper/${s}.py" ]; then
        echo "ERROR: paper/${s}.py not found (run from the repository root)." >&2
        fail=1
    fi
done
[ "$fail" -eq 0 ] || exit 1

# --- core selection ----------------------------------------------------------
# One logical CPU per *physical* core. Two hyperthreads of the same core are not
# two cores: pinning two single-threaded tools to siblings leaves them sharing
# the core's execution units, which is the very contention this script exists to
# remove. lscpu -p prints the columns in the order requested, so column 3 is CPU.
mapfile -t CORES < <(lscpu -p=CORE,SOCKET,CPU 2>/dev/null | grep -v '^#' \
                     | awk -F, '!seen[$1 "," $2]++ {print $3}')
if [ "${#CORES[@]}" -eq 0 ]; then
    # No usable lscpu: fall back to every logical CPU (a hyperthreaded box then
    # risks sibling sharing, so say so rather than pretend it is fine).
    mapfile -t CORES < <(seq 0 $(( $(nproc 2>/dev/null || echo 1) - 1 )))
    echo "WARNING: could not read the CPU topology; falling back to all" \
         "${#CORES[@]} logical CPUs." >&2
    echo "         On a hyperthreaded machine some pinned pairs may share a core." >&2
fi
if [ "$JOBS" -ge 1 ] && [ "$JOBS" -lt "${#CORES[@]}" ]; then
    CORES=("${CORES[@]:0:$JOBS}")
fi
NCORES="${#CORES[@]}"

# --- plan --------------------------------------------------------------------
echo "=============================================================================="
echo " Pinned single-format benchmark run"
echo "=============================================================================="
echo " repo root : $REPO_ROOT"
echo " logs      : $LOGDIR"
echo " physical cores available : $NCORES  (${CORES[*]})"
echo " scripts   : ${#SCRIPTS[@]}"
echo

# Waves: at most one script per physical core at a time.
WAVES=$(( (${#SCRIPTS[@]} + NCORES - 1) / NCORES ))
for ((w = 0; w < WAVES; w++)); do
    printf ' wave %d:' "$((w + 1))"
    for ((j = 0; j < NCORES; j++)); do
        idx=$(( w * NCORES + j ))
        [ "$idx" -lt "${#SCRIPTS[@]}" ] || break
        printf ' %s->cpu%s' "${SCRIPTS[$idx]}" "${CORES[$j]}"
    done
    echo
done
if [ "$WITH_MULTITHREAD" -eq 1 ]; then
    echo " then (unpinned, after the above): ${MULTITHREAD_SCRIPTS[*]}"
fi
echo

# A machine that is already busy will still distort the run; that is worth
# knowing before hours are spent, not after.
load1="$(awk '{print $1}' /proc/loadavg 2>/dev/null || echo 0)"
echo "current load average (1 min): $load1   (want it well under $NCORES)"
if awk -v l="$load1" -v n="$NCORES" 'BEGIN{exit !(l > n/2)}'; then
    echo "WARNING: something else is already using this machine; the timings may"
    echo "         still be distorted. Consider waiting for it to finish."
fi
echo

if [ "$DRY_RUN" -eq 1 ]; then
    echo "--dry-run: nothing executed."
    exit 0
fi

mkdir -p "$LOGDIR"

# Used by the report below to tell this run's JSONs from stale ones.
START_TS="$(date +%s)"

# --- clean up children on interruption ---------------------------------------
PIDS=()
cleanup() {
    trap - INT TERM
    echo
    echo "interrupted: stopping ${#PIDS[@]} running script(s)..."
    for pid in "${PIDS[@]}"; do kill "$pid" 2>/dev/null || true; done
    wait 2>/dev/null || true
    exit 130
}
trap cleanup INT TERM

# --- run ---------------------------------------------------------------------
run_one() {   # run_one <script> <cpu> [<pin:1|0>]
    local script="$1" cpu="$2" pin="${3:-1}"
    local log="$LOGDIR/${script}.log"
    echo "  -> $script  (log: $log)"
    if [ "$pin" -eq 1 ]; then
        taskset -c "$cpu" python3 "paper/${script}.py" >"$log" 2>&1
    else
        python3 "paper/${script}.py" >"$log" 2>&1
    fi
}

declare -A RC=()      # script -> exit status
declare -A FAILED=()  # script -> 1 when the run failed

for ((w = 0; w < WAVES; w++)); do
    PIDS=()
    WAVE_SCRIPTS=()
    echo "--- wave $((w + 1))/$WAVES ---"
    for ((j = 0; j < NCORES; j++)); do
        idx=$(( w * NCORES + j ))
        [ "$idx" -lt "${#SCRIPTS[@]}" ] || break
        script="${SCRIPTS[$idx]}"; cpu="${CORES[$j]}"
        WAVE_SCRIPTS+=("$script")
        run_one "$script" "$cpu" 1 &
        PIDS+=($!)
        printf '  core %-3s %s (pid %s)\n' "$cpu" "$script" "$!"
    done
    echo "  waiting..."
    for k in "${!PIDS[@]}"; do
        if wait "${PIDS[$k]}"; then RC["${WAVE_SCRIPTS[$k]}"]=0
        else RC["${WAVE_SCRIPTS[$k]}"]=$?; FAILED["${WAVE_SCRIPTS[$k]}"]=1; fi
    done
    echo
done

if [ "$WITH_MULTITHREAD" -eq 1 ]; then
    echo "--- multithread scalability (unpinned, sequential) ---"
    for script in "${MULTITHREAD_SCRIPTS[@]}"; do
        [ -f "paper/${script}.py" ] || { echo "  skip ${script} (not found)"; continue; }
        if run_one "$script" 0 0; then RC["$script"]=0
        else RC["$script"]=$?; FAILED["$script"]=1; fi
    done
    echo
fi

trap - INT TERM

# --- report ------------------------------------------------------------------
echo "=============================================================================="
echo " Results"
echo "=============================================================================="
for script in "${SCRIPTS[@]}" "${MULTITHREAD_SCRIPTS[@]}"; do
    [ -n "${RC[$script]+x}" ] || continue
    if [ "${FAILED[$script]:-0}" -eq 1 ]; then
        printf '  %-32s EXIT %s   <- see %s\n' "$script" "${RC[$script]}" \
               "$LOGDIR/${script}.log"
    else
        printf '  %-32s ok\n' "$script"
    fi
done
echo

# Per-row summary straight from the JSONs this run wrote: the table above shows
# each script's exit status, this shows whether the harness recorded the
# individual rows as successes and whether any run still paid for a disk read.
# Only files newer than the run start are read, so a stale JSON from an earlier
# session cannot be mistaken for this run's result.
# `|| rc=$?` rather than a bare `rc=$?`: under `set -e` a non-zero exit from the
# check would end the script before the summary below is printed.
rc=0
python3 - "$REPO_ROOT" "$START_TS" "${#SCRIPTS[@]}" <<'PY' || rc=$?
import glob, json, os, sys, time
root, start_ts, n_scripts = sys.argv[1], float(sys.argv[2]) - 2, int(sys.argv[3])

all_paths = sorted(glob.glob(os.path.join(root, "paper/results/benchmark_*.json")))
paths = [p for p in all_paths if os.path.getmtime(p) >= start_ts]
if not paths:
    print(f"no benchmark_*.json was written or updated by this run "
          f"({len(all_paths)} older file(s) present). Did the scripts get that far?")
    sys.exit(1)

# liftOver's -gff parser rejects the upstream RefSeq file on its 'Curated
# Genomic' source column. That failure is a property of the input, not of this
# harness, so it is reported but does not by itself fail the run.
def is_known_expected(r):
    return (r.get("tool") == "liftOver"
            and "RefSeq" in str(r.get("dataset_name", ""))
            and "Expecting number" in str(r.get("error_message", "")))

unexpected = 0
uncached_total = 0
expected = 0
for p in paths:
    try:
        d = json.load(open(p))
    except Exception as e:
        print(f"{os.path.basename(p)}: unreadable ({e})")
        unexpected += 1
        continue
    rows = d.get("results", [])
    bad = [r for r in rows if not r.get("success")]
    uncached = [(r.get("dataset_name") or r.get("tool"), r.get("all_disk_read_kb"))
                for r in rows
                if any((k or 0) > 1024 for k in (r.get("all_disk_read_kb") or []))]
    now_expected = sum(1 for r in bad if is_known_expected(r))
    expected += now_expected
    unexpected += len(bad) - now_expected
    uncached_total += len(uncached)
    flags = []
    if bad:      flags.append(f"{len(bad)} failed")
    if uncached: flags.append(f"{len(uncached)} read from disk")
    print(f"{os.path.basename(p)}: {len(rows)} rows"
          + (", " + ", ".join(flags) if flags else ""))
    for r in bad:
        tag = "expected" if is_known_expected(r) else "FAILED  "
        print(f"    {tag} {r.get('dataset_name')} / {r.get('tool')}: "
              f"{str(r.get('error_message')).splitlines()[0][:70]}")
    for name, reads in uncached:
        print(f"    UNCACHED {name}: all_disk_read_kb={reads}")

print()
if expected:
    print(f"{expected} known-expected failure(s): liftOver's -gff parser rejects the")
    print("upstream RefSeq file's 'Curated Genomic' source column. Not a defect here.")
print(f"rows checked in {len(paths)} JSON file(s) written by this run "
      f"({n_scripts} script(s) requested).")
if unexpected or uncached_total:
    print(f"{unexpected} unexpected failure(s), {uncached_total} uncached row(s).")
sys.exit(1 if (unexpected or uncached_total) else 0)
PY
echo
if [ "$rc" -eq 0 ]; then
    echo "OK: no unexpected failures and no uncached rows in what this run wrote."
else
    echo "PROBLEM: see the FAILED/UNCACHED lines above; do not quote those rows"
    echo "         before checking the matching log in $LOGDIR."
fi
echo "Logs: $LOGDIR"
exit "$rc"
