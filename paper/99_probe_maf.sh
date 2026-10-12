#!/usr/bin/env bash
#
# 99_probe_maf.sh —— 临时探针: MAF 的 -t 2 倒退, 到底是什么开销
#
# 这是一次性诊断脚本, 不是正式 benchmark。用完可以删。
#
# ---------------------------------------------------------------------------
# 要回答的问题
# ---------------------------------------------------------------------------
# 服务器结果里 maf 是唯一一个 -t 2 仍然变慢的格式 (t=1 1.900s → t=2 2.466s, 0.77x),
# 且 5 次重复全部复现 (2.52/2.47/2.46/2.42/2.46)。
#
# 但它的形态和当初的 wig/bigwig **不同**:
#
#   t    wall    cpu(cpu/wall)
#   1    1.902   1.894  (1.00)   <- 纯串行
#   2    2.466   3.484  (1.41)
#   4    1.736   3.300  (1.90)
#   8    1.414   3.282  (2.32)
#   16   1.288   3.742  (2.91)   <- 可并行上限约 3 核
#
# 两个特征:
#   (1) 从 t=1 跨到 t=2, 总 CPU 凭空多出 ~1.6s, 而 wall 也涨 —— 多出来的是开销不是功;
#       之后 t=4..16 的 CPU 稳定在 3.3-3.7, 说明这 ~1.5s 与线程数无关。
#   (2) 可并行占比很小 (cpu/wall 在 16 线程才到 2.9)。
#
# 关键事实: maf 的并发粒度是 maf.rs 里**硬编码**的 `const CHUNK_SIZE = 10_000`
#   (不像 wig/bigwig 吃 FCM_BATCH), 所以 500,001 条记录的 maf_synthetic_large 只触发
#   **约 50 次 `pool.install`** —— 而 wig 是 12,200 次。50 次 install 解释不了 1.5s,
#   所以 "install 过密" 这个 (在 wig 上成立的) 假设对 maf **不成立**。
#
# 于是还剩两种互斥的解释, 本探针用**现有 5 个尺寸**把它们分开:
#
#   A. 每次 install 的固定开销 —— 额外 CPU 随记录数 (≈ install 次数) 线性增长。
#      tiny→xlarge 记录数差 20 倍 (install 次数约 5→100), 额外 CPU 也应当涨约 20 倍。
#   B. 一次性固定成本 (线程池创建 / 预热 / 分配器首次多线程化) —— 额外 CPU 与记录数
#      无关, 是一条水平线。
#
# 判读: 看 pivot.tsv 的 "extra_cpu" 列 ( = 该尺寸 t 的 CPU 中位数 − 该尺寸 t=1 的 CPU
#   中位数)。随尺寸上涨 => A; 基本不变 => B。两者都不成立 (比如额外 CPU 不随尺寸涨
#   也不恒定) 也会在表里直接看出来。
#
# 注意: 本探针**不改源码**, 也不测 chunk 大小 —— CHUNK_SIZE 没有环境变量开关。如果结论
#   指向 A, 才值得再讨论「给它加一个临时 env 开关做 sweep」(那和当初 wig 的 FCM_BATCH
#   一样是决定性的), 或直接把 CHUNK_SIZE 调大复测。
#
# ---------------------------------------------------------------------------
# 用法 (仓库根目录, 建议在服务器上跑)
# ---------------------------------------------------------------------------
#   bash paper/99_probe_maf.sh
#
# 约 6 分钟 (5 个尺寸 × 4 个线程数 × 5 次重复, 每次 0.2-3.8s)。想更快:
#   THREADS="1 2" REPS=3 bash paper/99_probe_maf.sh              # 约 1 分钟
#   只跑部分尺寸: MAF_SIZES="large xlarge" bash paper/99_probe_maf.sh
#
# 环境变量 (未列出的与 20_benchmark_threads.sh 同名同义):
#   MAF_SIZES 要测的数据集 (默认: tiny small medium large xlarge)
#   THREADS REPS WARM SYNC KEEP OUT DATA FCM
#   MAF_BUILD MAF_REF (hg38.fa)
#
# 输出:
#   $OUT/timings.tsv  原始逐次数据 (size, records, threads, rep, wall, user, sys, rss, flush, ok)
#   $OUT/summary.tsv  每个 (size, threads) 的中位数 + cpu/wall 比
#   $OUT/pivot.tsv    焦点表: 每个 (size, threads) 一行, 含 extra_cpu (相对同尺寸 t=1)
#   $OUT/sha256.txt   一致性抽查 (同一 size, 不同 t 的产物哈希应完全相同)
#   $OUT/env.txt      运行环境
#
set -u

FCM="${FCM:-./target/release/fast-crossmap}"
DATA="${DATA:-paper/data}"
OUT="${OUT:-paper/results/probe-maf-$(date +%Y%m%d-%H%M%S)}"
THREADS="${THREADS:-1 2 4 8}"
REPS="${REPS:-5}"
WARM="${WARM:-1}"
SYNC="${SYNC:-1}"
KEEP="${KEEP:-0}"

MAF_SIZES="${MAF_SIZES:-tiny small medium large xlarge}"
MAF_BUILD="${MAF_BUILD:-GRCh38}"
CHAIN="$DATA/hg19ToHg38.over.chain.gz"
REF="${MAF_REF:-$DATA/hg38.fa}"

mkdir -p "$OUT"

[ -x "$FCM" ]   || { echo "ERROR: 找不到 $FCM (先 cargo build --release)" >&2; exit 1; }
[ -f "$CHAIN" ] || { echo "ERROR: 找不到 chain 文件 $CHAIN" >&2; exit 1; }
[ -f "$REF" ]   || { echo "ERROR: 找不到参考基因组 $REF" >&2; exit 1; }

# ---- 数据集 (只看存在的) ----------------------------------------------------
SIZES=()
declare -A IN_OF
for s in $MAF_SIZES; do
    f="$DATA/maf_synthetic_${s}.maf"
    if [ -f "$f" ]; then SIZES+=("$s"); IN_OF[$s]="$f"
    else echo "  skip $s: 无输入 $f" >&2; fi
done
[ "${#SIZES[@]}" -gt 0 ] || { echo "ERROR: 没有可跑的数据集" >&2; exit 1; }

# 记录数 (非 # 开头的行) 和它对应的 install 次数 —— 就是判读 A/B 的自变量。
# 每个尺寸只数一次: 输入可达 97 MB, 每次 measure_once 都 grep 一遍很浪费。
declare -A REC_OF
records_of() { echo "${REC_OF[$1]:-0}"; }
for s in "${SIZES[@]}"; do
    REC_OF[$s]=$(grep -vc '^#' "${IN_OF[$s]}" 2>/dev/null || echo 0)
done

# maf 的产物: 输出路径按原样使用 + 扩展名换成 maf.unmap。
artifacts_of() {  # artifacts_of <prefix>
    printf '%s\n' "$1" "${1%.*}.maf.unmap"
}
clean_outputs() {  # clean_outputs <prefix>
    local f; while IFS= read -r f; do rm -f "$f"; done < <(artifacts_of "$1")
}
prefix_of() { echo "$OUT/$1.out.maf"; }

FREE_KB=$(df -Pk "$OUT" 2>/dev/null | awk 'NR==2{print $4}')
if [ -n "${FREE_KB:-}" ] && [ "$FREE_KB" -lt $((2 * 1024 * 1024)) ]; then
    echo "WARNING: $OUT 剩余 $((FREE_KB / 1024 / 1024)) GiB, 可能不够 (xlarge 产物约 100 MB)" >&2
fi
case "$(df -T "$OUT" 2>/dev/null | awk 'NR==2{print $2}')" in
    9p|nfs|nfs4|cifs|fuse*) echo "WARNING: OUT 在网络文件系统上, wall 会被文件系统支配, 结果不可用" >&2 ;;
esac

if [ -x /usr/bin/time ]; then HAVE_GNU_TIME=1; else HAVE_GNU_TIME=0; fi

# ---- 展开 cell 列表: (size, threads), 每轮整体轮转一位 ----------------------
declare -a CELL_S CELL_T
for s in "${SIZES[@]}"; do
    for t in $THREADS; do CELL_S+=("$s"); CELL_T+=("$t"); done
done
NC=${#CELL_S[@]}

ENV_F="$OUT/env.txt"
{
    echo "date=$(date -Is)"
    echo "host=$(hostname)  nproc=$(nproc 2>/dev/null || echo '?')"
    echo "kernel=$(uname -sr)"
    echo "FCM=$FCM"
    "$FCM" --version 2>/dev/null | sed 's/^/version: /'
    echo "PROBE=maf_chunk  FORMAT=maf  MAF_BUILD=$MAF_BUILD"
    echo "MAF_SIZES=${SIZES[*]}  THREADS=$THREADS  REPS=$REPS  WARM=$WARM  SYNC=$SYNC"
    echo "gnu_time=$HAVE_GNU_TIME  cells=$NC"
    echo "OUT=$OUT (mount: $(df -T "$OUT" 2>/dev/null | awk 'NR==2{print $1" "$2}'))"
    echo "DATA=$DATA (mount: $(df -T "$DATA" 2>/dev/null | awk 'NR==2{print $1" "$2}'))"
    echo "dirty_ratio=$(cat /proc/sys/vm/dirty_ratio 2>/dev/null)  dirty_background_ratio=$(cat /proc/sys/vm/dirty_background_ratio 2>/dev/null)"
    echo "memtotal_kb=$(awk '/MemTotal/{print $2}' /proc/meminfo)"
    echo "loadavg_at_start=$(cat /proc/loadavg)"
    for s in "${SIZES[@]}"; do
        echo "input[$s]=${IN_OF[$s]} ($(stat -c%s "${IN_OF[$s]}" 2>/dev/null) bytes, $(records_of "$s") records)"
    done
} | tee "$ENV_F"

TIMINGS="$OUT/timings.tsv"; : > "$TIMINGS"
printf "size\trecords\tthreads\trep\twall_s\tuser_s\tsys_s\tmax_rss_kb\tflush_s\tok\n" >> "$TIMINGS"
: > "$OUT/sha256.txt"

# measure_once <size> <t> -> 回显 size<TAB>records<TAB>t<TAB>wall<TAB>user<TAB>sys<TAB>rss<TAB>ok
measure_once() {
    local s="$1" t="$2"
    local prefix; prefix=$(prefix_of "$s")
    clean_outputs "$prefix"
    local wall user sys rss rc
    if [ "$HAVE_GNU_TIME" = 1 ]; then
        local tf="$OUT/.time.$$"
        /usr/bin/time -f '%e\t%U\t%S\t%M' -o "$tf" \
            "$FCM" maf -b "$MAF_BUILD" "$CHAIN" "${IN_OF[$s]}" "$REF" "$prefix" -t "$t" \
            >/dev/null 2>&1
        rc=$?
        IFS=$'\t' read -r wall user sys rss < "$tf" || { wall=0; user=0; sys=0; rss=0; }
        rm -f "$tf"
    else
        local st e
        st=$(date +%s.%N)
        "$FCM" maf -b "$MAF_BUILD" "$CHAIN" "${IN_OF[$s]}" "$REF" "$prefix" -t "$t" >/dev/null 2>&1; rc=$?
        e=$(date +%s.%N)
        wall=$(awk -v a="$st" -v b="$e" 'BEGIN{printf "%.3f", b-a}')
        user=0; sys=0; rss=0
    fi
    local ok=0; [ "$rc" -eq 0 ] && ok=1
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$s" "$(records_of "$s")" "$t" "$wall" "$user" "$sys" "$rss" "$ok"
}

flush_dirty() {
    local s e; s=$(date +%s.%N); sync; e=$(date +%s.%N)
    awk -v a="$s" -v b="$e" 'BEGIN{printf "%.3f", b-a}'
}

echo
echo "=================== 探针: MAF 的并行开销 ==================="
echo "  数据集:  ${SIZES[*]}"
echo "  线程数:  $THREADS   重复: $REPS"
echo

# ============================ 热身 ============================
for s in "${SIZES[@]}"; do
    for _ in $(seq 1 "$WARM"); do
        clean_outputs "$(prefix_of "$s")"
        "$FCM" maf -b "$MAF_BUILD" "$CHAIN" "${IN_OF[$s]}" "$REF" "$(prefix_of "$s")" -t 1 >/dev/null 2>&1 || true
    done
    printf "  warm %-8s ok (%s records)\n" "$s" "$(records_of "$s")"
done
sync

# ============================ 轮转采样 ============================
echo
echo "=================== 采样 ==================="
for r in $(seq 1 "$REPS"); do
    for ci in $(seq 0 $((NC - 1))); do
        idx=$(( (ci + r - 1) % NC ))
        s="${CELL_S[$idx]}"; t="${CELL_T[$idx]}"

        [ "$SYNC" = 1 ] && sync
        row=$(measure_once "$s" "$t")
        if [ "$SYNC" = 1 ]; then flush=$(flush_dirty); else flush="-"; fi

        printf "%s\t%s\t%s\n" "$row" "$r" "$flush" \
            | awk -F'\t' '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$9,$4,$5,$6,$7,$10,$8}' \
            >> "$TIMINGS"

        printf "  r=%s  %-8s t=%-3s wall=%7ss  cpu=%5ss  rss=%6sKB  flush=%7ss  ok=%s\n" \
            "$r" "$s" "$t" \
            "$(printf '%s' "$row" | cut -f4)" \
            "$(awk -F'\t' '{printf "%.2f", $5+$6}' <<<"$row")" \
            "$(printf '%s' "$row" | cut -f7)" \
            "$flush" \
            "$(printf '%s' "$row" | cut -f8)"

        if [ "$r" -eq "$REPS" ]; then
            prefix=$(prefix_of "$s")
            while IFS= read -r f; do
                [ -f "$f" ] || continue
                sha256sum "$f" 2>/dev/null \
                    | awk -v s="$s" -v t="$t" -v a="$(basename "$f")" \
                          '{print s"\tt="t"\t"a"\t"$1}' >> "$OUT/sha256.txt"
            done < <(artifacts_of "$prefix")
            [ "$KEEP" = "1" ] || clean_outputs "$prefix"
        fi
    done
done

# ============================ 汇总 ============================
SUMMARY="$OUT/summary.tsv"
printf "size\trecords\tthreads\tmedian_wall\tmin_wall\tmax_wall\tmedian_cpu\tcpu_per_wall\tmax_rss_kb\tn\n" > "$SUMMARY"
awk -F'\t' '
    function median(src, m,   i) { return (m%2)?src[(m+1)/2]:(src[m/2]+src[m/2+1])/2 }
    function sort(arr, m,   i,j,t) {
        for(i=1;i<=m;i++) for(j=i+1;j<=m;j++) if(arr[j]<arr[i]){t=arr[i];arr[i]=arr[j];arr[j]=t}
    }
    NR>1 && $10==1 {
        k=$1"\t"$3; n[k]++; w[k,n[k]]=$5+0; c[k,n[k]]=($6+0)+($7+0); r[k,n[k]]=$8+0; rec[k]=$2
        if(!(k in lo)||$5+0<lo[k]) lo[k]=$5+0
        if(!(k in hi)||$5+0>hi[k]) hi[k]=$5+0
    }
    END{
        for(k in n){
            m=n[k]; delete a; delete b; delete rr
            for(i=1;i<=m;i++){ a[i]=w[k,i]; b[i]=c[k,i]; rr[i]=r[k,i] }
            sort(a,m); sort(b,m); sort(rr,m)
            mw=median(a,m); mc=median(b,m); mr=median(rr,m)
            split(k,p,"\t")
            printf "%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%d\t%d\n", \
                p[1], rec[k], p[2], mw, lo[k], hi[k], mc, (mw>0? mc/mw:0), mr, m
        }
    }' "$TIMINGS" | sort -k2,2n -k3,3n >> "$SUMMARY"

# 焦点表: 每个 (size, t) 一行, 外加 extra_cpu = 该 (size,t) 的 CPU 中位数 − 该 size 的 t=1 CPU 中位数。
PIVOT="$OUT/pivot.tsv"
printf "size\trecords\tthreads\tmedian_wall\tmedian_cpu\tcpu_per_wall\textra_cpu_vs_t1\n" > "$PIVOT"
awk -F'\t' '
    NR==1{next}
    { size=$1; rec=$2; t=$3+0; med[size"\t"t]=$4; cpu[size"\t"t]=$7; rec2[size]=rec }
    END{
        for(k in med){
            split(k,p,"\t"); size=p[1]; t=p[2]+0
            base=cpu[size"\t1"]
            extra=(t==1)?0:(base>0? cpu[k]-base : 0)
            printf "%s\t%s\t%d\t%.3f\t%.3f\t%.2f\t%+.3f\n", size, rec2[size], t, med[k], cpu[k], (med[k]>0?cpu[k]/med[k]:0), extra
        }
    }' "$SUMMARY" | sort -k2,2n -k3,3n >> "$PIVOT"

echo
echo "=================== 焦点表 ==================="
column -t "$PIVOT"
echo
echo "=================== 汇总 ==================="
column -t "$SUMMARY"
echo
echo "原始数据:   $TIMINGS"
echo "哈希抽查:   $OUT/sha256.txt  (同一 size 下不同 t 应完全相同)"
echo "环境:       $ENV_F"
echo
echo "怎么读:"
echo "  1) extra_cpu_vs_t1 列: 看它随 records 增长还是不变。"
echo "       随 records (≈install 次数) 线性增长 → 每次 install 的固定开销。"
echo "       t=2 起基本恒定 (一条水平线)   → 一次性固定成本 (池创建/预热/分配器)。"
echo "  2) cpu_per_wall 列: t=16 时的值就是可并行上限的核数 (若约 2.9, 说明只有 ~3 核的功)。"
echo "  3) 若 t=2 的 extra_cpu 在最小尺寸上就已接近大尺寸的值 => 与记录数无关 => 一次性成本。"
echo "  4) sha256.txt 里同一 size 各 t 的哈希若一致, 说明线程数不改输出字节。"
