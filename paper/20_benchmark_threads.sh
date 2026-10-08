#!/usr/bin/env bash
#
# 20_benchmark_threads.sh - FastCrossMap 多线程扩展性 benchmark
#
# 用法 (在仓库根目录):
#   bash paper/20_benchmark_threads.sh
#
# 环境变量:
#   FCM       fast-crossmap 可执行文件   (默认: ./target/release/fast-crossmap)
#   DATA      输入数据目录                (默认: paper/data)
#   CACHE     解压缓存目录                (默认: $DATA/decompressed)
#   OUT       输出目录                    (默认: paper/results/threads-<时间戳>)
#   THREADS   要测的线程数, 空格分隔       (默认: "1 2 4 8 16")
#   REPS      每个 (format,threads) 重复次数 (默认: 5)
#   FORMATS   只跑其中几个, 空格分隔       (默认: "bed vcf gff gvcf bam sam" = 全部)
#   WARM      热身次数 (不计时)            (默认: 1)
#   SYNC      1 = 每次计时前后都 sync      (默认: 1, 见下, 强烈建议保持)
#   KEEP      1 = 保留输出文件             (默认: 0, 跑完即删, 省磁盘)
#
# 输出:
#   $OUT/timings.tsv   原始耗时  (format, threads, rep, wall_s, user_s, sys_s, max_rss_kb, flush_s, ok)
#   $OUT/summary.tsv   汇总表    (median/min/max/spread + 加速比 + 并行效率 + 可靠性标记)
#   $OUT/sha256.txt    各线程输出哈希, 用于一致性抽查
#   $OUT/env.txt       运行环境 (nproc / 版本 / 内核 / 挂载 / dirty_ratio)
#
# ============================ 测量方法 ============================
#
# 目标只有一个: 让"同一个格式不同 -t"之间的对比是可信的。为此做了四件事。
#
#   1. 轮转采样 (最重要)。
#      所有 (format, threads) 组合被拆成 `REPS` 轮, 每轮把每个格式、每个线程数各跑一次;
#      并且**每轮把顺序轮转一位** (rotation)。
#
#      为什么不是"一个 t 连跑 5 次": 机器状态(其它进程、温度、写回)在几分钟内会漂移,
#      如果 t=1 的五次恰好都落在慢窗口, 后面每个 t 的加速比都被系统性放大。
#      为什么还要轮转: 如果顺序永远是 t=1,2,4,8,16, 而机器上存在**周期性**干扰, 那么
#      每一轮的第 1 个槽位都会撞上干扰 —— t=1 被稳定地惩罚, 这比随机漂移更坏, 因为它
#      不会在多次重复中平均掉。上一版 VCF 就出现了这个特征:
#         t=1: 56 / 74 / 22 / 265 / 266      t=4: 79 / 13 / 285 / 286 / 13
#         t=2: 13 / 13 / 221 / 13 / 13       t=8: 216 / 12 / 12 / 12 / 12
#      同一轮里 t=2 快而 t=4 慢、下一轮又反过来 —— 这不是线程数的效应, 是槽位效应。
#      轮转 + 跨轮平均把它摊平。
#
#   2. 每次计时前后都 sync (SYNC=1)。
#      未压缩的 VCF/GFF 输出有 1-11 GB。进程退出时大量脏页仍在异步回写, 会挤进**下一次**
#      计时, 表现为 wall 很长而 user+sys 很短。真正致命的是内核的 dirty page 限流:
#      当脏页超过 dirty_ratio, 正在写的进程会被按在 balance_dirty_pages() 里等,
#      wall 飙高、sys 却很低。上一版 t=8 的一次 216s 里 user 63s + sys 1.5s, 剩下 150s
#      就是这么等掉的 —— 它不是 CPU 或磁盘带宽的瓶颈, 是限流。
#      sync 放在计时之外, 让每次测量从同样的"脏页已清"状态起步。
#
#   3. 中位数 + 离散度, 而不是最小值。
#      min 只反映"最顺的一次", 对抖动大的格式给出乐观且不可复现的数字。汇总表同时给
#      median/min/max/spread_pct, 并在 spread 过大时打上标记 —— 那样的行不能用。
#
#   4. 记录 user/sys/maxRSS。
#      wall 远大于 user+sys → 卡在 I/O 或限流上, 加线程不会有用, 这一点必须在表里看得见。
#      另外给出并行效率 eff = 加速比 / 线程数, 一眼看出是否已经跑到收益递减。
#
#   5. 单独量出 flush_s (计时后那次 sync 的耗时)。
#      它不是被测对象, 而是一把尺子, 用来判定本次 wall 到底花在哪:
#        flush_s 大        → 本次 wall 主要就是"把输出真正写进磁盘", 加线程无用;
#        flush_s 小而 wall 大 → 卡在别处 (限流/回收/CPU 争抢), 与磁盘无关;
#        flush_s 忽大忽小    → 采样点落在不同的回写状态下, 该行的 wall 不可比。
#
#      上一版 VCF 的抖动就是这么解释的 (t=1: 22/56/74/265/267s, t=8: 12/12/12/216s):
#      VCF 输出 11 GB, 265s ≈ 11 GB / 265s = 42 MB/s —— 正好是这块 SATA HDD 的顺序写速;
#      而 12.3s 意味着 11 GB / 12.3s = 915 MB/s, 对 SATA 盘不可能, 说明那次输出还在
#      页缓存里没落盘。两种"快慢"的 user+sys 都 ≈ 64.5s, CPU 开销一模一样 ——
#      变的只是回写状态。所以 t=1 的 74s 中位数是被回写污染的, 不是真实基线。
#
# 一次跑完所有格式不会互相干扰, 因为:
#   - 每个测量点计时前后都 sync, 上一个格式的脏页不会渗进来;
#   - 每个格式跑完就删掉自己的产物 (KEEP=0), 所以**峰值磁盘 = 单个最大产物**,
#     不是所有格式之和; 各格式的输出文件从不同时存在;
#   - 加速比的分母是**同一格式自己的** t=1, 格式之间不共用基准。
#   唯一要单独跑的理由是机器中途要挪作他用时便于中断; 否则一次跑完更好 ——
#   轮转采样本来就依赖"把时间摊到所有格式上"。
#
# 关于 -t 的公平性: 所有格式都按 -t 1 出正式 benchmark 对比 CrossMap (CrossMap 无多线程)。
#   本脚本只是刻画 FCM 自身的扩展性曲线, 不改变 benchmark 的默认假设。
#
# ============================ 读表 ============================
#   spread_pct 很大 (比如 > 50%)  → 该行被干扰污染, 加 REPS 重跑, 不要引用
#   wall 远大于 cpu_s            → I/O / 限流受限, 加线程无用
#   flush_s 与 wall 同量级        → 该次 wall 主要花在落盘, 再看 timings.tsv 里同组的 flush_s
#   eff 在某个 t 之后掉下来      → 收益递减, 那个 t 就是推荐值
#   BAM 高 t 变慢是正常的: BGZF 解压/压缩池与转换池互相抢核, 且写回 BGZF 是串行装配。
#
# ============================ 数据说明 ============================
#   只有 VCF/GVCF 支持 .gz 输入; BED/GFF 需要明文, 脚本会自动解压到 $CACHE (只解一次)。
#   **输出目录请放本地盘**: 未压缩的 GFF/VCF 输出有 1-11 GB, 写到 NFS/9p 挂载上会让
#   wall 时间完全被网络 I/O 支配, 线程数怎么调都看不出来。
#   本地盘至少留 20 GB (脚本会预检并告警)。

set -u

FCM="${FCM:-./target/release/fast-crossmap}"
DATA="${DATA:-paper/data}"
CACHE="${CACHE:-$DATA/decompressed}"
OUT="${OUT:-paper/results/threads-$(date +%Y%m%d-%H%M%S)}"
THREADS="${THREADS:-1 2 4 8 16}"
REPS="${REPS:-5}"
FORMATS="${FORMATS:-bed vcf gff gvcf bam sam}"
WARM="${WARM:-1}"
SYNC="${SYNC:-1}"
KEEP="${KEEP:-0}"

CHAIN="$DATA/hg19ToHg38.over.chain.gz"
REF="$DATA/hg38.fa"

# ---- 输入文件 (与 paper/data 对应; 缺哪个就跳过哪个) ------------------------
BED_IN="$DATA/bed_ccre_grch38.bed"                       # 129 MB 明文
BED_GZ="$DATA/bed_k562_dnase.bed.gz"                     # 备选 (小, 14 MB)
VCF_IN="$DATA/1000g_chr22.vcf.gz"                        # 206 MB, gz 可直接读
GFF_IN="$DATA/gff_gencode_v44_basic_grch37.gff3.gz"      # 46 MB, 需解压
GVCF_IN="${GVCF_IN:-$DATA/gvcf_sample.gvcf}"             # paper/data 无, 默认跳过
BAM_IN="$DATA/bam_na12878_chr20.bam"                     # 312 MB
SAM_IN="$DATA/sam_na12878_chr20.sam"                     # 1.45 GB

mkdir -p "$OUT" "$CACHE"

[ -x "$FCM" ]     || { echo "ERROR: 找不到 $FCM (先 cargo build --release)" >&2; exit 1; }
[ -f "$CHAIN" ]   || { echo "ERROR: 找不到 chain 文件 $CHAIN" >&2; exit 1; }

# 磁盘预检: 峰值 = 单个最大产物 (各格式产物跑完即删, 从不同时存在), VCF/GFF 可达 11 GB。
FREE_KB=$(df -Pk "$OUT" 2>/dev/null | awk 'NR==2{print $4}')
if [ -n "${FREE_KB:-}" ] && [ "$FREE_KB" -lt $((20 * 1024 * 1024)) ]; then
    echo "WARNING: $OUT 所在盘剩余 $((FREE_KB / 1024 / 1024)) GiB; VCF/GFF 输出可达 11 GB, 可能不够" >&2
fi

if [ -x /usr/bin/time ]; then HAVE_GNU_TIME=1; else HAVE_GNU_TIME=0; fi

# plain <file>: 若格式不支持 gz, 就解压到 $CACHE 并返回明文路径
plain() {
    local f="$1"
    case "$f" in
        *.gz) ;;
        *) echo "$f"; return ;;
    esac
    local base; base="$(basename "$f" .gz)"
    local out="$CACHE/$base"
    if [ ! -s "$out" ]; then
        echo "  [decompress] $f -> $out" >&2
        gzip -dc "$f" > "$out" || { echo "  !! 解压失败" >&2; return 1; }
    fi
    echo "$out"
}

# ============================ 建立任务表 ============================
# 每个格式一条: 名称 + 传给 FCM 的参数 (最后一项是输出文件)。
# 参数用 0x1f 分隔存储, 避免路径里有空格时被拆错。
# 输入缺失的格式在这里就被丢掉, 后面的轮转循环只看到可跑的格式。
JOB=()                       # 格式名列表
declare -A ARGS_OF           # 格式名 -> 参数 (0x1f 分隔)
declare -A SKIP_REASON       # 格式名 -> 跳过原因

add_job() {  # add_job <name> <args...>
    local name="$1"; shift
    local joined="" a
    for a in "$@"; do
        [ -z "$joined" ] && joined="$a" || joined="$joined"$'\x1f'"$a"
    done
    JOB+=("$name")
    ARGS_OF["$name"]="$joined"
}

want() { case " $FORMATS " in *" $1 "*) return 0;; *) return 1;; esac; }

if want bed; then
    if   [ -f "$BED_IN" ]; then add_job bed bed "$CHAIN" "$BED_IN" "$OUT/bed.out.bed"
    elif [ -f "$BED_GZ" ]; then
        if p=$(plain "$BED_GZ"); then add_job bed bed "$CHAIN" "$p" "$OUT/bed.out.bed"; fi
    else SKIP_REASON[bed]="无输入"; fi
fi
if want vcf; then
    if [ -f "$VCF_IN" ]; then add_job vcf vcf "$CHAIN" "$VCF_IN" "$REF" "$OUT/vcf.out.vcf"
    else SKIP_REASON[vcf]="无输入"; fi
fi
if want gff; then
    if [ -f "$GFF_IN" ]; then
        if p=$(plain "$GFF_IN"); then add_job gff gff "$CHAIN" "$p" "$OUT/gff.out.gff3"; fi
    else SKIP_REASON[gff]="无输入"; fi
fi
if want gvcf; then
    if [ -f "$GVCF_IN" ]; then add_job gvcf gvcf "$CHAIN" "$GVCF_IN" "$REF" "$OUT/gvcf.out.gvcf"
    else SKIP_REASON[gvcf]="GVCF_IN 未设置 / 无输入"; fi
fi
if want bam; then
    if [ -f "$BAM_IN" ]; then add_job bam bam "$CHAIN" "$BAM_IN" "$OUT/bam.out.bam"
    else SKIP_REASON[bam]="无输入"; fi
fi
if want sam; then
    if [ -f "$SAM_IN" ]; then add_job sam bam "$CHAIN" "$SAM_IN" "$OUT/sam.out.sam"
    else SKIP_REASON[sam]="无输入"; fi
fi

if [ "${#JOB[@]}" -eq 0 ]; then
    echo "ERROR: 没有可跑的格式 (检查 FORMATS 和 $DATA 下的输入文件)" >&2
    exit 1
fi

# 取某格式的输出文件 = 参数表的最后一项
outfile_of() { local s="${ARGS_OF[$1]}"; echo "${s##*$'\x1f'}"; }

# 清掉某格式上一次的产物。不能简单地 `rm "$out"*`: unmap / bedGraph 等副产物的名字
# 是 output 的**扩展名被替换**得来的 (GFF 的 gff.out.gff3 -> gff.out.gff.unmap),
# 前缀通配符漏得掉。这里把每个格式的副产物显式列出来。
clean_outputs() {  # clean_outputs <name>
    local name="$1" out; out=$(outfile_of "$name")
    local base="${out%.*}"
    # 主产物
    rm -f "$out" "$out".bai
    # unmap: 名字是 output 的扩展名被替换 (gff.out.gff3 -> gff.out.gff.unmap);
    # 各格式的后缀不同, 已知的就显式列, 再用通配兜底。
    rm -f "$base".bed.unmap "$base".vcf.unmap "$base".gvcf.unmap \
          "$base".gff.unmap "$base".maf.unmap "$base".unmap
    # BigWig 的中间 bedGraph 与副产物
    rm -f "$base".bedGraph "$base".bw "$out".bedGraph "$out".bw
}

# 参数表 -> 调用 FCM。所有调用点都走这里, 保证参数不会被空格拆错。
# 用法: fcm_job <name> <t>   (stdout/stderr 由调用方决定怎么处理)
fcm_job() {
    local name="$1" t="$2"
    local IFS=$'\x1f'
    local args=(${ARGS_OF[$name]})
    unset IFS
    "$FCM" "${args[@]}" -t "$t"
}
run_job() { fcm_job "$1" "$2"; }

ENV_F="$OUT/env.txt"
{
    echo "date=$(date -Is)"
    echo "host=$(hostname)  nproc=$(nproc 2>/dev/null || echo '?')"
    echo "kernel=$(uname -sr)"
    echo "FCM=$FCM"
    "$FCM" --version 2>/dev/null | sed 's/^/version: /'
    echo "THREADS=$THREADS  REPS=$REPS  WARM=$WARM  SYNC=$SYNC"
    echo "gnu_time=$HAVE_GNU_TIME"
    echo "OUT=$OUT (mount: $(df -T "$OUT" 2>/dev/null | awk 'NR==2{print $1" "$2}'))"
    echo "DATA=$DATA (mount: $(df -T "$DATA" 2>/dev/null | awk 'NR==2{print $1" "$2}'))"
    # 脏页限流阈值: 大输出的 wall 抖动主要来自这里, 记下来便于事后解释。
    echo "dirty_ratio=$(cat /proc/sys/vm/dirty_ratio 2>/dev/null)  dirty_background_ratio=$(cat /proc/sys/vm/dirty_background_ratio 2>/dev/null)"
    echo "memtotal_kb=$(awk '/MemTotal/{print $2}' /proc/meminfo)"
    echo "loadavg_at_start=$(cat /proc/loadavg)"
    echo "jobs=${JOB[*]}"
} | tee "$ENV_F"

for k in "${!SKIP_REASON[@]}"; do
    echo "  skip $k: ${SKIP_REASON[$k]}"
done

TIMINGS="$OUT/timings.tsv"; : > "$TIMINGS"
printf "format\tthreads\trep\twall_s\tuser_s\tsys_s\tmax_rss_kb\tflush_s\tok\n" >> "$TIMINGS"
: > "$OUT/sha256.txt"

# measure_once <name> <t>  -> 回显 fmt<TAB>t<TAB>wall<TAB>user<TAB>sys<TAB>rss<TAB>ok
measure_once() {
    local name="$1" t="$2"
    local outfile; outfile=$(outfile_of "$name")
    clean_outputs "$name"
    local wall user sys rss rc
    if [ "$HAVE_GNU_TIME" = 1 ]; then
        # 用 GNU time 拿 user/sys/maxRSS。为了让 shell 函数能被 time 计时, 用一个
        # 干净的子 shell 调 FCM —— 参数仍从 ARGS_OF 取, 不受空格影响。
        local tf="$OUT/.time.$$"
        ARGS_JOINED="${ARGS_OF[$name]}" FCM_BIN="$FCM" T="$t" \
        /usr/bin/time -f '%e\t%U\t%S\t%M' -o "$tf" bash -c '
            IFS=$'"'"'\x1f'"'"'; args=($ARGS_JOINED); unset IFS
            exec "$FCM_BIN" "${args[@]}" -t "$T"
        ' >/dev/null 2>&1
        rc=$?
        IFS=$'\t' read -r wall user sys rss < "$tf" || { wall=0; user=0; sys=0; rss=0; }
        rm -f "$tf"
    else
        local s e
        s=$(date +%s.%N)
        run_job "$name" "$t" >/dev/null 2>&1; rc=$?
        e=$(date +%s.%N)
        wall=$(awk -v a="$s" -v b="$e" 'BEGIN{printf "%.3f", b-a}')
        user=0; sys=0; rss=0
    fi
    local ok=0; [ "$rc" -eq 0 ] && ok=1
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$name" "$t" "$wall" "$user" "$sys" "$rss" "$ok"
}

# 把本次留下的脏页刷进磁盘, 返回耗时 (秒)。
# 单独计时而不是混在 measure_once 里, 是为了让 user/sys 仍然只属于 FCM 进程。
flush_dirty() {
    local s e
    s=$(date +%s.%N)
    sync
    e=$(date +%s.%N)
    awk -v a="$s" -v b="$e" 'BEGIN{printf "%.3f", b-a}'
}

# ============================ 热身 ============================
# 把每个格式的输入读进页缓存, 建好输出目录项。不计时、不记录。
echo
for name in "${JOB[@]}"; do
    outfile=$(outfile_of "$name")
    for _ in $(seq 1 "$WARM"); do
        clean_outputs "$name"
        run_job "$name" "${THREADS%% *}" >/dev/null 2>&1 || true
    done
    printf "  warm %-6s ok\n" "$name"
done
sync

# ============================ 轮转采样 ============================
# 每轮: 所有格式 × 所有 t 各一次; 格式顺序和 t 顺序都按轮次轮转一位。
THREAD_ARR=($THREADS)
NT=${#THREAD_ARR[@]}
NF=${#JOB[@]}

echo
echo "=================== 采样 ($REPS 轮 × $NF 格式 × $NT 线程数) ==================="
for r in $(seq 1 "$REPS"); do
    for fi in $(seq 0 $((NF - 1))); do
        name="${JOB[$(((fi + r - 1) % NF))]}"
        for ti in $(seq 0 $((NT - 1))); do
            t="${THREAD_ARR[$(((ti + r - 1) % NT))]}"

            # 计时前 sync: 让本次测量从"脏页已清"起步 (见文件头的说明)。
            [ "$SYNC" = 1 ] && sync

            row=$(measure_once "$name" "$t")

            # 计时后立刻 sync, 并单独量出这次落盘花了多久。
            # 这一步的耗时不进 wall, 但它本身是一把尺子:
            #   flush 很大 → 本次的 wall 主要是把输出写进磁盘, 再加大线程也没用;
            #   flush 很小而 wall 很大 → 卡在别处 (内核回收/页限流/CPU 争抢)。
            # VCF/GFF 的输出有 11 GB, 各次 flush 之间的巨大差异正是 wall 忽大忽小的来源。
            if [ "$SYNC" = 1 ]; then flush=$(flush_dirty); else flush="-"; fi

            printf "%s\t%s\n" "$row" "$r" \
                | awk -F'\t' -v fl="$flush" '{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$8,$3,$4,$5,$6,fl,$7}' >> "$TIMINGS"
            printf "  r=%s  %-5s t=%-3s wall=%8ss  cpu=%5ss  rss=%6sKB  flush=%8ss  ok=%s\n" \
                "$r" "$name" "$t" \
                "$(printf '%s' "$row" | cut -f3)" \
                "$(awk -F'\t' '{printf "%.1f", $4+$5}' <<<"$row")" \
                "$(printf '%s' "$row" | cut -f6)" \
                "$flush" \
                "$(printf '%s' "$row" | cut -f7)"
        done
    done
done

# ============================ 一致性抽查 ============================
# 每个 (格式, t) 各跑一次, 比对哈希 —— 多线程不能改变输出。
echo
for name in "${JOB[@]}"; do
    outfile=$(outfile_of "$name")
    for t in "${THREAD_ARR[@]}"; do
        clean_outputs "$name"
        run_job "$name" "$t" >/dev/null 2>&1 || true
        sha256sum "$outfile" 2>/dev/null \
            | awk -v f="$name" -v t="$t" '{print f"\tt="t"\t"$1}' >> "$OUT/sha256.txt"
        [ "$KEEP" = "1" ] || clean_outputs "$name"
    done
done

# ============================ 汇总 ============================
# 每 (format,threads): 中位数为主, 附 min/max/spread; 加速比分母 = 同格式 t=1 的中位数;
# eff = 加速比 / 线程数; flag 标记不可信的行。
SUMMARY="$OUT/summary.tsv"
printf "format\tthreads\tmedian_s\tmin_s\tmax_s\tspread_pct\tspeedup_vs_t1\teff\tcpu_s\tn\tflag\n" > "$SUMMARY"
awk -F'\t' '
    function median(src, m,   i,j,t) {           # src 已排好序
        return (m%2) ? src[(m+1)/2] : (src[m/2]+src[m/2+1])/2
    }
    function sort(arr, m,   i,j,t) {
        for(i=1;i<=m;i++) for(j=i+1;j<=m;j++) if(arr[j]<arr[i]){t=arr[i];arr[i]=arr[j];arr[j]=t}
    }
    NR>1 && $9==1 {
        k=$1"\t"$2; n[k]++; v[k,n[k]]=($4+0); c[k,n[k]]=($5+0)+($6+0);
        if(!(k in lo)||$4+0<lo[k]) lo[k]=$4+0;
        if(!(k in hi)||$4+0>hi[k]) hi[k]=$4+0;
    }
    END{
        for(k in n){
            m=n[k];
            # 收集本组的 wall 与 cpu 到临时数组后排序 (awk 没有子数组, 只能按 k 重建)
            delete w; delete cc;
            for(i=1;i<=m;i++){ w[i]=v[k,i]; cc[i]=c[k,i] }
            sort(w,m); sort(cc,m);
            med=median(w,m); cmed=median(cc,m);
            split(k,p,"\t"); f=p[1];
            medf[k]=med;
            spread=(med>0)?100*(hi[k]-lo[k])/med:0;
            flag="";
            if(spread>50) flag="UNRELIABLE(spread)";
            else if(cmed>0 && med>2*cmed) flag="IO/throttle-bound";
            row[k]=sprintf("%s\t%s\t%.3f\t%.3f\t%.3f\t%.1f\t%s\t%s\t%.1f\t%d\t%s",
                f, p[2], med, lo[k], hi[k], spread, "-", "-", cmed, m, flag);
        }
        # t=1 的中位数作为加速比分母
        for(k in medf){ split(k,p,"\t"); if(p[2]=="1") base[p[1]]=medf[k]; }
        for(k in row){
            split(k,p,"\t"); f=p[1]; t=p[2]+0;
            sp=(base[f]>0 && medf[k]>0)? base[f]/medf[k] : 0;
            eff=(t>0)? sp/t : 0;
            n2=split(row[k],q,"\t"); q[7]=sprintf("%.2f", sp); q[8]=sprintf("%.2f", eff);
            printf "%s", q[1];
            for(i=2;i<=n2;i++) printf "\t%s", q[i];
            printf "\n";
        }
    }' "$TIMINGS" | sort -k1,1 -k2,2n >> "$SUMMARY"

echo; echo "=================== 汇总 ==================="; column -t "$SUMMARY"
echo
echo "逐次原始数据: $TIMINGS"
echo "哈希抽查:     $OUT/sha256.txt  (同一 format 下不同 t 应完全相同)"
echo "环境:         $ENV_F"
echo "有差异时:     bash paper/21_check_determinism.sh"
echo
echo "读表提示:"
echo "  flag=UNRELIABLE  → 该行抖动过大, 结论不可用, 用 FORMATS=<该格式> 单独加 REPS 重跑"
echo "  flag=IO/throttle → wall 远大于 cpu, 受 I/O 或脏页限流限制, 加线程无用"
echo "  eff 掉到 0.5 以下 → 已过收益递减点, 那之前的 t 就是推荐值"
echo "  只有 t=1..4 有稳定加速而 8/16 变慢, 是正常的: 转换池、BGZF 池和写回在抢同一批核"
