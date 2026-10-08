#!/usr/bin/env bash
#
# 21_check_determinism.sh - 校验 -t N 与 -t 1 的输出是否逐字节相同
#
# 用法 (在仓库根目录):
#   bash paper/21_check_determinism.sh
#
# 环境变量:
#   FCM     fast-crossmap 可执行文件   (默认: ./target/release/fast-crossmap)
#   DATA    输入数据目录                (默认: paper/data)
#   CACHE   解压缓存目录                (默认: paper/data/decompressed)
#   THREADS 要对比的线程数              (默认: "1 2 4 8")
#   OUT     中间文件目录                (默认: paper/results/determinism-<时间戳>)
#   KEEP    1 = 保留中间输出            (默认: 0)
#
# 退出码: 0 = 全部一致; 1 = 有差异

set -u

FCM="${FCM:-./target/release/fast-crossmap}"
DATA="${DATA:-paper/data}"
CACHE="${CACHE:-$DATA/decompressed}"
OUT="${OUT:-paper/results/determinism-$(date +%Y%m%d-%H%M%S)}"
THREADS="${THREADS:-1 2 4 8}"

CHAIN="$DATA/hg19ToHg38.over.chain.gz"
REF="$DATA/hg38.fa"

BED_IN="$DATA/bed_ccre_grch38.bed"
BED_GZ="$DATA/bed_k562_dnase.bed.gz"
VCF_IN="$DATA/1000g_chr22.vcf.gz"
GFF_IN="$DATA/gff_gencode_v44_basic_grch37.gff3.gz"
GVCF_IN="$DATA/gvcf_sample.gvcf"
BAM_IN="$DATA/bam_na12878_chr20.bam"
SAM_IN="$DATA/sam_na12878_chr20.sam"

mkdir -p "$OUT" "$CACHE"
FAIL=0

plain() {
    local f="$1"
    case "$f" in *.gz) ;; *) echo "$f"; return ;; esac
    local out="$CACHE/$(basename "$f" .gz)"
    if [ ! -s "$out" ]; then
        echo "  [decompress] $f -> $out" >&2
        gzip -dc "$f" > "$out" || { echo "  !! 解压失败" >&2; return 1; }
    fi
    echo "$out"
}

# check <format> <subcmd> <input> [ref]
check() {
    local fmt="$1" sub="$2" input="$3" ref="${4:-}"
    local base="$OUT/$fmt.t1"
    for t in $THREADS; do
        local o="$OUT/$fmt.t$t"; rm -f "$o" "$o".*
        local cmd=("$FCM" "$sub" "$CHAIN" "$input")
        [ -n "$ref" ] && cmd+=("$ref")
        cmd+=("$o" -t "$t")
        if ! "${cmd[@]}" >/dev/null 2>&1; then
            echo "  [$fmt] t=$t  !! 运行失败"; FAIL=1; continue
        fi
        if [ "$t" = "1" ]; then echo "  [$fmt] t=1  (基线)"; continue; fi
        local bad=0
        for suffix in "" ".unmap" ".gff.unmap" ".vcf.unmap" ".gvcf.unmap"; do
            if [ -f "$base$suffix" ] || [ -f "$o$suffix" ]; then
                if ! cmp -s "$base$suffix" "$o$suffix"; then
                    echo "  [$fmt] t=$t  DIFF: $fmt.t1$suffix  vs  $fmt.t$t$suffix"
                    diff "$base$suffix" "$o$suffix" 2>/dev/null | head -4 | sed 's/^/        /'
                    bad=1
                fi
            fi
        done
        [ "$bad" = 0 ] && echo "  [$fmt] t=$t  OK (byte-identical)"
    done
    [ "${KEEP:-0}" = "1" ] || rm -f "$OUT/$fmt".t* 2>/dev/null
}

echo "=== determinism check (out=$OUT) ==="

if   [ -f "$BED_IN" ]; then check bed bed "$BED_IN"
elif [ -f "$BED_GZ" ]; then p=$(plain "$BED_GZ") && check bed bed "$p"
else echo "  [bed] skip"; fi

[ -f "$VCF_IN" ]  && check vcf  vcf  "$VCF_IN"  "$REF" || echo "  [vcf] skip"
if [ -f "$GFF_IN" ]; then p=$(plain "$GFF_IN") && check gff gff "$p"; else echo "  [gff] skip"; fi
[ -f "$GVCF_IN" ] && check gvcf gvcf "$GVCF_IN" "$REF" || echo "  [gvcf] skip"
[ -f "$BAM_IN" ]  && check bam  bam  "$BAM_IN"         || echo "  [bam] skip"
[ -f "$SAM_IN" ]  && check sam  bam  "$SAM_IN"         || echo "  [sam] skip"

echo
[ "$FAIL" = 0 ] && echo "结果: 全部一致 ✅" || echo "结果: 存在差异 ❌ (见上方)"
exit "$FAIL"
