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
#   FORMATS 只校验其中几个, 空格分隔     (默认: "bed vcf gff gvcf bam sam wig bigwig maf")
#   MAF_BUILD  MAF 的 -b 目标 build      (默认: GRCh38)
#
# 覆盖的产物按格式不同 (artifacts_of 里显式列出):
#   wig    <前缀>.wig        + <前缀>.unmap.wig
#   bigwig <前缀>.bw         + <前缀>.unmap.bedGraph
#   maf    <输出>            + <输出扩展名换成 maf.unmap>
#   其余    <输出>           + <输出扩展名换成 <fmt>.unmap>
# 中间文件 (spool / bigtools 的 temp.bedGraph) 不参与比对。
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
WIG_IN="$DATA/wig_gm12878_h3k27ac.wig"
BIGWIG_IN="$DATA/bigwig_gm12878_h3k27ac.bigWig"
MAF_IN="$DATA/maf_synthetic_large.maf"

FORMATS="${FORMATS:-bed vcf gff gvcf bam sam wig bigwig maf}"
MAF_BUILD="${MAF_BUILD:-GRCh38}"
want() { case " $FORMATS " in *" $1 "*) return 0;; *) return 1;; esac; }

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

# 某格式在给定前缀下的全部产物, 一行一个。必须和 20_benchmark_threads.sh 里的
# artifacts_of 保持一致 —— 两边的产物命名规则来自各格式的实现:
#   wig/bigwig 是「前缀式」输出 (<前缀>.wig / <前缀>.bw), unmap 是 <前缀>.unmap.*;
#   其余格式的 unmap 是 output 的**扩展名被替换**得来的 (gff.t1 -> gff.gff.unmap),
#   所以要用 ${o%.*} 而不是 "$o"。
artifacts_of() {  # artifacts_of <fmt> <前缀>
    local fmt="$1" o="$2" base="${o%.*}"
    case "$fmt" in
        wig)    printf '%s\n' "$o.wig" "$o.unmap.wig" ;;
        bigwig) printf '%s\n' "$o.bw"  "$o.unmap.bedGraph" ;;
        maf)    printf '%s\n' "$o" "$base.maf.unmap" ;;
        *)      printf '%s\n' "$o" "$base.$fmt.unmap" ;;
    esac
}

# check <format> <subcmd> <input> [ref] [extra args...]
#   extra 直接跟在最后 (例如 maf 的 -b GRCh38); clap 允许选项出现在位置参数之后。
check() {
    local fmt="$1" sub="$2" input="$3" ref="${4:-}"
    local extra=("${@:5}")   # 可能为空 (仅当 $# == 4)
    local base="$OUT/$fmt.t1"
    for t in $THREADS; do
        local o="$OUT/$fmt.t$t"
        local f; while IFS= read -r f; do rm -f "$f"; done < <(artifacts_of "$fmt" "$o")
        local cmd=("$FCM" "$sub" "$CHAIN" "$input")
        [ -n "$ref" ] && cmd+=("$ref")
        cmd+=("$o" -t "$t")
        # 空数组在 set -u 下展开成未定义, 所以显式判空再加 (bash < 4.4 会报错)。
        [ "${#extra[@]}" -gt 0 ] && cmd+=("${extra[@]}")
        if ! "${cmd[@]}" >/dev/null 2>&1; then
            echo "  [$fmt] t=$t  !! 运行失败"; FAIL=1; continue
        fi
        if [ "$t" = "1" ]; then echo "  [$fmt] t=1  (基线)"; continue; fi
        local a1 a2 f1 f2 bad=0 i
        mapfile -t a1 < <(artifacts_of "$fmt" "$base")
        mapfile -t a2 < <(artifacts_of "$fmt" "$o")
        for i in "${!a1[@]}"; do
            f1="${a1[$i]}"; f2="${a2[$i]}"
            if [ -f "$f1" ] || [ -f "$f2" ]; then
                if ! cmp -s "$f1" "$f2"; then
                    echo "  [$fmt] t=$t  DIFF: $(basename "$f1")  vs  $(basename "$f2")"
                    diff "$f1" "$f2" 2>/dev/null | head -4 | sed 's/^/        /'
                    bad=1
                fi
            fi
        done
        [ "$bad" = 0 ] && echo "  [$fmt] t=$t  OK (byte-identical)"
    done
    [ "${KEEP:-0}" = "1" ] || {
        local f; while IFS= read -r f; do rm -f "$f"; done \
            < <(for t in $THREADS; do artifacts_of "$fmt" "$OUT/$fmt.t$t"; done)
    }
}

echo "=== determinism check (out=$OUT) ==="

if want bed; then
    if   [ -f "$BED_IN" ]; then check bed bed "$BED_IN"
    elif [ -f "$BED_GZ" ]; then p=$(plain "$BED_GZ") && check bed bed "$p"
    else echo "  [bed] skip"; fi
fi

if want vcf; then
    [ -f "$VCF_IN" ]  && check vcf  vcf  "$VCF_IN"  "$REF" || echo "  [vcf] skip"
fi
if want gff; then
    if [ -f "$GFF_IN" ]; then p=$(plain "$GFF_IN") && check gff gff "$p"; else echo "  [gff] skip"; fi
fi
if want gvcf; then
    [ -f "$GVCF_IN" ] && check gvcf gvcf "$GVCF_IN" "$REF" || echo "  [gvcf] skip"
fi
if want bam; then
    [ -f "$BAM_IN" ]  && check bam  bam  "$BAM_IN"         || echo "  [bam] skip"
fi
if want sam; then
    [ -f "$SAM_IN" ]  && check sam  bam  "$SAM_IN"         || echo "  [sam] skip"
fi

# --- 后三种格式 (wig / bigwig / maf) ---
# 输出是前缀 (wig/bigwig) 或需要参考基因组+目标 build (maf), 见 artifacts_of。
if want wig; then
    [ -f "$WIG_IN" ] && check wig wig "$WIG_IN" || echo "  [wig] skip"
fi
if want bigwig; then
    [ -f "$BIGWIG_IN" ] && check bigwig bigwig "$BIGWIG_IN" || echo "  [bigwig] skip"
fi
if want maf; then
    if   [ ! -f "$MAF_IN" ]; then echo "  [maf] skip (无输入 $MAF_IN)"
    elif [ ! -f "$REF" ];    then echo "  [maf] skip (缺参考基因组 $REF)"
    else check maf maf "$MAF_IN" "$REF" -b "$MAF_BUILD"
    fi
fi

echo
[ "$FAIL" = 0 ] && echo "结果: 全部一致 ✅" || echo "结果: 存在差异 ❌ (见上方)"
exit "$FAIL"
