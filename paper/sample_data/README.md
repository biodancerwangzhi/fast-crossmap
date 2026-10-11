# Sample Data

Small test files for the quick test and for checking a build.

- `input/` -- one real-data slice per supported format, 100 records each
- `output/` -- the FastCrossMap output committed for each `input/` file, i.e.
  the reference the quick test compares against
- `hg19ToHg38.over.chain.gz` -- UCSC chain file every conversion uses

## Input / output pairing

| Input | Committed output |
|-------|------------------|
| `input/sample.bed` | `output/sample_hg38.bed` (+ `.unmap`) |
| `input/sample.sam` | `output/sample_hg38.sam` |
| `input/sample.vcf` | `output/sample_hg38.vcf` (+ `.unmap`) |
| `input/sample.gff3` | `output/sample_hg38.gff3` (+ `.gff.unmap`) |
| `input/sample.wig` | `output/sample_hg38.wig.wig` (+ `.unmap.wig`) |
| `input/sample.bigWig` | `output/sample_hg38.bigWig.bw` (+ `.unmap.bedGraph`) |
| `input/sample.maf` | `output/sample_hg38.maf` (+ `.unmap`) |
| `input/sample.bam` | (no committed output; BAM shares the `bam` subcommand with SAM) |

## Quick Test

```bash
# From repository root:
./fast-crossmap bed \
    paper/sample_data/hg19ToHg38.over.chain.gz \
    paper/sample_data/input/sample.bed \
    output.bed

diff output.bed paper/sample_data/output/sample_hg38.bed   # prints nothing if OK
```

A correct build reports 100 records, all mapped, and both `diff`s are silent.
See `paper/README.md` for the full instructions.
