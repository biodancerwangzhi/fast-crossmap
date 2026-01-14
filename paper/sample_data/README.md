# Sample Data

Pre-generated test files for quick verification.

| File | Records | Description |
|------|---------|-------------|
| `sample.bed` | 10,000 | Synthetic BED records (hg19) |
| `sample.sam` | 1,000 | Synthetic SAM reads (hg19) |
| `hg19ToHg38.over.chain.gz` | - | UCSC chain file |

## Quick Test

```bash
# From repository root:
./fast-crossmap bed \
    paper/sample_data/hg19ToHg38.over.chain.gz \
    paper/sample_data/sample.bed \
    output.bed
```

See `paper/README.md` for full instructions.
