# madcats merge-pairs

Stitch overlapping paired-end reads into single fragment reads.

## Synopsis

```bash
madcats merge-pairs --fastq-dir <DIR> --out-dir <DIR> [OPTIONS]
```

## Description

For assays whose inserts are shorter than R1 + R2, `merge-pairs` collapses each pair into one
fragment before `normalize-reads` or `generate`:

- R2 is reverse-complemented and overlaps are scanned longest first.
- The merged read is R1's head, plus the overlap, plus R2's tail. In extended geometry that is the
  whole fragment. In engulfed (read-through) geometry, R1's adapter tail past the fragment end lies
  outside the overlap and never reaches the merged read, so no adapter trimming is needed. R2
  bases past R1's end are kept when the reads differ in length.
- Where overlap bases disagree, the higher-quality call wins.
- Mates are paired by read ID. Mismatched names or unequal read counts stop the run with an error.

Samples are independent and merge in parallel with `--jobs`. Each sample is written to
`OUT_DIR/<sample>.fastq.gz`, where `<sample>` is the filename with the mate suffix removed.

> **When to use:** only when mates overlap. For dual-anchored capture designs whose mates cover
> opposite halves of the fragment, pool both mates instead (see
> [Paired-End Reads and Tiling](../tutorials/paired-end.md)). MADCATS does not decide this
> for you.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--fastq-dir <DIR>` | Folder holding paired FASTQ(.gz) files |
| `--out-dir <DIR>` | Folder for merged per-sample FASTQs |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--r1-suffix <STR>` | `_R1.fastq.gz` | Filename suffix of the R1 mate |
| `--r2-suffix <STR>` | `_R2.fastq.gz` | Filename suffix of the R2 mate |
| `--min-overlap <INT>` | 15 | Shortest overlap, in bases, that lets a pair merge |
| `--max-mismatch-rate <FLOAT>` | 0.1 | Mismatch fraction tolerated inside the overlap |
| `--unmerged <MODE>` | `r1` | Unmergeable pairs: `r1` keeps the R1, `pool` keeps both mates, `drop` discards the pair |
| `--jobs <INT>` | up to 4 | Samples processed in parallel |

## Examples

### Standard Illumina naming

```bash
madcats merge-pairs \
  --fastq-dir raw/ \
  --out-dir merged/ \
  --jobs 8
```

### Custom mate suffixes

```bash
madcats merge-pairs \
  --fastq-dir raw/ \
  --r1-suffix _1.fq.gz \
  --r2-suffix _2.fq.gz \
  --out-dir merged/
```

## Accuracy

- Synthetic truth set with modeled adapter read-through: 100% of 513,715 pairs merged, 98.6%
  molecular count recovery, zero length errors.
- Compared with FLASH 1.2.11, vsearch and fastp on a four-regime exact-truth set: equal to the best
  (vsearch) in every regime, with zero false merges.

## See Also

- [madcats normalize-reads](normalize-reads.md): the usual next step
- [Paired-End Reads and Tiling](../tutorials/paired-end.md)
