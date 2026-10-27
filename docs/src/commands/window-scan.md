# madcats window-scan

Sweep candidate tag lengths on a sample subset and recommend one from the project's own reads.

## Synopsis

```bash
madcats window-scan --probe-design <FILE> --fastq-dir <DIR> --lengths <LIST> --out-dir <DIR> [OPTIONS]
```

## Description

Window length trades read retention against information. Longer windows carry more linked
variation, but fewer reads cover them fully and more reads carry errors. `window-scan` measures
three curves without external truth:

| Curve | Measures |
|-------|----------|
| **Retention** | Fraction of located reads that fully cover the best window. Predicted for a dense grid of lengths from one coverage pass |
| **Cleanliness** | Assigned reads not merged away by denoising or dropped for failing the evidence rules |
| **Information** | Evidence-passing novel haplotype rows per 100 markers |

Each candidate in `--lengths` runs the actual `normalize-reads` and `generate` code on a sample
subset.

### How the length is chosen

*Survival* = retention × assignment × cleanliness. Candidates within `--info-band` (default 10%) of
the best survival are eligible. Among them, the most informative length wins, with ties going to the
shorter length.

Information is not used on its own, because systematic errors (such as Ion homopolymer indels that
recur across samples) pass the reproducibility route and look like extra haplotypes. With this
rule, the scan picked the length that independent validation favored on every dataset tested:
41 bp on Ion oyster, 81 bp on GT-seq, and 61 bp on capture.

### Seed sensitivity

The recommendation can depend on which samples are scanned. `--scan-draws N` repeats the sweep on N
seeded draws, reports each draw's recommendation, and takes the most-voted length. If the
recommendation changes between draws, scan more samples or check the choice against external data.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--probe-design <FILE>` | Probe/marker design TSV or CSV |
| `--fastq-dir <DIR>` | Folder containing one FASTQ(.gz) per sample |
| `--lengths <LIST>` | Comma-separated candidate lengths to run in full, e.g. `41,61,81,110` |
| `--out-dir <DIR>` | Folder for per-length runs and the scan report |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--fastq-glob <GLOB>` | `*` | FASTQ filename glob within `--fastq-dir` |
| `--scan-samples <INT>` | 6 | Samples per candidate run: the first in sorted order, or a seeded draw |
| `--scan-seed <INT>` | none | Draw the subset at random with this seed |
| `--scan-draws <INT>` | 1 | Independent draws to sweep (seeds `--scan-seed`, +1, +2, ...) |
| `--grid-step <INT>` | 10 | Step for the dense predicted-retention grid |
| `--info-band <FLOAT>` | 0.10 | Survival tolerance for the information tiebreak |
| `--generate-args <STR>` | none | Extra `generate` flags for every candidate run, quoted as one string |
| `--group-routing <MODE>` | `nearest` | Passed to candidate `normalize-reads` runs (`covered` for merged tiled input) |
| `--anchor-length`, `--kmer-size`, `--str-penalty` | as `normalize-reads` | Passed to candidate `normalize-reads` runs |
| `--jobs <INT>` | up to 4 | Samples processed in parallel |

## Examples

### Basic scan

```bash
madcats window-scan \
  --probe-design markers.tsv \
  --fastq-dir raw/ \
  --lengths 31,41,61,81 \
  --out-dir scan/ \
  --jobs 8
```

### Check stability across three random subsets

```bash
madcats window-scan \
  --probe-design markers.tsv \
  --fastq-dir raw/ \
  --lengths 41,61,81,110 \
  --scan-seed 1 --scan-draws 3 \
  --generate-args "--min-supporting-freq 0.10" \
  --out-dir scan/
```

## Output

| File | Contents |
|------|----------|
| `OUT_DIR/window_scan.tsv` | One row per length: predicted and measured retention, assignment, cleanliness, information, survival, eligibility; then the dense retention grid |
| `OUT_DIR/window_scan.qc.json` | Parameters, per-length results and `recommended_tag_length`. With several draws, also `recommendation_votes` and `recommendation_stability` |
| `OUT_DIR/len_<length>/` | Per-candidate run directories (under `OUT_DIR/draw_N/` when `--scan-draws` > 1) |

If tiles cannot be given distinct windows at some length, that length is reported as failed and the
scan continues.

## Tips

1. Treat the result as a starting point. If you want more haplotype information and can check the
   result against external data or family structure, widen `--info-band`.
2. Pass your production evidence options with `--generate-args` so the scan measures the
   configuration you will run.

## See Also

- [madcats normalize-reads](normalize-reads.md)
