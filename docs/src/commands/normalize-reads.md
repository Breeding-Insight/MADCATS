# madcats normalize-reads

Normalize bidirectional, variable-start amplicon reads into pipeline-ready pseudo-tag FASTQs plus a
window panel.

## Synopsis

```bash
madcats normalize-reads --probe-design <FILE> --fastq-dir <DIR> --out-dir <DIR> [OPTIONS]
```

## Description

`generate` expects every read to start at a fixed anchor. Platforms such as AgriSeq (Ion Torrent)
sequence amplicons from both ends with scattered start positions. `normalize-reads` adapts them:

1. Each read is located against the probe-design targets by k-mer (`--kmer-size`) and oriented
   forward.
2. A fixed-length window (`--tag-length`) containing the variant is excised.
3. Per-sample FASTQs of excised windows are written to `--out-dir`, together with the matching
   window reference-tags CSV and a QC report.

`generate` then runs unchanged on the output.

Windows are placed from the reads, so a new run can place them differently. Every later run of the
same panel should therefore pass `--reuse-tags` with the first run's window CSV. That keeps
microhaplotype sequences comparable across runs. See
[AgriSeq and Other Non-Primer-Started Reads](../tutorials/agriseq.md).

Markers that share a physical locus (`LocusGroup` in the probe design) are tiles. Reads are located
to the group as a whole, then routed among member windows by `--group-routing`, and fresh placement
forces member windows to excise distinct sequences.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--probe-design <FILE>` | Probe/marker design TSV or CSV with bracketed variants |
| `--fastq-dir <DIR>` | Folder containing one raw FASTQ(.gz) per sample |
| `--out-dir <DIR>` | Folder for the normalized per-sample FASTQs |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--fastq-glob <GLOB>` | `*` | FASTQ filename glob within `--fastq-dir` |
| `--tag-length <INT>` | 41 | Excised window length. Shorter windows lose fewer reads to indels on indel-prone chemistries |
| `--window-placement <MODE>` | `auto` | `auto` learns each marker's window from observed read coverage (an extra pass); `center` needs no read data. With `--reuse-tags`, `auto` scores the frozen windows instead and `center` skips that |
| `--anchor-length <INT>` | 12 | Minimum assignment-anchor prefix length. Anchors grow per locus only where window prefixes collide; the variant is always kept out of the anchor |
| `--reuse-tags <FILE>` | none | Window reference-tags CSV from a previous run to reuse verbatim, freezing placement |
| `--group-routing <MODE>` | `nearest` | For reads covering several tiles of one `LocusGroup`: `nearest` feeds one tile, `covered` feeds every tile the read fully covers |
| `--str-penalty <FLOAT>` | 0.002 | Auto placement: coverage penalty per tandem-repeat bp (unit sizes 1–3) beyond 4 in a window; `0` disables |
| `--kmer-size <INT>` | 21 | k-mer length for locating reads against the design targets |
| `--jobs <INT>` | up to 4 | Samples processed in parallel |
| `--output-tags <FILE>` | `OUT_DIR/window_reference_tags.csv` | Window reference-tags CSV |
| `--report <FILE>` | `OUT_DIR/normalize_reads.qc.json` | QC/provenance JSON |

## Examples

### First run of a panel

```bash
madcats normalize-reads \
  --probe-design markers.tsv \
  --fastq-dir raw_fastqs/ \
  --out-dir normalized/ \
  --jobs 8
```

### Later runs: freeze windows

```bash
madcats normalize-reads \
  --probe-design markers.tsv \
  --fastq-dir raw_fastqs_run2/ \
  --out-dir normalized_run2/ \
  --reuse-tags normalized/window_reference_tags.csv \
  --jobs 8
```

### Merged paired-end fragments with tiles

```bash
madcats normalize-reads \
  --probe-design markers_with_tiles.tsv \
  --fastq-dir merged/ \
  --out-dir normalized/ \
  --tag-length 110 \
  --group-routing covered
```

## Output

| File | Contents |
|------|----------|
| `OUT_DIR/<sample>.fastq` | Excised, forward-oriented windows per sample |
| `OUT_DIR/window_reference_tags.csv` | Window panel for `generate --reference-tags` |
| `OUT_DIR/normalize_reads.qc.json` | Read fates, per-marker windows, `locus_groups`, and `freeze_health` when reusing tags |

## Tips

1. Pass the same `--tag-length` to `generate`. It must match the excised window length.
2. Use `--fastq-glob "*.fastq"` in `generate` so only normalized files are read.
3. With `--reuse-tags`, every row of the frozen panel is checked against the design. Rows that only
   partly match, and panels that contain two identical windows, are rejected.
4. With `--reuse-tags`, the QC (`freeze_health`) compares each frozen window's coverage on this run
   with the window fresh placement would pick, and lists markers that lose more than ten points.

## See Also

- [madcats window-scan](window-scan.md): choose `--tag-length`
- [madcats merge-pairs](merge-pairs.md): paired-end intake
- [AgriSeq tutorial](../tutorials/agriseq.md)
