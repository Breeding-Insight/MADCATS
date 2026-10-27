# AgriSeq and Other Non-Primer-Started Reads

AgriSeq reads (Ion Torrent) break the primer-started model three ways:

- amplicons are sequenced from both ends (about a 50/50 orientation split)
- primer digestion, and primers placed outside the submitted targets, scatter read start positions
- Ion Torrent's dominant error is the homopolymer indel

`madcats normalize-reads` adapts such reads. It locates each read against the probe-design targets
by k-mer, orients it forward, and excises a fixed-length window containing the variant. The output
is a folder of pipeline-ready per-sample FASTQs plus a matching window reference-tags CSV. After
that, `generate` runs unchanged.

## First Run of a Panel

```bash
madcats normalize-reads \
  --probe-design markers.tsv \
  --fastq-dir raw_fastqs/ \
  --out-dir normalized/ \
  --jobs 8

madcats generate \
  --fastq-dir normalized/ \
  --fastq-glob "*.fastq" \
  --reference-tags normalized/window_reference_tags.csv \
  --tag-length 41 \
  --min-supporting-freq 0.10 \
  --output project_MADC.csv
```

- `--tag-length 41` matches the default window of `normalize-reads`. The length comes from testing
  on honeybee, where longer windows lost coverage and indel-free reads faster than they added
  context. Use [`window-scan`](../commands/window-scan.md) to choose a length for a new panel.
- `--min-supporting-freq 0.10` is the recommended setting for AgriSeq. At AgriSeq depths (hundreds
  of reads per sample per marker), even sequencing errors pass the absolute read thresholds. On
  honeybee, default evidence admitted 1,752 novel rows, about 92% of them error haplotypes; the
  frequency floor brought that down to 147 while improving caller concordance.

### Window placement

With `--window-placement auto` (the default), each marker's window is placed from observed read
coverage, keeping the variant out of the assignment anchor and steering off tandem repeats
(`--str-penalty`). Repeats matter: when a variant allele lengthens an adjacent repeat, reads
carrying that allele are lost more often in a window that spans the repeat, and the marker can
look monomorphic.

## Every Later Run: Freeze the Windows

```bash
madcats normalize-reads \
  --probe-design markers.tsv \
  --fastq-dir raw_fastqs_run2/ \
  --out-dir normalized_run2/ \
  --reuse-tags normalized/window_reference_tags.csv \
  --jobs 8
```

Windows are placed from each run's own reads. Re-learning placement per run moved 42 of 85
windows between two honeybee chips, which made 58% of the second chip's haplotypes falsely appear
new. Freezing with `--reuse-tags` raised cross-chip haplotype sharing from 42% to
88%.

With `--reuse-tags`, `normalize-reads` still scores the frozen windows against the new run's reads
and reports `freeze_health` in the QC: the read-weighted coverage of the frozen panel against what
fresh placement would choose today, and the markers losing the most. One line on stderr summarizes
what freezing costs this run.

## Recommended Production Workflow

1. Run the first chip with auto placement.
2. Compare across runs and patch the few windows that coverage forced onto bad positions.
3. Freeze the curated panel.
4. Run every subsequent chip with `--reuse-tags`.

The curated frozen honeybee panel scored 93.4% / 93.9% vendor-caller concordance on two chips and
98.0% concordance across 30 technical replicate pairs.

## Large Panels

Assignment anchors start at `--anchor-length` (default 12) and grow per locus only where window
prefixes collide. A 17,500-marker genome-wide stress panel ran in about 820 MB with 17,291 markers
at the minimum and 209 grown to 13–32 bp.

## Other Platforms

The same route handles any amplicon or capture platform whose reads contain a panel window but do
not start at a fixed primer. For paired-end data, decide the geometry first; see
[Paired-End Reads and Tiling](paired-end.md).
