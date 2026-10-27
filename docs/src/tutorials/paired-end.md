# Paired-End Reads and Tiling

MADCATS does not detect paired-end data automatically. Check how the reads of your assay are laid
out and pick the matching route.

## Overlapping Pairs: Merge

When the insert is shorter than R1 + R2 (for example GT-seq PE100 over 130–200-bp amplicons, or
rhAmpSeq), merge first:

```bash
madcats merge-pairs \
  --fastq-dir raw/ \
  --out-dir merged/ \
  --jobs 8
```

The merger reverse-complements R2 and scans overlaps longest first. The merged read is always R1's
head, plus the overlap, plus R2's tail. In read-through geometry, R1's adapter tail past the
fragment end falls outside the overlap and never reaches the merged read, so no adapter trimming
is needed. Disagreeing overlap bases take the higher-quality call.

`--unmerged` decides what happens to pairs that do not merge:

| Value | Behaviour |
|-------|-----------|
| `r1` (default) | Keep the R1 |
| `pool` | Keep both mates as independent reads |
| `drop` | Discard the pair |

Mates are paired by read ID. Mismatched names or unequal read counts stop the run with an error.

On a synthetic paired-end truth set with modeled adapter read-through, 100% of 513,715 pairs merged
with 98.6% molecular count recovery.

## Non-Overlapping Pairs: Pool

Dual-anchored capture designs (Flex-Seq: R1 anchors the left probe, R2 the right) have mates
covering opposite halves of a fragment neither read spans. The mates cannot be merged, and pooling
them cannot count a window twice, because no window is covered by both mates. Concatenate both
mates into one FASTQ per sample and feed that to `normalize-reads`.

## Multi-Window Tiling

When a fragment is longer than one window, several markers can tile one physical locus. Add a
`LocusGroup` column to the probe design to declare them as tiles:

```
MarkerName	TargetSequence	LocusGroup
amp17_snp1	ACGT...[C/T]...ACGTACGT...	amp17
amp17_snp2	ACGT...C...ACG[A/G]TACGT...	amp17
```

Tiles must share an identical Ref expansion frame (validated at load), reads locate to the group,
and then they are routed to member windows with `--group-routing`:

| Mode | Behaviour | Use when |
|------|-----------|----------|
| `nearest` (default) | Each read feeds exactly one tile: the fully covered window whose center is nearest the read's center | Reads are short relative to tile spacing |
| `covered` | The read is excised at every member window it fully covers | Reads span several tiles, e.g. **merged fragments** |

Use `covered` with merged full-span fragments. They all share the amplicon's center, so with
`nearest` one tile gets almost every read: on the GT-seq benchmark one tile took effectively all
depth at 163 of 164 tiled loci. Each window still sees each fragment at most once in both modes.

Fresh placement forces tiles to excise distinct sequences, since twin windows would make every
read an exact tie. Nudged members are listed in the QC as `locus_groups.window_dedup_nudged`.

## Putting It Together

The recipe validated on GT-seq: merge, normalize against a frozen panel with tile rows, generate.

```bash
madcats merge-pairs --fastq-dir raw/ --out-dir merged/ --jobs 8

madcats normalize-reads \
  --probe-design markers_with_tiles.tsv \
  --fastq-dir merged/ \
  --out-dir normalized/ \
  --reuse-tags frozen_panel.csv \
  --group-routing covered \
  --tag-length 81 \
  --jobs 8

madcats generate \
  --fastq-dir normalized/ \
  --fastq-glob "*.fastq" \
  --reference-tags normalized/window_reference_tags.csv \
  --tag-length 81 \
  --output project_MADC.csv
```

Against an independent truth VCF, original markers reproduced the R1-only accuracy (98.80%) on the
same cells, while the added tiles contributed 54% more truth-validated genotype calls at 98.43%.

On the cranberry Flex-Seq capture benchmark, pooled mates with a three-tile panel and `covered`
routing gave 3.7x the callable genotypes of a single-window R1 run, at higher concordance.
