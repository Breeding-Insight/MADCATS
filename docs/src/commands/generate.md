# madcats generate

Generate a MADC microhaplotype count table from per-sample FASTQs and a panel.

## Synopsis

```bash
madcats generate --fastq-dir <DIR> --reference-tags <FILE> --output <FILE> [OPTIONS]
madcats generate --fastq-dir <DIR> --primers <FILE> --probe-design <FILE> --output <FILE> [OPTIONS]
```

## Description

Every workflow ends with `generate`. It runs in seven stages, and the option tables below follow
them:

1. Load the panel from a reference-tag CSV (or build it from DArT design files) and find one FASTQ
   per sample.
2. Assign each read to a marker by its left-primer/anchor prefix, keep the fixed-length tag, and
   accept it if it reaches `--min-identity` to the marker's Ref or Alt tag. Reads that fail for
   every marker get one more try with a gapped alignment under the same edit limit, so the counts
   include undesigned indels. Credit a read that fits several tiled markers of one amplicon equally
   to each of them, and drop other ties. Reads shorter than the tag count toward a full-length
   allele whose start they match. Cache results per sample.
3. Count exact tag sequences per sample.
4. Apply the evidence rules. Always write canonical Ref and Alt rows, and report any other sequence
   only if it passes an evidence route.
5. Denoise. Within each sample, merge rare sequences one edit away from a much more common one
   into it. A pooled pass then merges sequences that failed the evidence rules into a parent that
   carries the same target allele.
6. Classify each row as Ref, Alt, RefMatch, AltMatch or Other, and annotate it.
7. Write the MADC, the QC/provenance JSON and the variants sidecar.

The [Evidence Model](../reference/evidence-model.md) explains stages 2–6 in detail.

> **Process a cohort in one invocation.** The reproducibility route counts supporting samples
> within the run, so splitting a cohort into batches weakens discovery. `--cache-dir` makes
> re-running one large invocation cheap.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--fastq-dir <DIR>` | Folder containing one FASTQ(.gz) per sample |
| `--output <FILE>` | MADC CSV to write |
| `--reference-tags <FILE>` | Reference-tag CSV, or |
| `--primers <FILE>` + `--probe-design <FILE>` | DArT design files |

## Panel and Samples (stage 1)

| Option | Default | Description |
|--------|---------|-------------|
| `--fastq-glob <GLOB>` | `*` | FASTQ filename glob within `--fastq-dir` |
| `--target-metadata <FILE>...` | `targets_*.csv` in the FASTQ folder | Target-manifest CSV(s) for sample metadata |
| `--sample-name-field <FIELD>` | targetid / FASTQ stem | Unique metadata field used for sample columns |
| `--sample-lookup <FILE>` | none | Two-column CSV/TSV mapping FASTQ-derived ids to sample names. Every run sample must have an entry |
| `--tag-length <INT>` | 81 | Tag length |

## Read Assignment (stage 2)

| Option | Default | Description |
|--------|---------|-------------|
| `--min-identity <FLOAT>` | 0.80 | Minimum substitution-only identity to a Ref or Alt tag. A gapped alignment uses the same edit limit |
| `--leading-anchor <SEQ\|auto>` | none | Constant sequence before the amplicon in raw deliveries (for DArT, `GATGGATCCAGTGC`). Each read is cut where the anchor ends (it must end within 40 bases beyond its own length). `auto` detects it from the first FASTQ |
| `--primer-mismatches {0,1}` | 1 | Mismatches allowed while assigning a primer/anchor. With `0`, reads from Match haplotypes whose distinguishing SNP falls inside the anchor are lost |
| `--gapped-rescue {on,off}` | `on` | Give reads that fail substitution-only identity one gapped alignment, so the counts include undesigned indels |
| `--max-indel-length <INT>` | 10 | Widest gap (alignment band) a gapped alignment may use |
| `--tile-groups {auto,declared,off}` | `auto` | Link tiled markers (markers that read the same tag window of one amplicon), so a read is credited to every tile it fits instead of being dropped as a tie. `auto` links tiles declared with `LocusGroup` and detects undeclared ones, `declared` links only the declared ones, and `off` links none. See [Tiled Markers](../getting-started/input-formats.md#tiled-markers) |
| `--partial-read-mode {on,off}` | `on` | Add a read shorter than the tag to the one full-length allele (already passing the evidence rules) whose start it matches, only when it is longer than every prefix shared between that locus's alleles. `off` drops such reads |
| `--min-partial-length <INT>` | tag length − 9 | Minimum length of a credited partial read |
| `--cache-dir <DIR>` | none | Reusable per-sample count cache |
| `--jobs <INT>` | up to 4 | Samples processed in parallel (also parallelizes denoising) |

## Evidence (stage 4)

| Option | Default | Description |
|--------|---------|-------------|
| `--min-global-count <INT>` | 1 | Minimum project-wide reads for a de novo allele |
| `--reproducible-min-reads <INT>` | 10 | Reads required per supporting sample in the reproducibility route |
| `--reproducible-min-samples <INT>` | 2 | Supporting samples required by the reproducibility route. `0` disables the route |
| `--depth-min-reads <INT>` | 30 | Reads within one sample that alone establish an allele through the depth route. `0` disables |
| `--quality-min-clean-reads <INT>` | 0 | Project-wide clean reads that establish an allele through the quality route. `0` (default) disables |
| `--discovery-min-quality {10,20,30,40}` | 20 | Phred score every tag base must reach for a read to count as clean |
| `--min-supporting-freq <FLOAT>` | 0 | Minimum share of a sample's locus reads before that sample can support the reproducibility or depth route |

## Denoising (stage 5)

| Option | Default | Description |
|--------|---------|-------------|
| `--denoise-distance {0,1,2}` | 1 | Most differing bases allowed between a rare sequence and its parent in per-sample denoising. `0` disables |
| `--global-denoise-distance {0,1,2}` | 1 | Most differing bases allowed in the pooled same-allele merge. `0` disables |
| `--denoise-alpha <FLOAT>` | 2.0 | UNOISE-style abundance-skew exponent |

## Output (stages 6–7)

| Option | Default | Description |
|--------|---------|-------------|
| `--qc-output <FILE>` | `OUTPUT.qc.json` | QC/provenance JSON |
| `--annotations {on,off}` | `on` | Repurpose `PICSnp` as `Variants` and `AvgPIC` as `MarkerType`, and write `OUTPUT.variants.csv`. `off` keeps the exact vendor header |
| `--remove-primers {on,off}` | `off` | Write every row without its primers, and add together the rows that become identical. Needs `--annotations on`. See [Removing primers](#removing-primers) |
| `--number-match-ids` | off | Append stable per-file serials to Match IDs |
| `--dart-preamble` | off | Prepend the vendor-style seven-row starred preamble that HapApp expects |
| `--preamble-project <STR>` | output name before its first dot | Label for the first preamble row |

## Removing Primers

A row's `AlleleSequence` is the tag as it was read. It begins with the left primer. Where the
amplicon is shorter than the tag plus the right primer, it also ends with the first bases of the
right primer. Those bases copy the primer, not the sample. `--remove-primers on` writes the rows
without them:

- **Left primer.** Every row loses as many leading bases as the panel's `LeftPrimerSequence` is
  long. The read pass accepts no insertion or deletion in the left primer, so this cut is exact.
- **Right primer.** `generate` finds where the right primer starts from the panel alone. On a
  platform with one amplicon length, such as DArTag, the tags that run into their right primer all
  agree on that length. An insertion or deletion in a row moves the primer, so each row is cut
  where its own alignment to the Ref tag reaches the primer. Rows therefore differ in length.
- **Merged rows.** Rows of a marker that become the same sequence are written as one row with
  their reads added together. Most of them differed only inside a primer.
- **Annotations.** `ClusterConsensusSequence`, `Variants`, `MarkerType` and the variants sidecar
  describe the trimmed rows, and positions count from the first base after the left primer.

The option acts only when the table is written. Reads are assigned, counted, denoised and
classified on the full tags, as with `off`, and a row keeps the category it was given there. No
read is gained or lost.

What to know before using it:

- **The panel must list the real primers.** The cut uses `LeftPrimerSequence` and
  `RightPrimerSequence`. In a panel from `infer-design` or `normalize-reads`, or one built from a
  vendor MADC, `LeftPrimerSequence` is an anchor: a prefix of the tag chosen to assign reads. It
  can be shorter or longer than the primer, or hold sample sequence where reads do not start with
  a primer. With such a panel the option removes the anchor, so use it only where the panel comes
  from the primer file (`design`, or `--primers` with `--probe-design`).
- **The right primer is removed only where it is found.** The QC JSON gives the amplicon length
  and counts the markers by outcome (`primer_removal`). A panel without one amplicon length, or
  whose tags stop before the right primer, has no right primer to remove. Its rows lose the right
  primer only when a short product brings the whole primer, or 10 or more of its first bases, into
  the row.
- **The right cut is an estimate, not a guarantee.** Against the full-length reads of a pecan and a
  sweetpotato DArTag cohort, 99.4% and 99.2% of reads were in rows cut at exactly the right base.
  The rest were a few bases off, with primer bases left or sample bases removed. Nearly all of
  those rows show five or fewer bases of the primer, or a primer copy that differs from the
  primer, or a repeat in which two cuts fit equally well, so the tag itself cannot settle where
  the primer starts. An insertion or deletion within about 5 bases of the right primer is the
  usual cause, because in a fixed-length tag it looks like substitutions. A cut at the designed
  position in a row whose last 8 bases and primer bases are as designed was wrong in 0.2% to 0.3%
  of rows, and any other cut in 3% to 5%.
- **A row keeps its right primer** when the cut would remove the target base or leave fewer than
  10 bases. The QC JSON counts these rows. Its count of rows that lost right-primer bases
  (`rows_right_primer_removed`) includes the Ref and Alt rows of markers without reads.
- **Merged rows of different categories.** When rows of two categories become the same sequence,
  the row takes the category of the one with the most reads, or Ref or Alt if one of them is the
  canonical row. The QC JSON counts these merges (`rows_merged_across_allele_groups`). They are the
  only rows whose reads can move between the alleles of the target.
- **Sequences no longer match an untrimmed table.** Tools that compare sequences between runs need
  the same setting on both sides. `validate` against a vendor MADC is one. HapApp's allele
  database is the other, and the costly one: it matches whole sequences, so against a database
  built from tables with primers it reuses no fixed ID and issues a new one for every allele. It
  also writes the Ref and Alt rows from the database, with their primers, beside the other rows
  without them. Use one setting for every project that shares a database.
- **It cannot be combined with `--annotations off`.** The `Variants` column is the only record of
  how each row of such a table was read, and `madcats call` needs it. `generate` stops with an
  error if both are given. A step that drops the column later, such as HapApp's fixed-allele-ID
  output, leaves `call` with less to go on. See
  [Tables without primers](call.md#tables-without-primers).

## Examples

### DArTag with a reference-tag panel

```bash
madcats generate \
  --fastq-dir fastqs/ \
  --reference-tags reference_tags.csv \
  --cache-dir work/cache \
  --jobs 8 \
  --output results/project_MADC.csv
```

### Raw DArT delivery, HapApp-ready output

```bash
madcats generate \
  --fastq-dir raw_fastqs/ \
  --reference-tags reference_tags.csv \
  --leading-anchor auto \
  --dart-preamble \
  --sample-lookup sample_names.csv \
  --output results/project_MADC.csv
```

### High-depth AgriSeq after normalization

```bash
madcats generate \
  --fastq-dir normalized/ \
  --fastq-glob "*.fastq" \
  --reference-tags normalized/window_reference_tags.csv \
  --tag-length 41 \
  --min-supporting-freq 0.10 \
  --output results/project_MADC.csv
```

### Stricter discovery

```bash
madcats generate ... --reproducible-min-reads 20
```

### Exhaustive discovery

Setting every route's threshold to 0 emits every sequence with at least `--min-global-count` reads:

```bash
madcats generate ... --reproducible-min-reads 0 --depth-min-reads 0 --quality-min-clean-reads 0
```

## Output

| File | Contents |
|------|----------|
| `OUTPUT` | MADC table, described in [Output Formats](../reference/output-formats.md) |
| `OUTPUT.qc.json` | Parameters, provenance, per-sample read accounting |
| `OUTPUT.variants.csv` | One variant event per line (`--annotations on`) |

With `--remove-primers on`, all three describe the rows without their primers.

`OUTPUT.qc.json` and `OUTPUT.variants.csv` append to the full output name, so `--output MADC.csv`
writes `MADC.csv.qc.json` and `MADC.csv.variants.csv`.

## Tuning by Depth

Evidence thresholds are absolute read counts, so they do not scale with sequencing depth.

| Depth regime | Recommendation |
|--------------|----------------|
| Low (~60x) and standard DArTag depth | Defaults |
| Very high (above ~500x per sample per marker) | `--min-supporting-freq 0.10` |
| Unreplicated designs relying on the depth route | Consider `--min-supporting-freq 0.01` |
| Low-depth diploid, more sensitivity needed | `--reproducible-min-reads 5 --depth-min-reads 15 --min-supporting-freq 0.10`, confirming extra rows by segregation or replicates |

## See Also

- [Evidence Model](../reference/evidence-model.md)
- [Output Formats](../reference/output-formats.md)
- [madcats validate](validate.md)
