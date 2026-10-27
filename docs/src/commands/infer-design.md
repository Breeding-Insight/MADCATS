# madcats infer-design

Infer a reference-tag panel from the probe design plus FASTQ evidence, when the primer file is
unavailable.

## Synopsis

```bash
madcats infer-design --probe-design <FILE> --fastq-dir <DIR> --output <FILE> [OPTIONS]
```

## Description

The tag must be a tag-length window of the bracket-expanded target that contains the variant. So
inference only has to choose among a fixed set of candidate windows (every such window, on both
strands); it never assembles new sequences. Reads vote for windows in two tiers, in order:

1. The exact tier counts reads whose prefix exactly matches a candidate window. When a locus has any
   exact votes, they decide.
2. The fuzzy tier runs only for loci with no exact votes. It accepts reads that pass the
   `--min-identity` rule, and tags already resolved for other loci compete for those reads too, so
   one locus's reads are not counted for another.

A locus is resolved when its winning window has at least `--dominance-ratio` times the runner-up's
votes. The `LeftPrimerSequence` written to the CSV is the shortest tag prefix that is unique in the
panel (between `--anchor-min-length` and `--anchor-max-length`, and always ending before the
variant). Loci that cannot be resolved are left out, and the QC report gives the reason for each.

See [Panels Without a Primer File](../tutorials/primer-free.md) for a walkthrough and accuracy
figures.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--probe-design <FILE>` | Probe/marker design TSV or CSV |
| `--fastq-dir <DIR>` | Folder containing one FASTQ(.gz) per sample |
| `--output <FILE>` | Reference-tag CSV to write |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--fastq-glob <GLOB>` | `*` | FASTQ filename glob within `--fastq-dir` |
| `--leading-anchor <SEQ\|auto>` | none | Constant sequence preceding the amplicon in raw deliveries, or `auto` to detect it. Reads without it do not vote |
| `--tag-length <INT>` | 81 | Reported tag length of the panel (e.g. 54, 81, 109) |
| `--min-identity <FLOAT>` | 0.80 | Minimum identity for fuzzy-tier votes, as in `generate` |
| `--dominance-ratio <FLOAT>` | 10 | A winning window needs at least this many times the runner-up's votes |
| `--fuzzy-min-reads <INT>` | 30 | Fuzzy-tier votes needed to resolve a locus with no exact votes; `0` disables the fuzzy tier |
| `--fuzzy-min-samples <INT>` | 2 | Distinct samples that must contribute fuzzy votes |
| `--anchor-min-length <INT>` | 15 | Shortest pseudo-primer anchor to try |
| `--anchor-max-length <INT>` | 30 | Longest pseudo-primer anchor; also capped below each variant |
| `--jobs <INT>` | up to 4 | Samples processed in parallel |
| `--report <FILE>` | `OUTPUT.qc.json` | Inference QC JSON |

## Examples

### Basic usage

```bash
madcats infer-design \
  --probe-design probe-design.txt \
  --fastq-dir fastqs/ \
  --jobs 8 \
  --output inferred_reference_tags.csv
```

### 109-bp raw DArT delivery

```bash
madcats infer-design \
  --probe-design probe-design.txt \
  --fastq-dir raw_fastqs/ \
  --leading-anchor auto \
  --tag-length 109 \
  --output inferred_reference_tags.csv
```

### Exact evidence only

```bash
madcats infer-design \
  --probe-design probe-design.txt \
  --fastq-dir fastqs/ \
  --fuzzy-min-reads 0 \
  --output inferred_exact_only.csv
```

## Output

- `OUTPUT`: the reference-tag CSV, with one row per resolved locus.
- `OUTPUT.qc.json` (for example `inferred_reference_tags.csv.qc.json`): each locus's route
  (exact or fuzzy), support, margin, anchor length and flags, plus every excluded marker and the
  reason (no bracketed variant, flanks too short for the tag length, no reads within identity,
  support below threshold, ambiguous windows). Markers sharing a sequence-identical candidate
  window carry a paralog review flag.

## Tips

1. Check fuzzy-tier loci against a delivered MADC or a public panel resource when one is available.
2. Infer once, then reuse the CSV for every later order of the panel.
3. If you have a vendor MADC, use `scripts/madc_to_reference_tags.py` instead. It reconstructs the
   panel exactly.

## See Also

- [madcats design](design.md): use this instead when you have the primer file
- [Panels Without a Primer File](../tutorials/primer-free.md)
