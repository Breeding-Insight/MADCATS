# Output Formats

## MADC Table

`madcats generate --output` writes a CSV with one row per allele and one integer read-count column
per sample. It groups rows by marker and orders them Ref, Alt, RefMatch, AltMatch, Other.

```
AlleleID,CloneID,AlleleSequence,ClusterConsensusSequence,CallRate,...,RatioAvgCountRefAvgCountSnp,4404001,4404002,...
Chr01_000123|Ref,Chr01_000123,ACGT...,ACGT...,0.975,...,.,412,388,...
Chr01_000123|Alt,Chr01_000123,ACGT...,ACGT...,0.975,...,.,0,377,...
Chr01_000123|RefMatch,Chr01_000123,ACGT...,ACGT...,0,...,.,57,0,...
```

### Columns

| Column | Contents |
|--------|----------|
| `AlleleID` | Marker name and category joined by a pipe, e.g. `Chr01_000123\|AltMatch`. With `--number-match-ids`, Match/Other rows get a stable serial, such as `Chr01_000123\|RefMatch_0001` |
| `CloneID` | Marker name |
| `AlleleSequence` | The row's tag sequence |
| `ClusterConsensusSequence` | Canonical rows carry their own sequence. Match/Other rows carry the Ref tag, following vendor convention |
| `CallRate` | For canonical rows, the fraction of samples with canonical Ref or Alt reads at the locus. For other rows, `0` |
| `OneRatioRef` … `PICRef` | Vendor downstream statistics, not FASTQ evidence, so emitted as `.` |
| `PICSnp` → **`Variants`** | With `--annotations on` (default), the row's events against the Ref tag. Otherwise `.` |
| `AvgPIC` → **`MarkerType`** | With `--annotations on` (default), what the row carries beyond its canonical. Otherwise `.` |
| `AvgCountRef`, `AvgCountSnp` | For canonical rows, nonzero means of the Ref and Alt counts. For other rows, the row's own nonzero mean, in both columns |
| `RatioAvgCountRefAvgCountSnp` | `.` |
| one column per sample | Integer read counts |

MADCATS names sample columns by the FASTQ stem (or a unique `--sample-name-field` from target
manifests), or `--sample-lookup` renames them.

### Annotations

`--annotations on` (the default) repurposes two statistics columns the pipeline never computes. The
column *positions* are unchanged, so positional consumers see the same 16-column layout.
`--annotations off` restores the exact vendor header and writes no sidecar.

**`MarkerType`** records what the row carries beyond its canonical:

| Value | Meaning |
|-------|---------|
| `.` | The Ref row, or a row that only realises an ambiguity code |
| `snp` | Substitutions only (the Alt row of a designed SNP is `snp`) |
| `indel` | Gaps only (the Alt row of a designed indel is `indel`) |
| `both` | Substitutions and gaps |

To keep SNP-only microhaplotypes, drop `indel` and `both` rows by this label.

**`Variants`** lists the row's events against the Ref tag, `;`-separated, in 1-based coordinates on
the tag as printed:

| Token | Meaning |
|-------|---------|
| `91:A>T` | Substitution |
| `55:C>CAAAA` | Insertion, VCF-anchored on the preceding reference base, left-normalised |
| `70:TAG>T` | Deletion |
| `clip=4` | An insertion pushed the last 4 Ref bases out of the window |
| `overhang=3` | A deletion pulled 3 downstream bases into the window |

An example value is `55:C>CAAAA;91:A>T;clip=4`. The Ref row is `.`, and the Alt row of a designed SNP
is exactly the design (e.g. `46:C>G`).

## Variants Sidecar

`OUTPUT.variants.csv` (with `--annotations on`) lists the same events one per line, keyed by
`CloneID` and sequence so it survives allele-ID renumbering downstream:

| Column | Contents |
|--------|----------|
| `CloneID`, `AlleleID`, `AlleleSequence` | Row identity |
| `Position` | 1-based position on the Ref tag |
| `Ref`, `Alt` | Event alleles |
| `Type` | `snp`, `ins`, `del`, `clip` (Ref holds the unobserved bases) or `overhang` (Alt holds the pulled-in bases) |

## QC / Provenance JSON

`OUTPUT.qc.json` records:
- the pipeline version and every parameter
- input paths
- per-sample read accounting of records read, reads assigned, gapped and partial assignments, reads
  in emitted rows, reads moved by denoising, reads dropped and why
- the tile families (`tile_families`) and the extra credits they add (`tile_extra_credits_*`, described in
  [Read Accounting](evidence-model.md#read-accounting))
- row counts by category, by evidence route, and by `marker_types`
- the variants sidecar path

`normalize-reads`, `infer-design`, `window-scan` and `filter-paralogs` write their own QC JSONs. See each
command's page.

## Genome Labels (filter-paralogs)

[`madcats filter-paralogs`](../commands/filter-paralogs.md) writes a filtered copy of any MADC. It
copies kept rows unchanged. In flag mode (`--paralogs flag` or `--ambiguous flag`) it repurposes a
third statistics column, keeping the 16-column layout:

| Column | Contents |
|--------|----------|
| `PICRef` → **`GenomeLabel`** | `on_target`, `paralog`, `ambiguous`, `unmapped` or `unplaced`, depending on where the reference genome places the row's haplotype |

It also writes `OUTPUT.genome_labels.csv` (one line per input row, dropped rows included, with
the row's score at its marker and its best score elsewhere) and `OUTPUT.genome_markers.csv` (where
it placed each marker). The command page describes both.

## Genotype VCFs

`madcats call` writes up to three VCFs from a MADC, chosen with `--vcf`. They are `PREFIX.targets.vcf`
(the designed variant of each marker), `PREFIX.all.vcf` (targets plus other base substitutions in
the tags) and `PREFIX.haplotypes.vcf` (whole tag sequences as alleles). It also writes a QC JSON and,
for the haplotype VCF, a table tying each allele number to its MADC rows. Coordinates are
marker-local, so `CHROM` is the marker and `POS` a position on its Ref tag. The fields and the model behind them are
on the [command page](../commands/call.md#outputs).

## Vendor-Style Preamble (HapApp)

`--dart-preamble` prepends the seven starred rows that vendor raw MADCs carry above the header, so
preamble-expecting consumers such as HapApp and HapApp_utils accept the file:

1. service/project label (`--preamble-project`, by default the output filename before its first dot)
2. plate barcode
3. plate number
4. plate row
5. plate column
6. plate name
7. `targetid`, the FASTQ-derived sample ids, directly above the header's sample names

Every empty preamble cell is `*`, as in vendor output. The file's first cell holds a full
provenance string (pipeline, version, and every generation parameter as `key=value`), prefixed with
`* ` so consumers that skip starred lines still skip it. The file includes no timestamp, so output stays
byte-reproducible.

Combine with `--sample-lookup` to rename sample columns while the targetid row keeps the FASTQ ids,
matching the raw vendor layout.

## Downstream Processing

A MADCATS table is a raw MADC. Its allele IDs are per-project (`|RefMatch`, `|AltMatch`, ...) and
are not yet matched to the panel's fixed allele IDs. Run it through HapApp or HapApp_utils
(generate with `--dart-preamble`) before using it in BIGr, updog or other downstream tools. HapApp
assigns the fixed allele IDs and maintains the panel's cross-project allele database, and MADCATS does
neither. [`madcats call`](../commands/call.md) is the exception. It calls genotypes from the raw
MADC directly, identifying haplotypes by their sequences.

polyRAD's `readDArTag` is one example. In polyRAD 2.0.1 it rejects raw MADCs, from MADCATS or the
vendor. Match rows share allele IDs, and it reads the statistics columns as samples. It does read
the processed layout (unique allele IDs, then `AlleleID`, `CloneID`, `AlleleSequence` and the sample
columns).

The tested path is MADCATS → HapApp_utils (fixed allele IDs) → `BIGr::madc2vcf_all()` → updog.
Genotypes from this path agreed with independent reference genotypes at 98.42% on GT-seq and
95.13% on the 192-sample capture cohort.

HapApp's fixed-allele-ID step keeps only columns 1–3 and the sample columns, so the `Variants` and
`MarkerType` annotations do not carry through it. Use the `OUTPUT.variants.csv` sidecar if you need
them later.

## Format References

- [MADC file format](https://breeding-insight.github.io/BIGr_documentation/posts/madc-file-format/) from Breeding Insight
- [MADC data processing](https://breeding-insight.github.io/BIGr_documentation/posts/madc-data-processing/) from Breeding Insight
