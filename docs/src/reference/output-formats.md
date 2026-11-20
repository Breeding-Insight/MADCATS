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
| `AlleleSequence` | The row's tag sequence. With `--remove-primers on`, the tag without its primers |
| `ClusterConsensusSequence` | Canonical rows carry their own sequence. Match/Other rows carry the Ref tag (the Ref row's sequence), following vendor convention |
| `CallRate` | For canonical rows, the fraction of samples with canonical Ref or Alt reads at the locus. For other rows, `0` |
| `OneRatioRef` … `PICRef` | Vendor downstream statistics, not FASTQ evidence, so emitted as `.` |
| `PICSnp` → `Variants` | With `--annotations on` (default), the row's events against the Ref tag. Otherwise `.` |
| `AvgPIC` → `MarkerType` | With `--annotations on` (default), what the row carries beyond its canonical. Otherwise `.` |
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
| `unresolved` | The row's end is not read (`unresolved=k` in `Variants`), whatever the row carries before it. A gap there is likely but unverified |

To keep SNP-only microhaplotypes, drop `indel`, `both` and `unresolved` rows by this label.

**`Variants`** lists the row's events against the Ref tag, `;`-separated, in 1-based coordinates on
the tag as printed. They come from the alignment the row's category was read from, so the two show
the same base on the target. One exception: in a table without primers, a row of another length
than the Ref row is read again as printed, which can differ from the reading its category came
from. The row still holds the base its category was read from.

| Token | Meaning |
|-------|---------|
| `91:A>T` | Substitution |
| `55:C>CAAAA` | Insertion, VCF-anchored on the preceding reference base, left-normalised |
| `70:TAG>T` | Deletion |
| `clip=4` | An insertion pushed the last 4 Ref bases out of the window |
| `overhang=3` | A deletion pulled 3 downstream bases into the window |
| `unresolved=9` | The last 9 Ref bases are not read in this row. Its end has two explanations, a gap or substitutions, with too few bases after the gap to tell them apart ([Classification](evidence-model.md#classification)). Always the last token |
| `1:.>GG`, `1:TT>.` | An insertion before the very first base, or a deletion that starts on it. No base precedes it, so it is written without an anchor base, and `.` stands for no bases. Here `GG` is inserted before base 1, and bases 1 and 2 are deleted |
| `31:G>A;32:CACG>.`, `31:G>A;32:.>TT` | An insertion or deletion right after a substituted base. The substitution uses that base, so the gap is written without an anchor base as well, from its own first base. Here bases 32 to 35 are deleted, or `TT` is inserted before base 32. It appears on a row that carries the Alt base with a gap directly after it |

An example value is `55:C>CAAAA;91:A>T;clip=4`. The Ref row is `.`, and the Alt row of a designed SNP
is exactly the design (e.g. `46:C>G`).

Each `position:REF>ALT` token reads "replace these Ref bases, starting at this position, by these
bases". `clip=k` drops the last `k` Ref bases, and `overhang=k` appends the row's own last `k`
bases. `unresolved=k` drops the last `k` Ref bases too, and then appends whatever of the row has
not been built yet. No two tokens of a row use the same Ref base, so applying them in order
rebuilds the row from the Ref row.

On a marker whose own design is an insertion or deletion, every row is written through gaps
against the Ref row, and the unread end is counted on the Ref row as everywhere else. The designed
insertion or deletion does not settle another gap near the row's end.

A row with `unresolved=k` says nothing about Ref positions from `length - k + 1` on: no
substitution, no insertion or deletion, and no agreement with the Ref row either. A consumer that
counts reads per position should leave the row's reads out there, as it does for `clip`.
A reader that does not know a token should not guess the row from its sequence: the token states
something about the row that the sequence alone does not show.

With `--remove-primers on` the frame is the Ref row as printed, so position 1 is the first base
after the left primer. Whether a row is read through an insertion or deletion is decided on its
full tag, as without the option, and a row whose length differs from the Ref row's is always read
through one. A row cut at its right primer, on a marker whose Ref row also ends at its primer, is
compared with the Ref row end to end. A difference at its end is then an insertion or a deletion.
Every other row is a fixed-length window and keeps `clip` and `overhang`. A row that differs from
the Ref row in length by more than `--max-indel-length` is compared base by base from its start,
and the rest is `clip` or `overhang`.

An unresolved end is decided on the full tag as well, and carried: the row keeps `unresolved=k`,
counted on the Ref row as printed, whatever length it was cut to. A right primer found at another
place shows that the row gained or lost bases, not which reading of its end was right. Only an
unread stretch that lay wholly in the primer cut from the Ref row leaves nothing to mark.

## Variants Sidecar

`OUTPUT.variants.csv` (with `--annotations on`) lists the same events one per line, keyed by
`CloneID` and sequence so it survives allele-ID renumbering downstream:

| Column | Contents |
|--------|----------|
| `CloneID`, `AlleleID`, `AlleleSequence` | Row identity |
| `Position` | 1-based position on the Ref tag |
| `Ref`, `Alt` | Event alleles |
| `Type` | `snp`, `ins`, `del`, `clip` (Ref holds the unobserved bases), `overhang` (Alt holds the pulled-in bases) or `unresolved` (Position is the first Ref base not read, Ref the Ref bases from it on, Alt the row's bases that stand there) |

## QC / Provenance JSON

`OUTPUT.qc.json` records:
- the pipeline version and every parameter
- input paths
- per-sample read accounting of records read, reads assigned, gapped and partial assignments, reads
  in emitted rows, reads moved by denoising, reads dropped and why
- the tile families (`tile_families`) and the extra credits they add (`tile_extra_credits_*`,
  described in [Read Accounting](evidence-model.md#read-accounting))
- row counts by category, by evidence route, and by `marker_types`
- `unresolved_ends`: the rows whose end is not read and their reads, and how many of them have the
  target in that end (those rows are `Other`). Counted on the full tags, so it is the same with and
  without `--remove-primers on` and with `--annotations off`
- `read_continuation`: the number of bases after the tag that reads were compared over
  (`bases_after_tag`, 0 when reads end at the tag), and the rows and reads whose unread target
  they settled ([Classification](evidence-model.md#classification))
- the variants sidecar path
- with `--remove-primers on`, `primer_removal`: the panel's amplicon length, the markers by whether
  their right primer was located (and if not, why), the rows that were merged, the rows that lost
  right-primer bases, and the rows that kept their right primer

`normalize-reads`, `infer-design`, `window-scan` and `filter-paralogs` write their own QC JSONs. See
each command's page.

## Genome Labels (filter-paralogs)

[`madcats filter-paralogs`](../commands/filter-paralogs.md) writes a filtered copy of any MADC. It
copies kept rows unchanged. In flag mode (`--paralogs flag` or `--ambiguous flag`) it repurposes a
third statistics column, keeping the 16-column layout:

| Column | Contents |
|--------|----------|
| `PICRef` → `GenomeLabel` | `on_target`, `paralog`, `ambiguous`, `unmapped` or `unplaced`, depending on where the reference genome places the row's haplotype |

It also writes `OUTPUT.genome_labels.csv` (one line per input row, dropped rows included, with
the row's score at its marker and its best score elsewhere) and `OUTPUT.genome_markers.csv` (where
it placed each marker). The command page describes both.

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
neither.

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
