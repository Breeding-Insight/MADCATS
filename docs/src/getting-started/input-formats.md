# Input Formats

MADCATS reads CSV/TSV design tables and FASTQ files. It finds columns by their header names, so it
ignores extra columns.

## FASTQs

One FASTQ per sample, optionally gzipped, in a single folder (`--fastq-dir`). Restrict which files
MADCATS reads with `--fastq-glob` (default `*`).

- `generate` expects **primer-started single-end reads**, so every read must begin at the
  assignment anchor. `--leading-anchor` handles raw DArT deliveries with a barcode and common
  sequence in front.
- `normalize-reads` accepts any read that contains a panel window, in either orientation.
- `merge-pairs` expects mate files named by a shared stem plus `--r1-suffix` / `--r2-suffix`
  (defaults `_R1.fastq.gz` / `_R2.fastq.gz`), with mates in the same order.

MADCATS names sample columns by the unique FASTQ filename stem unless a `--sample-lookup` table
renames them. It picks up target manifests (`targets_*.csv`) in the FASTQ folder automatically for
sample metadata. They are optional.

## Primer CSV (DArT)

Used by `design` and by `generate --primers`.

| Column | Description |
|--------|-------------|
| `MarkerName` | Marker identifier, matching the probe design |
| `LeftPrimerSequence` | The primer every read starts with, which may contain IUPAC codes for degenerate oligos |
| `RightPrimerSequence` | Reverse primer |
| `VariantRef`, `VariantAlt` | Alleles. Indels use dash notation, with `-` on the empty side |
| `VariantPositionInTagSequence` | Position of the variant (or the start of the indel span) in the tag |

MADCATS reconstructs each tag as the primer sequence verbatim plus the genomic continuation,
matching the vendor's representation of degenerate primers.

## Probe Design (TSV/CSV)

Used by `design`, `infer-design`, `normalize-reads` and `window-scan`.

| Column | Description |
|--------|-------------|
| `MarkerName` | Marker identifier |
| `TargetSequence` | Target region with exactly one bracketed variant |
| `LocusGroup` | *Optional.* Declares markers as tiles of one physical locus. See [Tiled Markers](#tiled-markers) |

The bracketed variant is a SNP such as `[C/T]` or an anchored-base indel such as `[A/ACTTG]` or
`[TTTTTTG/T]`. The first allele is the reference. Any DArT marker table with `MarkerName` and a
bracket-notation `TargetSequence` serves, extra markers and columns included.

```
MarkerName	TargetSequence
Chr01_000123	TGCAGGTCAAGTTCAAGGT[C/T]GATCCATTGACAAGG...
```

Panel construction rejects:
- designs whose variant lies inside the left primer (Alt reads could never carry the primer)
- identical alleles
- SNP designs whose constructed Alt tag is not among the bracket expansions (a primer-table /
  probe-design disagreement)

For `normalize-reads` and `window-scan`, members of a `LocusGroup` must present the same target
region (identical Ref expansions differing only in which variant is bracketed). Violations are
load-time errors. `design` and `generate --primers --probe-design` copy `LocusGroup` onto the panel
without that requirement.

## Reference Tags (CSV)

This is the platform-neutral panel definition that `generate --reference-tags` consumes. `design`,
`infer-design` and `normalize-reads` write it, or you can curate it by hand from any source.

| Column | Required | Description |
|--------|----------|-------------|
| `MarkerName` | yes | Marker identifier |
| `LeftPrimerSequence` | yes | Assignment anchor, the prefix MADCATS assigns reads by |
| `RefSequence` | yes | Full-length Ref tag |
| `AltSequence` | yes | Full-length Alt tag |
| `VariantPositionInTagSequence` | when Ref and Alt differ at more than one base | Validated to point at a real Ref/Alt difference |
| `VariantRef`, `VariantAlt` | no | Validated for IUPAC compatibility with the tag bases |
| `RightPrimerSequence` | no | |
| `LocusGroup` | no | Tile grouping, written only when some locus uses it. See [Tiled Markers](#tiled-markers) |

MADCATS rejects identical Ref/Alt sequences. Curated tags may contain IUPAC codes, and the
[evidence model](../reference/evidence-model.md#iupac-codes-and-n-bases) explains what that means
for counts.

## Tiled Markers

DArTag panels often tile several target SNPs inside one amplicon. The tiles share a left primer
(or one primer is the other plus a base or two), so every read from the amplicon starts with both,
and each tile's tag carries its neighbours' SNPs as IUPAC codes:

```
MarkerName,LeftPrimerSequence,RefSequence,...,VariantPositionInTagSequence,LocusGroup
Chr1_121673941,GGAAGACATGCTGAGGGYAC,GGAAGACATGCTGAGGGYACGTGTCATCRATTTCGG...,23,amp_Chr1_121673
Chr1_121673947,GGAAGACATGCTGAGGGYAC,GGAAGACATGCTGAGGGYACGTKTCATCGATTTCGG...,29,amp_Chr1_121673
```

A read from this amplicon fits both tiles' tags equally well. `generate` drops such ties as
ambiguous, so without tile handling both markers would get zero reads. Instead, it credits the read
to every tile it fits, as DArT's own pipeline does (see the
[evidence model](../reference/evidence-model.md#tiled-markers)).

Two markers are tiles of one amplicon when they **read the same tag window**. Their tags have the
same length, one marker's left primer matches the start of the other's tag, and their Ref tags
agree, up to IUPAC codes, at every position except the two markers' own variant positions. A tile
family is a set of markers linked this way. `generate --tile-groups` chooses which links it uses:

| `--tile-groups` | Linked markers |
|-----------------|----------------|
| `auto` (default) | Markers sharing a `LocusGroup`, plus undeclared markers that pass the same-window test. Markers declared in different groups are never linked |
| `declared` | Only markers sharing a `LocusGroup`, and only when they pass the same-window test |
| `off` | None, so every tie is dropped as ambiguous |

Declare `LocusGroup` in the reference-tag CSV, or in the probe design when you build the panel with
`--primers`/`--probe-design`. `design` copies it into the reference tags it writes. Detection uses
only the primers and tags. A reference-tag CSV has no chromosome coordinates, and two markers that
read the same tag window tie wherever their designs place them. `generate` names every family
on stderr and in the QC JSON (`tile_families`). With `declared` or `off` it warns about markers
that pass the same-window test but stay unlinked.

`normalize-reads` excises tiles of one `LocusGroup` as different windows of the locus. They
never read the same tag window, so `generate` keeps giving each of their reads to one marker.

## Sample Lookup (CSV/TSV)

`generate --sample-lookup` takes a two-column table with the FASTQ-derived sample id, then the
desired sample name. It detects and skips a header line.

```
fastq_id,sample_name
4404075,Beauregard_rep1
4404076,Covington_rep1
```

MADCATS checks the table strictly, because a missing or wrong entry would put the wrong name on a
count column. Every run sample must be mapped, ids that match no sample are errors, and duplicate names
are errors.

## Expected MADC

`validate --expected` reads any MADC, with or without the vendor's starred preamble rows, because it
finds the header by its `AlleleID` cell.
