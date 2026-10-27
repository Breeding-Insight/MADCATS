# Input Formats

MADCATS reads plain CSV/TSV design tables and FASTQ files. Column names are matched by header, so
extra columns are ignored.

## FASTQs

One FASTQ per sample, optionally gzipped, in a single folder (`--fastq-dir`). Restrict which files
are read with `--fastq-glob` (default `*`).

- `generate` expects **primer-started single-end reads**: every read begins at the assignment
  anchor. Raw DArT deliveries with a barcode and common sequence in front are handled with
  `--leading-anchor`.
- `normalize-reads` accepts any read that contains a panel window, in either orientation.
- `merge-pairs` expects mate files named by a shared stem plus `--r1-suffix` / `--r2-suffix`
  (defaults `_R1.fastq.gz` / `_R2.fastq.gz`), with mates in the same order.

Sample columns are named by the unique FASTQ filename stem unless a `--sample-lookup` table renames
them. Target manifests (`targets_*.csv`) in the FASTQ folder are picked up automatically for sample
metadata; they are optional.

## Primer CSV (DArT)

Used by `design` and by `generate --primers`.

| Column | Description |
|--------|-------------|
| `MarkerName` | Marker identifier, matching the probe design |
| `LeftPrimerSequence` | The primer every read starts with; may contain IUPAC codes for degenerate oligos |
| `RightPrimerSequence` | Reverse primer |
| `VariantRef`, `VariantAlt` | Alleles. Indels use dash notation: `-` on the empty side |
| `VariantPositionInTagSequence` | Position of the variant (or the start of the indel span) in the tag |

Tags are reconstructed as the primer sequence verbatim plus the genomic continuation, matching the
vendor's representation of degenerate primers.

## Probe Design (TSV/CSV)

Used by `design`, `infer-design`, `normalize-reads` and `window-scan`.

| Column | Description |
|--------|-------------|
| `MarkerName` | Marker identifier |
| `TargetSequence` | Target region with exactly one bracketed variant |
| `LocusGroup` | *Optional.* Declares markers as tiles of one physical locus |

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

Members of a `LocusGroup` must present the same target region (identical Ref expansions differing
only in which variant is bracketed). Violations are load-time errors.

## Reference Tags (CSV)

The platform-neutral panel definition consumed by `generate --reference-tags`. Written by `design`,
`infer-design` and `normalize-reads`, or curated by hand from any source.

| Column | Required | Description |
|--------|----------|-------------|
| `MarkerName` | yes | Marker identifier |
| `LeftPrimerSequence` | yes | Assignment anchor: the prefix reads are assigned by |
| `RefSequence` | yes | Full-length Ref tag |
| `AltSequence` | yes | Full-length Alt tag |
| `VariantPositionInTagSequence` | when Ref and Alt differ at more than one base | Validated to point at a real Ref/Alt difference |
| `VariantRef`, `VariantAlt` | no | Validated for IUPAC compatibility with the tag bases |
| `RightPrimerSequence` | no | |
| `LocusGroup` | no | Tile grouping; written only when some locus uses it |

Identical Ref/Alt sequences are rejected. Curated tags may contain IUPAC codes; see the
[evidence model](../reference/evidence-model.md#iupac-codes-and-n-bases) for what that means for
counts.

## Sample Lookup (CSV/TSV)

`generate --sample-lookup` takes a two-column table: FASTQ-derived sample id, then the desired
sample name. A header line is detected and skipped.

```
fastq_id,sample_name
4404075,Beauregard_rep1
4404076,Covington_rep1
```

The table is checked strictly, because a missing or wrong entry would put the wrong name on a count
column: every run sample must be mapped, ids that match no sample are errors, and duplicate names
are errors.

## Expected MADC

`validate --expected` reads any MADC, with or without the vendor's starred preamble rows: the header
is found by its `AlleleID` cell.
