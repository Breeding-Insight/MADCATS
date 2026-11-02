# madcats filter-paralogs

Drop or flag MADC haplotypes that a reference genome places at another locus. Optional; needs the
`mappy` aligner.

## Synopsis

```bash
madcats filter-paralogs --madc <FILE> --genome <FASTA> [--genome-index <MMI>] --output <FILE> [OPTIONS]
madcats filter-paralogs --madc <FILE> --genome-index <MMI> --output <FILE> [OPTIONS]
```

## Description

In some panels, paralogous copies of a target co-amplify with it, and their reads land in the
marker's rows. A paralog that differs from the Ref tag only at non-target bases becomes a RefMatch
or AltMatch row. When counts are summed per target SNP (as `BIGr::madc2vcf_all` does), those rows
push homozygous samples towards false heterozygous calls. `filter-paralogs` uses a reference genome
to decide where each haplotype comes from:

1. **Place each marker.** Its Ref tag is aligned to the genome. A CloneID of the form `CHROM_POS`
   (zero padding allowed, e.g. `Chr1_000333995`) is a cross-check. When the Ref tag has several
   equal-best sites, the CloneID position picks one; if it cannot, the marker is reported as
   multi-copy. When the Ref tag does not align, a window around the CloneID position is used.
2. **Align every distinct haplotype.** Each haplotype gets two scores: its best alignment at its
   marker's site, and its best alignment anywhere else. A haplotype the genome-wide search found
   nowhere near its marker is realigned to the marker sites alone with short seeds, so a real allele
   carrying several SNPs is still scored at home.
3. **Label it**:

   | Label | Meaning | Default action |
   |-------|---------|----------------|
   | `on_target` | The marker's site scores strictly best | kept |
   | `paralog` | Another site scores strictly better, or the haplotype aligns only elsewhere | dropped (`--paralogs`) |
   | `ambiguous` | The marker's site ties with another site, or the marker's own Ref tag is multi-copy and unresolved | dropped (`--ambiguous`) |
   | `unmapped` | No alignment anywhere | kept |
   | `unplaced` | The marker itself could not be placed | kept |

4. **Write the cleaned MADC.** Canonical Ref and Alt rows are never dropped, because downstream
   dosage tools need both designed alleles. A canonical row whose label would have dropped it is
   counted as `protected`. Kept rows are copied unchanged, and the output keeps the input's layout:
   the vendor preamble, the statistics columns, HapApp's `Code_version` column and the column order.

The filter is opt-in, meant for panels known to co-amplify paralogs. Check the share of reads in
`paralog` rows in the QC JSON: in the lettuce 3K panel it was 63% of all reads in the vendor MADC.
On panels without the problem the filter removes little and is not needed.

It reads three layouts: a MADC from `generate`, a vendor DArT MADC, and the fixed-allele-ID MADC
that HapApp writes (`Code_version`, `AlleleID`, `CloneID`, `AlleleSequence`, then the samples, with
IDs such as `|Ref_0001`). It can therefore run on the raw MADC before HapApp, or on HapApp's output.

Scores use minimap2's short-read scoring (match +2, mismatch −8, gap open −12 and extension −2,
ambiguous base −1), recomputed from each alignment, so ties are exact.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--madc <FILE>` | MADC to clean: from `generate`, a vendor DArT MADC (with or without its starred preamble), or a HapApp fixed-allele-ID MADC. Gzipped input is read directly |
| `--output <FILE>` | Cleaned MADC to write (must differ from `--madc`) |

Give `--genome`, `--genome-index`, or both.

## Genome

| Option | Default | Description |
|--------|---------|-------------|
| `--genome <FASTA>` | none | Reference genome FASTA, optionally gzipped |
| `--genome-index <MMI>` | none | minimap2 index. Loaded when it exists; otherwise built from `--genome` and saved here for the next run |

Indexing a large genome is slow and memory-heavy, so build the index once and reuse it. The first
run with both options writes it:

```bash
madcats filter-paralogs --madc MADC.csv --genome genome.fa --genome-index genome.sr.mmi --output MADC_filtered.csv
```

Later runs need only `--genome-index`. An index from the minimap2 command line works too, built with
the short-read preset (`minimap2 -x sr -d genome.sr.mmi genome.fa`). On the 2.4 Gb lettuce genome
indexing took 50 s and 10.6 GB of memory, and the index file is 5.7 GB.

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--paralogs {drop,flag}` | `drop` | `paralog` rows: drop them, or keep them labelled in the `GenomeLabel` column |
| `--ambiguous {drop,flag,keep}` | `drop` | `ambiguous` rows: drop, keep labelled, or keep unlabelled |
| `--qc-output <FILE>` | `OUTPUT.qc.json` | QC JSON |
| `--jobs <INT>` | up to 4 | Alignment threads |

Dropping ambiguous rows is the default because, on lettuce, keeping them left false heterozygous
calls at 3.4% of inbred genotypes, against 0.87% when they were dropped.

## Flag Mode

With either option set to `flag`, the `PICRef` column is renamed `GenomeLabel` and every row
carries its label there: `on_target`, `paralog`, `ambiguous`, `unmapped` or `unplaced`. Flagged
rows are kept. Column positions do not change, so positional consumers still see the 16-column
layout. In a vendor MADC this overwrites the vendor's `PICRef` statistic. Flag mode needs the
`PICRef` (or an earlier `GenomeLabel`) column. HapApp's fixed-allele-ID MADC has none, so on that
layout use drop mode and read the labels sidecar.

A flagged file is for review. Tools that check the vendor's column names, such as HapApp-online's
upload check, reject `GenomeLabel` in place of `PICRef`; give them the drop-mode output.

## Outputs

| File | Contents |
|------|----------|
| `OUTPUT` | The cleaned MADC |
| `OUTPUT.genome_labels.csv` | One line per input row, dropped rows included |
| `OUTPUT.genome_markers.csv` | One line per marker: where it was placed and how |
| `OUTPUT.qc.json` | Parameters, aligner settings, marker placement summary, rows and reads per label, per outcome and per allele category, and `count_overflow` (see below) |

Read totals are summed as stored. If the input has negative counts, the command warns and
`count_overflow` records how many cells and markers carry them. A negative count with a maximum of
exactly 127 points to signed 8-bit overflow in a vendor file; see the
[FAQ](../reference/faq.md#my-vendor-madc-has-negative-read-counts).

When the input carries a MADCATS provenance cell (`generate --dart-preamble`), the filter's
settings are appended to it.

### Labels sidecar

| Column | Contents |
|--------|----------|
| `CloneID`, `AlleleID`, `AlleleSequence` | Row identity |
| `Label` | `on_target`, `paralog`, `ambiguous`, `unmapped` or `unplaced` |
| `Action` | `kept`, `flagged`, `dropped`, or `protected` (a canonical row kept despite its label) |
| `TargetScore` | Best alignment score at the marker's site (`.` if none) |
| `TargetFrom` | `genome` (found by the genome-wide search), `site` (by the realignment to the marker sites) or `.` |
| `OffTargetScore` | Best alignment score anywhere else (`.` if none) |
| `OffTargetSites` | Distinct sites elsewhere at that score |
| `OffTargetHit` | The best site elsewhere, 1-based: `Chr4:61708465-61708545(+)` |

A perfect 81-base match scores 162 and each mismatch costs 10, so `TargetScore` and
`OffTargetScore` show how many mismatches separate a haplotype from each copy.

### Markers sidecar

| Column | Contents |
|--------|----------|
| `CloneID` | Marker |
| `Placement` | `ref_tag`, `ref_tag+clone_id` (multi-copy Ref tag resolved by the CloneID position), `ref_tag_multicopy` (unresolved), `clone_id` (Ref tag did not align) or `unplaced` |
| `RefTagSites` | Equal-best genome sites of the Ref tag (0 when it did not align) |
| `CloneIDCheck` | `agree`, `disagree`, `no_position` (CloneID is not `CHROM_POS`), `contig_absent` or `fallback` |
| `Sites` | The marker's site(s), 1-based, `;`-separated |

When the alignment and the CloneID disagree, the alignment wins, since CloneID coordinates may come
from another assembly version. The command warns when more than 10% of markers disagree, which
usually means the genome is not the assembly the panel was designed on.

## Example

```bash
madcats filter-paralogs \
  --madc DLet22-7611_MADC.csv \
  --genome-index ref/Lsativa_467_v8_Chr.sr.mmi \
  --output DLet22-7611_MADC_filtered.csv \
  --jobs 8
```

On the lettuce 3K DArT MADC (732,534 rows, 376 samples) this took 2 min 15 s and 9.4 GB of memory.
It kept 113,655 rows, with 35.7% of the reads. See the
[paralog filter tutorial](../tutorials/paralog-filter.md) for the whole workflow.

## Limitations

- The genome must be the assembly the panel's tags come from, or one close to it. Markers whose
  region is missing or collapsed in the assembly cannot be judged.
- A copy that exists in your germplasm but not in the reference assembly cannot be seen; its
  haplotypes look like diverged alleles of the target.
- minimap2 finds fewer distant copies than bwa mem, so a haplotype far from every genome copy is
  judged on the hits it has. Labels agreed with a bwa-mem classification for 98.2% of rows of a
  MADCATS MADC and 94.9% of a vendor MADC, whose rows are more diverged.
- Only the tag is aligned (81 bases by default), not the full read.

## See Also

- [Paralog filter tutorial](../tutorials/paralog-filter.md)
- [madcats generate](generate.md)
- [Output Formats](../reference/output-formats.md)
