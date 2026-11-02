# Filtering Paralogs with a Reference Genome

Some panels co-amplify paralogous copies of their targets. A paralog read that carries the target
base lands in the marker's RefMatch or AltMatch rows. Tools that sum every row at the target SNP
(`BIGr::madc2vcf_all` does) then count those reads as alleles of the target, and homozygous samples
look heterozygous. This tutorial cleans such a MADC with
[`madcats filter-paralogs`](../commands/filter-paralogs.md), which uses a reference genome to find
where each haplotype comes from.

The filter is opt-in: use it on panels that show the problem.

## Is My Panel Affected?

Signs of co-amplified paralogs:

- Inbred or doubled-haploid lines that should be homozygous show many heterozygous calls.
- Many markers have RefMatch/AltMatch rows that carry a large share of the marker's reads.
- Heterozygosity drops sharply when `generate --min-identity` is raised from 0.80 to 0.95 or 0.97.
  Stricter identity rejects the more diverged copies at read assignment.

The filter's QC JSON measures it: `rows.by_label.paralog.reads` is the number of
reads in rows the genome places elsewhere. In the lettuce 3K DArTag panel that was 63% of the reads
in the vendor MADC and 16% of a MADCATS MADC generated at `--min-identity 0.95`.

## 1. Install the Aligner

```bash
pip install 'madcats[genome]'
```

This adds `mappy`, minimap2's Python binding. No other command needs it.

## 2. Index the Genome Once

Use the assembly the panel was designed on, or one close to it. The first run with both
`--genome` and `--genome-index` builds the index and saves it. Later runs load it in seconds:

```bash
madcats filter-paralogs \
  --madc project_MADC.csv \
  --genome ref/genome.fa \
  --genome-index ref/genome.sr.mmi \
  --output project_MADC_filtered.csv
```

On the 2.4 Gb lettuce genome indexing took 50 s and 10.6 GB of memory. Run it on a machine with
enough memory, then reuse the `.mmi` anywhere.

## 3a. A MADCATS MADC

Generate as usual. On a paralog-heavy panel, a stricter identity gate removes the most diverged
copies before the genome filter sees them:

```bash
madcats generate --fastq-dir fastqs/ --reference-tags tags.csv \
  --min-identity 0.95 --jobs 8 --output project_MADC.csv

madcats filter-paralogs \
  --madc project_MADC.csv \
  --genome-index ref/genome.sr.mmi \
  --output project_MADC_filtered.csv \
  --jobs 8
```

On lettuce this pair (identity 0.95, then the genome filter dropping paralog and ambiguous rows)
matched the published realignment-based genotypes for PCA, kinship and linkage mapping. That was
measured with the bwa-mem prototype of this filter; this command keeps the same rows for 99.7% of
the reads.

## 3b. A Vendor MADC

The filter reads vendor DArT MADCs as they are delivered, with or without the starred preamble, and
keeps the preamble, the statistics columns and the column order:

```bash
madcats filter-paralogs \
  --madc DLet22-7611_MADC.csv \
  --genome-index ref/genome.sr.mmi \
  --output DLet22-7611_MADC_filtered.csv \
  --jobs 8
```

The lettuce vendor MADC (732,534 rows, 376 samples) took 2 min 15 s and 9.4 GB of memory. 113,655
rows were kept, with 36% of the reads. Genotype agreement with the published realignment VCF rose
from 89.9% to 98.5%. That figure was measured with the bwa-mem prototype of this filter; this
command keeps the same rows for 99.4% of the reads.

## 4. Read the QC

`OUTPUT.qc.json` summarises the run. On the lettuce MADCATS MADC (17,847 rows):

| Label | Rows | Reads |
|-------|-----:|------:|
| `on_target` | 12,687 | 59.8 M (82.4%) |
| `paralog` | 4,135 | 11.6 M (16.0%) |
| `ambiguous` | 1,025 | 1.1 M (1.6%) |

Check:

- `markers.placement` and `markers.clone_id_check`: nearly every marker should be placed by its Ref
  tag, and CloneID positions should agree. Many disagreements mean the genome is not the panel's
  assembly; the command warns when more than 10% disagree.
- `markers.multicopy_ref_tags`: markers whose own Ref tag has several equal-best sites. Their
  canonical rows are kept, but their counts mix copies; consider dropping these markers downstream.
- `rows.by_outcome.protected`: canonical Ref/Alt rows that the genome would have dropped. They are
  always kept, so dosage tools still see both designed alleles.
- `rows.by_category`: which allele categories lost reads. Paralogs concentrate in RefMatch and
  AltMatch rows.

## 5. Review Before Dropping (Optional)

To see the labels without removing anything, flag instead of dropping:

```bash
madcats filter-paralogs --madc project_MADC.csv --genome-index ref/genome.sr.mmi \
  --paralogs flag --ambiguous flag --output project_MADC_flagged.csv
```

Every row then carries its label in the `GenomeLabel` column (the `PICRef` position). The
`OUTPUT.genome_labels.csv` sidecar lists, for every row, its score at the marker and its best score
elsewhere. A perfect 81-base match scores 162 and each mismatch costs 10, so a row with
`TargetScore` 132 and `OffTargetScore` 162 is three mismatches from the target and a perfect copy
of another locus.

`--ambiguous keep` keeps ambiguous rows without a label. On lettuce, keeping them left false
heterozygous calls at 3.4% of inbred genotypes, against 0.87% when they were dropped, which is why
`drop` is the default.

## 6. Continue Downstream

The cleaned file keeps the layout it came in. A cleaned raw MADC continues as any MADCATS or vendor
MADC does: HapApp or HapApp_utils for fixed allele IDs, then `BIGr::madc2vcf_all()` and dosage
calling. See [Output Formats](../reference/output-formats.md#downstream-processing).

The filter also reads HapApp's fixed-allele-ID output, so it can run after HapApp instead:

```bash
madcats filter-paralogs \
  --madc project_MADC_rename_v1.csv \
  --genome-index ref/genome.sr.mmi \
  --output project_MADC_rename_v1_filtered.csv
```

The fixed allele IDs and the column layout are kept, so the cleaned file goes wherever HapApp's
output would have gone. On the lettuce MADCATS MADC, rewritten in HapApp's layout, every row got
the same label as in the raw MADC.

## How the Labels Are Decided

Each marker is placed by aligning its Ref tag. The CloneID's `CHROM_POS` checks the placement and
breaks ties between equal-best sites. Every distinct haplotype is then aligned genome-wide. A
haplotype with no hit near its marker is realigned to the marker sites alone with short seeds:
real alleles that differ from the reference by several SNPs are too diverged for the genome-wide
seeds, and would otherwise be judged only by their paralogs. The haplotype's best score at its
marker is then compared with its best score anywhere else: higher at the marker is `on_target`,
higher elsewhere is `paralog`, and a tie is `ambiguous`.

Labels were checked against an independent bwa-mem classification of the same lettuce haplotypes.
The drop or keep decision agreed for 98.2% of the rows of the MADCATS MADC and 94.9% of the vendor
MADC. The vendor MADC's rows are more diverged, and bwa mem finds more distant copies than
minimap2. Weighted by reads, the filtered vendor MADC agreed with the bwa-based prototype for 99.4%
of reads.
