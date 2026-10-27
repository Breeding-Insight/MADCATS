# DArTag from Design Files

This tutorial walks through a complete DArTag project using the vendor's primer and probe-design
files. The numbers quoted come from the 80-sample hexaploid sweetpotato project MADCATS was first
developed on (3,120 markers, 81-bp tags, 154.6M reads).

## Inputs

```
project/
├── fastqs/                                   # one FASTQ(.gz) per sample
│   ├── 4404001.fastq.gz
│   ├── 4404002.fastq.gz
│   └── ...
├── SweetPotato_DArTag_primer_sequence.csv    # primer CSV
└── Sweetpotato_3K_SNPs_DArTag-probe-design.txt
```

## 1. Reconstruct the Panel

```bash
madcats design \
  --primers SweetPotato_DArTag_primer_sequence.csv \
  --probe-design Sweetpotato_3K_SNPs_DArTag-probe-design.txt \
  --output sweetpotato_reference_tags.csv
```

For each marker, `design` locates the left primer in the target sequence (either orientation),
extracts the primer-started tag that carries the bracketed variant, verifies the variant position,
and builds the Alt tag. If the design files disagree with each other, the run stops with an error.

Open the CSV and check it: one row per marker with `RefSequence` and `AltSequence` of the tag
length. This file is now the panel. You can hand it to collaborators or use it on later orders
without the vendor files.

## 2. Generate the MADC

```bash
madcats generate \
  --fastq-dir fastqs/ \
  --reference-tags sweetpotato_reference_tags.csv \
  --cache-dir work/cache \
  --jobs 8 \
  --output results/sweetpotato_MADC.csv
```

What happens, in order:

1. Each read is assigned to a marker by its left-primer prefix (one mismatch allowed by default).
2. The fixed-length tag is kept and accepted if it is within 80% identity of the marker's Ref or
   Alt tag. Reads that fail get one more try with a gapped alignment, so undesigned indels are
   still counted.
3. Reads shorter than the tag are added to the count of a full-length allele whose start they
   match exactly. They cannot create a new allele.
4. Canonical Ref and Alt rows are always written. Any other sequence is reported only if it passes
   one of the [evidence routes](../reference/evidence-model.md).
5. Rare sequences that differ by one base from a much more common sequence in the same sample are
   treated as sequencing errors and merged into it (abundance-skew denoising).

## 3. Read the QC

`results/sweetpotato_MADC.csv.qc.json` accounts for every read. Look first at the per-sample totals:
samples with very few assignable tags (on the sweetpotato project, twelve samples had fewer than
100,000; the lowest had 23) usually contain mostly short or primer-less reads. Check them for
missing data before using them downstream.

## 4. Optional: Compare with the Vendor MADC

```bash
madcats validate \
  --generated results/sweetpotato_MADC.csv \
  --expected Report-DSp25-10607/DSp25-10607_MADC.csv \
  --report results/validation.json
```

The report gives row precision and recall by allele category plus cell-level count agreement. On
this project every shared row had identical classification, and shared-cell count correlation was
0.9997.

Do not expect 100% row recall against a vendor table: 27.3% of the vendor's positive-count rows on
this project had fewer than 10 reads in fewer than 2 samples of the vendor's *own* counts, below
any within-project evidence threshold.

## Tuning

The defaults are the recommendation for diploid and polyploid DArTag projects alike.

| Goal | Change |
|------|--------|
| Stricter, higher-precision discovery on a high-depth project | `--reproducible-min-reads 20` |
| Very high depth (above ~500x per sample per marker) | `--min-supporting-freq 0.10` |
| Count only full-length reads | `--partial-read-mode off` |
| Exact vendor header, no annotation columns | `--annotations off` |
| HapApp input | `--dart-preamble` (and `--sample-lookup` to rename columns) |

## Raw DArT Deliveries

Older DArT orders deliver reads with an inline barcode and the library's common sequence
`GATGGATCCAGTGC` in front of the amplicon, so nothing is primer-started. MADCATS checks the first
FASTQ before any read pass. When almost no read starts with a panel tag but most carry one after a
constant prefix, the run stops and names that prefix. Re-run with:

```bash
madcats generate ... --leading-anchor auto
# or explicitly
madcats generate ... --leading-anchor GATGGATCCAGTGC
```

Reads without the anchor are counted as `reads_no_leading_anchor` in the QC.
