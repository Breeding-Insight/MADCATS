# MADCATS <img src="https://raw.githubusercontent.com/Breeding-Insight/MADCATS/main/docs/src/images/madcats-logo.png" align="right" width="250" alt="MADCATS logo" />

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI](https://img.shields.io/pypi/v/madcats.svg)](https://pypi.org/project/madcats/)
[![Tests](https://github.com/Breeding-Insight/MADCATS/actions/workflows/tests.yml/badge.svg)](https://github.com/Breeding-Insight/MADCATS/actions/workflows/tests.yml)
[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue.svg)](https://madcats.readthedocs.io)

**Marker-Anchored Discovery and Counting of Alleles from Targeted Sequencing.**

`madcats` is a command-line tool that turns targeted-sequencing FASTQs into
[MADC](https://softwarehub.breedinginsight.org/learn/what-is-a-madc-file.html) microhaplotype count
tables. A MADC (Missing Allele Discovery Counts) table records the read count of every
microhaplotype in every sample. A microhaplotype is the short sequence a marker amplifies, so a
marker can have more than two alleles, and the variants beside its target SNP are kept. `madcats`
counts the alleles a panel was designed for and finds additional ones in the reads.

The output is a raw MADC table. It goes through HapApp, which assigns fixed allele IDs for building
a microhaplotype database and comparing sequencing runs. For existing downstream pipelines, a MADC
converts to a VCF of SNP genotype calls, with
[BIGr and updog](https://softwarehub.breedinginsight.org/learn/madc-to-vcf.html) after HapApp or
with `madcats call` (experimental) straight from the raw table.

The same code handles DArTag, AgriSeq (Ion Torrent), GT-seq, rhAmpSeq and Flex-Seq capture data.
It needs only the Python 3 standard library. Every analysis threshold is a command-line option, and
each run writes a JSON file recording its parameters and read counts.

> **A prior MADC is not an input.** Generation uses only the FASTQs plus the assay design, either the
> DArTag primer/probe files, a probe design alone, or a platform-neutral Ref/Alt tag manifest.

## Highlights

- Closely matches vendor output. Read counts correlate at r = 0.9997 with the vendor MADC on a
  sweetpotato DArTag panel, and dosage calls agree 95–99% across five DArTag panels.
- Builds the panel from vendor design files (`madcats design`), from the probe design and the reads
  when there is no primer file (`madcats infer-design`), or from your own Ref/Alt tag CSV.
- Reports a new microhaplotype only when several samples carry it, one sample has many reads of it,
  or (optionally) it has enough high-quality reads.
- Merges sequencing errors back into the alleles they came from, without changing any sample's
  Ref/Alt dosage.
- Handles indels, including ones the panel was not designed for.
- Works with AgriSeq-style reads that don't start at a primer (`madcats normalize-reads`) and with
  paired-end reads (`madcats merge-pairs`).
- Calls diploid and polyploid genotypes from the MADC (`madcats call`, experimental) and writes
  VCFs of the designed targets, of all tag variants, or of whole microhaplotypes.
- Optional genome-based paralog filter (`madcats filter-paralogs`, experimental, needs
  `pip install 'madcats[genome]'`) for panels where co-amplified paralogous copies inflate
  heterozygosity. It works on vendor and HapApp fixed-allele-ID MADCs too.
- Runs samples in parallel and caches per-sample results, so reruns are fast.

## Installation

Requires Python 3.9 or newer.

```bash
pip install madcats
madcats --help
```

Each [release](https://github.com/Breeding-Insight/MADCATS/releases) also has standalone executables
for Linux and macOS that need no Python. See
[Installation](https://madcats.readthedocs.io/en/latest/getting-started/installation.html) for
those, for `pipx`, for running from a clone, and for the optional aligner that `filter-paralogs`
needs.

## Quick Start

To try MADCATS on a small simulated project first, see [Try the Example](https://madcats.readthedocs.io/en/latest/getting-started/example.html).

```bash
# 1) Reconstruct each marker's Ref/Alt tags from the DArTag design files
madcats design \
  --primers panel_primer_sequence.csv \
  --probe-design panel_probe-design.txt \
  --output reference_tags.csv

# 2) Generate the MADC count table from one FASTQ per sample
madcats generate \
  --fastq-dir fastqs/ \
  --reference-tags reference_tags.csv \
  --cache-dir work/cache \
  --jobs 8 \
  --output results/project_MADC.csv

# 3) Optional: compare against a vendor-delivered MADC
madcats validate \
  --generated results/project_MADC.csv \
  --expected vendor_MADC.csv \
  --report results/validation.json

# 4) Optional: call genotypes and write a VCF
madcats call \
  --madc results/project_MADC.csv \
  --ploidy 4 \
  --output-prefix results/project
```

The tutorials cover the other routes, which are panels
[without a primer file](https://madcats.readthedocs.io/en/latest/tutorials/primer-free.html),
[AgriSeq and other reads that do not start at a primer](https://madcats.readthedocs.io/en/latest/tutorials/agriseq.html)
and [paired-end reads](https://madcats.readthedocs.io/en/latest/tutorials/paired-end.html).

## Documentation

Full documentation is available at [madcats.readthedocs.io](https://madcats.readthedocs.io), including:

- [Getting Started](https://madcats.readthedocs.io/en/latest/getting-started/installation.html) for installation and first steps
- [Tutorials](https://madcats.readthedocs.io/en/latest/tutorials/dartag.html) with step-by-step examples for DArTag, primer-free, AgriSeq and paired-end data, the paralog filter and genotype calling
- [Command Reference](https://madcats.readthedocs.io/en/latest/commands/overview.html) for every subcommand and option
- [Evidence Model](https://madcats.readthedocs.io/en/latest/reference/evidence-model.html) on how MADCATS counts, discovers and denoises alleles

The [development record](https://github.com/Breeding-Insight/MADCATS/blob/main/docs/DEVELOPMENT.md)
gives the design rationale and version-by-version measurements behind every default.

## Commands

| Command | Description |
|---------|-------------|
| `madcats design` | Reconstruct Ref/Alt tags from vendor primer + probe-design files |
| `madcats infer-design` | Infer a reference-tag panel from the probe design plus FASTQ evidence |
| `madcats merge-pairs` | Merge overlapping paired-end reads into single fragments |
| `madcats normalize-reads` | Orient and window non-primer-started amplicon reads (e.g. AgriSeq) |
| `madcats window-scan` | Sweep candidate tag lengths and recommend one from the project's own reads |
| `madcats generate` | Generate the MADC count table |
| `madcats validate` | Compare a generated MADC against a reference MADC |
| `madcats filter-paralogs` | Drop or flag MADC rows that a reference genome places at another locus (experimental, needs `mappy`) |
| `madcats call` | Call genotypes from a MADC and write target, all-variant or haplotype VCFs (experimental) |

Run `madcats <command> --help` for full options. `filter-paralogs` and `call` are new in 0.14.0 and
experimental, so their options, defaults and outputs may change.

## Input and Output

`generate` reads one FASTQ per sample (optionally gzipped) and a panel definition, which is a
reference-tag CSV from `design`, `infer-design` or `normalize-reads`, or one curated by hand. It
writes the MADC table, a QC JSON with the parameters and per-sample read accounting, and a table
of each row's variants. See
[Input Formats](https://madcats.readthedocs.io/en/latest/getting-started/input-formats.html) and
[Output Formats](https://madcats.readthedocs.io/en/latest/reference/output-formats.html). Breeding
Insight's [MADC file format overview](https://breeding-insight.github.io/BIGr_documentation/posts/madc-file-format/)
describes the table's columns.

## Validation

Each platform was checked against the best independent comparison available for its dataset.
"Genotype agreement" means updog dosage calls made from the MADCATS table (through HapApp and BIGr)
compared with the reference genotypes.

| Platform | Accuracy | Validated against |
|----------|----------|-------------------|
| DArTag (sweetpotato, red clover, blackberry, alfalfa, pecan) | count r = 0.9997, dosage agreement 95.4–98.8% | vendor MADC, vendor/updog dosage, families, skim-seq VCF |
| AgriSeq, Ion Torrent (honeybee, oyster) | 94.9% vs Torrent VariantCaller, 98.0% replicate concordance | vendor caller + replicates |
| GT-seq (delta smelt) | 99.49% | independent BWA/bcftools truth VCF |
| Flex-Seq capture (cranberry) | 95.13% genotype agreement | published population VCF + 20 trios |
| rhAmpSeq (8 grape families) | 99.92–99.99% family transmission | Mendelian constraints + same-reads head-to-head |

See [Validation & Accuracy](https://madcats.readthedocs.io/en/latest/reference/validation.html) for details.

## Development

```bash
python3 -m unittest discover -s tests -v
```

See [Contributing](https://madcats.readthedocs.io/en/latest/development/contributing.html) for the
docs build and release process.

## Citation

If MADCATS is useful in your work, please cite:

> Sandercock, A.M. (2026). MADCATS: Marker-Anchored Discovery and Counting of Alleles from Targeted
> Sequencing. GitHub repository: https://github.com/Breeding-Insight/MADCATS

*(A formal article/DOI is planned, and this section will be updated when available.)*

**AI assistance.** MADCATS was developed in collaboration with Claude Opus 5.5 and Fable 5.1
(Anthropic), which helped write its code and documentation. The authors directed the work, reviewed
and tested the results, and are ultimately responsible for the application.

## License

GPL-3.0-or-later
