<img src="images/madcats-logo.png" align="right" width="180" alt="MADCATS logo">

# Welcome to MADCATS!

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI](https://img.shields.io/pypi/v/madcats.svg)](https://pypi.org/project/madcats/)
[![Tests](https://github.com/Breeding-Insight/MADCATS/actions/workflows/tests.yml/badge.svg)](https://github.com/Breeding-Insight/MADCATS/actions/workflows/tests.yml)

**MADCATS** (Marker-Anchored Discovery and Counting of Alleles from Targeted Sequencing) is a
platform-independent command-line tool that turns targeted-sequencing FASTQs into MADC
microhaplotype count tables. It recovers the alleles a panel was designed to detect, discovers the
additional microhaplotypes present in the reads, and writes a raw MADC table. That table goes
through HapApp (or HapApp_utils) before dosage calling and VCF conversion with BIGr and related
tools, or straight to [`madcats call`](commands/call.md) for genotypes and VCFs.

The name also reflects the tool's origins as a FASTQ-to-MADC workflow, where MADC is the
*Missing Allele Discovery Counts* format.

## What can MADCATS do?

| Command | Description |
|---------|-------------|
| [`madcats design`](commands/design.md) | Reconstruct Ref/Alt tags from vendor primer + probe-design files |
| [`madcats infer-design`](commands/infer-design.md) | Infer the panel from the probe design plus reads when the primer file is withheld |
| [`madcats merge-pairs`](commands/merge-pairs.md) | Merge overlapping paired-end reads into single fragments |
| [`madcats normalize-reads`](commands/normalize-reads.md) | Orient and window non-primer-started reads (e.g. AgriSeq) |
| [`madcats window-scan`](commands/window-scan.md) | Choose the tag/window length from the project's own reads |
| [`madcats generate`](commands/generate.md) | Generate the MADC count table |
| [`madcats validate`](commands/validate.md) | Compare a generated MADC against a reference MADC |
| [`madcats filter-paralogs`](commands/filter-paralogs.md) | Drop or flag MADC rows that a reference genome places at another locus |
| [`madcats call`](commands/call.md) | Call genotypes from a MADC and write target, all-variant or haplotype VCFs |

## Key Features

- **Open and deterministic.** Every threshold is on the command line, and MADCATS writes a
  QC/provenance JSON beside every result.
- **Platform-independent.** The same code has been validated on DArTag, AgriSeq (Ion Torrent),
  GT-seq, rhAmpSeq and Flex-Seq capture data.
- **Evidence-based discovery.** MADCATS reports a new microhaplotype only when the reads support it
  (see the [evidence model](reference/evidence-model.md)), and denoising never changes a sample's
  Ref/Alt dosage.
- **Indel-aware.** MADCATS fully supports designed indel markers and recovers undesigned indels
  with a gapped alignment that caps the number of substitutions and gaps.
- **No dependencies.** MADCATS uses only the Python 3 standard library, with parallel workers and a
  reusable per-sample cache. The optional genome paralog filter needs `mappy`.
- **Optional paralog filter.** For panels where co-amplified paralogs inflate heterozygosity,
  [`filter-paralogs`](commands/filter-paralogs.md) drops or flags haplotypes a reference genome
  places at another locus.
- **Genotype calling.** [`call`](commands/call.md) estimates diploid and polyploid genotypes from
  the MADC and writes VCFs of the designed targets, of all tag variants, or of whole
  microhaplotypes.

## A prior MADC is not an input

Generation uses only the FASTQs plus the assay design, either the DArTag primer and probe files, the
probe design alone, or a platform-neutral Ref/Alt tag manifest. You use a vendor MADC, when you have
one, only *after* generation, as an external check with [`madcats validate`](commands/validate.md).

## Quick Example

```bash
# Reconstruct the panel from the vendor design files
madcats design \
  --primers panel_primer_sequence.csv \
  --probe-design panel_probe-design.txt \
  --output reference_tags.csv

# Count alleles in every sample
madcats generate \
  --fastq-dir fastqs/ \
  --reference-tags reference_tags.csv \
  --jobs 8 \
  --output project_MADC.csv
```

## Getting Started

New to MADCATS? Begin with these pages:

1. [Installation](getting-started/installation.md) with pip or a downloaded binary
2. [Quick Start](getting-started/quickstart.md) for your first MADC in two commands
3. [Try the Example](getting-started/example.md) to run a small simulated project and check it against the truth
4. [Choosing a Workflow](getting-started/workflows.md) to find which commands your platform and inputs need

## Scope

MADCATS is an independent implementation of the MADC format and does not aim to reproduce any
vendor's table exactly. [Validation & Accuracy](reference/validation.md) records how closely it
agrees, and against which comparators.

## Getting Help

- Browse the [FAQ](reference/faq.md)
- Report issues on [GitHub](https://github.com/Breeding-Insight/MADCATS/issues)
- Run `madcats <command> --help` for any command's full option list
