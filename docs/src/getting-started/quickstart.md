# Quick Start

This page takes a DArTag project from FASTQs to a MADC table. Other platforms add a step or two
before `generate` (see [Choosing a Workflow](workflows.md)), and all of them finish with it.

## What You Need

- One FASTQ file per sample (`.fastq` or `.fastq.gz`), all in one folder
- The assay design. For DArTag, that is the vendor's primer-sequence CSV and probe-design file.

## Step 1: Build the Panel

`madcats design` reconstructs each marker's canonical Ref and Alt tags from the design files. No
reads are involved; the result is a platform-neutral reference-tag CSV you can inspect and reuse.

```bash
madcats design \
  --primers panel_primer_sequence.csv \
  --probe-design panel_probe-design.txt \
  --output reference_tags.csv
```

The default tag length is 81 bp. Pass `--tag-length` for panels reported at another length (for
example 54 or 109).

> **No primer file?** Use [`madcats infer-design`](../commands/infer-design.md) to recover the same
> panel from the probe design plus the reads.

## Step 2: Generate the MADC

```bash
madcats generate \
  --fastq-dir fastqs/ \
  --reference-tags reference_tags.csv \
  --cache-dir work/cache \
  --jobs 8 \
  --output results/project_MADC.csv
```

`generate` can also read the design files directly (`--primers` + `--probe-design` in place of
`--reference-tags`). The defaults are the recommended configuration; no evidence flags are needed.

This writes:

| File | Contents |
|------|----------|
| `results/project_MADC.csv` | The MADC table: one row per allele, one read-count column per sample |
| `results/project_MADC.csv.qc.json` | Parameters, provenance and per-sample read accounting |
| `results/project_MADC.csv.variants.csv` | One line per variant event per row |

`--cache-dir` stores per-sample assignment results, so re-running with different evidence or
denoising settings skips the expensive read pass. As a reference point, an 80-sample, 155M-read
DArTag project processes from scratch in about 4.5 minutes with 8 workers and regenerates from a
warm cache in about a minute.

> **Raw DArT deliveries.** Older DArT orders keep an inline barcode and the common sequence
> `GATGGATCCAGTGC` in front of every amplicon. MADCATS detects this and stops with a message
> naming the prefix; re-run with `--leading-anchor auto` (or pass the sequence explicitly).

## Step 3 (Optional): Compare Against a Vendor MADC

If the provider delivered its own MADC, compare the two row-by-row and cell-by-cell:

```bash
madcats validate \
  --generated results/project_MADC.csv \
  --expected vendor_MADC.csv \
  --report results/validation.json
```

## Process the Whole Cohort at Once

Run **one** `generate` over every sample whenever novel-allele discovery matters. The
reproducibility evidence route counts supporting samples within the run, so splitting a cohort
into batches weakens discovery in proportion to the split. Expensive upstream stages such as
`normalize-reads` can be batched freely; feed all of their outputs to a single `generate`.

## Next: HapApp

The MADC that `generate` writes is a raw table. Run it through HapApp (or HapApp_utils), which
assigns the panel's fixed allele IDs, before using it in BIGr, updog or other downstream tools.
HapApp expects the vendor-style preamble rows above the header:

```bash
madcats generate ... --dart-preamble --sample-lookup sample_names.csv
```

See [Output Formats](../reference/output-formats.md) for details.

## Next Steps

- [Try the Example](example.md): run these steps on a small simulated project and check the output
- [Choosing a Workflow](workflows.md): routes for other platforms and inputs
- [madcats generate](../commands/generate.md): every option, grouped by pipeline stage
- [Evidence Model](../reference/evidence-model.md): when a new sequence is reported as a row
