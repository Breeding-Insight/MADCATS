# madcats design

Reconstruct each marker's canonical Ref and Alt tags from vendor design files.

## Synopsis

```bash
madcats design --primers <FILE> --probe-design <FILE> --output <FILE> [OPTIONS]
```

## Description

`design` turns a vendor's primer-sequence CSV and bracket-notation probe-design file into a
platform-neutral reference-tag CSV. It uses no sequencing data and gives the same result every
time. The other panel routes were validated against it.

For each marker it:

1. locates the left primer inside the probe-design target, in either orientation;
2. keeps only primer-started tag windows that contain the bracketed variant (a primer occurrence
   in a duplicated flank that never reaches the variant is not a candidate);
3. builds the Ref and Alt tags as the primer sequence verbatim plus the genomic continuation, so
   degenerate (IUPAC) primer bases persist exactly as the vendor reports them;
4. checks the result and stops with an error if the files disagree.

SNPs, insertions, deletions and homopolymer-shifted indels are handled by one rule: the variant
position is where the two full allele expansions first differ.

> **When to use:** whenever you have both design files. `generate` can also read them directly
> (`--primers` + `--probe-design`), but running `design` first gives you a reviewable panel you
> can reuse without the vendor files.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--primers <FILE>` | Primer-sequence CSV |
| `--probe-design <FILE>` | Probe/marker design TSV or CSV |
| `--output <FILE>` | Reference-tag CSV to write |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--tag-length <INT>` | 81 | Tag length to reconstruct (e.g. 54, 81, 109) |

## Examples

### Standard 81-bp DArTag panel

```bash
madcats design \
  --primers SweetPotato_DArTag_primer_sequence.csv \
  --probe-design Sweetpotato_3K_SNPs_DArTag-probe-design.txt \
  --output sweetpotato_reference_tags.csv
```

### 109-bp panel

```bash
madcats design \
  --primers primers.csv \
  --probe-design probes.txt \
  --tag-length 109 \
  --output reference_tags_109.csv
```

## Output Format

A reference-tag CSV with one row per marker; see
[Input Formats](../getting-started/input-formats.md#reference-tags-csv). The number of definitions
written is printed to stderr.

## Errors

The run stops, naming the marker, when:
- the variant lies inside the left primer (Alt reads could never carry the primer)
- both alleles are identical
- a SNP's constructed Alt tag is not among the bracket expansions (the primer table and probe design
  disagree)

## Accuracy

| Panel | Design tags vs vendor canonical rows |
|-------|--------------------------------------|
| Blackberry (3,015 markers, 10 indels) | 2,881 identical; all 134 others differ only by vendor IUPAC codes |
| Pecan 3K | 3,084 identical of 3,100; 16 differ only by vendor IUPAC codes |
| Sweetpotato | All 6,238 81-bp prefixes of the public PanelHub 109-bp design sequences agree |

## See Also

- [madcats infer-design](infer-design.md): the same panel without the primer file
- [madcats generate](generate.md): use the panel
