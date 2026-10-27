# madcats validate

Compare a generated MADC against a reference MADC, row by row and cell by cell.

## Synopsis

```bash
madcats validate --generated <FILE> --expected <FILE> [--report <FILE>]
```

## Description

`validate` measures how a generated table agrees with another MADC, typically the vendor's delivery
for the same FASTQs. It reports row precision and recall by allele category, classification
agreement on shared rows, and count agreement across shared row × sample cells.

The reference MADC is used only for this comparison. It has no effect on `generate`.

Either file may carry the vendor's starred preamble rows; the header is located by its `AlleleID`
cell.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--generated <FILE>` | The MADC produced by `madcats generate` |
| `--expected <FILE>` | The reference MADC to compare against |

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--report <FILE>` | none | Also write the JSON report to this file |

The report is always printed to stdout as JSON.

## Example

```bash
madcats validate \
  --generated results/project_MADC.csv \
  --expected Report-DSp25-10607/DSp25-10607_MADC.csv \
  --report results/validation.json
```

## Report Contents

| Section | Contents |
|---------|----------|
| `sample_overlap` | Samples in each file and in common |
| `allele_rows` | Rows in common, precision, recall, generated-only and expected-only counts, `classification_agreement` on shared rows, category counts, and `by_category` precision/recall |
| `common_row_count_cells` | Over shared rows × shared samples: exact-match fraction, mean absolute error, Pearson r, and count totals with their ratio |
| `nonzero_union_count_cells` | The same agreement statistics over cells nonzero in either file |

## Interpreting the Results

Expect high classification and count agreement but well under 100% row recall against a vendor
table. Vendor row lists include many alleles with too little support in the project's own reads to
pass any within-project evidence route. On the red clover validation project, 36.7% of the vendor's
positive-count rows failed every route *on the vendor's own counts*. Against vendor rows with enough
read support to pass an evidence route, recall was 94.1%.

See [Validation & Accuracy](../reference/validation.md) for results across projects.

## See Also

- [madcats generate](generate.md)
- [Validation & Accuracy](../reference/validation.md)
