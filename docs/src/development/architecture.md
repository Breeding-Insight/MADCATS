# Architecture

MADCATS is a single module, `madcats.py`, that uses only the Python standard library. There is
nothing to compile and no dependency to pin, and the same file runs as a script (`python3
madcats.py`), as the installed `madcats` command, or frozen into a standalone binary.

## Repository Layout

| Path | Contents |
|------|----------|
| `madcats.py` | The whole tool: library functions and the command line |
| `tests/` | Unit and end-to-end tests (`unittest`) |
| `scripts/` | Diagnostic and conversion scripts (e.g. `madc_to_reference_tags.py`, `validate_inferred_design.py`) |
| `docs/` | This documentation (mdBook), plus `DEVELOPMENT.md`, the long-form development record |
| `examples/madcats-demo/` | Simulated example project and the script that generates it |
| `assets/` | Full-size logo (2028 × 2028 PNG); the docs use a 600 px copy in `docs/src/images/` |
| `packaging/` | PyInstaller entry point for the standalone binaries |
| `pyproject.toml` | Package metadata for PyPI |
| `.github/workflows/` | Tests on every push; PyPI and binary releases on version tags |

## Reading Map of `madcats.py`

The module is laid out in the order a read travels. Each section opens with a banner comment:

| Section | Contents |
|---------|----------|
| 0. Constants, errors and data types | `Locus`, `Sample`, `SampleResult`, MADC column definitions |
| 1. Sequence and file utilities | IUPAC matching, mismatch counts, CSV/FASTQ helpers |
| 2. Alignment engine | Myers bit-vector prescreen, banded DP, `accepted_gapped_reading` / `row_reading` |
| 3. Panel definition | Design-file reconstruction (`build_panel`), reference-tag tables, primer-free inference |
| 4. Platform intake | Pair merging, window location and freezing, tiling and routing, window-scan |
| 5. Per-sample read assignment | `PrimerIndex`, `select_marker`, `process_sample` and its cache |
| 6. Evidence, denoising and partial-read salvage | `EvidenceRule`, `collect_evidence`, `finalize_counts`, `apply_partial_reads` |
| 7. Classification and annotation | `classify_sequence`, `variant_events`, `marker_type_label` |
| 8. MADC input/output and QC | `read_madc`, `write_madc`, `validate_madc` |
| 9. Commands and command line | One `run_*` per subcommand, `build_parser`, `main` |

## The `generate` Pipeline

```
FASTQs ──► 5. assign reads ──► (per-sample cache)
              │  anchor prefix → identity check → gapped alignment
              ▼
           6. partial reads
              ▼
           6. evidence (reproducibility / depth / quality routes)
              ▼
           6. denoise (per-sample skew, then pooled same-allele merge)
              ▼
           7. classify + annotate each row
              ▼
           8. write MADC + QC JSON + variants sidecar
```

Rows hold the exact sequences observed in reads, with no consensus correction, so the output can be
checked against the raw FASTQs.

Decisions that several stages need are made in one place. The same `EvidenceRule` decides
discovery, which sequences denoising protects, and partial-read crediting, and the same
`row_reading` is used for classification and annotation, so those stages always agree.

## Parallelism and Caching

- Per-sample read assignment, pair merging, normalization and inference run across a
  `ProcessPoolExecutor` sized by `--jobs`.
- Since 0.12.1, `generate`'s finalisation also runs across the workers: per-sample denoising per
  sample, and the pooled pass per locus. The result is identical to the sequential passes.
- The per-sample cache (`--cache-dir`) is keyed by the panel and every read-assignment parameter,
  so downstream thresholds never invalidate it. `CACHE_FORMAT_VERSION` is bumped whenever the
  serialized cache or raw-read processing changes.

## Vocabulary

| Term | Meaning |
|------|---------|
| *tag* | The fixed-length sequence a read is measured as |
| *window* | The position of that sequence on a longer read (intake code) |
| *marker* / *locus* | The same panel entry (`Locus` is the record) |
| *canonical* | The design's Ref and Alt tags |
| *category* | Ref, Alt, RefMatch, AltMatch or Other, read at the target base |
| *anchor* | The prefix a read is assigned by: the left primer or a pseudo-primer |
