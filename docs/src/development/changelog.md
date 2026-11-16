# Changelog

Notable changes by version. Measurements behind each change are in the
[development record](https://github.com/Breeding-Insight/MADCATS/blob/main/docs/DEVELOPMENT.md) on GitHub.

## 0.14.0

- New command [`call`](../commands/call.md) (experimental), in its own module, `madcats_call.py`.
  It calls genotypes from a MADC and writes one VCF per type chosen with `--vcf`, which is
  `target` (the default), `all` or `haplotype`. Ploidy can differ between samples. The model is
  adapted from polyRAD, with the allelic bias and the `norm` prior adapted from updog, and it uses
  only the standard library.
- `call` takes its prior from `--prior` (`hwe`, `norm`, `structure` or `uniform`), with
  `--selfing-rate` for inbred and partly selfing samples and `--allele-bias` (`off`, `fit` or
  `auto`) for markers that read one allele more readily. It estimates the overdispersion and a
  background-read rate from the table unless they are given.
- New optional command [`filter-paralogs`](../commands/filter-paralogs.md) (experimental), in its
  own module, `madcats_paralogs.py`. It drops or flags the MADC rows that a reference genome
  places at another locus, and writes a labels table, a markers table and a QC JSON. It needs the
  `mappy` aligner (`pip install 'madcats[genome]'`), so MADCATS still has no required
  dependencies. `--clone-positions ignore` is for a genome other than the one the CloneIDs refer
  to.
- For tiled markers, which read the same tag window, `generate` credits a read that fits several
  tiles of one amplicon to each of them. Before, it dropped the read as a tie. `--tile-groups`
  (`auto`, `declared` or `off`, default `auto`) controls it, and `off` gives the output of 0.13.0.
  The tile families are part of the per-sample cache key, so panels without tiles keep their
  caches.
- The MADC readers accept HapApp's fixed-allele-ID layout, whose header starts with
  `Code_version`, and pass over a `readCountSum` column.
- `validate` and `filter-paralogs` warn when a MADC has negative read counts, which some vendor
  MADCs carry from counts stored in a signed byte, and record it in their QC reports.

## 0.13.0

- Removed the `madc_pipeline` compatibility module. Import `madcats` instead, as the bundled
  scripts already do.
- `--partial-read-mode` takes `on`/`off` like the other switches, and the old value `assign` still
  works as `on`, with a warning.
- Command-line help is clearer. Every command has a description, every option has help text, and
  usage lines name the program `madcats`.
- The docs site no longer publishes the development record (`docs/DEVELOPMENT.md`), which stays in
  the repository.
- A small simulated example project, with its true haplotypes and genotypes, is included in
  `examples/madcats-demo` and attached to each release as `madcats-demo.tar.gz`.
- The test suite checks that the command reference pages match each command's `--help`.

## 0.12.1

- Packaged for PyPI, so `pip install madcats` installs a `madcats` command.
- Each release comes with standalone binaries for Linux and macOS.
- New documentation site. The former README is now the development record, `docs/DEVELOPMENT.md`.
- Licensed under GPL-3.0-or-later.
- `generate`'s finalisation (per-sample denoising and the pooled pass) runs across the `--jobs`
  workers. The selection is identical to the sequential passes. On a 44-run rhAmpSeq family the
  pooled pass alone had taken 1,123 s in one process.

## 0.12.0

- For raw DArT deliveries, `generate` and `infer-design` take `--leading-anchor SEQ|auto`. They
  cut reads where the common sequence ends, and runs on an undeclared constant prefix stop and name
  it.
- On the design route, a reconstructed tag must contain the bracketed variant, so primer occurrences
  in duplicated flanks are no longer candidates.
- The denoiser scales better. At distance 1, it finds parents by enumerating single-substitution
  neighbours once centroids outnumber them, fixing quadratic time at highly paralogous loci.
  Output is identical.

## 0.11.0

- With partial reads (`--partial-read-mode assign`, now `on`, the default), reads shorter than the
  tag count toward an allele that already passed the evidence rules and whose start they match,
  only when this cannot change an allele ratio.
  `--min-partial-length` defaults to tag length − 9. The experimental `extend` mode was removed.
- Classification and annotation share one reading per row (`row_reading`), so they cannot
  disagree.
- For freeze health, `normalize-reads --reuse-tags` scores frozen windows against the new run's
  reads.
- `window-scan --scan-draws` reports recommendation stability across seeded sample subsets.

## 0.10.0

- For undesigned indel recovery, reads that fail substitution-only identity get one gapped
  alignment under the same edit limit (`--gapped-rescue`, `--max-indel-length`), with substitution
  and gap limits that keep paralogs out and a bit-vector prescreen for speed.
- Target base read through the gapped alignment when that reading is cheaper and verified.
- Annotations (`--annotations on`, default) add `Variants` and `MarkerType` columns and the
  `OUTPUT.variants.csv` sidecar.
- The pooled denoising pass merges only into a same-allele parent, and credits or drops N-bearing
  tags, never emitting them.
- `window-scan` gains `--scan-seed` and `--generate-args`.
- Cache format 0.5.0 (caches rebuild once).

## 0.9.0

- Paired-end intake, with the `merge-pairs` command and pooled mates for dual-anchored designs.
- Multi-window tiling, with the `LocusGroup` column, `--group-routing nearest|covered`, and the
  check that tiles excise distinct windows.
- The `window-scan` command.
- `--primer-mismatches` default changed from 0 to 1.

## 0.8.3

- `normalize-reads` assigns anchors per locus, growing past `--anchor-length` only where prefixes
  collide.

## 0.8.2

- Uncomputed statistics columns emit `.` instead of blank cells.
- `--dart-preamble`, `--preamble-project` and `--sample-lookup` for HapApp-compatible output.

## 0.8.1

- `normalize-reads --reuse-tags` freezes window placement across runs.
- Auto placement penalizes tandem repeats (`--str-penalty`).

## 0.8.0

- The `normalize-reads` command for bidirectional, variable-start platforms such as AgriSeq.

## 0.7.0

- Indel markers and degenerate (IUPAC) primers, with one general window rule for SNPs, insertions,
  deletions and homopolymer shifts, and tags built as primer overlay plus genomic suffix.

## 0.5.0

- Evidence model redesign, with reproducibility (primary), depth and optional quality routes,
  combined by OR. Flags renamed so old commands stop with an error:
  - `--min-quality-support` → `--quality-min-clean-reads` (route now off by default)
  - `--novel-min-reads-per-sample` / `--novel-min-samples` → `--reproducible-min-reads` /
    `--reproducible-min-samples`
  - `--depth-min-reads` added
  - `--allele-catalogue` / `--known-only` removed (now downstream HapApp responsibilities)
- Per-sample processing 7.1x faster through memoized assignment decisions.

## Earlier

- `infer-design` for panels without a primer file.
- 0.3.x–0.4.x brought quality-gated discovery and the exploratory audits that motivated 0.5.0.
- The tool was renamed from `madc_pipeline` to MADCATS, with `madc_pipeline.py` kept as a
  compatibility shim until 0.13.0.
