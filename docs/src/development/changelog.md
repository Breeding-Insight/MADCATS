# Changelog

Notable changes by version. Measurements behind each change are in the
[development record](https://github.com/Breeding-Insight/MADCATS/blob/main/docs/DEVELOPMENT.md) on GitHub.

## Unreleased

- `generate --remove-primers on` writes every row without its primers. The left primer is cut at
  its length. The right primer is cut where it begins in each row, found from the panel's own
  amplicon length and the row's alignment to the Ref tag, with no extra input. Rows that become
  identical are added together. `ClusterConsensusSequence`, `Variants`, `MarkerType` and the
  variants sidecar describe the trimmed rows, with positions counted from the first base after
  the left primer. Reads are still assigned, counted, denoised and classified on the full tags,
  so a row keeps its category. The QC JSON gains `primer_removal`. The default is `off`, and the
  table written with `off` is unchanged.
- `Variants` writes an insertion or deletion on the very first base of a tag without an anchor
  base (`1:.>GG`, `1:TT>.`). Before, it borrowed the base after it, which another event could use
  as well, and the row then could not be rebuilt from its `Variants`. A tag that starts with a
  primer never has such an event, so tables written without `--remove-primers on` are unchanged.

## 0.14.0

- New optional command `filter-paralogs`, in its own module, `madcats_paralogs.py`. It filters a
  MADC (from `generate`, a vendor DArT MADC with or without its preamble, or HapApp's
  fixed-allele-ID output) against a reference genome. It places each marker by its Ref tag
  and labels every haplotype `on_target`, `paralog` or `ambiguous`. Paralog and ambiguous rows
  are dropped (default) or flagged in a `GenomeLabel` column, while canonical rows are always kept.
  It writes a labels sidecar, a markers sidecar and a QC JSON. It needs `mappy`
  (`pip install 'madcats[genome]'`), imported only by this command, so MADCATS still has no
  required dependencies. On the lettuce 3K DArTag panel the filtered vendor MADC agreed with a
  bwa-mem-based prototype for 99.4% of reads.
- `filter-paralogs --clone-positions ignore` keeps CloneID coordinates out of marker placement, for
  a genome other than the one the CloneIDs refer to.
- The MADC readers find a header row whose `AlleleID` is not the first cell, as in HapApp's
  fixed-allele-ID output, which starts with `Code_version`.
- `validate` and `filter-paralogs` warn when a MADC has negative read counts. Some vendor MADCs
  store counts in a signed byte, so counts of 128–255 appear as count − 256, and counts of 256 or
  more wrap silently. Both commands report a negative count with a maximum of exactly 127 as this
  overflow and record it in their QC reports. Neither command corrects the counts, since the
  wrapped ones cannot be recovered.
  The lettuce 3K vendor MADC has 90,913 such cells across 1,108 markers.
- For tiled markers, `generate` credits a read that fits several tiles of one amplicon to each of
  them, instead of dropping it as a tie. Tiles are markers that read the same tag window, with a
  shared or nested left primer and Ref tags that agree except at the markers' variant positions.
  `--tile-groups auto|declared|off` (default `auto`) links `LocusGroup` members and detects
  undeclared tiles. The design route now also reads `LocusGroup` from the probe design.
  Ties across families are still dropped. The QC read accounting counts extra credits as
  `tile_extra_credits_*`, so `reads_assigned_full_length` still counts each read once. On the
  Lettuce 3K DArTag panel (376 samples), 13 markers in 6 amplicons now get, in one full-panel
  run, the same table cells as running each tile in a panel of its own. 11 had zero reads before
  (now 6,562–103,852 each), and one pair with identical tags had split its reads 12,299 / 80,306
  (now about 112,400 each). All other markers' rows are unchanged.
- The tile families are part of the per-sample cache key. Panels without tiles keep their caches.

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
