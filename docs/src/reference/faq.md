# FAQ

## General

### Do I need a vendor MADC to run MADCATS?

No. Generation uses only the FASTQs plus the assay design. A vendor MADC can be used afterwards
with `madcats validate`, or to bootstrap a panel when the design files are private
(`scripts/madc_to_reference_tags.py`).

### Does MADCATS reproduce the vendor's MADC exactly?

No, and it does not try to. MADCATS is an independent implementation. Canonical-row counts agree
closely (count r = 0.9997 on sweetpotato), while the novel-row lists differ: many vendor rows have
too little support in the project's own reads to pass any within-project evidence route. See
[Validation & Accuracy](validation.md).

### What dependencies does it need?

None beyond Python 3.9+. The standalone binaries need nothing at all.

## Running

### My run stopped saying reads carry a constant prefix. What do I do?

You have a raw DArT delivery: an inline barcode and the common sequence `GATGGATCCAGTGC` precede
every amplicon. Re-run with `--leading-anchor auto`, or pass the sequence the message names.

### Should I trim adapters first?

No. Reads are located from the middle out, windows stay inside the amplicon, and `merge-pairs` drops
adapter read-through as part of merging. Trimming changed concordance by 0.05 points on the capture
benchmark while costing 2 points of extraction.

### Can I split a large cohort into batches?

Batch `normalize-reads` freely, but run **one** `generate` over every sample. The reproducibility
route counts supporting samples within the run, so splitting a cohort into batches finds fewer new
alleles. Use `--cache-dir` to
make re-runs cheap.

### How many workers should I use?

`--jobs` sets the number of samples processed in parallel (default: up to 4). Each worker processes
one sample at a time, so memory scales with `--jobs`. The finalisation (denoising) stage also runs
across the workers.

### Re-running is slow. Can I reuse earlier work?

Pass `--cache-dir`. Per-sample assignment results are cached and keyed by the read-assignment
settings, so changing evidence, denoising or output options regenerates from cache in a fraction of
the time. Changing an assignment setting (identity, anchor, gapped alignment, partial-read length, tag
length) rebuilds the affected caches.

## Panels

### Which panel route should I use?

Use the design files if you have them. Otherwise, build the panel from a vendor MADC, and use
primer-free inference only when neither is available. See
[Choosing a Workflow](../getting-started/workflows.md).

### Why do some canonical rows have zero counts?

Two common reasons: the locus differs from the design in your germplasm, or the vendor's curated
tag contains an IUPAC ambiguity code. Rows and counts use exact sequence equality, so an
ambiguity-coded canonical never matches a concrete read; the reads appear in RefMatch/AltMatch rows
instead.

### What is `--reuse-tags` for, and when do I need it?

For `normalize-reads` panels only. Window placement is learned from coverage, so it can shift
between runs and make identical biology look like new haplotypes. Pass the first run's window CSV
to every later run of the same panel.

## Discovery

### My high-depth project has thousands of novel rows. Are they real?

Probably not all. Evidence thresholds are absolute read counts, and at several hundred reads per
sample per marker even sequencing errors can pass them. Add `--min-supporting-freq 0.10`.

### How do I get only SNP microhaplotypes?

Filter on the `MarkerType` column: drop `indel` and `both` rows. Designed indel markers stay
intact.

### How do I get the exact vendor header?

`--annotations off` restores the vendor column names and writes no variants sidecar.

## Output

### HapApp won't read my MADC.

Generate with `--dart-preamble`. HapApp expects the vendor's starred preamble rows above the header.

### Where did a sample's reads go?

`OUTPUT.qc.json` accounts for every read per sample: assigned, ambiguous, without the leading
anchor, too short, moved by denoising, dropped with low-evidence sequences, and in emitted rows.
