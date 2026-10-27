# Panels Without a Primer File

Some providers do not share the primer file. The only thing MADCATS needs from it is where each tag
starts inside the probe-design target, and that can be recovered from the reads.

## Option A: Bootstrap From a Delivered MADC

If any vendor MADC from the same panel exists, it reconstructs the panel exactly. The canonical
Ref/Alt rows carry the tags, and their first difference gives the variant position.

```bash
python3 scripts/madc_to_reference_tags.py \
  --madc vendor_MADC.csv \
  --output reference_tags.csv
```

Markers that cannot become loci are announced with a `WARNING:` and listed with their reason in
`reference_tags.csv.skipped.tsv`. The most common reason is Ref and Alt rows carrying the same
sequence, which means the variant lies one base past the window or the design is presence/absence.

Prefer this route when a vendor MADC is available. It reproduces the vendor's tags, including
IUPAC curation.

## Option B: Infer From Reads

With only the probe-design file:

```bash
madcats infer-design \
  --probe-design probe-design.txt \
  --fastq-dir fastqs/ \
  --jobs 8 \
  --output inferred_reference_tags.csv
```

### How it works

The tag must be a tag-length window of the bracket-expanded target that contains the variant. So
inference only has to choose among a fixed set of candidate windows (about 4x tag length per locus);
it never assembles new sequences. Because reads are primer-started, exact read prefixes almost
always point clearly to one window. Ref/Alt identity comes from bracket order.

Votes are counted in two tiers, in order:

1. The exact tier counts reads whose prefix exactly matches a candidate window. When a locus has any
   exact votes, they decide.
2. The fuzzy tier runs only for loci with no exact votes. It accepts reads that pass the same 80%
   identity rule as `generate`.

A winning window needs at least `--dominance-ratio` (default 10) times the runner-up's votes. The
`LeftPrimerSequence` written to the CSV is the shortest tag prefix that is unique in the panel,
which may be shorter or longer than the real primer. The reads do not show where the real primer
ends, and assignment does not need it.

### Check the result

The `inferred_reference_tags.csv.qc.json` report records each locus's route, support, margin,
anchor length and flags, plus the loci excluded and why.

Review the fuzzy-tier loci. When the fuzzy tier picks the wrong window, it is usually because a
paralog supplied the reads, and the wrong window then wins by a large margin, so the QC statistics
cannot flag it. Cross-check `fuzzy_anchor` loci against any delivered MADC's
canonical rows or a public panel resource.

If you have a ground-truth panel for comparison:

```bash
python3 scripts/validate_inferred_design.py \
  --inferred inferred_reference_tags.csv \
  --truth design_reference_tags.csv \
  --inference-report inferred_reference_tags.csv.qc.json
```

### Then generate as usual

```bash
madcats generate \
  --fastq-dir fastqs/ \
  --reference-tags inferred_reference_tags.csv \
  --output project_MADC.csv
```

Inference is a one-time operation per panel: freeze the CSV and reuse it for every later order.

## How Accurate Is Inference?

| Project | Result |
|---------|--------|
| Sweetpotato DArTag (blinded) | 3,068/3,120 loci inferred: exact tier 100.0% correct (2,742), fuzzy tier 97.9% (319/326); zero Ref/Alt swaps. `generate` on the inferred panel reproduced 99.1% of rows with 100% classification agreement. |
| Delta smelt GT-seq | 368/368 loci exactly as the physical primers define them |
| Alfalfa DArTag (raw reads) | 2,999/3,000 loci in the exact tier, all identical to the design tags |
| Blackberry DArTag | 2,943/3,015 inferred; 18 wrong windows (0.6%), all at read-starved or short-amplicon loci |

## Raw Deliveries and Tag Length

`infer-design` takes the same `--leading-anchor` option as `generate` for raw DArT deliveries.
`--tag-length` sets the panel's reported tag length (54, 81, 109, ...). Reads must reach it, and the
design flanks must contain a full window; both are validated per locus.
