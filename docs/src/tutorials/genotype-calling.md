# Calling Genotypes from a MADC

A MADC holds read counts. Most analyses need genotypes, the number of copies of each allele a
sample carries. In a diploid the counts alone nearly show the genotype. In a tetraploid or a
hexaploid, neighbouring dosages differ by a small share of the reads, and read counts scatter.
[`madcats call`](../commands/call.md) estimates each sample's genotype from the counts and writes
VCFs. It reads the MADC as it is, from `generate`, `filter-paralogs`, the vendor or HapApp, and
needs no other tool.

The steps below use the [example project](../getting-started/example.md), which is diploid. They
are the same for polyploids.

## 1. Call the Targets

```bash
madcats call --madc results/demo_MADC.csv --output-prefix results/demo
```

```text
note: too few read counts to estimate the read model (96 sample-by-marker cells with reads; 100 needed).
Using background-read rate: 0.002
Using overdispersion: 9
results/demo.targets.vcf: 8 records, 96 of 96 genotypes called (100.0%)
Wrote QC to results/demo.qc.json
```

`call` called every genotype of the eight markers and twelve samples. All 96 equal the simulated
genotypes in `truth_genotypes.csv`.

The note is about the size of the example. On a real table `call` measures from the
counts themselves how far read fractions scatter (the overdispersion) and how often a read shows
an allele it did not come from (the background-read rate). Eight markers are too few for that, so
it uses fixed values and prints them.

## 2. Read the VCF

One record of `results/demo.targets.vcf`, with two of its twelve samples:

```text
#CHROM           POS  REF  ALT  INFO                            FORMAT                       Line_01                              Line_02
Chr01_000412883  34   G    C    TARGET;NS=12;AF=0.5;PMC=0.0001  GT:AD:DP:XD:GQ:MPP:DS:GP:PL  0/0:76,0:76:0:48:1:0:1,0,0:0,51,142  0/1:35,35:70:0:40:0.9999:1:0.0001,0.9999,0.0001:40,0,40
```

`Line_01` has 76 reads of G and none of C, so it is `0/0`. `Line_02` has 35 of each and is `0/1`.

| Field | Meaning |
|-------|---------|
| `GT` | The genotype, with one allele per chromosome copy, such as `0/0/0/1` in a tetraploid |
| `AD`, `DP` | Reads of each allele, and their sum |
| `GQ` | The model's confidence in `GT`, Phred-scaled, so 20 is one chance in 100 of a wrong genotype and 30 one in 1,000 |
| `MPP` | The same confidence as a probability, which is the largest value of `GP` |
| `DS` | The expected number of alternate copies, for analyses that take a dosage |
| `GP` | The probability of each genotype |
| `PMC`, in INFO | The share of the marker's samples expected to have a wrong genotype, which is one minus the mean of `MPP` |

**Coordinates are marker-local.** `CHROM` is the marker and `POS` is a position in its tag, not on
a genome. A `##madcatsCoordinates` line in the header marks the coordinates as marker-local.

`call` writes a genotype as missing (`./.`) when the sample has fewer than 10 reads or its `GQ` is
below 10. `--min-dp` and `--min-gq` change that. The read counts stay in the record either way.

## 3. Give the Ploidy

`call` takes every sample to be diploid unless told otherwise:

```bash
madcats call --madc project_MADC.csv --ploidy 4 --output-prefix results/project
```

When ploidy differs between samples, list the exceptions in a two-column table and pass it with
`--ploidy-file`. Samples it does not list take `--ploidy`:

```text
sample	ploidy
Clone_017	2
Clone_044	6
```

## 4. Choose the Settings for Your Samples

Start with the defaults for diploids, and add `--allele-bias auto` for polyploids:

```bash
madcats call --madc project_MADC.csv --ploidy 4 --allele-bias auto --output-prefix results/project
```

`call` estimates how the read counts scatter from the table itself. With `auto` it also fits an
allelic bias at every marker, and uses it only where the bias gains enough and stays within the
search limits. Markers left without a bias get the same calls as without the option.

Change this when one of these describes your samples:

| Your samples | Change | Why |
|--------------|--------|-----|
| Inbred lines, doubled haploids, selfing species | `--selfing-rate` of `0.95` to `0.99` for inbred lines, or the species' own rate for partly selfing ones. Leave the bias off for inbred lines | The default prior expects the heterozygotes of random mating, and a bias has nothing to correct where there are none |
| Families (F1, half-sibs) | Add `--prior norm` | The default prior assumes random mating. It takes two to five times as long |
| Hexaploids, or fewer than about 50 reads, on a panel known to be biased | `--allele-bias fit` | `auto` keeps a bias only where it gains enough, which takes depth, while `fit` uses every bias |
| Subpopulations at 25 reads or fewer, on a panel without bias | `--prior structure` | Each sample's prior comes from the samples most like it. You cannot combine it with a bias |
| A mix of these, or you do not know | Add `--prior norm` | It fits the spread of dosages at each marker and assumes no mating pattern. It gains in families and partly inbred samples, and costs little in one random-mating population |

Read depth changes how much these settings matter:

- **The prior matters most at low depth**, where a sample's reads leave its dosage uncertain. In
  simulated tetraploid families `--prior norm` gained 2.0 points at 10 reads, 0.5 at 60 and 0.1
  at 200.
- **A bias takes depth to find.** In simulated random-mating tetraploids on a biased panel, `auto`
  kept 39% of what `fit` gained at 10 reads, 73% at 25 and 96% at 60. Where the panel had no
  bias, `fit` lost 0.9 points at 10 reads and 0.2 at 60.

[Which settings to use](../commands/call.md#which-settings-to-use) has the measurements behind
each row.

## 5. Write Other Record Types

`--vcf` takes one or more record types, and each goes to its own file:

```bash
madcats call --madc results/demo_MADC.csv --output-prefix results/demo --vcf target all haplotype
```

| `--vcf` | File | A record is |
|---------|------|-------------|
| `target` | `PREFIX.targets.vcf` | The variant a marker was designed for |
| `all` | `PREFIX.all.vcf` | A target, or any other base substitution seen in the tags |
| `haplotype` | `PREFIX.haplotypes.vcf` | A marker, with whole tag sequences as its alleles |

Records with more than two alleles, which is every haplotype record, carry `MPC` in place of `GQ`
and `MPP`, and have no `PMC`. `MPC` is the lowest probability among the allele dosages called,
and `--min-marginal-posterior` (default 0.9) is its threshold.

## 6. Check the Run

`call` prints what to look at, and `PREFIX.qc.json` keeps it:

- **How many genotypes were called.** A low share at the default thresholds usually means low
  depth.
- **Markers without a target record.** `call` does not guess a target when a marker's Ref and Alt
  tags differ in more than one place. It reports how many markers that left out. Pass the panel's
  reference tags with `--reference-tags` to settle them.
- **The overdispersion estimate.** A warning that the estimate is the lowest value tried means the
  table has no heterozygous samples to measure scatter from, as with inbred lines.
- **The allelic bias**, when one is fitted. The `allele_bias` block counts the records called with
  their bias and those called without it, by reason. A warning that most fitted biases ended at a
  limit of the search means the samples have few heterozygotes, so call them without a bias.

## 7. Filter Markers Afterwards

No genotype threshold removes a marker that does not behave as one locus. On two potato panels
with SNP-array genotypes for every marker:

- **Filter markers on call rate first.** The share of samples that pass the default thresholds
  separated the markers that disagreed with the array from the rest better than depth, `GQ`, the
  fitted bias or allele frequency. Removing markers with a call rate under 0.9 raised agreement
  with the array from 91.3% to 94.1% on one panel and from 93.3% to 95.0% on the other.
- **`PMC` finds fewer of them.** Removing markers with a `PMC` of 0.05 or more raised agreement
  to 92.7% and 94.2%, and every marker it removed had a call rate under 0.9 as well. `PMC` also
  rises with ploidy and falls with depth, so a cut that suits one data set removes good markers
  in another.
- **Raise `--min-gq` to 20 at most.** Beyond that, most of the genotypes lost were right ones.
- **`INFO/BIAS` is a weak filter.** Fitting the bias gives a strongly biased marker the right
  calls. Removing markers with a bias outside 0.5 to 2 took out 62 good markers among 162 on one
  panel and 121 among 179 on the other.

These figures are from runs with the bias fitted.

## How `call` Relates to polyRAD and updog

`call` adapts its model from two R packages. Most of it follows
[polyRAD](https://github.com/lvclark/polyRAD) (Clark et al. 2019), which gives it the likelihood
of the read counts, the Hardy-Weinberg prior with its selfing rate, and the repeated estimate of
the allele frequencies. Two parts follow [updog](https://github.com/dcgerard/updog) (Gerard et al.
2018, Gerard and Ferrão 2020). They are the allelic bias (`--allele-bias`) and the `norm` prior
(`--prior norm`).

`call` is not a port of either package. It estimates one overdispersion and one background-read
rate for the whole table. polyRAD takes a fixed overdispersion and contamination rate, and updog
estimates an overdispersion and a sequencing-error rate at every marker. The search for the bias,
the `auto` rule and the `structure` prior are MADCATS's own.

If you have used one of the packages, these settings come closest to it. The last column is the
share of genotypes that were the same in four simulated tables (one diploid, two tetraploid and
one hexaploid, with 18,000 genotypes each), with the package at its defaults:

| Package and model | Closest settings in `call` | Same genotypes |
|-------------------|----------------------------|----------------|
| polyRAD, `IterateHWE` | `--prior hwe --overdispersion 9 --contamination-rate 0.001 --background-rate 0` | All of them. These are polyRAD's settings. `call`'s defaults differ, because they estimate the overdispersion and the background-read rate |
| updog, `norm` | `--prior norm --allele-bias fit` | 99% or more in the diploid and the tetraploids, and 95% in the hexaploid |
| updog, `hw` | `--prior hwe --allele-bias fit` | 99% or more in the diploid and the tetraploids, and 96% in the hexaploid |
| polyRAD, `IteratePopStruct` | `--prior structure` | Not compared. The purpose is the same and the method is different |

`MPP` and `PMC` are updog's `maxpostprob` and `prop_mis`, and `INFO/BIAS` is updog's bias. They
carry the names BIGr's `updog2vcf` gives them, so BIGr's `filterVCF` filters on all three.

## Limits

- **Allopolyploids.** `call` takes every chromosome copy to be able to carry any allele. A marker
  that amplifies all subgenomes but varies on one is outside that model. A marker that amplifies
  one subgenome is a diploid marker, so call it with `--ploidy 2`.
- **Pedigrees.** There is no prior for a cross with known parents, and the model does not use
  linkage between markers. polyRAD has both.
- **Haplotype records.** Their accuracy has been measured on simulated tables only.

## Next Steps

- [madcats call](../commands/call.md) for every option, the model and its measurements
- [Output Formats](../reference/output-formats.md) for the VCF fields and the QC JSON
- [Validation & Accuracy](../reference/validation.md) for what `call` was checked against
