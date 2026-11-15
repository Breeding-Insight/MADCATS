# madcats call

Call genotypes from a MADC and write them as VCFs.

## Synopsis

```bash
madcats call --madc <FILE> --output-prefix <PREFIX> [--vcf {target,all,haplotype} ...] [OPTIONS]
```

## Description

`call` estimates each sample's genotype from the read counts in a MADC. It needs no other tool.
The MADC can come straight from `generate` or `filter-paralogs`, from a vendor delivery (with or
without the starred preamble), or from HapApp (the fixed-allele-ID layout). Ploidy can differ
between samples.

`--vcf` chooses what a record is. `call` writes each type to its own file, and one run can write
several:

| `--vcf` | File | One record per | Alleles |
|---------|------|----------------|---------|
| `target` (default) | `PREFIX.targets.vcf` | marker | The designed variant. Ref and RefMatch reads count for the reference allele, and Alt and AltMatch reads for the alternate allele. Reads of Other rows do not count |
| `all` | `PREFIX.all.vcf` | variant | The target records, plus one record for every other tag position where a row with reads shows a base substitution |
| `haplotype` | `PREFIX.haplotypes.vcf` | marker | Whole tag sequences (microhaplotypes), first the Ref tag and then every other sequence with reads |

The RefMatch and AltMatch rows are part of a target genotype's evidence. On a potato vendor MADC
of diploid and tetraploid samples, target genotypes from every row had the dosage of a SNP array
for 89.1% of 202,143 genotypes. From the exact Ref and Alt rows alone, 82.1% did, with 2% of the
genotypes left without reads.

Coordinates are marker-local. `CHROM` is the marker (its CloneID) and `POS` is a 1-based
position on its Ref tag, so these files are not placed on a genome. The header declares every
marker as a contig, and a `##madcatsCoordinates` line records that the coordinates are
marker-local.

### Which settings to use

The defaults suit one outbred population or a diversity panel, diploid or autopolyploid. They use
the Hardy-Weinberg prior, estimate the overdispersion and background-read rate from the table, and
fit no allelic bias. Two kinds of samples do better with a change:

| Your samples | Use | What it changed |
|--------------|-----|-----------------|
| Polyploids on an amplicon panel such as DArTag | `--allele-bias auto` | In 85 tetraploid potato clones, the share of dosages matching a SNP array rose from 92.3% to 96.4% |
| Inbred lines, doubled haploids, selfing species | `--selfing-rate 0.95` to `0.99` | In simulated inbred tetraploids, the share of correct dosages rose from 97.4% to 99.6–99.9% |

Two other priors helped under particular conditions and not otherwise:

| Your samples | Consider | Where it helped | Where it did not |
|--------------|----------|-----------------|------------------|
| Families (F1, half-sibs) | `--prior norm` | In simulated tetraploid F1 families at 20 to 60 reads, over five draws, it gained 1.0 point (0.6 to 1.5), and 2.5 with `--allele-bias auto` on a biased panel. With it, `auto` kept no bias that was not there. In two F1 populations at about 10 reads, the error was 0.186 against 0.200, and 0.385 against 0.409 | In alfalfa and potato families at high read depth, it ranged from 0.8 points below `hwe` to 0.2 above. In one population it was up to 0.7 below |
| Subpopulations at 60 reads or fewer, on a panel without bias | `--prior structure` | In simulated tetraploids it gained 0.8 to 5.5 points at 8 reads, 0.5 to 2.5 at 25, and 0.6 at 20 to 60. In potato clones of ten families thinned to 10 and 25 reads, it gained 0.9 and 0.7 points | The same potato clones at 300 reads showed no gain. In one F1 family it matched `hwe`. `call` refuses to combine it with a fitted bias or a selfing rate, and on a biased panel `norm` with `auto` was 1.5 points ahead of it |

Three things to know before changing them:

- **A bias matters only where there are heterozygotes.** It changes only the read fraction that
  the model expects a heterozygous genotype to show. Diploids, with one heterozygous class, moved
  by −0.1 to +0.6 points with it, and tetraploid potato gained 4 to 5. Inbred lines gain nothing.
  `fit` miscalls a quarter of their genotypes, and `auto` leaves them as without a bias and warns
  (see [Allelic bias](#allelic-bias)).
- **The other priors and the bias cost time.** `norm` and `structure` take two to five times as
  long as `hwe`, and a fitted bias seven to eight times.
- **The QC JSON counts the records called with the bias and those at a search limit.**
  `allele_bias.records_called_with_the_bias` was half to two thirds of the records on the potato
  data, 13 in 100 on diploid smelt and 3 in 100 or fewer on polyploids pooled from reads without
  bias. `allele_bias.at_search_limit` was 3 in 100 or fewer on every real data set tested, and
  every record on simulated inbred lines.

The measurements behind this section are in
[Validation & Accuracy](../reference/validation.md#genotype-calling-call).

### The genotype model

`call` adapts the model of [polyRAD](https://github.com/lvclark/polyRAD) (Clark, Lipka and Sacks
2019), with the `norm` prior of [updog](https://github.com/dcgerard/updog) (Gerard and Ferrão
2020) and a prior for structured populations as options. Within a record the model considers each
allele separately, against all the other alleles together:

1. **Likelihood.** In a sample of ploidy *P* an allele can be present in 0 to *P* copies (its
   dosage). The probability of the allele's read count, out of the sample's reads in the record,
   is beta-binomial. The expected read fraction of a dosage is dosage / *P*, changed by the
   allelic bias if one is fitted (`--allele-bias`) and then by the share of reads that show the
   wrong allele (`--background-rate`). `--overdispersion` sets how far read fractions may stray
   from the expected fraction.
2. **Prior.** `--prior uniform` makes every dosage equally probable. `--prior hwe` (the default)
   uses Hardy-Weinberg proportions for the allele's frequency, estimated over all the samples
   together. `--prior norm` gives the dosages a normal shape, with a mean and a spread fitted to
   the samples. `--prior structure` uses Hardy-Weinberg proportions at each sample's own allele
   frequency, taken from the samples most like it.
3. **Posterior.** The posterior is the prior times the likelihood, normalised over the dosages.
4. **Iteration** (`hwe`, `norm` and `structure`). `call` re-estimates the allele frequencies, and
   with `norm` the prior's mean and spread, from the posteriors and repeats steps 1–3 until they
   change by no more than `--tol` or it reaches `--max-iterations`. Each record converges
   separately.

All four priors apply to every VCF type.

**`hwe`** assumes random mating and that every chromosome copy is drawn independently, which holds
for diploids and autopolyploids. It is not the expectation under disomic inheritance in an
allopolyploid, in a biparental family, or when samples of different ploidy come from different
populations.

For inbred lines and partly selfing species, `--selfing-rate` adjusts it. A share of each
generation comes from self-fertilisation, which moves prior probability from the intermediate
dosages to the two homozygous ones (the equilibrium proportions of de Silva et al. 2005, as in
polyRAD). At 1 only the two homozygous dosages have prior probability, so `call` cannot call a
heterozygote. Use a value slightly below 1, such as 0.99, to allow residual heterozygosity. A rate
below 1 applies to samples of even ploidy, while samples of odd ploidy keep random-mating
proportions. A rate of 1 applies to every ploidy, so samples of odd ploidy also get only the two
homozygous dosages.

**`norm`** makes no assumption about mating. Its dosage probabilities are proportional to a normal
density at the dosages 0 to *P*, so it can be narrower than Hardy-Weinberg proportions (one
family) or wider (subpopulations), and when the dosages have two peaks (inbred lines) it flattens
out. It is updog's `norm` genotype distribution, fitted the way updog fits it. With the same
likelihood, `call` and updog reach the same prior and the same genotypes. `call` fits it separately
for each ploidy among the samples, since a shape over the dosages of one ploidy says nothing about
another. A ploidy with fewer than 10 samples uses Hardy-Weinberg proportions instead, and `call`
prints a note and lists those ploidies in the QC JSON. In a haplotype record every allele gets its
own fit.

**`structure`** is for samples that are not one randomly mating population, such as
subpopulations, a cline, or a set of families. It keeps Hardy-Weinberg proportions but gives each
sample its own allele frequency in every record. That frequency is the mean over the samples most
like it (its neighbours), with the record's frequency over all samples counted as one more of
them. A sample's own reads never enter its own prior. `call` finds the neighbours once for the
whole table, before it calls any record:

1. `call` calls the target counts of the markers (at most 1,000, at even spacing) with the `hwe`
   prior, which gives every sample an expected share of alternate copies at each marker.
2. Two samples are alike when those shares differ little over the markers. This ranks pairs of
   samples as a genomic relationship matrix would.
3. `call` chooses the number of neighbours from 0, 5, 10, 20, 40 and 80 by prediction. The
   neighbours found on one half of the markers predict each sample's shares at the other half, and
   `call` uses the number with the smallest error. If none predicts better than the table's own
   frequencies, the number is 0 and the prior is the `hwe` prior. `call` prints a note.

`call` then calls records one at a time as with the other priors, each first with the `hwe` prior
and then with the neighbours. The number of neighbours used and the error of each number tried are
in the QC JSON. The method is MADCATS's own, not a port of polyRAD's.

- It needs at least 20 samples and 100 markers with reads to compare samples on. Below that it is
  the `hwe` prior, with a note.
- `call` refuses to combine it with `--allele-bias fit` or `auto`, or with `--selfing-rate`. These
  combinations have not been tested.
- It takes two to three times as long as `--prior hwe`.

### What the model can and cannot say

With two alleles, the dosage of the alternate allele fixes the genotype, so its posterior is
the genotype posterior. These records carry `GQ`, `MPP`, `GP` and `PL`, and `PMC` in INFO. `call`
computes `PL` from the likelihoods and `GP` from the posteriors, so they differ by the prior.

With more than two alleles, there is one posterior per allele, not one over whole genotypes.
The most probable dosages of the alleles need not add up to the ploidy. When they do not, `call`
adds or removes copies one at a time where that costs the least probability (polyRAD's method).
These records carry `MPC`, the lowest posterior among the called dosages. It is an upper bound on
the probability that the whole genotype is right, not an estimate of it, which is why it is not
called `GQ`. `call` writes every haplotype record this way, and every `all` record with three or
more bases.

The model assumes that every copy of a sample's chromosomes can carry any allele (autopolyploid
inheritance), so a hexaploid's dosage runs from 0 to 6. In an allopolyploid whose marker
amplifies all subgenomes but varies on one, the dosage can only be 0 to 2 of 6, and `call` has no
model for that. On the simulated allohexaploid of the polyRAD paper it got 82% of genotypes right
with `hwe` and 87% with `norm`, where polyRAD, told the three subgenomes, did better. A marker that
amplifies one subgenome only is a diploid marker, and `call` calls it correctly with `--ploidy 2`.

### Missing genotypes

`call` writes a genotype as missing (`./.`, with one dot per chromosome copy) when:

- the sample has fewer than `--min-dp` reads in the record,
- two dosages are equally probable,
- its `GQ` is below `--min-gq` (two-allele records) or its `MPC` is below
  `--min-marginal-posterior` (other records), or
- more than `--max-other-fraction` of its reads were left out of the call (when that option is
  given).

`call` always keeps `AD` and `DP`. A sample with reads but below `--min-dp` keeps its probabilities
too. A sample with no reads has `.` for them, since they would only repeat the prior.

### Finding the target

The target record needs the designed variant's position and alleles. `call` takes them from the
marker's Ref and Alt tags (the Ref and Alt rows of the MADC, or `--reference-tags` when given):

- If the two tags differ by a single substitution, insertion or deletion, that is the target.
  `call` writes a designed indel the VCF way, from the base before it.
- If they differ in more places, the target is the substitution at the reference tags'
  `VariantPositionInTagSequence`. Without reference tags the marker gets no target record.
- A marker without an Alt row gets no target record.

`call` never guesses the target. The QC JSON counts skipped markers by reason, and they still get
haplotype records and off-target records.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--madc <FILE>` | MADC to call from, written by `generate` or `filter-paralogs`, or a vendor DArT MADC, or a HapApp fixed-allele-ID MADC. `call` reads gzipped input directly |
| `--output-prefix <PREFIX>` | Path prefix of the outputs |

## Output Selection

| Option | Default | Description |
|--------|---------|-------------|
| `--vcf {target,all,haplotype}` | `target` | VCF types to write, one file each. Give one or more, for example `--vcf target haplotype` |
| `--qc-output <FILE>` | `PREFIX.qc.json` | QC JSON |

## Samples and Panel

| Option | Default | Description |
|--------|---------|-------------|
| `--ploidy <INT>` | 2 | Ploidy of every sample not listed in `--ploidy-file` |
| `--ploidy-file <FILE>` | none | Two columns, sample name and ploidy, separated by a tab or a comma, with or without a header line. Names are the MADC's sample column names |
| `--reference-tags <CSV>` | none | The panel's reference-tag CSV (from `design` or `infer-design`). Its tags and `VariantPositionInTagSequence` identify each marker's target |
| `--max-indel-length <INT>` | 10 | Widest insertion or deletion that `call` considers when it compares a row with the Ref tag. Use the value given to `generate` |

## Model

| Option | Default | Description |
|--------|---------|-------------|
| `--prior {uniform,hwe,norm,structure}` | `hwe` | Genotype prior (see [the genotype model](#the-genotype-model)) |
| `--selfing-rate <FLOAT>` | 0 | With `--prior hwe`, the share of individuals that come from self-fertilisation, where 0 is random mating and 1 is fully inbred lines |
| `--allele-bias {off,fit,auto}` | `off` | Allelic bias, where a copy of one allele yields more reads than a copy of the other. `fit` estimates it for each record with two alleles. `auto` estimates it likewise and calls with it only the records that pass two rules (see below). Not with `--prior structure` |
| `--overdispersion <auto\|FLOAT>` | `auto` | Beta-binomial overdispersion parameter, as in polyRAD. Lower values allow read fractions to stray further from the expected fraction. `auto` estimates it from the read counts (see below), and a number sets it |
| `--background-rate <auto\|FLOAT>` | `auto` | Probability that a read shows one allele although it came from another. `auto` estimates it from the read counts (see below), and a number sets it |
| `--contamination-rate <FLOAT>` | 0 | polyRAD's contamination, the share of a sample's reads that come from other samples and show an allele at its frequency among them. `--background-rate` covers the same reads and more, so this is off unless set. polyRAD uses 0.001 |
| `--tol <FLOAT>` | 1e-5 | With `--prior hwe`, `norm` or `structure`, stop iterating a record once its allele frequencies change by no more than this (with `norm`, also its prior probabilities, and with `structure`, also the samples' own allele frequencies) |
| `--max-iterations <INT>` | 100 | With `--prior hwe`, `norm` or `structure`, the most iterations a record may take, in each of its two stages with `structure` |

### Overdispersion and background reads

Two numbers describe how read counts stray from what a genotype predicts, and `call` estimates
both from the table unless you set them.

- **Overdispersion** is how far read fractions scatter around the expected fraction. It matters
  most in polyploids, where neighbouring dosages differ by a small read fraction. polyRAD chose
  its default of 9 for RAD-seq. On two tetraploid alfalfa DArTag populations it left 59% of target
  genotypes below the confidence thresholds, against 8–11% with estimated values.
- **The background-read rate** is the share of reads that show an allele they did not come from.
  It covers sequencing and amplification errors, index hopping and contamination together, as they
  remain after the MADC was built. If it is set too low, the model can explain such reads only as
  scatter, and the overdispersion estimate drops with it.

`call` chooses the pair that makes the read counts most probable. It tries background rates from
0.0002 to 0.02 and overdispersions from 2 to 500, on each marker's reference and alternate read
counts with the genotypes summed out, using up to 20,000 sample-by-marker cells from markers at
even spacing through the table. `call` can tell the two apart because they show in different
samples. Stray reads show in samples that do not carry an allele, and scatter in samples that
carry both. The values used, and the log-likelihood of every value tried, are in the QC JSON, and
the VCF header states them.

- The estimate does not use `--prior`. It gives every dosage its own free probability at each
  marker, because the estimate is only as good as the genotype proportions it assumes.
- A table with no heterozygous samples and no stray reads has no scatter to measure. The
  overdispersion estimate then ends at the lowest value tried, and `call` warns and uses 9. Inbred
  lines whose reads do include stray reads of the other allele can be estimated from those.
- An overdispersion of 500, the highest value tried, means the counts scatter no more than random
  sampling of reads would.
- A table with fewer than 100 sample-by-marker cells with reads is too small to estimate from.
  `call` prints a note and uses an overdispersion of 9 and a background rate of 0.002.
- With `--allele-bias fit` or `auto` there is a second estimate, which fits every record's bias.
  Without the bias, biased markers look like scatter and the overdispersion comes out lower (30
  and 40 against 150 on the two alfalfa populations). `fit` calls every record with the second
  estimate. `auto` calls with it the records that keep their bias, and the others with the first.
  The QC JSON has the first under `without_bias` in the `overdispersion` and `background_rate`
  blocks.

To run polyRAD's model as polyRAD runs it, use
`--prior hwe --overdispersion 9 --contamination-rate 0.001 --background-rate 0`.

### Allelic bias

In amplicon data a copy of one allele often yields more reads than a copy of the other, so a
sample with two copies of each allele in four does not show half its reads for each. Called
without regard to this, such a marker gets shifted dosages. `--allele-bias fit` estimates the bias
of every record with two alleles (targets, and off-target sites with two bases) and calls with it.
It calls haplotype records and sites with three or more bases without a bias.

`--allele-bias auto` fits the bias in the same way and then calls with it only the records that
pass two rules. It calls the other records without a bias:

- **The bias must gain enough.** It has to raise the log-likelihood of the read counts, after the
  prior on the bias, by 10 or more over no bias. A bias fitted to a marker that has none follows
  the noise in its counts and costs accuracy, most of all in hexaploids.
- **The bias must not end at a limit of the search** (0.14 or 7.4). Samples with few
  heterozygotes, such as inbred lines, drive the fit to a limit, where homozygous samples pass for
  heterozygous ones. `call` warns when it rejects a bias for this, marks the record in the bias
  table and lists it in the QC JSON. `--allele-bias fit` uses such a bias.

`auto` calls a record without a bias the same way `--allele-bias off` calls it, with the
overdispersion and background rate estimated without a bias. Only the records that keep their bias
use the read model estimated with the bias fitted, which is the one `fit` uses for every record.
The VCF header and the QC JSON give both. Both rules are MADCATS's own, while updog always uses the
bias it fits.

Every two-allele record carries its fitted bias in `INFO/BIAS` with `fit` and with `auto` alike,
so you can filter the records on it afterwards whichever was used. With `auto`, a record whose
genotypes were called without its bias also carries the flag `INFO/BIASUNUSED`.

The bias is updog's (Gerard et al. 2018). With bias *h*, an expected read fraction *f* of the
reference allele becomes *f* / (*h*(1 − *f*) + *f*), so below 1 the reference allele is
over-represented and above 1 under-represented. Like updog, `call` pulls the log of the bias
towards 0 with a normal prior of standard deviation 0.7, starts from updog's five starting values
(0.37 to 2.7) and keeps the fit with the highest log-likelihood after that prior. With the other
settings equal, the fitted bias and the genotypes match updog's.

- `call` finds the bias by calling the record at a series of bias values, 17 or more, so a run with
  a bias takes seven to eight times as long. The estimate of the read model does not grow with the
  table, and giving `--overdispersion` and `--background-rate` as numbers skips it.
- `call` applies it to the sample's own allele copies and adds background reads afterwards.
  updog applies its sequencing error first. In this order, samples that do not carry an allele say
  nothing about its bias, which matters because the background rate is one value for the table.
- The search covers biases from 0.14 to 7.4. The QC JSON counts a record whose fit ends at either
  limit, and the bias table marks it. Such a marker is usually not behaving as one locus.
- The VCF carries the fitted bias as `INFO/BIAS`, and `PREFIX.bias.csv` lists every fitted record
  with what its bias gained over no bias (`Gain`), whether its genotypes were called with it
  (`Used`) and if not why (`Rejected`). The QC JSON counts the records called with their bias and
  those called without it, by reason.

When to use it:

- **Use `auto` for polyploids on an amplicon panel.** It needs no knowledge of the panel. On
  potato it raised agreement with a SNP array from 92.3% to 96.4%, where `fit` reached 96.6%. On
  polyploids pooled from reads without bias it stayed within 0.2 points of `off`, where `fit` lost
  1 point in tetraploids and 7 in hexaploids.
- **Use `fit` for a panel known to be biased, in hexaploids or at fewer than about 50 reads.**
  `auto` keeps a bias only where it gains enough, which takes depth. Of what `fit` gained on
  simulated biased tetraploids it kept 98% at 150 reads, 78% at 25, and 64% at 20 to 60 reads with
  smaller biases. In hexaploids it kept 83% at 150 reads and 19% at 20 to 60.
- **Use neither for diploids.** A bias moved accuracy by −0.1 to +0.6 points.
- **Use neither for inbred lines.** A bias has nothing to correct in them. `fit` miscalled a
  quarter of their genotypes on simulated tables. `auto` rejects their biases, warns, and has the
  genotypes of `off`. Give `--selfing-rate`.
- **Use it with `--selfing-rate` or `--prior norm` for partly inbred samples.** On simulated
  half-inbred tetraploids without bias, `auto` with the `hwe` prior alone kept a bias at 12 markers
  in 100 and lost 0.8 points. With `--selfing-rate 0.5` it lost 0.2, and with `norm` 0.1.
- **Use it with `--prior norm` for families.** On simulated F1 families without bias, `auto` with
  `hwe` kept a bias at 2 markers in 100 and lost 0.7 points. With `norm` it kept none. With bias,
  `norm` was 2.5 points ahead of `hwe` under `auto` and 4.4 under `fit`. On real alfalfa and
  potato families at high depth the two priors were within 0.2 points with a bias fitted.

`auto` is not the default because it takes as long as `fit`, and because on diploid smelt it
stayed 0.1 points under `off`.

## Filters

| Option | Default | Description |
|--------|---------|-------------|
| `--min-dp <INT>` | 10 | Write a genotype as missing when the sample has fewer reads than this in the record |
| `--min-gq <INT>` | 10 | For two-allele records, write a genotype as missing when its `GQ` is below this |
| `--min-marginal-posterior <FLOAT>` | 0.9 | For haplotype records and other records with more than two alleles, write a genotype as missing when its `MPC` is below this |
| `--max-other-fraction <FLOAT>` | none | For `target` and `all` records, write a genotype as missing when more than this share of the sample's reads was left out of the call |

## Performance

| Option | Default | Description |
|--------|---------|-------------|
| `--jobs <INT>` | up to 4 | Worker processes. `call` calls each record independently, and the output does not depend on this |

## Examples

```bash
# Target genotypes of a diploid project
madcats call --madc results/project_MADC.csv --output-prefix results/project

# Target and haplotype VCFs of a mixed-ploidy project
madcats call --madc results/project_MADC.csv --output-prefix results/project \
  --vcf target haplotype --ploidy 4 --ploidy-file ploidy.tsv

# A tetraploid diversity panel at low read depth: each sample's prior follows the samples most like it
madcats call --madc results/panel_MADC.csv --output-prefix results/panel \
  --ploidy 4 --prior structure
```

## Outputs

| File | Contents |
|------|----------|
| `PREFIX.targets.vcf`, `PREFIX.all.vcf`, `PREFIX.haplotypes.vcf` | The VCFs chosen with `--vcf` (VCF 4.3, uncompressed) |
| `PREFIX.haplotypes.alleles.csv` | With `--vcf haplotype`, for each marker and allele number, the MADC AlleleIDs it was counted from, its sequence and its total reads |
| `PREFIX.bias.csv` | With `--allele-bias fit` or `auto`, for each fitted record, the fitted bias of the reference allele, the log-likelihood at that bias (after the prior on the bias), what the bias gained over no bias (`Gain`), whether the genotypes were called with it (`Used`, always with `fit`) and if not why (`Rejected`, either `gain_below_margin` or `at_search_limit`), the starting value the search continued from, the number of bias values tried, and whether the fit converged or ended at a search limit |
| `PREFIX.qc.json` | Parameters, the model in words, the overdispersion and background-rate estimates (with a bias fitted, also the estimates made without one), the number of neighbours chosen for the structure prior, a summary of the fitted biases (how many records were called with their bias, how many without and why, and which were rejected at a search limit), samples per ploidy, and for each VCF the records, genotypes called and missing by reason, reads used and excluded, records that did not converge, and targets skipped by reason |

### VCF fields

| Field | In | Meaning |
|-------|----|---------|
| `GT` | all | Genotype. Allele copies are unphased, so `/` joins them |
| `AD` | all | Reads of each allele |
| `DP` | all | Reads used for the call, equal to the sum of `AD` |
| `XD` | `target`, `all` | Reads left out of `AD` and `DP`. At a target these are Other rows, and elsewhere rows that do not show the position (a deletion or an insertion removed it from the tag) |
| `GQ` | two-allele records | Phred-scaled posterior probability that `GT` is wrong, capped at 99 |
| `MPP` | two-allele records | Posterior probability of the most probable genotype. It is the largest value of `GP`, and what `GQ` gives on the Phred scale. updog's `maxpostprob` |
| `DS` | all | Posterior mean number of copies of each alternate allele |
| `GP` | two-allele records | Posterior probability of each genotype, from 0 to *P* alternate copies |
| `PL` | two-allele records | Phred-scaled likelihood of each genotype, before the prior |
| `MPC` | records with more than two alleles, and every haplotype record | Lowest marginal posterior among the called dosages |
| `INFO/TARGET` | `target`, `all` | The record is a marker's designed variant |
| `INFO/BIAS` | two-allele records, with `--allele-bias fit` or `auto` | Fitted allelic bias of the reference allele, as updog reports it for reference read counts |
| `INFO/BIASUNUSED` | two-allele records, with `--allele-bias auto` | The record's genotypes were called without its fitted bias, which gained less than the margin or ended at a limit of the search |
| `INFO/NS` | all | Samples with a called genotype |
| `INFO/AF` | all | Estimated frequency of each alternate allele |
| `INFO/PMC` | two-allele records | Estimated proportion of samples with a wrong genotype, which is one minus the mean of `MPP` over the samples with reads. updog's `prop_mis` |

In the haplotype VCF, `call` adds together rows with the same sequence, and an alternate sequence
with no reads in any sample is not an allele. It writes an indel haplotype as the tag that was
read, which ends a few bases away from where the Ref tag ends.

In the `all` VCF, the three record kinds come from the same reads, so agreement between them is a
check on consistency, not independent evidence. `call` does not yet write off-target insertions
and deletions as records. A substitution carried by several haplotypes gets one record with their
reads added together, and gets the same calls as a two-allele marker with those reads.

A position where the Ref tag holds an ambiguity code gets no off-target record. DArTag panels
write a second designed variant inside a tag that way, and it has a marker of its own.

Tiled markers (`generate --tile-groups`) share the reads of their amplicon, and `call` calls each
from all of them. Each tile's target is its own variant. Where a tile's tag holds the reference
base at another tile's variant, not an ambiguity code, `call` also writes that variant as an
off-target record of this tile, with the same genotypes as the other tile's target record. The
tiles' haplotype records describe the same sequences. Count such a variant once downstream.

### MPP and PMC

`MPP` and `PMC` are updog's `maxpostprob` and `prop_mis`, under the names BIGr's `updog2vcf` gives
them. BIGr's `filterVCF` filters on these names (`filter.MPP`, `filter.PMC`) and on `INFO/BIAS`.

- `MPP` is the probability the model gives a sample's most probable genotype. `GQ` is the same
  number on the Phred scale, so an `MPP` of 0.9 is a `GQ` of 10, and 0.99 is 20.
- `PMC` is one minus the mean of `MPP` over a record's samples, which is the share of them
  expected to have a wrong genotype.

`call` takes `PMC` over the samples with reads, and `--min-dp` and `--min-gq` do not change it.
updog averages over every sample it is given a depth for, and a sample with a depth of 0 enters at
the largest probability of the prior, which raises the value. `call` writes what updog gives when
you pass those samples to it as missing.

Both say how sure the model is, and the model is only as right as it fits the samples:

- Where it fits, `PMC` averaged within 0.4 points of the share of wrong genotypes. This held for
  simulated random mating in diploids, tetraploids and hexaploids, with `--allele-bias off` or
  `auto`.
- Where it does not, `PMC` is too low, by 2 to 11 points on biased polyploid panels called
  without a bias, and by 9 to 15 on a biased F1 family under `--prior hwe`. updog's own value
  was 7 points too low on that family.
- `PMC` rises with ploidy and falls with depth. Of simulated tetraploid markers at 20 to 60
  reads, one in six had a `PMC` under 0.05, the cut updog's documentation gives as an example,
  and nothing was wrong with the others. In hexaploids at that depth almost none did.
- Against the potato SNP array `PMC` averaged 1.4% where 3.6% of the genotypes disagreed (1.1% and
  3.3% for updog). The reads and the array also differ for reasons outside the model.

[Validation & Accuracy](../reference/validation.md#mpp-and-pmc) has the comparison with updog and
the measurements.

## See Also

- [Calling Genotypes from a MADC](../tutorials/genotype-calling.md) for a tutorial on the example
  project
- [Validation & Accuracy](../reference/validation.md#genotype-calling-call) for what `call` was
  checked against
- [madcats generate](generate.md)
- [madcats filter-paralogs](filter-paralogs.md)
- [Output Formats](../reference/output-formats.md)
