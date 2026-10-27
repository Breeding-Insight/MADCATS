# MADCATS example: a simulated DArTag-style project

Eight markers, twelve diploid samples, and the true haplotypes and genotypes behind every read.
The data are simulated by `make_demo_data.py`; rerunning it reproduces these files exactly.

```bash
madcats design --primers primers.csv --probe-design probe_design.tsv \
  --output results/reference_tags.csv

madcats generate --fastq-dir fastq --reference-tags results/reference_tags.csv \
  --sample-lookup sample_names.csv --output results/demo_MADC.csv
```

Then compare `results/demo_MADC.csv` with `truth_haplotypes.csv` and `truth_genotypes.csv`.

The walkthrough, including what was simulated and what to look for in the output, is at
https://madcats.readthedocs.io/en/latest/getting-started/example.html
