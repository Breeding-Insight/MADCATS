"""Write the MADCATS example dataset: a small simulated DArTag-style project.

Everything in this folder except README.md is produced by this script, so the
data can be regenerated exactly and every planted feature is visible here:

    python make_demo_data.py

The panel has 8 markers with 81-bp tags and the project has 12 diploid
samples. On top of each marker's designed Ref and Alt alleles the reads carry:

- Chr02 and Chr04: an additional microhaplotype (a second SNP on the Ref or
  Alt background) that `generate` should discover;
- Chr06: an undesigned 1-bp insertion haplotype, recovered by gapped alignment;
- Chr08: no sample carries the Alt allele, so its Alt row is still written but
  holds only a stray sequencing-error read;
- Chr03 and Chr07: primers on the reverse strand of the probe-design target;
- sequencing errors (0.2% per base), reads shorter than the tag, and
  off-target reads, so denoising, partial-read crediting and the QC read
  accounting all have something to do.

truth_haplotypes.csv and truth_genotypes.csv record what was simulated.
"""
import csv
import gzip
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEED = 2026
TAG_LENGTH = 81
PRIMER_LENGTH = 20
GENOME_LENGTH = 140
PRIMER_OFFSET = 10          # bases of target sequence before the left primer
AMPLICON_LENGTH = 110       # bases sequenced from the start of the primer
N_SAMPLES = 12
READS_PER_COPY = (25, 45)   # reads per haplotype copy per sample
ERROR_RATE = 0.002
PARTIAL_FRACTION = 0.08
OFF_TARGET_READS = (10, 20)

COMPLEMENT = str.maketrans("ACGT", "TGCA")

# (marker, primer on the reverse strand of the target, variant position in the tag)
MARKERS = [
    ("Chr01_000412883", False, 34),
    ("Chr02_001257730", False, 48),
    ("Chr03_000086214", True, 41),
    ("Chr04_002318006", False, 57),
    ("Chr05_000970341", False, 29),
    ("Chr06_001644492", False, 52),
    ("Chr07_000203517", True, 63),
    ("Chr08_001105968", False, 45),
]


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def other_base(rng, base):
    return rng.choice([candidate for candidate in "ACGT" if candidate != base])


def substitute(sequence, position, base):
    return sequence[:position] + base + sequence[position + 1:]


def build_marker(rng, index, name, reverse, variant_position):
    """Return the design rows and the simulated haplotypes (full amplicons) of one marker."""
    genome = "".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH))
    snp = PRIMER_OFFSET + variant_position - 1
    ref_base = genome[snp]
    alt_base = other_base(rng, ref_base)
    left_primer = genome[PRIMER_OFFSET:PRIMER_OFFSET + PRIMER_LENGTH]
    right_primer = reverse_complement(genome[-PRIMER_OFFSET - PRIMER_LENGTH:-PRIMER_OFFSET])

    if reverse:
        target = reverse_complement(genome)
        bracket = GENOME_LENGTH - 1 - snp
        bracketed = f"[{reverse_complement(ref_base)}/{reverse_complement(alt_base)}]"
    else:
        target = genome
        bracket = snp
        bracketed = f"[{ref_base}/{alt_base}]"
    target_sequence = target[:bracket] + bracketed + target[bracket + 1:]

    amplicon = genome[PRIMER_OFFSET:PRIMER_OFFSET + AMPLICON_LENGTH]
    tag_snp = variant_position - 1
    haplotypes = {
        "Ref": (amplicon, "designed Ref allele"),
        "Alt": (substitute(amplicon, tag_snp, alt_base), "designed Alt allele"),
    }
    if index == 1:
        second = 66
        haplotypes["RefMatch"] = (
            substitute(amplicon, second, other_base(rng, amplicon[second])),
            f"Ref allele plus a second SNP at tag position {second + 1}",
        )
    if index == 3:
        second = 24
        alt = haplotypes["Alt"][0]
        haplotypes["AltMatch"] = (
            substitute(alt, second, other_base(rng, alt[second])),
            f"Alt allele plus a second SNP at tag position {second + 1}",
        )
    if index == 5:
        insertion = 30
        haplotypes["Insertion"] = (
            amplicon[:insertion] + other_base(rng, amplicon[insertion]) + amplicon[insertion:-1],
            f"Ref allele with an undesigned 1-bp insertion near tag position {insertion}",
        )

    primer_row = [name, left_primer, right_primer, ref_base, alt_base, variant_position]
    design_row = [name, target_sequence]
    return primer_row, design_row, haplotypes


def genotypes_for(index, haplotypes):
    """Fixed diploid genotypes: every allele is carried by at least three samples."""
    dosage_cycle = [0, 1, 2, 1, 0, 1, 2, 2, 1, 0, 1, 1]
    calls = []
    for sample in range(N_SAMPLES):
        alt_copies = dosage_cycle[(sample + index) % N_SAMPLES]
        if index == 7:
            alt_copies = 0
        pair = ["Alt"] * alt_copies + ["Ref"] * (2 - alt_copies)
        extra = [name for name in haplotypes if name not in ("Ref", "Alt")]
        if extra and sample % 3 == 0:
            pair[sample % 2] = extra[0]
        calls.append(sorted(pair))
    return calls


def add_errors(rng, sequence):
    return "".join(other_base(rng, base) if rng.random() < ERROR_RATE else base for base in sequence)


def sample_reads(rng, markers, sample):
    reads = []
    for (_, haplotypes, calls) in markers:
        for copy in calls[sample]:
            amplicon = haplotypes[copy][0]
            for _ in range(rng.randint(*READS_PER_COPY)):
                read = add_errors(rng, amplicon)
                if rng.random() < PARTIAL_FRACTION:
                    read = read[:rng.randint(TAG_LENGTH - 9, TAG_LENGTH - 1)]
                reads.append(read)
    for _ in range(rng.randint(*OFF_TARGET_READS)):
        reads.append("".join(rng.choice("ACGT") for _ in range(AMPLICON_LENGTH)))
    rng.shuffle(reads)
    return reads


def write_fastq(path, sample_name, reads):
    # mtime=0 and no stored filename keep regenerated files byte-identical
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
        for number, read in enumerate(reads, start=1):
            handle.write(f"@{sample_name}:{number}\n{read}\n+\n{'I' * len(read)}\n".encode())


def main():
    rng = random.Random(SEED)
    markers = []
    primer_rows, design_rows = [], []
    for index, (name, reverse, variant_position) in enumerate(MARKERS):
        primer_row, design_row, haplotypes = build_marker(rng, index, name, reverse, variant_position)
        primer_rows.append(primer_row)
        design_rows.append(design_row)
        markers.append((name, haplotypes, genotypes_for(index, haplotypes)))

    with (HERE / "primers.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["MarkerName", "LeftPrimerSequence", "RightPrimerSequence",
                         "VariantRef", "VariantAlt", "VariantPositionInTagSequence"])
        writer.writerows(primer_rows)
    with (HERE / "probe_design.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["MarkerName", "TargetSequence"])
        writer.writerows(design_rows)

    samples = [f"S{number:02d}" for number in range(1, N_SAMPLES + 1)]
    with (HERE / "sample_names.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["fastq_id", "sample_name"])
        writer.writerows([sample, f"Line_{sample[1:]}"] for sample in samples)

    with (HERE / "truth_haplotypes.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["MarkerName", "Haplotype", "TagSequence", "Description"])
        for name, haplotypes, _ in markers:
            for haplotype, (amplicon, description) in haplotypes.items():
                writer.writerow([name, haplotype, amplicon[:TAG_LENGTH], description])
    with (HERE / "truth_genotypes.csv").open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Sample", "MarkerName", "Genotype"])
        for sample_index, sample in enumerate(samples):
            for name, _, calls in markers:
                writer.writerow([f"Line_{sample[1:]}", name, "/".join(calls[sample_index])])

    fastq_dir = HERE / "fastq"
    fastq_dir.mkdir(exist_ok=True)
    for sample_index, sample in enumerate(samples):
        write_fastq(fastq_dir / f"{sample}.fastq.gz", sample, sample_reads(rng, markers, sample_index))


if __name__ == "__main__":
    main()
