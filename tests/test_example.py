"""Run the bundled example (examples/madcats-demo) end to end against its simulated truth."""
import contextlib
import csv
import io
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

import madcats

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "madcats-demo"


def read_csv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


@unittest.skipUnless(EXAMPLE.is_dir(), "example data not present")
class ExampleProjectTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        out = Path(cls.temporary.name)
        cls.tags = out / "reference_tags.csv"
        cls.madc = out / "demo_MADC.csv"
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            madcats.main([
                "design", "--primers", str(EXAMPLE / "primers.csv"),
                "--probe-design", str(EXAMPLE / "probe_design.tsv"), "--output", str(cls.tags),
            ])
            madcats.main([
                "generate", "--fastq-dir", str(EXAMPLE / "fastq"), "--reference-tags", str(cls.tags),
                "--sample-lookup", str(EXAMPLE / "sample_names.csv"), "--jobs", "1", "--output", str(cls.madc),
            ])
        cls.truth = read_csv(EXAMPLE / "truth_haplotypes.csv")
        cls.rows = read_csv(cls.madc)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_design_reconstructs_the_simulated_tags(self):
        tags = {row["MarkerName"]: row for row in read_csv(self.tags)}
        for haplotype in self.truth:
            if haplotype["Haplotype"] in ("Ref", "Alt"):
                column = "RefSequence" if haplotype["Haplotype"] == "Ref" else "AltSequence"
                self.assertEqual(tags[haplotype["MarkerName"]][column], haplotype["TagSequence"])

    def test_rows_are_exactly_the_simulated_haplotypes(self):
        emitted = {(row["CloneID"], row["AlleleSequence"]) for row in self.rows}
        simulated = {(row["MarkerName"], row["TagSequence"]) for row in self.truth}
        self.assertEqual(emitted, simulated)

    def test_undesigned_insertion_is_annotated(self):
        insertion = next(row for row in self.truth if row["Haplotype"] == "Insertion")
        row = next(row for row in self.rows if row["AlleleSequence"] == insertion["TagSequence"])
        self.assertEqual(row["MarkerType"], "indel")
        self.assertRegex(row["Variants"], r"\d+:[ACGT]>[ACGT]{2}")

    def test_counts_follow_the_simulated_genotypes(self):
        samples = [column for column in self.rows[0] if column.startswith("Line_")]
        self.assertEqual(len(samples), 12)
        names = {(row["MarkerName"], row["TagSequence"]): row["Haplotype"] for row in self.truth}
        counts = defaultdict(lambda: defaultdict(int))
        for row in self.rows:
            haplotype = names[(row["CloneID"], row["AlleleSequence"])]
            for sample in samples:
                counts[(sample, row["CloneID"])][haplotype] += int(row[sample])
        for call in read_csv(EXAMPLE / "truth_genotypes.csv"):
            observed = counts[(call["Sample"], call["MarkerName"])]
            total = sum(observed.values())
            carried = call["Genotype"].split("/")
            for haplotype, reads in observed.items():
                share = reads / total
                with self.subTest(sample=call["Sample"], marker=call["MarkerName"], haplotype=haplotype):
                    if haplotype not in carried:
                        self.assertLess(share, 0.05)
                    elif carried.count(haplotype) == 1:
                        self.assertTrue(0.2 < share < 0.8, share)
                    else:
                        self.assertGreater(share, 0.95)


if __name__ == "__main__":
    unittest.main()
