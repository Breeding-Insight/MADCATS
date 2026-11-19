"""Tests of the unresolved end of a row.

A gap near the end of a tag can explain a row with fewer edits than
substitutions and still have too few bases after it to be believed.  Such a
row has two readings of its end.  ``generate`` takes neither: the row is read
up to the first base the two readings disagree on, its Variants close with
``unresolved=k``, and a target in that stretch is ``Other``.  These tests pin
that one decision where it is used: the category, the Variants and their
replay, the table with and without primers, the pooled denoiser, and ``call``.
"""
import contextlib
import csv
import io
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))   # the repository's modules, never an installed madcats
sys.path.insert(0, str(Path(__file__).resolve().parent))   # the helpers of test_remove_primers

import madcats as mp            # noqa: E402
import madcats_call as call     # noqa: E402

BAND = 10
#          0         1         2         3         4         5         6         7         8
#          012345678901234567890123456789012345678901234567890123456789012345678901234567890
REF = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTG"
LEFT = 20   # bases of REF that are the left primer


def locus_at(target, alt_base=None, index=0, marker="m1"):
    """A marker on REF with its target at the 0-based ``target``."""
    alt_base = alt_base or next(base for base in "ACGT" if base != REF[target])
    return mp.Locus(index, marker, REF[:LEFT], "", REF[target], alt_base, target, REF, REF[:target] + alt_base + REF[target + 1 :])


def cell(locus, row, gap_band=BAND):
    return mp.format_variants(mp.variant_events(locus, row, gap_band))


def replay(ref, text, row):
    """Rebuild a row from the Ref row and its Variants cell, by the rule the documentation gives.

    ``pos:REF>ALT`` replaces REF at POS by ALT (``.`` is no bases) and may not start on a Ref
    base an earlier token used.  ``clip=k`` and ``unresolved=k`` drop the last k Ref bases;
    ``overhang=k`` then appends the row's last k bases, ``unresolved`` the rest of the row.
    """
    out, cursor, drop, tail = "", 0, 0, None
    for token in (token for token in text.split(";") if token not in ("", ".")):
        if token.startswith("clip="):
            drop = int(token[5:])
        elif token.startswith("overhang="):
            tail = int(token[9:])
        elif token.startswith("unresolved="):
            drop, tail = int(token[11:]), "rest"
        else:
            position, alleles = token.split(":")
            removed, added = (allele.strip(".") for allele in alleles.split(">"))
            start = int(position) - 1
            assert start >= cursor, f"two tokens use Ref base {position}: {text}"
            assert ref[start : start + len(removed)] == removed, f"{token} does not name the Ref bases: {text}"
            out += ref[cursor:start] + added
            cursor = start + len(removed)
    limit = len(ref) - drop
    assert cursor <= limit, f"a token reaches into the bases the last one drops: {text}"
    out += ref[cursor:limit]
    if tail == "rest":
        out += row[len(out) :]
    elif tail:
        out += row[len(row) - tail :]
    return out


# A Ref haplotype that lost base 71 (0-based 70), with one more base of its amplicon at the end.
ONE_BASE_DELETION = REF[:70] + REF[71:] + "C"
# One that lost bases 73 and 74 (0-based 72 and 73).
TWO_BASE_DELETION = REF[:72] + REF[74:] + "GT"
# A verified two-base deletion after base 25, then one inserted base near the end.
DELETION_THEN_LATE_GAP = "GCTAAAGACAATTACATAACATACATCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTCAGGTAAGTGTGT"


class DecisionTests(unittest.TestCase):
    def test_the_modules_under_test_are_the_repository_s(self):
        self.assertEqual(Path(mp.__file__).resolve().parent, REPO)
        self.assertEqual(Path(call.__file__).resolve().parent, REPO)

    def test_target_before_the_unread_end_keeps_its_category(self):
        locus = locus_at(40)
        # Nine bases follow the deletion, three fewer than a gapped reading needs.  The gap costs
        # two edits and substitutions six, so the row's end has two readings and is not read.
        self.assertEqual(mp.tag_mismatches(TWO_BASE_DELETION, REF), 6)
        decision = mp.row_decision(locus, TWO_BASE_DELETION, BAND)
        self.assertIsNone(decision.reading)
        self.assertEqual(decision.unresolved, (72, 72))
        # The target, base 41, is read where both readings put it.
        self.assertEqual(mp.classify_sequence(locus, TWO_BASE_DELETION, BAND), "RefMatch")
        events = mp.variant_events(locus, TWO_BASE_DELETION, BAND)
        self.assertEqual(events, [("unresolved", 73, REF[72:], TWO_BASE_DELETION[72:])])
        self.assertEqual(mp.format_variants(events), "unresolved=9")
        self.assertEqual(mp.marker_type_label(locus, TWO_BASE_DELETION, BAND, events), "unresolved")

    def test_target_in_the_unread_end_is_other(self):
        # The target is base 77.  The deletion moves base 78, a T, onto it, and T is the Alt allele:
        # read at the fixed offset this Ref haplotype is an Alt allele, read through the gap it is
        # a Ref allele.  Ten bases follow the gap, so neither reading is taken.
        locus = locus_at(76, alt_base="T")
        self.assertEqual((REF[76], REF[77]), ("G", "T"))
        self.assertEqual(ONE_BASE_DELETION[76], locus.alt_base)
        self.assertEqual(mp.row_decision(locus, ONE_BASE_DELETION, BAND).unresolved, (71, 71))
        self.assertEqual(mp.classify_sequence(locus, ONE_BASE_DELETION, BAND), "Other")
        self.assertEqual(cell(locus, ONE_BASE_DELETION), "unresolved=10")
        # What each reading alone would say: without gapped readings, the Alt allele ...
        self.assertEqual(mp.classify_sequence(locus, ONE_BASE_DELETION, 0), "AltMatch")
        # ... and with every cheaper gapped reading believed, the Ref allele.
        with mock.patch.object(mp, "GAP_VERIFIED_SUFFIX", 0):
            self.assertEqual(mp.classify_sequence(locus, ONE_BASE_DELETION, BAND), "RefMatch")
            self.assertEqual(cell(locus, ONE_BASE_DELETION), "70:TA>T;overhang=1")

    def test_a_gap_with_enough_bases_after_it_is_read_as_before(self):
        locus = locus_at(76, alt_base="T")
        row = locus.alt_tag[:30] + locus.alt_tag[32:] + "GT"   # an Alt haplotype without bases 31 and 32
        decision = mp.row_decision(locus, row, BAND)
        self.assertIsNotNone(decision.reading)
        self.assertIsNone(decision.unresolved)
        self.assertEqual(mp.classify_sequence(locus, row, BAND), "AltMatch")
        self.assertEqual(cell(locus, row), "30:AGC>A;77:G>T;overhang=2")

    def test_substitutions_near_the_end_stay_substitutions(self):
        locus = locus_at(40)
        row = REF[:78] + "A" + REF[79:]
        self.assertEqual(mp.row_decision(locus, row, BAND), mp.RowDecision(None, None))
        self.assertEqual(cell(locus, row), "79:G>A")
        self.assertEqual(mp.marker_type_label(locus, row, BAND), "snp")

    def test_a_verified_gap_before_the_unread_end_is_kept(self):
        # Against the Alt tag the row is the two-base deletion alone, verified; against the Ref tag
        # a second gap near the end is cheaper.  The two agree up to base 70: the deletion is
        # written, the rest is not read, and the row's 11 last bases stand for the Ref row's 11.
        locus = locus_at(70)
        decision = mp.row_decision(locus, DELETION_THEN_LATE_GAP, BAND)
        self.assertIsNotNone(decision.reading)
        self.assertEqual(decision.unresolved, (70, 68))
        self.assertEqual(mp.classify_sequence(locus, DELETION_THEN_LATE_GAP, BAND), "Other")
        events = mp.variant_events(locus, DELETION_THEN_LATE_GAP, BAND)
        self.assertEqual(
            events, [("del", 25, "ACG", "A"), ("unresolved", 71, REF[70:], DELETION_THEN_LATE_GAP[68:])]
        )
        self.assertEqual(mp.format_variants(events), "25:ACG>A;unresolved=11")
        self.assertEqual(replay(REF, "25:ACG>A;unresolved=11", DELETION_THEN_LATE_GAP), DELETION_THEN_LATE_GAP)
        # A target before the stretch is read through the verified deletion.
        earlier = locus_at(40)
        self.assertEqual(mp.classify_sequence(earlier, DELETION_THEN_LATE_GAP, BAND), "RefMatch")

    def test_canonical_rows_and_rescue_off_have_nothing_unresolved(self):
        locus = locus_at(76, alt_base="T")
        for row in (locus.ref_tag, locus.alt_tag):
            self.assertEqual(mp.row_decision(locus, row, BAND), mp.RowDecision(None, None))
        self.assertEqual(mp.row_decision(locus, ONE_BASE_DELETION, 0), mp.RowDecision(None, None))
        self.assertNotIn("unresolved", cell(locus, ONE_BASE_DELETION, 0))

    def test_marker_designed_as_an_indel_marks_the_unread_end_in_ref_coordinates(self):
        # The Alt tag lacks base 41, so every row is written through gaps against the Ref tag.
        # A designed deletion does not settle a second gap near the row's end: that end is
        # unread on such a marker too, counted on the Ref tag whichever tag the row was read from.
        alt = REF[:40] + REF[41:] + "A"
        locus = mp.Locus(0, "m1", REF[:LEFT], "", REF[40], alt[40], 40, REF, alt)
        self.assertEqual(mp.format_variants(mp._designed_events(REF, alt, BAND)), "39:CT>C;overhang=1")
        # A Ref haplotype without bases 73 and 74.
        self.assertEqual(mp.unresolved_end(locus, TWO_BASE_DELETION, BAND), (72, 72))
        self.assertEqual(cell(locus, TWO_BASE_DELETION), "unresolved=9")
        # The same two bases missing from an Alt haplotype: Alt base 72 is Ref base 73.
        on_alt = alt[:71] + alt[73:] + "GT"
        self.assertEqual(mp.row_decision(locus, on_alt, BAND).unresolved, (71, 71))
        self.assertEqual(mp.unresolved_end(locus, on_alt, BAND), (72, 71))
        self.assertEqual(cell(locus, on_alt), "39:CT>C;unresolved=9")
        self.assertEqual(mp.classify_sequence(locus, on_alt, BAND), "AltMatch")   # the target lies before it
        for row in (TWO_BASE_DELETION, on_alt):
            events = mp.variant_events(locus, row, BAND)
            self.assertEqual(replay(REF, mp.format_variants(events), row), row)
            self.assertEqual(mp.parse_variants(mp.format_variants(events), REF, row), events)
            # ``call`` is shown no base in the unread end.
            reading = call.RowReading(REF, events)
            self.assertIsNone(reading.base_at(77))
            self.assertEqual(reading.base_at(60), REF[59])


class VariantsCellTests(unittest.TestCase):
    def rows(self):
        """Ref and Alt haplotypes with zero or one gap well inside, then an edit in the last 15 bases."""
        rng = random.Random(5)
        locus = locus_at(70)
        rows = set()
        while len(rows) < 600:
            row = rng.choice([locus.ref_tag, locus.alt_tag])
            if rng.random() < 0.5:
                cut, lost = rng.randrange(22, 50), rng.randrange(1, 4)
                row = row[:cut] + row[cut + lost :] + "".join(rng.choice("ACGT") for _ in range(lost))
            for _ in range(rng.randrange(1, 3)):
                at, kind = rng.randrange(66, 81), rng.choice("sdi")
                if kind == "s":
                    row = row[:at] + rng.choice("ACGT") + row[at + 1 :]
                elif kind == "d":
                    row = row[:at] + row[at + 1 :] + rng.choice("ACGT")
                else:
                    row = (row[:at] + rng.choice("ACGT") + row[at:])[:81]
            if row not in (locus.ref_tag, locus.alt_tag):
                rows.add(row)
        return locus, sorted(rows)

    def test_every_row_is_rebuilt_from_its_cell_and_the_cell_parses_back(self):
        locus, rows = self.rows()
        unresolved = after_a_gap = 0
        for row in rows:
            events = mp.variant_events(locus, row, BAND)
            text = mp.format_variants(events)
            self.assertEqual(replay(REF, text, row), row, text)
            self.assertEqual(mp.parse_variants(text, REF, row), events, text)
            if "unresolved=" in text:
                unresolved += 1
                after_a_gap += len(text.split(";")) > 1 and mp.row_decision(locus, row, BAND).reading is not None
                self.assertTrue(text.split(";")[-1].startswith("unresolved="), text)
                self.assertEqual(mp.marker_type_label(locus, row, BAND, events), "unresolved")
        # The rows do exercise the outcome, with and without a verified gap before it.
        self.assertGreater(unresolved, 100)
        self.assertGreater(after_a_gap, 10)

    def test_category_and_cell_agree_about_the_target(self):
        # A target at or after the first Ref base not read is Other; one before it is read.
        _, rows = self.rows()
        for target in (40, 70, 76):
            locus = locus_at(target)
            for row in rows:
                stretch = mp.row_decision(locus, row, BAND).unresolved
                if stretch is not None and target >= stretch[0]:
                    self.assertEqual(mp.classify_sequence(locus, row, BAND), "Other")

    def test_the_token_counts_ref_bases_from_the_end(self):
        events = mp.parse_variants("25:ACG>A;unresolved=11", REF, DELETION_THEN_LATE_GAP)
        self.assertEqual(events[-1], ("unresolved", 71, REF[70:], DELETION_THEN_LATE_GAP[68:]))
        self.assertEqual(
            mp.parse_variants("unresolved=9", REF, TWO_BASE_DELETION), [("unresolved", 73, REF[72:], TWO_BASE_DELETION[72:])]
        )

    def test_a_token_this_version_does_not_know_is_an_error(self):
        for text in ("79:G>A;someday=3", "unresolved=200", "unresolved=x"):
            with self.assertRaises(ValueError, msg=text):
                mp.parse_variants(text, REF, TWO_BASE_DELETION)


class TableTests(unittest.TestCase):
    """The table ``generate`` writes: two markers on REF, one row with an unread end at each."""

    def setUp(self):
        self.before = locus_at(40, index=0, marker="before")                 # target before the unread end
        self.inside = locus_at(76, alt_base="T", index=1, marker="inside")   # target inside it
        self.loci = [self.before, self.inside]
        self.rows = {
            0: {self.before.ref_tag: (40, 0), self.before.alt_tag: (0, 30), TWO_BASE_DELETION: (9, 11)},
            1: {self.inside.ref_tag: (40, 0), self.inside.alt_tag: (0, 30), ONE_BASE_DELETION: (7, 6)},
        }

    def write(self, directory, **options):
        results = []
        for sample in range(2):
            counts = {(index, tag): reads[sample] for index, rows in self.rows.items() for tag, reads in rows.items() if reads[sample]}
            results.append(mp.SampleResult(mp.Sample(f"s{sample}", f"s{sample}", "unused.fastq", {}), counts, {}, {}, {}, {}))
        output = Path(directory) / "madc.csv"
        summary = mp.write_madc(
            output, self.loci, results, mp.EvidenceRule(1, 20, 0, 0, 0, 5, 0.0),
            mp.OutputOptions(variants_path=output.with_suffix(".variants.csv"), **options), BAND,
        )
        with open(output, newline="") as handle:
            table = {(row["CloneID"], row["AlleleSequence"]): row for row in csv.DictReader(handle)}
        sidecar = []
        if output.with_suffix(".variants.csv").exists():
            with open(output.with_suffix(".variants.csv"), newline="") as handle:
                sidecar = [row for row in csv.DictReader(handle) if row["Type"] == "unresolved"]
        return summary, table, sidecar

    TALLY = {"rows": 2, "reads": 33, "rows_target_not_read": 1, "reads_target_not_read": 13}

    def test_table_with_primers(self):
        with tempfile.TemporaryDirectory() as directory:
            summary, table, sidecar = self.write(directory)
        before, inside = table[("before", TWO_BASE_DELETION)], table[("inside", ONE_BASE_DELETION)]
        self.assertEqual((before["AlleleID"], before["Variants"], before["MarkerType"]), ("before|RefMatch", "unresolved=9", "unresolved"))
        self.assertEqual((inside["AlleleID"], inside["Variants"], inside["MarkerType"]), ("inside|Other", "unresolved=10", "unresolved"))
        self.assertEqual((before["s0"], before["s1"], inside["s0"], inside["s1"]), ("9", "11", "7", "6"))
        self.assertEqual(
            [(row["CloneID"], row["Position"], row["Ref"], row["Alt"]) for row in sidecar],
            [("before", "73", REF[72:], TWO_BASE_DELETION[72:]), ("inside", "72", REF[71:], ONE_BASE_DELETION[71:])],
        )
        self.assertEqual(summary["unresolved_ends"], self.TALLY)
        self.assertEqual(summary["variants"]["events"], {"snp": 2, "unresolved": 2})
        self.assertEqual(summary["marker_types"]["unresolved"], 2)
        self.assertEqual(summary["categories"], {"Ref": 2, "Alt": 2, "RefMatch": 1, "Other": 1})

    def test_table_without_annotations_has_the_same_categories_and_tally(self):
        with tempfile.TemporaryDirectory() as directory:
            summary, table, sidecar = self.write(directory, annotations=False)
        self.assertEqual(table[("before", TWO_BASE_DELETION)]["AlleleID"], "before|RefMatch")
        self.assertEqual(table[("inside", ONE_BASE_DELETION)]["AlleleID"], "inside|Other")
        self.assertEqual(summary["unresolved_ends"], self.TALLY)
        self.assertEqual(sidecar, [])

    def test_table_without_primers_carries_the_decision(self):
        # The same decision in the frame of the rows as printed: the categories and the token are
        # those of the full tags, positions count from the first base after the left primer.
        with tempfile.TemporaryDirectory() as directory:
            summary, table, sidecar = self.write(directory, remove_primers=True)
        before, inside = table[("before", TWO_BASE_DELETION[LEFT:])], table[("inside", ONE_BASE_DELETION[LEFT:])]
        self.assertEqual((before["AlleleID"], before["Variants"], before["MarkerType"]), ("before|RefMatch", "unresolved=9", "unresolved"))
        self.assertEqual((inside["AlleleID"], inside["Variants"], inside["MarkerType"]), ("inside|Other", "unresolved=10", "unresolved"))
        self.assertEqual(
            [(row["CloneID"], row["Position"], row["Ref"], row["Alt"]) for row in sidecar],
            [("before", "53", REF[72:], TWO_BASE_DELETION[72:]), ("inside", "52", REF[71:], ONE_BASE_DELETION[71:])],
        )
        for (marker, sequence), row in table.items():
            ref_row = next(shown for (other, shown), entry in table.items() if other == marker and entry["AlleleID"].endswith("|Ref"))
            self.assertEqual(replay(ref_row, row["Variants"], sequence), sequence)
        self.assertEqual(summary["unresolved_ends"], self.TALLY)

    def test_gapped_rescue_off_writes_no_unresolved_row(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "madc.csv"
            results = [mp.SampleResult(mp.Sample("s0", "s0", "unused.fastq", {}), {(1, ONE_BASE_DELETION): 9}, {}, {}, {}, {})]
            summary = mp.write_madc(output, self.loci, results, mp.EvidenceRule(1, 20, 0, 0, 0, 5, 0.0), mp.OutputOptions(), 0)
        self.assertEqual(summary["unresolved_ends"]["rows"], 0)
        self.assertNotIn("unresolved", summary["variants"]["events"])


class PrimerCutTests(unittest.TestCase):
    """Rows of a panel whose tags end with the first bases of the right primer."""

    def trimmed_row(self, target):
        """A marker with its target at this insert base, and a row without insert bases 30 and 31.

        One insert base and ten primer bases follow the gap on the full tag, so the row's end is
        unread from insert base 30.  The row's primer is found two bases early and cut off, which
        leaves the row two bases shorter than the Ref row."""
        import test_remove_primers as primers

        rng = random.Random(11)
        left, insert, right = primers.dna(rng, primers.LEFT), primers.dna(rng, primers.INSERT), primers.dna(rng, primers.RIGHT)
        alt_insert = primers.substituted(insert, target)
        tag = lambda body: (left + body + right + primers.ADAPTER)[: primers.TAG]
        locus = mp.Locus(0, "m", left, right, insert[target], alt_insert[target], primers.LEFT + target, tag(insert), tag(alt_insert))
        row = tag(insert[:29] + insert[31:])
        self.assertEqual(mp.unresolved_end(locus, row, BAND), (49, 49))
        category = mp.classify_sequence(locus, row, BAND)
        shown, _, _, readings = mp.trim_locus_rows(
            locus, [(1, "Ref", locus.ref_tag), (2, "Alt", locus.alt_tag), (3, category, row)],
            {locus.ref_tag: [5], locus.alt_tag: [5], row: [5]}, (primers.LEFT + primers.INSERT, right), BAND, mock.MagicMock(),
        )
        trimmed = insert[:29] + insert[31:]
        self.assertEqual((len(shown.ref_tag), len(trimmed)), (32, 30))
        self.assertIn(trimmed, readings)
        text = mp.format_variants(mp._describe(shown.ref_tag, trimmed, BAND, *readings[trimmed]))
        self.assertEqual(replay(shown.ref_tag, text, trimmed), trimmed)
        return category, readings[trimmed], text

    def test_row_cut_at_its_primer_to_another_length_keeps_its_unread_end(self):
        # The target, insert base 11, is read: the category stands.  The row's other length shows
        # that something is missing, but not which reading of its end was right, so the end stays
        # unread in the row as printed: the last three Ref bases, where one base of the row stands.
        category, reading, text = self.trimmed_row(target=10)
        self.assertEqual(category, "RefMatch")
        self.assertEqual(reading[2], (29, 29))
        self.assertEqual(text, "unresolved=3")

    def test_row_that_is_other_for_an_unread_target_never_shows_a_reading_of_it(self):
        # The target is the last insert base, inside the unread end: the row is Other.  Cut at
        # its primer, the row must not be written through a deletion that puts a base on it.
        category, reading, text = self.trimmed_row(target=31)
        self.assertEqual(category, "Other")
        self.assertEqual(text, "unresolved=3")


class PooledDenoiserTests(unittest.TestCase):
    def test_row_with_an_unread_target_is_not_folded_into_an_allele(self):
        # The row sits one base closer to the Alt tag than to the Ref tag, because the shifted
        # base at the target is the Alt allele.  Read at the fixed offset it was an Alt-allele row
        # and the pooled pass could fold its reads into the Alt row.  Its target is not read now,
        # so it belongs to neither allele and stays where it is.
        locus = locus_at(76, alt_base="T")
        self.assertLess(mp.tag_mismatches(ONE_BASE_DELETION, locus.alt_tag), mp.tag_mismatches(ONE_BASE_DELETION, locus.ref_tag))
        counts = {locus.ref_tag: 10**6, locus.alt_tag: 10**6, ONE_BASE_DELETION: 3}
        distance = mp.tag_mismatches(ONE_BASE_DELETION, locus.alt_tag)
        self.assertEqual(distance, 8)
        self.assertEqual(mp._pooled_moves(locus, counts, (), distance, 2.0, BAND), {})
        # With gapped readings off the same row is an Alt-allele row and is folded.
        self.assertEqual(mp._pooled_moves(locus, counts, (), distance, 2.0, 0), {ONE_BASE_DELETION: locus.alt_tag})


class CallTests(unittest.TestCase):
    """``call`` on a table with an unresolved row: the row gives no reads to sites in its unread end."""

    SAMPLES = ["s1", "s2", "s3", "s4"]

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.dir = Path(self.temporary.name)
        self.locus = locus_at(76, alt_base="T")
        self.site_row = REF[:78] + "A" + REF[79:]         # Ref at the target, G>A at base 79
        self.rows = [
            ("Ref", REF, [30, 12, 0, 20], "."),
            ("Alt", self.locus.alt_tag, [0, 14, 25, 0], "77:G>T"),
            ("RefMatch", self.site_row, [0, 0, 0, 16], "79:G>A"),
            ("Other", ONE_BASE_DELETION, [0, 0, 0, 18], "unresolved=10"),
        ]

    def tearDown(self):
        self.temporary.cleanup()

    def table(self, variants=True, replace=None):
        header = [mp.ANNOTATION_COLUMNS.get(column, column) if variants else column for column in mp.MADC_METADATA_COLUMNS]
        path = self.dir / "madc.csv"
        with open(path, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header + self.SAMPLES)
            for category, sequence, counts, text in self.rows:
                cells = dict.fromkeys(header, ".")
                cells.update(AlleleID=f"m1|{category}", CloneID="m1", AlleleSequence=sequence)
                if variants:
                    cells["Variants"] = (replace or {}).get(sequence, text)
                writer.writerow([cells[column] for column in header] + counts)
        return path

    def records(self, path):
        prefix = self.dir / "calls"
        captured = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(captured):
            mp.main(["call", "--madc", str(path), "--output-prefix", str(prefix), "--jobs", "1", "--vcf", "all"])
        self.stderr = captured.getvalue()
        with open(str(prefix) + ".qc.json") as handle:
            self.qc = json.load(handle)
        records = {}
        with open(str(prefix) + call.VCF_FILE_SUFFIX["all"]) as handle:
            for line in handle:
                if not line.startswith("#"):
                    fields = line.rstrip("\n").split("\t")
                    keys = fields[8].split(":")
                    samples = [dict(zip(keys, column.split(":"))) for column in fields[9:]]
                    records[int(fields[1])] = (fields[3], fields[4], [(sample["AD"], sample["XD"]) for sample in samples])
        return records

    def test_the_cells_this_table_states_are_what_generate_writes(self):
        for category, sequence, _, text in self.rows:
            self.assertEqual(cell(self.locus, sequence), text)
            if category not in ("Ref", "Alt"):
                self.assertEqual(mp.classify_sequence(self.locus, sequence, BAND), category)

    def test_unread_end_gives_no_reads_to_the_target_or_to_sites_in_it(self):
        records = self.records(self.table())
        # The target, base 77: the unresolved row is Other, so its 18 reads are for neither allele.
        self.assertEqual(records[77], ("G", "T", [("30,0", "0"), ("12,14", "0"), ("0,25", "0"), ("36,0", "18")]))
        # The site at base 79 lies in the unread end: the 18 reads show no base there.
        self.assertEqual(records[79], ("G", "A", [("30,0", "0"), ("26,0", "0"), ("25,0", "0"), ("20,16", "18")]))
        self.assertEqual(sorted(records), [77, 79])
        self.assertEqual(self.qc["markers"]["rows_with_reads_whose_variants_were_not_understood"], 0)

    def test_table_without_a_variants_column_is_read_the_same_way(self):
        stated = self.records(self.table())
        self.assertEqual(self.records(self.table(variants=False)), stated)

    def test_fixed_offset_reading_would_have_invented_sites(self):
        # What the unread end replaces: read base by base, the row's shifted end is a run of
        # substitutions, each a site nobody carries.
        events = mp._describe(REF, ONE_BASE_DELETION, BAND, gapped=False)
        self.assertEqual(mp.format_variants(events), "72:A>G;73:G>T;74:T>A;76:A>G;77:G>T;78:T>G;79:G>T;80:T>G;81:G>C")
        self.assertEqual(cell(self.locus, ONE_BASE_DELETION, 0), mp.format_variants(events))

    def test_cell_with_a_token_call_does_not_know_leaves_the_row_unread(self):
        # A later ``generate`` may write a token this version has never seen.  The row is not read
        # from its sequence in the cell's place: its reads show no base anywhere, and it is counted.
        records = self.records(self.table(replace={self.site_row: "79:G>A;someday=3"}))
        self.assertEqual(sorted(records), [77])            # the only row that showed the site at 79 is unread
        self.assertEqual(self.qc["markers"]["rows_with_reads_whose_variants_were_not_understood"], 1)
        self.assertIn("1 rows with reads have a Variants cell this version of call cannot read", self.stderr)
        row = call.MadcRow("m1|RefMatch", "RefMatch", self.site_row, (1,), "79:G>A;someday=3")
        self.assertIs(call.table_events(row, REF), call.NOT_UNDERSTOOD)
        reading = call.read_row(self.locus, row, BAND)
        self.assertFalse(reading.understood)
        self.assertIsNone(reading.base_at(79))
        self.assertIsNone(reading.base_at(10))


if __name__ == "__main__":
    unittest.main()
