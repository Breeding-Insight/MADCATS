"""A row that the Ref tag and the Alt tag both explain through gaps, and how it is written.

``row_decision`` takes the reading with fewer edits, then the one with fewer gap openings
(``_reading_cost``), and the Ref tag's only when the two are equal in that as well.  Before, the Ref
tag's reading stood on every tie: an Alt haplotype with one insertion or deletion beside its target
was read as a Ref haplotype with two, or as a row whose target was deleted.

A row read from the Alt tag is also written from that reading (``_read_from_alt``): its Variants are
its alignment to the Alt tag in the Ref tag's positions.  Before, they were always its alignment to
the Ref tag, which could show the target unchanged or deleted on an AltMatch row.
"""

import csv
import random
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import madcats as mp            # noqa: E402
import madcats_call as call     # noqa: E402

BAND = 10
#          0         1         2         3         4         5         6         7         8
#          012345678901234567890123456789012345678901234567890123456789012345678901234567890
REF = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGATGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTG"
LEFT = 20   # bases of REF that are the left primer


def locus_at(target, alt_base, ref=REF):
    """A marker on ``ref`` with its target at the 0-based ``target``."""
    return mp.Locus(0, "m1", ref[:LEFT], "", ref[target], alt_base, target, ref, ref[:target] + alt_base + ref[target + 1 :])


def readings(locus, row):
    """The gapped reading of ``row`` against the Ref tag and against the Alt tag, as row_decision compares them."""
    bound = min(mp.tag_mismatches(row, locus.ref_tag), mp.tag_mismatches(row, locus.alt_tag))
    return tuple(mp._unverified_gapped_reading(row, tag, BAND, bound) for tag in (locus.ref_tag, locus.alt_tag))


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


def shown_at(events, position):
    """What a row's events put on a 1-based Ref position: the Ref base ("ref"), another base, or None (no base)."""
    for kind, start, removed, added in events:
        if kind == "snp" and start == position:
            return added
        if kind == "del":
            first = start + 1 if added.strip(".") else start     # an anchored deletion keeps its first base
            if first <= position < start + len(removed.strip(".")):
                return None
        if kind in ("clip", "unresolved") and start <= position:
            return None
    return "ref"


def side_of(category, alt_base):
    """What a category says stands on the target, in shown_at's terms; None for a category that names no allele."""
    return {"Ref": "ref", "RefMatch": "ref", "Alt": alt_base, "AltMatch": alt_base}.get(category)


def haplotypes(locus, count, seed):
    """Ref and Alt haplotypes of a marker with one insertion or deletion, some with one more substitution.

    Half of the gaps lie within three bases of the target, where the two tags' readings tie most often;
    inserted bases are copied from the neighbourhood half of the time, which makes repeats."""
    rng = random.Random(seed)
    rows = set()
    while len(rows) < count:
        base = rng.choice([locus.ref_tag, locus.alt_tag])
        near = locus.variant_pos + rng.randint(-3, 4)
        at = near if rng.random() < 0.5 else rng.randrange(LEFT + 2, len(base) - 14)
        size = rng.randint(1, 6)
        if rng.random() < 0.5:
            row = base[:at] + base[at + size :] + "".join(rng.choice("ACGT") for _ in range(size))
        else:
            piece = base[at : at + size] if rng.random() < 0.5 else "".join(rng.choice("ACGT") for _ in range(size))
            row = (base[:at] + piece + base[at:])[: len(base)]
        if rng.random() < 0.3:
            index = rng.randrange(LEFT, len(row))
            row = row[:index] + rng.choice([b for b in "ACGT" if b != row[index]]) + row[index + 1 :]
        if row not in (locus.ref_tag, locus.alt_tag):
            rows.add(row)
    return sorted(rows)


class ReadingTieTests(unittest.TestCase):
    def test_the_module_under_test_is_the_repository_s(self):
        self.assertEqual(Path(mp.__file__).resolve().parent, REPO)

    def test_one_insertion_beside_the_alt_base_is_not_read_as_two_around_a_ref_base(self):
        # The target is base 42, G in Ref and A in Alt, inside TGATG.  An Alt haplotype (TAATG) that
        # gained ATGATT after it reads TAATGATGATT: the Alt tag plus one insertion, or the Ref tag
        # with TAA inserted before the T and ATT after the last G.  Six inserted bases either way.
        locus = locus_at(41, "A")
        self.assertEqual((locus.ref_tag[40:45], locus.alt_tag[40:45]), ("TGATG", "TAATG"))
        row = (locus.alt_tag[:45] + "ATGATT" + locus.alt_tag[45:])[: len(REF)]
        self.assertEqual(row, (REF[:40] + "TAA" + REF[40:45] + "ATT" + REF[45:])[: len(REF)])
        on_ref, on_alt = readings(locus, row)
        self.assertEqual((on_ref.edits, on_ref.gap_events), (6, 2))
        self.assertEqual((on_alt.edits, on_alt.gap_events), (6, 1))
        decision = mp.row_decision(locus, row, BAND)
        self.assertIsNone(decision.unresolved)
        self.assertEqual(decision.reading.gap_events, 1)
        self.assertEqual(decision.reading.ref_to_read[41], 41)   # the target is read where it stands
        self.assertEqual(row[41], locus.alt_base)
        self.assertEqual(mp.classify_sequence(locus, row, BAND), "AltMatch")

    def test_the_same_insertion_in_a_ref_haplotype_stays_with_the_ref_tag(self):
        # The mirror marker: TAATG is the Ref tag and G the Alt base.  The Ref tag's reading has the
        # one insertion now, the Alt tag's the two.
        ref = REF[:41] + "A" + REF[42:]
        locus = locus_at(41, "G", ref)
        row = (ref[:45] + "ATGATT" + ref[45:])[: len(ref)]
        on_ref, on_alt = readings(locus, row)
        self.assertEqual((on_ref.edits, on_ref.gap_events), (6, 1))
        self.assertEqual((on_alt.edits, on_alt.gap_events), (6, 2))
        self.assertEqual(mp.row_decision(locus, row, BAND).reading.gap_events, 1)
        self.assertEqual(mp.classify_sequence(locus, row, BAND), "RefMatch")

    def test_one_deletion_after_the_alt_base_is_not_read_as_two_that_delete_the_target(self):
        # The target is base 31, G in Ref and A in Alt, followed by CACG.  An Alt haplotype without
        # those four bases shows its A where the Ref tag has G C A: against the Ref tag the cheapest
        # reading deletes G and C, keeps the A of CACG, and deletes C and G.  Four deleted bases
        # either way, in one piece against the Alt tag.
        locus = locus_at(30, "A")
        self.assertEqual(REF[30:35], "GCACG")
        row = locus.alt_tag[:31] + locus.alt_tag[35:] + "ACGT"
        on_ref, on_alt = readings(locus, row)
        self.assertEqual((on_ref.edits, on_ref.gap_events), (4, 2))
        self.assertEqual((on_alt.edits, on_alt.gap_events), (4, 1))
        self.assertNotIn(30, on_ref.ref_to_read)                 # the Ref tag's reading deletes the target
        decision = mp.row_decision(locus, row, BAND)
        self.assertEqual(decision.reading.ref_to_read[30], 30)
        self.assertEqual(mp.classify_sequence(locus, row, BAND), "AltMatch")

    def test_readings_equal_in_edits_and_openings_leave_the_ref_tag_s(self):
        # With five bases gone (CACGA) the A that follows them is the Alt base again: the row is the
        # Alt haplotype without CACGA and equally the Ref haplotype without GCACG.  One deletion
        # either way; no reading is better, and the Ref tag's has no base on the target.
        locus = locus_at(30, "A")
        row = locus.alt_tag[:31] + locus.alt_tag[36:] + "ACGTA"
        self.assertEqual(row, REF[:30] + REF[35:] + "ACGTA")
        on_ref, on_alt = readings(locus, row)
        self.assertEqual((on_ref.edits, on_ref.gap_events), (5, 1))
        self.assertEqual((on_alt.edits, on_alt.gap_events), (5, 1))
        self.assertEqual(mp.row_decision(locus, row, BAND).reading.ops, on_ref.ops)
        self.assertEqual(mp.classify_sequence(locus, row, BAND), "Other")

    def test_fewer_edits_come_before_fewer_openings(self):
        def reading(edits, openings):
            return mp.GappedReading(edits, 0, edits, openings, (), 0, 0, 0, {})

        self.assertLess(mp._reading_cost(reading(5, 3)), mp._reading_cost(reading(6, 1)))
        self.assertLess(mp._reading_cost(reading(6, 1)), mp._reading_cost(reading(6, 2)))
        self.assertEqual(mp._reading_cost(reading(6, 2)), mp._reading_cost(reading(6, 2)))

    def test_rescue_off_reads_the_alt_base_at_the_fixed_offset(self):
        # Without gapped readings there is nothing to tie: both rows show the Alt base on the target.
        locus = locus_at(41, "A")
        row = (locus.alt_tag[:45] + "ATGATT" + locus.alt_tag[45:])[: len(REF)]
        self.assertEqual(mp.row_decision(locus, row, 0), mp.RowDecision(None, None))
        self.assertEqual(mp.classify_sequence(locus, row, 0), "AltMatch")


class VariantsFollowTheReadingTests(unittest.TestCase):
    def test_insertion_beside_the_alt_base_is_written_as_one_insertion(self):
        locus = locus_at(41, "A")
        row = (locus.alt_tag[:45] + "ATGATT" + locus.alt_tag[45:])[: len(REF)]
        self.assertTrue(mp._read_from_alt(locus, row, BAND))
        events = mp.variant_events(locus, row, BAND)
        self.assertEqual(events, [("snp", 42, "G", "A"), ("ins", 45, "G", "GATGATT"), ("clip", 76, REF[75:], ".")])
        self.assertEqual(mp.format_variants(events), "42:G>A;45:G>GATGATT;clip=6")
        self.assertEqual(replay(REF, "42:G>A;45:G>GATGATT;clip=6", row), row)
        self.assertEqual(mp.marker_type_label(locus, row, BAND, events), "indel")
        # The row's own alignment to the Ref tag, which was written before: two insertions and no Alt base.
        self.assertEqual(mp.format_variants(mp._describe(REF, row, BAND, True)), "40:T>TTAA;45:G>GATT;clip=6")

    def test_gap_right_after_the_alt_base_is_written_without_an_anchor(self):
        # The base before the gap is the target, and the row shows the Alt base there: the
        # substitution uses that Ref base, so the gap starts on its own first base.
        locus = locus_at(30, "A")
        deleted = locus.alt_tag[:31] + locus.alt_tag[35:] + "ACGT"
        events = mp.variant_events(locus, deleted, BAND)
        self.assertEqual(events, [("snp", 31, "G", "A"), ("del", 32, "CACG", "."), ("overhang", 82, ".", "ACGT")])
        self.assertEqual(mp.format_variants(events), "31:G>A;32:CACG>.;overhang=4")
        inserted = (locus.alt_tag[:31] + "TTT" + locus.alt_tag[31:])[: len(REF)]
        self.assertEqual(mp.classify_sequence(locus, inserted, BAND), "AltMatch")
        self.assertEqual(cell(locus, inserted), "31:G>A;32:.>TTT;clip=3")
        for row in (deleted, inserted):
            text = cell(locus, row)
            self.assertEqual(replay(REF, text, row), row)
            self.assertEqual(mp.parse_variants(text, REF, row), mp.variant_events(locus, row, BAND))
            self.assertEqual(mp.marker_type_label(locus, row, BAND), "indel")

    def test_call_reads_a_deletion_without_an_anchor_from_its_first_base(self):
        locus = locus_at(30, "A")
        row = locus.alt_tag[:31] + locus.alt_tag[35:] + "ACGT"
        reading = call.RowReading(REF, mp.parse_variants("31:G>A;32:CACG>.;overhang=4", REF, row))
        self.assertEqual(reading.base_at(31), "A")
        self.assertEqual([reading.base_at(position) for position in range(32, 36)], [None] * 4)
        self.assertEqual(reading.base_at(36), REF[35])
        # An anchored deletion still shows the base it is written from.
        anchored = call.RowReading(REF, [("del", 30, "AGCAC", "A")])
        self.assertEqual(anchored.base_at(30), "A")
        self.assertEqual([anchored.base_at(position) for position in range(31, 35)], [None] * 4)

    def test_a_row_read_from_the_ref_tag_is_written_as_before(self):
        ref = REF[:41] + "A" + REF[42:]
        locus = locus_at(41, "G", ref)
        row = (ref[:45] + "ATGATT" + ref[45:])[: len(ref)]
        self.assertFalse(mp._read_from_alt(locus, row, BAND))
        self.assertEqual(mp.variant_events(locus, row, BAND), mp._describe(ref, row, BAND, True))
        self.assertEqual(cell(locus, row), "45:G>GATGATT;clip=6")

    def test_marker_designed_as_an_indel_keeps_its_alignment_to_the_ref_tag(self):
        alt = REF[:40] + REF[41:] + "A"
        locus = mp.Locus(0, "m1", REF[:LEFT], "", REF[40], alt[40], 40, REF, alt)
        row = alt[:60] + alt[62:] + "GT"
        self.assertFalse(mp._read_from_alt(locus, row, BAND))
        self.assertEqual(mp.variant_events(locus, row, BAND), mp._describe(REF, row, BAND, True))

    def test_variants_and_category_show_the_same_base_on_the_target(self):
        # Every row is rebuilt from its Variants with no Ref base used twice, and its Variants put
        # on the target what its category says: the Ref base, the Alt base, or no allele.
        for target, alt_base, seed in ((41, "A", 1), (30, "A", 2), (55, "C", 3)):
            locus = locus_at(target, alt_base)
            read_from_alt = 0
            for row in haplotypes(locus, 500, seed):
                category = mp.classify_sequence(locus, row, BAND)
                events = mp.variant_events(locus, row, BAND)
                text = mp.format_variants(events)
                self.assertEqual(replay(REF, text, row), row, text)
                self.assertEqual(mp.parse_variants(text, REF, row), events, text)
                shown = shown_at(events, target + 1)
                if side_of(category, alt_base) is not None:
                    self.assertEqual(shown, side_of(category, alt_base), (category, text))
                else:
                    self.assertNotIn(shown, ("ref", alt_base), (category, text))
                read_from_alt += mp._read_from_alt(locus, row, BAND)
            self.assertGreater(read_from_alt, 50)


class TableTests(unittest.TestCase):
    """The table ``generate`` writes for Alt haplotypes with a gap beside the target, with and without primers."""

    def setUp(self):
        self.loci = [
            mp.Locus(0, "inserted", REF[:LEFT], "", REF[41], "A", 41, REF, REF[:41] + "A" + REF[42:]),
            mp.Locus(1, "deleted", REF[:LEFT], "", REF[30], "A", 30, REF, REF[:30] + "A" + REF[31:]),
        ]
        ins, dele = self.loci
        self.row_inserted = (ins.alt_tag[:45] + "ATGATT" + ins.alt_tag[45:])[: len(REF)]
        self.row_deleted = dele.alt_tag[:31] + dele.alt_tag[35:] + "ACGT"
        self.rows = {
            0: {ins.ref_tag: 40, ins.alt_tag: 30, self.row_inserted: 12},
            1: {dele.ref_tag: 40, dele.alt_tag: 30, self.row_deleted: 9},
        }
        for index, locus in enumerate(self.loci):
            for row in haplotypes(locus, 150, 10 + index):
                self.rows[index].setdefault(row, 6)

    def write(self, directory, **options):
        counts = {(index, tag): reads for index, rows in self.rows.items() for tag, reads in rows.items()}
        results = [mp.SampleResult(mp.Sample("s0", "s0", "unused.fastq", {}), counts, {}, {}, {}, {})]
        output = Path(directory) / "madc.csv"
        mp.write_madc(
            output, self.loci, results, mp.EvidenceRule(1, 20, 0, 0, 0, 5, 0.0),
            mp.OutputOptions(variants_path=output.with_suffix(".variants.csv"), **options), BAND,
        )
        with open(output, newline="") as handle:
            table = {(row["CloneID"], row["AlleleSequence"]): row for row in csv.DictReader(handle)}
        with open(output.with_suffix(".variants.csv"), newline="") as handle:
            sidecar = list(csv.DictReader(handle))
        return table, sidecar

    def check(self, table, shift):
        """Every row is rebuilt from the marker's Ref row and its Variants, which show on the target what its category says."""
        checked = 0
        for locus in self.loci:
            ref_row = next(seq for (marker, seq), row in table.items() if marker == locus.marker and row["AlleleID"].endswith("|Ref"))
            for (marker, sequence), row in table.items():
                if marker != locus.marker:
                    continue
                self.assertEqual(replay(ref_row, row["Variants"], sequence), sequence, row["Variants"])
                category = row["AlleleID"].split("|")[1].rstrip("_0123456789")
                shown = shown_at(mp.parse_variants(row["Variants"], ref_row, sequence), locus.variant_pos - shift + 1)
                if side_of(category, locus.alt_base) is not None:
                    self.assertEqual(shown, side_of(category, locus.alt_base), (category, row["Variants"]))
                else:
                    self.assertNotIn(shown, ("ref", locus.alt_base), (category, row["Variants"]))
                checked += 1
        return checked

    def test_table_with_primers(self):
        with tempfile.TemporaryDirectory() as directory:
            table, sidecar = self.write(directory)
        inserted, deleted = table[("inserted", self.row_inserted)], table[("deleted", self.row_deleted)]
        self.assertEqual((inserted["AlleleID"].split("_")[0], inserted["Variants"], inserted["MarkerType"]), ("inserted|AltMatch", "42:G>A;45:G>GATGATT;clip=6", "indel"))
        self.assertEqual((deleted["AlleleID"].split("_")[0], deleted["Variants"], deleted["MarkerType"]), ("deleted|AltMatch", "31:G>A;32:CACG>.;overhang=4", "indel"))
        gap = [row for row in sidecar if row["AlleleSequence"] == self.row_deleted and row["Type"] == "del"]
        self.assertEqual([(row["Position"], row["Ref"], row["Alt"]) for row in gap], [("32", "CACG", ".")])
        self.assertGreater(self.check(table, 0), 300)

    def test_table_without_primers(self):
        # The same reading in the rows as printed: positions count from the first base after the left primer.
        with tempfile.TemporaryDirectory() as directory:
            table, _ = self.write(directory, remove_primers=True)
        inserted, deleted = table[("inserted", self.row_inserted[LEFT:])], table[("deleted", self.row_deleted[LEFT:])]
        self.assertEqual((inserted["AlleleID"].split("_")[0], inserted["Variants"]), ("inserted|AltMatch", "22:G>A;25:G>GATGATT;clip=6"))
        self.assertEqual((deleted["AlleleID"].split("_")[0], deleted["Variants"]), ("deleted|AltMatch", "11:G>A;12:CACG>.;overhang=4"))
        self.assertGreater(self.check(table, LEFT), 300)


class PrimerCutTests(unittest.TestCase):
    """Rows of a panel whose tags end with the first bases of the right primer, cut at that primer."""

    def test_row_cut_at_its_primer_is_written_from_the_alt_tag_s_reading(self):
        import test_remove_primers as primers
        from unittest import mock

        checked = 0
        for seed in range(40):
            rng = random.Random(seed)
            left, insert, right = primers.dna(rng, primers.LEFT), primers.dna(rng, primers.INSERT), primers.dna(rng, primers.RIGHT)
            target = primers.TARGET
            alt_insert = primers.substituted(insert, target)
            tag = lambda body: (left + body + right + primers.ADAPTER)[: primers.TAG]
            locus = mp.Locus(0, "m", left, right, insert[target], alt_insert[target], primers.LEFT + target, tag(insert), tag(alt_insert))
            # An Alt haplotype without the three bases right after its target.
            body = alt_insert[: target + 1] + alt_insert[target + 4 :]
            row = tag(body)
            category = mp.classify_sequence(locus, row, BAND)
            if category != "AltMatch" or not mp._read_from_alt(locus, row, BAND):
                continue
            shown, rows, _, described = mp.trim_locus_rows(
                locus, [(1, "Ref", locus.ref_tag), (2, "Alt", locus.alt_tag), (3, category, row)],
                {locus.ref_tag: [5], locus.alt_tag: [5], row: [5]}, (primers.LEFT + primers.INSERT, right), BAND, mock.MagicMock(),
            )
            self.assertEqual((shown.ref_tag, shown.alt_tag), (insert, alt_insert))
            self.assertIn(body, described)            # cut at its primer: three bases shorter than the Ref row
            self.assertEqual(described[body][3], alt_insert)
            events = mp._describe(shown.ref_tag, body, BAND, *described[body])
            text = mp.format_variants(events)
            self.assertEqual(replay(insert, text, body), body, text)
            self.assertEqual(shown_at(events, target + 1), locus.alt_base, text)
            self.assertEqual(events[0], ("snp", target + 1, insert[target], alt_insert[target]), text)
            self.assertEqual(events[1], ("del", target + 2, insert[target + 1 : target + 4], "."), text)
            checked += 1
        self.assertGreater(checked, 20)


if __name__ == "__main__":
    unittest.main()
