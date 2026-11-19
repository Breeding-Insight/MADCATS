"""A row that the Ref tag and the Alt tag both explain through gaps, in the same number of edits.

``row_decision`` then takes the reading with fewer gap openings (``_reading_cost``), and the Ref tag's
only when the two are equal in that as well.  Before, the Ref tag's reading stood on every tie: an
Alt haplotype with one insertion or deletion beside its target was read as a Ref haplotype with two,
or as a row whose target was deleted.
"""

import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import madcats as mp            # noqa: E402

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


if __name__ == "__main__":
    unittest.main()
