"""Tests of read continuation.

A row whose target lies in an unresolved end has two readings and too few
bases to choose.  Reads longer than the tag show what follows it.  ``generate``
keeps the commonest of those bases, pools them over the samples, and decides
such a row again on the longer window.  These tests pin that: the row settled
through its gap, the row settled as substitutions, the rows that stay
unresolved, the guard on the added bases, and a cohort whose reads end at the
tag, where nothing changes.
"""
import contextlib
import csv
import gzip
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))   # the repository's modules, never an installed madcats
sys.path.insert(0, str(Path(__file__).resolve().parent))   # the helpers of test_unresolved

import madcats as mp            # noqa: E402
from test_unresolved import BAND, LEFT, REF, locus_at, replay   # noqa: E402

# The amplicon goes on after the 81-base tag.
AFTER = "CATTGACCTGAAGTCCATCGGTAAGCTTGA"
SEEN = 28   # bases of it a read shows

# The target is the tag's last base, G; the Alt allele is C.
LOCUS = locus_at(80, alt_base="C")
# The tag of a Ref haplotype that lost base 79: its last three bases are those of REF moved up by
# one.  It is also the tag of an Alt haplotype with two substitutions before the target.  The tag
# cannot tell the two apart; the bases after it can.
LATE = REF[:78] + REF[79:] + AFTER[0]
LOST_A_BASE = AFTER[1 : SEEN + 1]    # what follows the tag of the haplotype that lost a base
SUBSTITUTED = AFTER[:SEEN]           # what follows the tag of the one with substitutions


# A Ref haplotype with three bases inserted after base 76.  They push the last three Ref bases,
# the target among them, out of the tag.  On the tag alone the cheapest alignment is another one:
# two substitutions and a deletion that the haplotype does not have.
INSERTION = (REF[:76] + "ACG" + REF[76:])[:81]
AFTER_INSERTION = (REF[78:] + AFTER)[:SEEN]


def continuations(after_row, after_ref=AFTER[:SEEN], after_alt=AFTER[:SEEN], locus=LOCUS, row=LATE):
    """What row_decision is handed for one marker: {tag: the bases its reads carry after it}."""
    found = {locus.ref_tag: after_ref, locus.alt_tag: after_alt, row: after_row}
    return {tag: bases for tag, bases in found.items() if bases}


def cell(locus, row, found):
    return mp.format_variants(mp.variant_events(locus, row, BAND, found))


def substituted(bases, count):
    """``bases`` with ``count`` bases changed, spread out."""
    changed = list(bases)
    for index in (5, 11, 17, 23)[:count]:
        changed[index] = "A" if changed[index] != "A" else "C"
    return "".join(changed)


class DecisionTests(unittest.TestCase):
    def test_the_module_under_test_is_the_repository_s(self):
        self.assertEqual(Path(mp.__file__).resolve().parent, REPO)

    def test_the_tag_alone_leaves_the_target_unread(self):
        self.assertEqual(mp.row_decision(LOCUS, LATE, BAND), mp.RowDecision(None, (78, 78)))
        self.assertEqual(mp.classify_sequence(LOCUS, LATE, BAND), "Other")
        self.assertEqual(cell(LOCUS, LATE, None), "unresolved=3")

    def test_a_row_settled_through_its_gap(self):
        found = continuations(LOST_A_BASE)
        reading, unresolved = mp.row_decision(LOCUS, LATE, BAND, found)
        self.assertIsNone(unresolved)
        self.assertEqual((reading.gap_bases, reading.end_read, reading.end_ref), (1, 80, 81))
        # Read through the deletion, the target shows the Ref base.
        self.assertEqual(mp.classify_sequence(LOCUS, LATE, BAND, found), "RefMatch")
        events = mp.variant_events(LOCUS, LATE, BAND, found)
        text = mp.format_variants(events)
        self.assertEqual(text, "78:TG>T;overhang=1")
        self.assertEqual(mp.marker_type_label(LOCUS, LATE, BAND, events), "indel")
        self.assertEqual(replay(REF, text, LATE), LATE)
        self.assertEqual(mp.parse_variants(text, REF, LATE), events)

    def test_a_row_settled_as_substitutions(self):
        found = continuations(SUBSTITUTED)
        self.assertEqual(mp.row_decision(LOCUS, LATE, BAND, found), mp.RowDecision(None, None))
        # Read at the fixed offsets, the target shows the Alt base.
        self.assertEqual(mp.classify_sequence(LOCUS, LATE, BAND, found), "AltMatch")
        events = mp.variant_events(LOCUS, LATE, BAND, found)
        self.assertEqual(mp.format_variants(events), "79:G>T;80:T>G;81:G>C")
        self.assertEqual(mp.marker_type_label(LOCUS, LATE, BAND, events), "snp")

    def test_a_row_stays_unresolved_without_enough_bases_after_the_tag(self):
        alone = mp.row_decision(LOCUS, LATE, BAND)
        short = mp.GAP_VERIFIED_SUFFIX - 1
        cases = {
            "no continuation of the row": continuations(""),
            "no continuation of the Ref tag": continuations(LOST_A_BASE, after_ref=""),
            "the row's is too short": continuations(LOST_A_BASE[:short]),
            "the Ref tag's is too short": continuations(LOST_A_BASE, after_ref=AFTER[:short]),
            "none at all": {},
        }
        for name, found in cases.items():
            with self.subTest(name):
                self.assertEqual(mp.row_decision(LOCUS, LATE, BAND, found), alone)
                self.assertEqual(mp.classify_sequence(LOCUS, LATE, BAND, found), "Other")
                self.assertEqual(cell(LOCUS, LATE, found), "unresolved=3")

    def test_the_shortest_window_that_settles_adds_twelve_bases(self):
        enough = mp.GAP_VERIFIED_SUFFIX
        found = continuations(LOST_A_BASE[:enough], after_ref=AFTER[:enough], after_alt=AFTER[:enough])
        self.assertEqual(mp.classify_sequence(LOCUS, LATE, BAND, found), "RefMatch")

    def test_the_alt_tag_borrows_the_ref_tag_s_continuation(self):
        found = continuations(LOST_A_BASE, after_alt="")
        self.assertEqual(mp.classify_sequence(LOCUS, LATE, BAND, found), "RefMatch")

    def test_the_added_bases_must_be_read_cleanly(self):
        # Two substitutions in the bases after the tag are sequence variation; three are a row
        # that does not follow this amplicon, and the longer window then settles nothing.
        allowed = mp.CONTINUATION_MAX_MISMATCHES
        for after_row, settled in ((LOST_A_BASE, "RefMatch"), (SUBSTITUTED, "AltMatch")):
            with self.subTest(settled):
                found = continuations(substituted(after_row, allowed))
                self.assertEqual(mp.classify_sequence(LOCUS, LATE, BAND, found), settled)
                found = continuations(substituted(after_row, allowed + 1))
                self.assertEqual(mp.row_decision(LOCUS, LATE, BAND, found), mp.row_decision(LOCUS, LATE, BAND))
                self.assertEqual(cell(LOCUS, LATE, found), "unresolved=3")

    def test_a_target_pushed_past_the_tag_stays_other(self):
        # Three inserted bases push the last three Ref bases out of the tag.  The longer window
        # reads the insertion and finds the target after the tag.  The row does not hold it.
        row = (REF[:74] + "ACG" + REF[74:])[:81]
        after_row = (REF[78:] + AFTER)[:SEEN]
        pushed_out = locus_at(79)
        found = continuations(after_row, locus=pushed_out, row=row)
        self.assertEqual(mp.row_decision(pushed_out, row, BAND).unresolved, (75, 75))
        self.assertIsNone(mp.row_decision(pushed_out, row, BAND, found).unresolved)
        self.assertEqual(mp.classify_sequence(pushed_out, row, BAND, found), "Other")
        self.assertEqual(cell(pushed_out, row, found), "74:T>TACG;clip=3")
        self.assertEqual(replay(REF, "74:T>TACG;clip=3", row), row)
        # A target the insertion left in the tag is read from the same row.
        kept_in = locus_at(76, alt_base="T")
        found = continuations(after_row, locus=kept_in, row=row)
        self.assertEqual(mp.classify_sequence(kept_in, row, BAND), "Other")
        self.assertEqual(mp.classify_sequence(kept_in, row, BAND, found), "RefMatch")

    def test_a_settled_row_is_written_from_the_longer_window_s_alignment(self):
        found = continuations(AFTER_INSERTION, row=INSERTION)
        self.assertEqual(cell(LOCUS, INSERTION, None), "77:G>A;78:T>C;unresolved=2")
        # The alignment the tag alone would be given, were its gap believed:
        self.assertEqual(mp.format_variants(mp._describe(REF, INSERTION, BAND, True)), "77:G>A;78:T>C;79:GT>G;overhang=1")
        # The longer window reads the insertion.  The target is past the tag, so the row stays Other.
        self.assertEqual(cell(LOCUS, INSERTION, found), "76:A>AACG;clip=3")
        self.assertEqual(replay(REF, "76:A>AACG;clip=3", INSERTION), INSERTION)
        self.assertEqual(mp.classify_sequence(LOCUS, INSERTION, BAND, found), "Other")

    def test_a_row_whose_target_the_tag_reads_is_never_asked(self):
        before = locus_at(40)
        found = continuations(LOST_A_BASE, locus=before)
        self.assertEqual(mp.row_decision(before, LATE, BAND, found), mp.row_decision(before, LATE, BAND))
        self.assertEqual(mp.classify_sequence(before, LATE, BAND, found), "RefMatch")
        self.assertEqual(cell(before, LATE, found), "unresolved=3")

    def test_canonical_rows_and_rescue_off_are_never_asked(self):
        found = continuations(LOST_A_BASE)
        for row in (LOCUS.ref_tag, LOCUS.alt_tag):
            self.assertEqual(mp.row_decision(LOCUS, row, BAND, found), mp.RowDecision(None, None))
        self.assertEqual(mp.row_decision(LOCUS, LATE, 0, found), mp.RowDecision(None, None))


class PoolingTests(unittest.TestCase):
    def pooled(self, *samples):
        results = [mp.SampleResult(mp.Sample("s", "s", "unused.fastq", {}), {}, {}, {}, {}, {}, continuations=kept) for kept in samples]
        return mp.pooled_continuations(results)

    def test_samples_add_up_and_the_commonest_is_taken(self):
        key = (0, LATE)
        found, length = self.pooled({key: [(LOST_A_BASE, 2), (SUBSTITUTED, 1)]}, {key: [(LOST_A_BASE, 5)]})
        self.assertEqual((found, length), ({0: {LATE: LOST_A_BASE}}, SEEN))

    def test_continuations_are_cut_to_the_length_most_reads_reach(self):
        # Nine reads in ten show 16 bases, one shows 28: all are compared over 16.
        key = (0, LATE)
        found, length = self.pooled({key: [(LOST_A_BASE[:16], 9), (LOST_A_BASE, 1)]})
        self.assertEqual((found, length), ({0: {LATE: LOST_A_BASE[:16]}}, 16))
        # Reads that stop short of that length are left out.
        found, length = self.pooled({key: [(LOST_A_BASE, 9), (SUBSTITUTED[:16], 1)]})
        self.assertEqual((found, length), ({0: {LATE: LOST_A_BASE}}, SEEN))

    def test_a_tag_followed_by_two_sequences_has_no_continuation(self):
        key = (0, LATE)
        self.assertEqual(self.pooled({key: [(LOST_A_BASE, 5), (SUBSTITUTED, 5)]})[0], {})
        self.assertEqual(self.pooled({key: [(LOST_A_BASE, 6), (SUBSTITUTED, 4)]})[0], {0: {LATE: LOST_A_BASE}})

    def test_a_continuation_needs_three_reads(self):
        key = (0, LATE)
        self.assertEqual(self.pooled({key: [(LOST_A_BASE, 2)]})[0], {})
        self.assertEqual(self.pooled({key: [(LOST_A_BASE, 3)]})[0], {0: {LATE: LOST_A_BASE}})

    def test_no_continuation_anywhere(self):
        self.assertEqual(self.pooled({}, {}), ({}, 0))


class GenerateTests(unittest.TestCase):
    """``generate`` from FASTQ files: one marker, two samples, reads of the Ref, Alt and LATE tags."""

    PLAIN = REF[:50] + "A" + REF[51:]   # a Ref haplotype with one substitution: the tag decides it

    def write_reads(self, directory, name, after_tag, after_late):
        """Two samples' reads, each tag followed by ``after_tag`` bases; LATE by those of ``after_late``."""
        reads = Path(directory) / f"reads_{name}"
        reads.mkdir()
        rows = [(LOCUS.ref_tag, AFTER, 40), (LOCUS.alt_tag, AFTER, 30), (self.PLAIN, AFTER, 12), (LATE, after_late, 12)]
        for sample in ("s1", "s2"):
            with gzip.open(reads / f"{sample}.fastq.gz", "wt") as handle:
                number = 0
                for tag, following, copies in rows:
                    sequence = tag + following[:after_tag]
                    for _ in range(copies):
                        number += 1
                        handle.write(f"@{sample}_{number}\n{sequence}\n+\n{'I' * len(sequence)}\n")
        return reads

    def generate(self, directory, reads, cache=None):
        """Run generate on ``reads``: (table path, {sequence: row}, the QC file's output summary)."""
        panel = Path(directory) / "tags.csv"
        mp.write_reference_tags(panel, [LOCUS])
        output = Path(directory) / f"{reads.name}.csv"
        arguments = [
            "generate", "--fastq-dir", str(reads), "--reference-tags", str(panel), "--tag-length", "81",
            "--jobs", "1", "--output", str(output),
        ]
        if cache is not None:
            arguments += ["--cache-dir", str(cache)]
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(mp.main(arguments), 0)
        with open(output, newline="") as handle:
            table = {row["AlleleSequence"]: row for row in csv.DictReader(handle)}
        with open(str(output) + ".qc.json") as handle:
            summary = json.load(handle)["output_summary"]
        return output, table, summary

    def run_generate(self, directory, name, after_tag, after_late):
        return self.generate(directory, self.write_reads(directory, name, after_tag, after_late))

    def test_reads_that_end_at_the_tag_change_nothing(self):
        # Reads that stop at the tag, or two bases after it, give the same table: the row's
        # target stays unread and nothing is kept from after the tag.
        with tempfile.TemporaryDirectory() as directory:
            at_tag, table, summary = self.run_generate(directory, "at_tag", 0, LOST_A_BASE)
            two_more, _, summary_two = self.run_generate(directory, "two_more", 2, LOST_A_BASE)
            self.assertEqual(at_tag.read_bytes(), two_more.read_bytes())
        for found in (summary, summary_two):
            self.assertEqual(found["read_continuation"], {"bases_after_tag": 0, "rows_settled": 0, "reads_settled": 0})
            self.assertEqual(found["unresolved_ends"]["rows_target_not_read"], 1)
        self.assertEqual((table[LATE]["AlleleID"], table[LATE]["Variants"]), ("m1|Other", "unresolved=3"))

    def test_longer_reads_settle_the_row_and_only_that_row(self):
        with tempfile.TemporaryDirectory() as directory:
            _, before, _ = self.run_generate(directory, "at_tag", 0, LOST_A_BASE)
            cache = Path(directory) / "cache"
            reads = self.write_reads(directory, "longer", SEEN, LOST_A_BASE)
            output, table, summary = self.generate(directory, reads, cache)
            self.assertEqual(
                (table[LATE]["AlleleID"], table[LATE]["Variants"], table[LATE]["MarkerType"]),
                ("m1|RefMatch", "78:TG>T;overhang=1", "indel"),
            )
            self.assertEqual(summary["read_continuation"], {"bases_after_tag": SEEN, "rows_settled": 1, "reads_settled": 24})
            self.assertEqual(summary["unresolved_ends"]["rows"], 0)
            self.assertEqual(set(table), set(before))
            for sequence in set(table) - {LATE}:
                self.assertEqual(table[sequence], before[sequence])
            self.assertEqual((table[LATE]["s1"], table[LATE]["s2"]), (before[LATE]["s1"], before[LATE]["s2"]))
            # The read pass kept what follows the two canonical tags and the unread row, each from
            # its third read on, and nothing of the row the tag decides.
            with gzip.open(sorted(cache.glob("s1.*.json.gz"))[0], "rt") as handle:
                kept = {tag: found for _, tag, _, *found in json.load(handle)["counts"]}
            self.assertEqual(
                kept,
                {LOCUS.ref_tag: [[SUBSTITUTED, 38]], LOCUS.alt_tag: [[SUBSTITUTED, 28]], LATE: [[LOST_A_BASE, 10]], self.PLAIN: []},
            )
            # A second run takes them from the cache and writes the same table.
            first = output.read_bytes()
            with mock.patch.object(mp, "_write_cache", side_effect=AssertionError("the cache was not used")):
                again, _, _ = self.generate(directory, reads, cache)
            self.assertEqual(again.read_bytes(), first)

    def test_the_same_tag_settles_the_other_way_when_its_reads_say_so(self):
        with tempfile.TemporaryDirectory() as directory:
            _, table, _ = self.run_generate(directory, "substituted", SEEN, SUBSTITUTED)
        self.assertEqual(
            (table[LATE]["AlleleID"], table[LATE]["Variants"], table[LATE]["MarkerType"]),
            ("m1|AltMatch", "79:G>T;80:T>G;81:G>C", "snp"),
        )

    def test_table_without_primers_carries_the_settled_reading(self):
        # The rows without their primers are not aligned again: each settled row is written from
        # the alignment it was settled with, in the positions of the Ref row as printed.
        with tempfile.TemporaryDirectory() as directory:
            results = []
            for sample in ("s1", "s2"):
                counts = {(0, LOCUS.ref_tag): 40, (0, LOCUS.alt_tag): 30, (0, LATE): 12, (0, INSERTION): 12}
                kept = {
                    (0, LOCUS.ref_tag): [(SUBSTITUTED, 40)], (0, LATE): [(LOST_A_BASE, 12)], (0, INSERTION): [(AFTER_INSERTION, 12)],
                }
                results.append(mp.SampleResult(mp.Sample(sample, sample, "unused.fastq", {}), counts, {}, {}, {}, {}, continuations=kept))
            output = Path(directory) / "madc.csv"
            mp.write_madc(
                output, [LOCUS], results, mp.EvidenceRule(1, 20, 0, 0, 0, 5, 0.0), mp.OutputOptions(remove_primers=True), BAND,
            )
            with open(output, newline="") as handle:
                table = {row["AlleleSequence"]: row for row in csv.DictReader(handle)}
        row = table[LATE[LEFT:]]
        self.assertEqual((row["AlleleID"], row["Variants"], row["MarkerType"]), ("m1|RefMatch", "58:TG>T;overhang=1", "indel"))
        self.assertEqual(replay(REF[LEFT:], row["Variants"], LATE[LEFT:]), LATE[LEFT:])
        row = table[INSERTION[LEFT:]]
        self.assertEqual((row["AlleleID"], row["Variants"], row["MarkerType"]), ("m1|Other", "56:A>AACG;clip=3", "indel"))
        self.assertEqual(replay(REF[LEFT:], row["Variants"], INSERTION[LEFT:]), INSERTION[LEFT:])

    def test_a_settled_row_with_a_gap_in_its_left_primer_is_described_from_the_bases_shown(self):
        # Two bases inserted in the left primer shift the row against the Ref row from there on, so
        # the two rows as printed start at different bases of the alignment.  It is not carried
        # then: the printed rows are aligned as they stand, and the cell still rebuilds the row.
        whole = REF[:8] + "GG" + REF[8:75] + "ACG" + REF[75:] + AFTER
        row, after_row = whole[:81], whole[81 : 81 + SEEN]
        found = continuations(after_row, row=row)
        self.assertEqual(cell(LOCUS, row, None), "unresolved=73")   # the two readings part at the first inserted base
        self.assertEqual(cell(LOCUS, row, found), "8:A>AGG;75:A>AACG;clip=5")
        with tempfile.TemporaryDirectory() as directory:
            results = []
            for sample in ("s1", "s2"):
                counts = {(0, LOCUS.ref_tag): 40, (0, LOCUS.alt_tag): 30, (0, row): 12}
                kept = {(0, LOCUS.ref_tag): [(SUBSTITUTED, 40)], (0, row): [(after_row, 12)]}
                results.append(mp.SampleResult(mp.Sample(sample, sample, "unused.fastq", {}), counts, {}, {}, {}, {}, continuations=kept))
            output = Path(directory) / "madc.csv"
            mp.write_madc(
                output, [LOCUS], results, mp.EvidenceRule(1, 20, 0, 0, 0, 5, 0.0), mp.OutputOptions(remove_primers=True), BAND,
            )
            with open(output, newline="") as handle:
                table = {entry["AlleleSequence"]: entry for entry in csv.DictReader(handle)}
        shown = table[row[LEFT:]]
        self.assertEqual((shown["AlleleID"], shown["Variants"]), ("m1|Other", "1:.>AC;56:A>AC;58:T>A;clip=3"))
        self.assertEqual(replay(REF[LEFT:], shown["Variants"], row[LEFT:]), row[LEFT:])


if __name__ == "__main__":
    unittest.main()
