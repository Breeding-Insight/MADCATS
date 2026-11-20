"""Tests of ``generate --remove-primers on``: where the primers are, the cut in each row, and the table.

The panels here are built the way a fixed-amplicon-length platform makes them:
every amplicon is left primer + insert + right primer, and a tag is its first
bases, so tags longer than left primer + insert end with the start of the
right primer.
"""
import csv
import random
import contextlib
import dataclasses
import io
import json
import tempfile
import unittest
from collections import Counter, defaultdict
from pathlib import Path
from unittest import mock

import madcats as mp

LEFT, INSERT, RIGHT = 20, 32, 18   # lengths of the three parts of an amplicon
TAG = 60                           # so a Ref tag holds the first 8 bases of its right primer
TARGET = 10                        # index of the target base in the insert
ADAPTER = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"   # what follows the right primer in a read
GAP_BAND = 10
EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "madcats-demo"


def dna(rng, length):
    return "".join(rng.choice("ACGT") for _ in range(length))


def other_base(base):
    return "ACGT"[("ACGT".index(base) + 1) % 4]


def substituted(sequence, index):
    return sequence[:index] + other_base(sequence[index]) + sequence[index + 1 :]


def make_panel(markers=30, tag=TAG, seed=7, reverse=False):
    """Loci plus, per marker, the parts its tags were cut from: (left primer, insert, right primer)."""
    rng = random.Random(seed)
    loci, parts = [], {}
    for index in range(markers):
        left, insert, right = dna(rng, LEFT), dna(rng, INSERT), dna(rng, RIGHT)
        alt_insert = substituted(insert, TARGET)
        marker = f"m{index:02d}"
        parts[marker] = (left, insert, right)
        loci.append(
            mp.Locus(
                index, marker, left, mp.reverse_complement(right) if reverse else right,
                insert[TARGET], alt_insert[TARGET], LEFT + TARGET,
                (left + insert + right + ADAPTER)[:tag], (left + alt_insert + right + ADAPTER)[:tag],
            )
        )
    return loci, parts


def tag_of(parts, insert, tag=TAG, left=None, right=None):
    """The tag a read of this insert gives: a fixed-length window of its amplicon."""
    return ((left or parts[0]) + insert + (right or parts[2]) + ADAPTER)[:tag]


def rebuild(ref, events):
    """Apply a row's events to its Ref row; every event reads "replace REF at POS by ALT".

    Plain replay, left to right: an event that starts on a Ref base an earlier event already
    used is an error, because a reader of the Variants column could not apply it either."""
    out, done = [], 0
    for kind, position, ref_allele, alt_allele in events:
        start = position - 1
        if start < done:
            raise AssertionError(f"two events use Ref base {position}: {mp.format_variants(events)}")
        if kind == "clip":
            out.append(ref[done:start])
            done = len(ref)
        elif kind == "overhang":
            out.append(ref[done:] + alt_allele)
            done = len(ref)
        elif kind == "unresolved":   # the Ref bases from here on are not read: the row's own bases stand there
            out.append(ref[done:start] + alt_allele)
            done = len(ref)
        else:
            out.append(ref[done:start] + alt_allele.replace(".", ""))
            done = start + len(ref_allele.replace(".", ""))
    return "".join(out) + ref[done:]


class GeometryTests(unittest.TestCase):
    def test_amplicon_length_is_read_from_the_tags(self):
        loci, _ = make_panel()
        self.assertEqual(mp.amplicon_geometry(loci), (LEFT + INSERT + RIGHT, False))
        site, reason = mp.right_primer_site(loci[0], mp.amplicon_geometry(loci))
        self.assertEqual((site, reason), ((LEFT + INSERT, loci[0].right_primer), "located"))

    def test_primers_listed_reverse_complemented(self):
        loci, parts = make_panel(reverse=True)
        geometry = mp.amplicon_geometry(loci)
        self.assertEqual(geometry, (LEFT + INSERT + RIGHT, True))
        self.assertEqual(mp.right_primer_site(loci[3], geometry)[0], (LEFT + INSERT, parts["m03"][2]))

    def test_no_amplicon_length_without_enough_tags_that_reach_the_primer(self):
        self.assertIsNone(mp.amplicon_geometry(make_panel(tag=LEFT + INSERT)[0]))        # tags stop at the primer
        self.assertIsNone(mp.amplicon_geometry(make_panel(markers=mp.AMPLICON_MIN_MARKERS - 1)[0]))
        no_primers = [dataclasses.replace(locus, right_primer="") for locus in make_panel()[0]]
        self.assertIsNone(mp.amplicon_geometry(no_primers))
        self.assertEqual(mp.right_primer_site(no_primers[0], None), (None, "no_right_primer"))

    def test_primer_past_the_tag_end_is_still_located(self):
        loci, _ = make_panel()
        geometry = mp.amplicon_geometry(loci)
        short = dataclasses.replace(loci[0], ref_tag=loci[0].ref_tag[:45], alt_tag=loci[0].alt_tag[:45])
        self.assertEqual(mp.right_primer_site(short, geometry)[0][0], LEFT + INSERT)

    def test_marker_that_cannot_be_cut_is_not_located(self):
        loci, parts = make_panel()
        geometry = mp.amplicon_geometry(loci)
        left, insert, right = parts["m00"]
        # The tag's last bases are not the start of the primer the panel lists.
        clash = dataclasses.replace(loci[0], right_primer=substituted(right, 2))
        self.assertEqual(mp.right_primer_site(clash, geometry), (None, "disagrees_with_amplicon_length"))
        # The target lies where the primer starts.
        in_primer = dataclasses.replace(loci[0], variant_pos=LEFT + INSERT)
        self.assertEqual(mp.right_primer_site(in_primer, geometry), (None, "target_in_right_primer"))
        # Without an amplicon length, only a primer that lies whole in the Ref tag is located.
        self.assertEqual(mp.right_primer_site(loci[0], None), (None, "no_amplicon_length"))
        whole = dataclasses.replace(loci[0], ref_tag=(left + insert[:22] + right + ADAPTER)[:TAG])
        self.assertEqual(mp.right_primer_site(whole, None), ((LEFT + 22, right), "located"))
        # With one, such a tag does not read as the primer where the panel's length puts it: the
        # marker's amplicon is another length, and its primer is the whole copy.
        self.assertEqual(mp.right_primer_site(whole, geometry), ((LEFT + 22, right), "located"))

    def test_copy_of_the_primer_in_the_insert_is_not_the_primer(self):
        loci, parts = make_panel()
        geometry = mp.amplicon_geometry(loci)
        left, insert, right = parts["m00"]
        # The insert repeats the marker's right primer after the target.  The tag still reads as the
        # primer where the panel's amplicon length puts it, so that is the primer.
        repeat = insert[: TARGET + 2] + right + insert[TARGET + 2 + RIGHT :]
        tags = lambda length: {
            "ref_tag": (left + repeat + right + ADAPTER)[:length],
            "alt_tag": (left + substituted(repeat, TARGET) + right + ADAPTER)[:length],
        }
        for what, length in (("tag with 8 bases of the primer", TAG), ("tag that stops at the primer", LEFT + INSERT)):
            with self.subTest(tag=what):
                locus = dataclasses.replace(loci[0], **tags(length))
                site, reason = mp.right_primer_site(locus, geometry)
                self.assertEqual((site, reason), ((LEFT + INSERT, right), "located"))
                # No row is cut at the copy: every row keeps its insert.
                for row_insert in (repeat, substituted(repeat, TARGET), substituted(repeat, 5)):
                    row = (left + row_insert + right + ADAPTER)[:length]
                    cut, kept = mp.right_primer_cut(locus, row, site, GAP_BAND)
                    self.assertEqual((row[LEFT:][:cut], kept), (row_insert, None))


class CutTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.loci, cls.parts = make_panel()
        cls.geometry = mp.amplicon_geometry(cls.loci)

    def trimmed(self, locus, row, site="panel"):
        if site == "panel":
            site = mp.right_primer_site(locus, self.geometry)[0]
        cut, kept = mp.right_primer_cut(locus, row, site, GAP_BAND)
        return row[len(locus.left_primer) :][:cut], kept

    def test_rows_are_cut_where_their_primer_begins(self):
        for locus in self.loci[:10]:
            parts = self.parts[locus.marker]
            left, insert, right = parts
            inserts = {
                "the Ref tag": insert,
                "a substitution": substituted(insert, 25),
                "a deletion, which brings more of the primer in": insert[:18] + insert[21:],
                "an insertion, which pushes the primer toward the end": insert[:18] + "TTG" + insert[18:],
                "an insertion that pushes the primer out of the row": insert[:14] + "TTGACGTCA" + insert[14:],
            }
            for what, row_insert in inserts.items():
                with self.subTest(marker=locus.marker, row=what):
                    expected = row_insert[: TAG - LEFT]   # the last case fills the row
                    self.assertEqual(self.trimmed(locus, tag_of(parts, row_insert)), (expected, None))
            with self.subTest(marker=locus.marker, row="substitutions in the primer copies"):
                row = tag_of(parts, insert, left=substituted(left, 4), right=substituted(right, 1))
                self.assertEqual(self.trimmed(locus, row), (insert, None))

    def test_indel_next_to_the_primer_is_read_from_the_primer(self):
        # Four primer bases in the tag, and a primer whose third and fourth bases repeat its first two:
        # an indel of two bases right before it costs the same as two substitutions.
        left, insert, right = dna(random.Random(3), LEFT), "GCTTAGGCTAACCGGTTACTTAGCAGTTCAGT", "ACACGGTCTGAACTGGAT"
        parts = (left, insert, right)
        tag = LEFT + INSERT + 4
        locus = mp.Locus(0, "m", left, right, insert[TARGET], "T", LEFT + TARGET, tag_of(parts, insert, tag),
                         tag_of(parts, substituted(insert, TARGET), tag))
        site = (LEFT + INSERT, right)
        deletion = tag_of(parts, insert[:-2], tag)
        self.assertEqual(self.trimmed(locus, deletion, site), (insert[:-2], None))
        insertion = tag_of(parts, insert + "TG", tag)
        self.assertEqual(self.trimmed(locus, insertion, site), (insert + "TG", None))
        self.assertEqual(self.trimmed(locus, locus.alt_tag, site), (substituted(insert, TARGET), None))

    def test_short_product_is_found_without_a_site(self):
        locus = self.loci[0]
        parts = self.parts[locus.marker]
        left, insert, right = parts
        # The whole primer lies in the row, 17 bases earlier than the indel limit could explain.
        self.assertEqual(self.trimmed(locus, tag_of(parts, insert[:15]), None), (insert[:15], None))
        self.assertEqual(self.trimmed(locus, tag_of(parts, insert[:15])), (insert[:15], None))
        # The same with the primer listed on the other strand.
        flipped = dataclasses.replace(locus, right_primer=mp.reverse_complement(right))
        self.assertEqual(self.trimmed(flipped, tag_of(parts, insert[:15]), None), (insert[:15], None))
        # A row that ends in the first 10 bases of the primer, with no site to say so.
        tail = (left + insert + right)[: LEFT + INSERT + mp.PRIMER_TAIL_UNLOCATED]
        self.assertEqual(self.trimmed(locus, tail, None), (insert, None))
        # Nine bases could be chance there.  With a site, eight are enough.
        self.assertEqual(self.trimmed(locus, tail[:-1], None), (tail[LEFT:-1], None))
        elsewhere = (LEFT + INSERT + 14, right)   # a site too far for the alignment to reach
        self.assertEqual(self.trimmed(locus, tail[:-2], elsewhere), (insert, None))
        self.assertEqual(self.trimmed(locus, tail[:-3], elsewhere), (tail[LEFT:-3], None))

    def test_degenerate_right_primer_is_found_whole(self):
        # One marker, so no amplicon length: the primer is looked for itself, and the panel lists
        # it with an N where the tag holds a base.
        locus = self.loci[0]
        left, insert, right = self.parts[locus.marker]
        degenerate = dataclasses.replace(
            locus, right_primer=right[:5] + "N" + right[6:],
            ref_tag=(left + insert + right + ADAPTER)[:80], alt_tag=(left + substituted(insert, TARGET) + right + ADAPTER)[:80],
        )
        site = mp.right_primer_site(degenerate, mp.amplicon_geometry([degenerate]))
        self.assertEqual(site, ((LEFT + INSERT, degenerate.right_primer), "located"))
        self.assertEqual(self.trimmed(degenerate, degenerate.ref_tag, site[0]), (insert, None))
        self.assertEqual(self.trimmed(degenerate, degenerate.ref_tag, None), (insert, None))
        flipped = dataclasses.replace(degenerate, right_primer=mp.reverse_complement(degenerate.right_primer))
        self.assertEqual(self.trimmed(flipped, flipped.ref_tag, None), (insert, None))

    def test_canonical_rows_never_merge(self):
        locus = self.loci[0]
        rows = [(0, "Ref", locus.ref_tag), (1, "Alt", locus.alt_tag)]
        counts = {locus.ref_tag: [3], locus.alt_tag: [4]}
        tally = Counter()
        with mock.patch.object(mp, "right_primer_cut", return_value=(TARGET - 2, None)):   # a cut that loses the target
            shown, trimmed, reads, _ = mp.trim_locus_rows(locus, rows, counts, None, GAP_BAND, tally)
        self.assertEqual(tally["markers_right_primer_kept"], 1)
        self.assertEqual((shown.ref_tag, shown.alt_tag), (locus.ref_tag[LEFT:], locus.alt_tag[LEFT:]))
        self.assertEqual([row[1:] for row in trimmed], [("Ref", locus.ref_tag[LEFT:]), ("Alt", locus.alt_tag[LEFT:])])
        self.assertEqual(reads, {locus.ref_tag[LEFT:]: [3], locus.alt_tag[LEFT:]: [4]})

    def test_row_keeps_its_primer_rather_than_lose_the_target_or_most_of_itself(self):
        locus = self.loci[0]
        parts = self.parts[locus.marker]
        whole = lambda row: row[LEFT:]
        before_target = tag_of(parts, parts[1][: TARGET - 2])
        self.assertEqual(self.trimmed(locus, before_target), (whole(before_target), "target"))
        early_target = dataclasses.replace(locus, variant_pos=LEFT + 2)
        few_bases = tag_of(parts, parts[1][: mp.TRIM_MIN_BASES - 1])
        self.assertEqual(self.trimmed(early_target, few_bases), (whole(few_bases), "short"))

    def test_row_keeps_the_base_its_category_was_read_from(self):
        # A diverged haplotype: seven substitutions, and two bases missing shortly before the
        # target, which is the last insert base but one.  Seven substitutions are too many for a
        # gapped reading, so the row is read at the fixed offsets.  There its first primer base
        # stands on the target, and that base is the Alt base.  The alignment that finds the
        # primer reads the deletion and places the target two bases earlier, so it would cut
        # the primer off together with the base the category was read from.
        rng = random.Random(2)
        parts = left, insert, right = dna(rng, LEFT), dna(rng, INSERT), dna(rng, RIGHT)
        target = INSERT - 2
        alt_insert = insert[:target] + right[0] + insert[target + 1 :]
        locus = mp.Locus(0, "m", left, right, insert[target], right[0], LEFT + target, tag_of(parts, insert), tag_of(parts, alt_insert))
        diverged = insert
        for index in range(1, 20, 3):
            diverged = substituted(diverged, index)
        row = tag_of(parts, diverged[:26] + diverged[28:])
        site = (LEFT + INSERT, right)
        self.assertIsNone(mp.row_reading(locus, row, GAP_BAND))
        self.assertEqual(mp.classify_sequence(locus, row, GAP_BAND), "AltMatch")
        self.assertTrue(right.startswith(row[LEFT + target :]))   # the primer begins on the base that was read
        self.assertEqual(self.trimmed(locus, row, site), (row[LEFT:], "target"))
        # On a marker with its target elsewhere, the same row is cut where its primer begins.
        early = dataclasses.replace(locus, variant_pos=LEFT + TARGET, alt_tag=tag_of(parts, substituted(insert, TARGET)))
        self.assertEqual(self.trimmed(early, row, site), (row[LEFT : LEFT + INSERT - 2], None))


class DescriptionTests(unittest.TestCase):
    """Events of rows of another length than their Ref row, as a table without primers has them."""

    REF = "TTACGGATCCTTAGCATGCAAGTCGGATTACA"

    def events(self, row, whole=False):
        locus = mp.Locus(0, "m", "", "", "G", "C", 5, self.REF, substituted(self.REF, 5))
        events = mp._describe(self.REF, row, GAP_BAND, mp._read_through_gaps(locus, row, GAP_BAND), whole)
        if not whole:
            self.assertEqual(events, mp.variant_events(locus, row, GAP_BAND))
        self.assertEqual(rebuild(self.REF, events), row)
        self.assertEqual(mp.parse_variants(mp.format_variants(events), self.REF, row), events)
        return mp.format_variants(events)

    def test_rows_of_another_length_are_aligned(self):
        for whole in (False, True):
            self.assertEqual(self.events(self.REF[:12] + self.REF[15:], whole), "12:TAGC>T")
            self.assertEqual(self.events(self.REF[:12] + "GGA" + self.REF[12:], whole), "12:T>TGGA")
            self.assertEqual(self.events(self.REF[:12] + self.REF[15:20] + "T" + self.REF[21:], whole), "12:TAGC>T;21:A>T")

    def test_end_of_a_row_cut_at_its_primer_is_an_indel_and_of_any_other_row_is_not_shown(self):
        self.assertEqual(self.events(self.REF[:-2], whole=True), "29:TAC>T")   # the last two bases, left-normalised
        self.assertEqual(self.events(self.REF + "GT", whole=True), "32:A>AGT")
        self.assertEqual(self.events(self.REF[:-2]), "clip=2")
        self.assertEqual(self.events(self.REF + "GT"), "overhang=2")
        # A deletion that brought in bases past the end of the Ref row.
        self.assertEqual(self.events(self.REF[:12] + self.REF[15:] + "G"), "12:TAGC>T;overhang=1")

    def test_shared_end_is_never_inferred_from_the_lengths(self):
        # Two tags of different lengths need not end at the same place: the shorter may be cut at
        # its primer while an insertion pushed the primer out of the longer one.
        self.assertEqual(mp._designed_events(self.REF, self.REF[:12] + self.REF[15:], GAP_BAND), (("del", 12, "TAGC", "T"),))
        self.assertEqual(mp.format_variants(mp._designed_events(self.REF, self.REF[:-2], GAP_BAND)), "clip=2")
        pushed_out = self.REF[:12] + "GGACCTTGAC" + self.REF[12:-2]
        self.assertEqual(mp.format_variants(mp._designed_events(self.REF, pushed_out, GAP_BAND)), "12:T>TGGACCTTGAC;clip=2")

    def test_gap_on_the_first_base_is_written_without_an_anchor(self):
        for whole in (False, True):
            self.assertEqual(self.events(self.REF[1:], whole), "1:T>.")
            self.assertEqual(self.events("GG" + self.REF, whole), "1:.>GG")
            # A change on the base after the gap is its own event.
            self.assertEqual(self.events("C" + self.REF[2:], whole), "1:T>.;2:T>C")

    def test_no_two_events_use_the_same_ref_base(self):
        # Each of these put two events on base 1, or on the base after a first-base deletion, when
        # the gap was written from the base after it; replaying them repeated that base.
        ref = self.REF
        cases = {
            "an insertion before base 1 and one after it": ("G" + ref[0] + "G" + ref[1:], "1:.>G;1:T>TG"),
            "a longer insertion before base 1 and one after it": ("CAAC" + ref[0] + "GCC" + ref[1:], "1:.>CAAC;1:T>TGCC"),
            "bases 1-2 deleted, then a deletion after base 3": (ref[2] + ref[4:], "1:TT>.;3:AC>A"),
        }
        for what, (row, expected) in cases.items():
            for whole in (False, True):
                with self.subTest(case=what, whole=whole):
                    self.assertEqual(self.events(row, whole), expected)   # events() replays them too

    def test_rows_too_different_in_length_are_compared_from_the_start(self):
        short = self.REF[:14]
        self.assertEqual(self.events(short), "clip=18")
        self.assertEqual(self.events(substituted(short, 3)), f"4:C>{other_base('C')};clip=18")
        self.assertEqual(self.events(self.REF + "ACGTACGTACGT"), "overhang=12")

    def test_rows_of_the_ref_length_are_described_as_before(self):
        self.assertEqual(self.events(substituted(self.REF, 20)), "21:A>C")
        self.assertEqual(self.events(self.REF[:8] + self.REF[11:] + "GGT"), "7:ATCC>A;overhang=3")


class TableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.loci, cls.parts = make_panel()

    def write(self, rows_by_marker, remove_primers, directory, gap_band=GAP_BAND, **options):
        """Write a two-sample table; ``rows_by_marker`` maps a marker to {tag: (reads, reads)}."""
        results = []
        index_of = {locus.marker: locus.index for locus in self.loci}
        for sample in range(2):
            counts = {
                (index_of[marker], tag): reads[sample]
                for marker, rows in rows_by_marker.items()
                for tag, reads in rows.items()
                if reads[sample]
            }
            results.append(mp.SampleResult(mp.Sample(f"s{sample}", f"s{sample}", "unused.fastq", {}), counts, {}, {}, {}, {}))
        output = Path(directory) / f"madc_{'on' if remove_primers else 'off'}.csv"
        rule = mp.EvidenceRule(1, 20, 0, 0, 0, 5, 0.0)   # a row needs 5 reads in one sample
        summary = mp.write_madc(
            output, self.loci, results, rule,
            mp.OutputOptions(variants_path=output.with_suffix(".variants.csv"), remove_primers=remove_primers, **options),
            gap_band,
        )
        with open(output, newline="") as handle:
            table = list(csv.DictReader(handle))
        return summary, table

    def test_table_without_primers(self):
        locus = self.loci[0]
        parts = self.parts[locus.marker]
        left, insert, right = parts
        alt_insert = substituted(insert, TARGET)
        rows = {
            locus.ref_tag: (40, 0),
            locus.alt_tag: (0, 30),
            tag_of(parts, insert, left=substituted(left, 3)): (7, 0),          # a miscopied left primer
            tag_of(parts, alt_insert, right=substituted(right, 2)): (0, 6),   # a miscopied right primer
            tag_of(parts, substituted(insert, 25)): (9, 8),
            tag_of(parts, insert[:18] + insert[21:]): (5, 11),
        }
        with tempfile.TemporaryDirectory() as directory:
            summary, table = self.write({locus.marker: rows}, True, directory)
            _, untrimmed = self.write({locus.marker: rows}, False, directory)
        marker_rows = [row for row in table if row["CloneID"] == locus.marker]
        self.assertEqual(len([row for row in untrimmed if row["CloneID"] == locus.marker]), 6)
        shown = [
            (row["AlleleID"].split("|")[1], row["AlleleSequence"], row["Variants"], row["MarkerType"], row["s0"], row["s1"])
            for row in marker_rows
        ]
        deleted = insert[:18] + insert[21:]
        self.assertEqual(
            sorted(shown),
            sorted([
                ("Ref", insert, ".", ".", "47", "0"),
                ("Alt", alt_insert, f"{TARGET + 1}:{insert[TARGET]}>{alt_insert[TARGET]}", "snp", "0", "36"),
                ("RefMatch", substituted(insert, 25), f"26:{insert[25]}>{other_base(insert[25])}", "snp", "9", "8"),
                ("RefMatch", deleted, mp.format_variants(mp._describe(insert, deleted, GAP_BAND, True)), "indel", "5", "11"),
            ]),
        )
        for row in marker_rows:
            canonical = row["AlleleID"].endswith(("|Ref", "|Alt"))
            self.assertEqual(row["ClusterConsensusSequence"], row["AlleleSequence"] if canonical else insert)
        # Every marker's rows lose their primers, with or without reads.
        for row in table:
            if row["CloneID"] != locus.marker:
                insert = self.parts[row["CloneID"]][1]
                self.assertIn(row["AlleleSequence"], (insert, substituted(insert, TARGET)))
        self.assertEqual(
            summary["primer_removal"],
            {
                "amplicon_length": LEFT + INSERT + RIGHT,
                "right_primers_reverse_complemented": False,
                "markers_by_right_primer": {"located": len(self.loci)},
                "rows_merged": 2,
                "rows_merged_across_allele_groups": 0,
                "rows_right_primer_removed": 2 * len(self.loci) + 4,
                "rows_right_primer_kept_target": 0,
                "rows_right_primer_kept_short": 0,
                "markers_right_primer_kept": 0,
            },
        )
        self.assertEqual(summary["rows"], 2 * len(self.loci) + 2)
        self.assertEqual(summary["reads_in_emitted_rows"], sum(sum(reads) for reads in rows.values()))

    def variants_of(self, table, sequence):
        return next(row["Variants"] for row in table if row["AlleleSequence"] == sequence)

    def test_known_target_alleles_and_events_before_and_after_primer_removal(self):
        locus = self.loci[0]
        parts = self.parts[locus.marker]
        insert = parts[1]
        # Fixed simulated truth: the target is G/T at full position 31, insert
        # position 11. Expected events are literal Ref coordinates, independent
        # of the aligner and annotation helpers being tested.
        self.assertEqual(insert, "ATACACGTCAGCACGAAACTTGTTGGCCCAGT")
        alt = substituted(insert, TARGET)
        cases = [
            (insert[:4] + insert[6:], "RefMatch", "22:TAC>T;overhang=2", "2:TAC>T", "indel"),
            (alt[:4] + alt[6:], "AltMatch", "22:TAC>T;31:G>T;overhang=2", "2:TAC>T;11:G>T", "indel"),
            (insert[:4] + "AG" + insert[4:], "RefMatch", "24:C>CAG;clip=2", "4:C>CAG", "indel"),
            (alt[:4] + "AG" + alt[4:], "AltMatch", "24:C>CAG;31:G>T;clip=2", "4:C>CAG;11:G>T", "indel"),
            (insert[:10] + insert[11:], "Other", "30:AG>A;overhang=1", "10:AG>A", "indel"),
            (alt[:17] + alt[19:], "AltMatch", "31:G>T;37:AAC>A;overhang=2", "11:G>T;17:AAC>A", "indel"),
            (substituted(insert, 25), "RefMatch", "46:G>T", "26:G>T", "snp"),
        ]
        counts = {locus.ref_tag: (40, 0), locus.alt_tag: (0, 30)}
        counts.update({tag_of(parts, row): (9, 11) for row, *_ in cases})
        with tempfile.TemporaryDirectory() as directory:
            for trimmed in (False, True):
                with self.subTest(remove_primers=trimmed):
                    summary, table = self.write({locus.marker: counts}, trimmed, directory)
                    by_sequence = {row["AlleleSequence"]: row for row in table if row["CloneID"] == locus.marker}
                    for insert_row, category, full_events, trimmed_events, marker_type in cases:
                        sequence = insert_row if trimmed else tag_of(parts, insert_row)
                        row = by_sequence[sequence]
                        self.assertEqual(row["AlleleID"], f"{locus.marker}|{category}")
                        self.assertEqual(row["Variants"], trimmed_events if trimmed else full_events)
                        self.assertEqual(row["MarkerType"], marker_type)
                        self.assertEqual((row["s0"], row["s1"]), ("9", "11"))
                    self.assertEqual(summary["reads_in_emitted_rows"], 70 + 20 * len(cases))
                    if trimmed:
                        self.assertEqual(summary["primer_removal"]["rows_merged_across_allele_groups"], 0)

    def test_rescue_off_can_misclassify_an_indel_that_passes_assignment(self):
        # A known Alt haplotype has an insertion upstream of its target. The
        # base at the old offset happens to be Ref, despite passing the default
        # 80% identity filter. Primer removal cannot repair this fixed-offset
        # interpretation. Exercise all four switch combinations through FASTQ.
        ref = "CCGTAATGCCTTTCCCTAACAGAGTTTTTCGAACTCGTGTTGTCGAGCGACGGAATTAGATCAGTTAAATGGCAGAAAACT"
        target = 65
        alt = ref[:target] + "A" + ref[target + 1:]
        indel = (alt[:60] + "C" + alt[60:])[:len(ref)]
        self.assertEqual((ref[target], indel[target], indel[target + 1]), ("T", "T", "A"))
        locus = mp.Locus(0, "late", ref[:20], "", "T", "A", target, ref, alt)
        self.assertEqual(mp.select_marker(indel, [0], [locus], 16, gap_band=0), (0, "assigned"))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tags = root / "tags.csv"
            mp.write_reference_tags(tags, [locus])
            fastqs = root / "fastq"
            fastqs.mkdir()
            with (fastqs / "sample.fastq").open("w") as handle:
                for number, sequence in enumerate([ref] * 40 + [alt] * 40 + [indel] * 40):
                    handle.write(f"@r{number}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            for rescue in ("on", "off"):
                for trimmed in ("on", "off"):
                    with self.subTest(gapped_rescue=rescue, remove_primers=trimmed):
                        output = root / "madc.csv"
                        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                            mp.main([
                                "generate", "--fastq-dir", str(fastqs), "--reference-tags", str(tags),
                                "--output", str(output), "--jobs", "1", "--gapped-rescue", rescue,
                                "--remove-primers", trimmed,
                            ])
                        with output.open(newline="") as handle:
                            rows = list(csv.DictReader(handle))
                        row = next(row for row in rows if row["AlleleSequence"] == (indel[20:] if trimmed == "on" else indel))
                        self.assertEqual(row["AlleleID"], "late|AltMatch" if rescue == "on" else "late|RefMatch")
                        self.assertEqual(row["sample"], "40")
                        self.assertEqual(sum(int(row["sample"]) for row in rows), 120)
                        if rescue == "on":
                            self.assertEqual(row["Variants"], "40:A>AC;46:T>A;clip=1" if trimmed == "on" else "60:A>AC;66:T>A;clip=1")
                            self.assertEqual(row["MarkerType"], "indel")
                        else:
                            self.assertEqual(row["MarkerType"], "snp")
                            self.assertNotIn("clip=", row["Variants"])
                            self.assertNotIn("46:T>A" if trimmed == "on" else "66:T>A", row["Variants"])

    def test_gapped_reading_of_the_full_tag_is_kept_after_the_primer_that_verified_it_is_cut(self):
        # A deletion and an insertion of three bases each: the row keeps the Ref length.  On the full
        # tag the primer bases after the insertion verify the reading; the trimmed row has too few
        # bases left to verify it again, so the decision made on the full tag is the one used.
        locus = self.loci[0]
        parts = self.parts[locus.marker]
        insert = parts[1]
        row_insert = insert[:12] + insert[15:21] + "TGA" + insert[21:]
        tag = tag_of(parts, row_insert)
        rows = {locus.marker: {locus.ref_tag: (40, 0), locus.alt_tag: (0, 30), tag: (9, 9)}}
        with tempfile.TemporaryDirectory() as directory:
            _, on = self.write(rows, True, directory)
            _, off = self.write(rows, False, directory)
        full = mp.variant_events(locus, tag, GAP_BAND)
        self.assertEqual([event[0] for event in full], ["del", "ins"])
        self.assertEqual(self.variants_of(off, tag), mp.format_variants(full))
        moved = [(kind, position - LEFT, ref_allele, alt_allele) for kind, position, ref_allele, alt_allele in full]
        self.assertEqual(self.variants_of(on, row_insert), mp.format_variants(moved))
        self.assertEqual(len(row_insert), INSERT)

    def test_designed_insertion_that_pushes_the_alt_primer_out_of_the_tag(self):
        locus = self.loci[0]
        parts = self.parts[locus.marker]
        insert = parts[1]
        alt_insert = insert[:TARGET] + "ACGTTGCAAT" + insert[TARGET:]
        design = dataclasses.replace(locus, alt_tag=tag_of(parts, alt_insert), alt_base=alt_insert[TARGET])
        self.assertNotEqual(design.ref_base, design.alt_base)
        loci = [design] + self.loci[1:]
        results = [
            mp.SampleResult(mp.Sample(name, name, "unused.fastq", {}), counts, {}, {}, {}, {})
            for name, counts in (("s0", {(0, design.ref_tag): 40}), ("s1", {(0, design.alt_tag): 30}))
        ]
        with tempfile.TemporaryDirectory() as directory:
            madc = Path(directory) / "madc.csv"
            mp.write_madc(madc, loci, results, mp.EvidenceRule(1, 20, 0, 0, 0, 5, 0.0), mp.OutputOptions(remove_primers=True), GAP_BAND)
            with open(madc, newline="") as handle:
                ref_row, alt_row = [row for row in csv.DictReader(handle) if row["CloneID"] == design.marker]
            # The Ref row is cut at its primer; the Alt row fills the tag and stops two bases short of its own.
            self.assertEqual(ref_row["AlleleSequence"], insert)
            self.assertEqual(alt_row["AlleleSequence"], alt_insert[: TAG - LEFT])
            expected = mp.format_variants(mp._describe(insert, alt_insert[: TAG - LEFT], GAP_BAND, True))
            self.assertEqual((alt_row["Variants"], alt_row["MarkerType"]), (expected, "indel"))
            self.assertRegex(alt_row["Variants"], r"^\d+:[ACGT]>[ACGT]{11};clip=2$")
            # call finds the designed insertion as the target, from the table alone and with the panel.
            tags = Path(directory) / "tags.csv"
            mp.write_reference_tags(tags, loci)
            for options in ([], ["--reference-tags", str(tags)]):
                prefix = Path(directory) / "calls"
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    mp.main(["call", "--madc", str(madc), "--output-prefix", str(prefix), "--jobs", "1", *options])
                with open(str(prefix) + ".targets.vcf") as handle:
                    records = [line.split("\t") for line in handle if line.startswith(design.marker + "\t")]
                self.assertEqual(len(records), 1, options)
                position, anchor, inserted = alt_row["Variants"].split(";")[0].replace(">", ":").split(":")
                self.assertEqual(records[0][1:5], [position, design.marker, anchor, inserted])
            # Without the Variants column the two tags' lengths alone do not say where they end:
            # no target record, and the QC file says why.
            plain = Path(directory) / "plain.csv"
            mp.write_madc(plain, loci, results, mp.EvidenceRule(1, 20, 0, 0, 0, 5, 0.0),
                          mp.OutputOptions(remove_primers=True, annotations=False), GAP_BAND)
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                mp.main(["call", "--madc", str(plain), "--output-prefix", str(prefix), "--jobs", "1"])
            with open(str(prefix) + ".targets.vcf") as handle:
                self.assertEqual([line for line in handle if line.startswith(design.marker + "\t")], [])
            with open(str(prefix) + ".qc.json") as handle:
                self.assertEqual(json.load(handle)["vcfs"]["target"]["targets_skipped_tags_differ_in_length"], 1)

    def test_off_is_the_default_and_adds_nothing_to_the_table_or_its_first_line(self):
        arguments = ["generate", "--fastq-dir", "fastq", "--output", "madc.csv", "--reference-tags", "tags.csv"]
        args = mp.build_parser().parse_args(arguments)
        self.assertEqual(args.remove_primers, "off")
        mp._resolve_generate_defaults(args)
        self.assertNotIn("remove_primers", mp._generate_provenance(args))
        args.remove_primers = "on"
        self.assertIn(" remove_primers=on ", mp._generate_provenance(args))
        self.assertFalse(mp.OutputOptions().remove_primers)
        locus = self.loci[0]
        rows = {locus.marker: {locus.ref_tag: (40, 0), locus.alt_tag: (0, 30)}}
        with tempfile.TemporaryDirectory() as directory:
            summary, table = self.write(rows, False, directory)
        self.assertNotIn("primer_removal", summary)
        self.assertEqual(table[0]["AlleleSequence"], locus.ref_tag)

    @unittest.skipUnless(EXAMPLE.is_dir(), "example data not present")
    def test_generate_refuses_a_table_without_primers_and_without_variants(self):
        with tempfile.TemporaryDirectory() as directory:
            tags, table = Path(directory) / "tags.csv", Path(directory) / "madc.csv"
            arguments = ["generate", "--fastq-dir", str(EXAMPLE / "fastq"), "--reference-tags", str(tags), "--jobs", "1",
                         "--output", str(table), "--remove-primers", "on"]
            captured = io.StringIO()
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(captured):
                mp.main(["design", "--primers", str(EXAMPLE / "primers.csv"), "--probe-design",
                         str(EXAMPLE / "probe_design.tsv"), "--output", str(tags)])
                with self.assertRaises(SystemExit):
                    mp.main(arguments + ["--annotations", "off"])
                self.assertFalse(table.exists())
                self.assertIn("cannot be combined with --annotations off", captured.getvalue())
                mp.main(arguments)                       # annotations are on by default
            with open(table, newline="") as handle:
                self.assertIn("Variants", next(csv.reader(handle)))
        self.assertIn("Removed primers:", captured.getvalue())

    def test_every_row_rebuilds_from_its_ref_row_and_variants_and_no_read_changes_allele(self):
        rng = random.Random(11)
        rows_by_marker = {}
        for locus in self.loci[:12]:
            parts = self.parts[locus.marker]
            rows = {locus.ref_tag: (rng.randint(5, 40), rng.randint(0, 40)), locus.alt_tag: (rng.randint(0, 40), rng.randint(5, 40))}
            for _ in range(25):
                insert = substituted(parts[1], TARGET) if rng.random() < 0.5 else parts[1]
                for _ in range(rng.randint(0, 2)):
                    insert = substituted(insert, rng.randrange(INSERT))
                kind = rng.choice(("none", "none", "deletion", "insertion"))
                at = rng.randrange(0, INSERT - 3)
                size = rng.randint(1, 4)
                if kind == "deletion":
                    insert = insert[:at] + insert[at + size :]
                elif kind == "insertion":
                    insert = insert[:at] + dna(rng, size) + insert[at:]
                left = substituted(parts[0], rng.randrange(LEFT)) if rng.random() < 0.2 else parts[0]
                rows.setdefault(tag_of(parts, insert, left=left), (rng.randint(5, 30), rng.randint(0, 30)))
            rows_by_marker[locus.marker] = rows
        with tempfile.TemporaryDirectory() as directory:
            summary, on = self.write(rows_by_marker, True, directory)
            _, off = self.write(rows_by_marker, False, directory)
            with open(Path(directory) / "madc_on.variants.csv", newline="") as handle:
                sidecar = defaultdict(list)
                for line in csv.DictReader(handle):
                    event = (line["Type"], int(line["Position"]), line["Ref"], line["Alt"])
                    sidecar[(line["CloneID"], line["AlleleSequence"])].append(event)
        ref_rows = {row["CloneID"]: row["AlleleSequence"] for row in on if row["AlleleID"].endswith("|Ref")}
        for row in on:
            events = sidecar[(row["CloneID"], row["AlleleSequence"])]
            self.assertEqual(rebuild(ref_rows[row["CloneID"]], events), row["AlleleSequence"], row["Variants"])
            self.assertEqual(mp.format_variants(events), row["Variants"])
            # The Variants cell alone gives the events back.
            self.assertEqual(mp.parse_variants(row["Variants"], ref_rows[row["CloneID"]], row["AlleleSequence"]), events)

        def reads_by_allele(table):
            totals = defaultdict(int)
            for row in table:
                group = mp.ALLELE_GROUP[mp.allele_type(row["AlleleID"])]
                for sample in ("s0", "s1"):
                    totals[(row["CloneID"], group, sample)] += int(row[sample])
            return totals

        self.assertEqual(summary["primer_removal"]["rows_merged_across_allele_groups"], 0)
        self.assertEqual(reads_by_allele(on), reads_by_allele(off))
        self.assertGreater(summary["primer_removal"]["rows_merged"], 5)
        self.assertLess(len(on), len(off))


if __name__ == "__main__":
    unittest.main()
