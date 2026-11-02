"""``madcats filter-paralogs``: clean a finished MADC with a reference genome.

In panels where paralogous copies co-amplify with a target, a paralog read
that differs from the Ref tag only at non-target bases becomes a RefMatch or
AltMatch row of that marker.  This module places each marker on the genome
by its Ref tag, aligns every haplotype, and drops or flags the rows whose
best alignment lies at another locus (paralogs) or ties between loci
(ambiguous).

It is kept apart from ``madcats.py`` because it is the only command with a
dependency outside the standard library: the ``mappy`` aligner (minimap2),
imported on first use (``pip install 'madcats[genome]'``), so ``madcats``
runs without it.  The MADC reading helpers it shares with ``validate``
(header detection, sample columns, allele categories, the count-overflow
check) stay in ``madcats``.

Reading map: scoring and alignment (``hit_score``, ``GenomeAligner``,
``SiteIndex``, ``map_sequences``); marker placement and labels
(``place_marker``, ``genome_label``); the two-pass table filter
(``filter_paralogs``); the command line (``add_filter_paralogs_command``,
``run_filter_paralogs``).
"""

from __future__ import annotations

import argparse
import bisect
import concurrent.futures
import csv
import datetime as dt
import re
import sys
import threading
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, NamedTuple, Optional, Sequence, Set, Tuple

from madcats import (
    ALLELE_TYPE_ORDER,
    CountOverflowCheck,
    MadcError,
    __version__,
    _default_jobs,
    _write_json,
    allele_type,
    is_madc_header,
    madc_identity_columns,
    madc_sample_start,
    open_text,
    positive_int,
)

GENOME_LABEL_COLUMN = "PICRef"     # the vendor statistics column the GenomeLabel flag replaces
GENOME_LABEL_HEADER = "GenomeLabel"
GENOME_LABELS = ("on_target", "paralog", "ambiguous", "unmapped", "unplaced")
GENOME_ROW_OUTCOMES = ("kept", "flagged", "protected", "dropped")
# minimap2 -x sr scoring: match, mismatch, gap open/extend, long-gap open/extend, ambiguous base.
SR_SCORING = (2, 8, 12, 2, 24, 1, 1)
# Map-time sensitivity on top of -x sr (the index is unchanged): single-seed chains and up to 50
# hits, so an 81-bp allele with several SNPs still finds its copies.  On lettuce this lifted
# drop/keep agreement with bwa mem from 97.3% to 98.2%.
SR_MAP_OPTIONS = {"min_cnt": 1, "min_chain_score": 20, "best_n": 50}
GENOME_SITES_KEPT = 8              # distinct genome sites remembered per sequence, best first
CLONE_POSITION_SLACK = 10          # a CloneID position this far outside the placed Ref tag still agrees
CLONE_POSITION_RE = re.compile(r"(.+)_0*(\d+)")
SITE_PAD = 20                      # genome bases either side of a marker site in the site realignment
SITE_INDEX_K, SITE_INDEX_W = 11, 3 # the site realignment's seeds: short, so diverged alleles still seed
GENOME_LABEL_SIDECAR_COLUMNS = [
    "CloneID", "AlleleID", "AlleleSequence", "Label", "Action",
    "TargetScore", "TargetFrom", "OffTargetScore", "OffTargetSites", "OffTargetHit",
]
GENOME_MARKER_SIDECAR_COLUMNS = ["CloneID", "Placement", "RefTagSites", "CloneIDCheck", "Sites"]


class GenomeHit(NamedTuple):
    """One alignment to the genome: a 0-based, half-open reference interval and its score."""

    contig: str
    start: int
    end: int
    strand: int
    score: int


class GenomePlacement(NamedTuple):
    """A sequence's distinct genome sites, best score first (hits overlapping a better one are merged into it)."""

    sites: Tuple[GenomeHit, ...]

    @property
    def best(self) -> Tuple[GenomeHit, ...]:
        """The equal-best sites."""
        return tuple(site for site in self.sites if site.score == self.sites[0].score)


class MarkerPlacement(NamedTuple):
    """Where a marker's own tag lies in the genome.

    ``source``: ``ref_tag`` (the Ref tag has one best site), ``ref_tag+clone_id``
    (several equal-best sites; the CloneID position picks one),
    ``ref_tag_multicopy`` (several, none picked: every site is the marker's),
    ``clone_id`` (the Ref tag did not align; a window around the CloneID
    position stands in) or ``unplaced``.  ``clone_check`` compares the CloneID
    position with the Ref-tag sites: ``agree``, ``disagree``, ``no_position``
    (the CloneID is not CHROM_POS), ``contig_absent``, ``fallback``, or
    ``ignored`` (``--clone-positions ignore``).
    """

    sites: Tuple[GenomeHit, ...]
    source: str
    ref_tag_sites: int
    clone_check: str


class GenomeLabel(NamedTuple):
    """A haplotype's label and the two scores that decided it."""

    label: str
    target_score: Optional[int]       # best alignment at the marker's own site
    target_from: str                  # genome (the genome-wide search found it), site (realigned there) or '.'
    off_target_score: Optional[int]   # best alignment anywhere else
    off_target_sites: int             # distinct sites elsewhere at that score
    off_target_hit: Optional[GenomeHit]


def hit_score(
    cigar: Sequence[Sequence[int]], edit_distance: int, matches: int, block_length: int, scoring: Sequence[int] = SR_SCORING
) -> int:
    """An alignment's DP score recomputed from what mappy reports (it gives neither AS nor XS).

    minimap2's block length counts aligned columns and gap bases but not
    ambiguous columns, and NM counts all three kinds of difference, so
    mismatches = block - matches - gap bases and ambiguous = NM - (block -
    matches).  Gaps cost the cheaper of minimap2's two affine pieces.
    """
    match, mismatch, gap_open, gap_extend, long_open, long_extend, ambiguous = scoring
    gap_bases = gap_cost = 0
    for length, operation in cigar:
        if operation in (1, 2):  # I, D
            gap_bases += length
            gap_cost += min(gap_open + gap_extend * length, long_open + long_extend * length)
    differences = block_length - matches
    mismatches = max(0, differences - gap_bases)
    ambiguous_bases = max(0, edit_distance - differences)
    return match * matches - mismatch * mismatches - ambiguous * ambiguous_bases - gap_cost


def summarize_hits(hits: Iterable[GenomeHit], keep: int = GENOME_SITES_KEPT) -> GenomePlacement:
    """A sequence's distinct sites, best first.  A hit overlapping a better-scoring
    one on the genome is the same site (a tandem repeat inside the tag can yield two)."""
    sites: List[GenomeHit] = []
    for hit in sorted(hits, key=lambda hit: (-hit.score, hit.contig, hit.start)):
        if any(site.contig == hit.contig and hit.start < site.end and site.start < hit.end for site in sites):
            continue
        sites.append(hit)
        if len(sites) == keep:
            break
    return GenomePlacement(tuple(sites))


def _load_mappy():
    try:
        import mappy
    except ImportError as exc:
        raise MadcError(
            "filter-paralogs needs the optional minimap2 aligner: pip install 'madcats[genome]' (or pip install mappy)"
        ) from exc
    return mappy


class GenomeAligner:
    """minimap2's short-read preset, through mappy, with the scoring pinned to SR_SCORING."""

    def __init__(self, aligner, mappy_module) -> None:
        self._aligner = aligner
        self._mappy = mappy_module
        self.contigs = frozenset(aligner.seq_names)
        self.description = {
            "mappy": getattr(mappy_module, "__version__", None),
            "preset": "sr",
            "k": aligner.k,
            "w": aligner.w,
            "scoring": list(SR_SCORING),
            "map_options": dict(SR_MAP_OPTIONS),
            "site_realignment": {"k": SITE_INDEX_K, "w": SITE_INDEX_W, "pad": SITE_PAD},
        }

    @classmethod
    def open(cls, genome: Optional[Path], index: Optional[Path], threads: int) -> "GenomeAligner":
        """Load ``index`` when it exists; otherwise index ``genome`` (saving to ``index`` when given)."""
        mappy = _load_mappy()
        if index is not None and index.exists():
            source, save = index, None
        elif genome is not None:
            source, save = genome, index
        else:
            raise MadcError(f"Genome index not found: {index} (pass --genome to build it)")
        if not source.exists():
            raise MadcError(f"Genome file not found: {source}")
        if save is not None:
            save.parent.mkdir(parents=True, exist_ok=True)
        verb = "Loading index" if source == index else "Indexing genome"
        print(f"{verb} {source}" + (f" (saving {save})" if save else ""), file=sys.stderr)
        aligner = mappy.Aligner(
            str(source), preset="sr", n_threads=threads, fn_idx_out=str(save) if save else None,
            scoring=list(SR_SCORING), **SR_MAP_OPTIONS,
        )
        if not aligner:
            raise MadcError(f"Could not load or build a minimap2 index from {source}")
        if (aligner.k, aligner.w) != (21, 11):
            print(
                f"warning: {source} was indexed with k={aligner.k}, w={aligner.w}; "
                "short-read indexes (minimap2 -x sr) use k=21, w=11",
                file=sys.stderr,
            )
        return cls(aligner, mappy)

    def new_buffer(self):
        return self._mappy.ThreadBuffer()

    def hits(self, sequence: str, buffer) -> List[GenomeHit]:
        return [
            GenomeHit(hit.ctg, hit.r_st, hit.r_en, hit.strand, hit_score(hit.cigar, hit.NM, hit.mlen, hit.blen))
            for hit in self._aligner.map(sequence, buf=buffer)
        ]

    def fetch(self, contig: str, start: int, end: int) -> str:
        return self._aligner.seq(contig, start, end) or ""

    def site_index(self, segments: Sequence[str]) -> "SiteIndex":
        return SiteIndex(self._mappy, segments)


class SiteIndex:
    """A sensitive index of the marker sites alone, so an allele too diverged
    to seed genome-wide (k=21) is still scored at its own marker.  Segments
    are joined with N spacers; hits come back as (segment number, score)."""

    def __init__(self, mappy_module, segments: Sequence[str]) -> None:
        spacer = "N" * 100
        self._mappy = mappy_module
        self._starts: List[int] = []
        offset = 0
        for segment in segments:
            self._starts.append(offset)
            offset += len(segment) + len(spacer)
        self._aligner = mappy_module.Aligner(
            seq=spacer.join(segments) or "N", preset="sr", k=SITE_INDEX_K, w=SITE_INDEX_W, scoring=list(SR_SCORING)
        )

    def new_buffer(self):
        return self._mappy.ThreadBuffer()

    def hits(self, sequence: str, buffer) -> List[Tuple[int, int]]:
        found = []
        for hit in self._aligner.map(sequence, buf=buffer):
            segment = bisect.bisect_right(self._starts, hit.r_st) - 1
            found.append((segment, hit_score(hit.cigar, hit.NM, hit.mlen, hit.blen)))
        return found


def map_sequences(aligner, sequences: Sequence[str], jobs: int = 1, summarize=summarize_hits, batch: int = 2000) -> Dict[str, object]:
    """Every sequence's hits, summarized.  Worker threads share the index; mappy
    releases the GIL while it maps, and each thread keeps its own buffer."""
    local = threading.local()

    def place(chunk: Sequence[str]) -> List[object]:
        buffer = getattr(local, "buffer", None)
        if buffer is None:
            buffer = local.buffer = aligner.new_buffer()
        return [summarize(aligner.hits(sequence, buffer)) for sequence in chunk]

    chunks = [sequences[start : start + batch] for start in range(0, len(sequences), batch)]
    results: Dict[str, object] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        for chunk, summaries in zip(chunks, pool.map(place, chunks)):
            results.update(zip(chunk, summaries))
    return results


def clone_position(clone_id: str, contigs: Set[str]) -> Tuple[Optional[Tuple[str, int]], str]:
    """The CHROM_POS a CloneID names (1-based; zero padding allowed), or why it names none."""
    match = CLONE_POSITION_RE.fullmatch(clone_id)
    if not match:
        return None, "no_position"
    if match.group(1) not in contigs:
        return None, "contig_absent"
    return (match.group(1), int(match.group(2))), "agree"


def _site_covers(site: GenomeHit, contig: str, position: int) -> bool:
    return site.contig == contig and site.start - CLONE_POSITION_SLACK < position <= site.end + CLONE_POSITION_SLACK


def place_marker(
    clone_id: str,
    ref_placement: Optional[GenomePlacement],
    contigs: Set[str],
    tag_length: int,
    use_clone_position: bool = True,
) -> MarkerPlacement:
    """Locate a marker's tag: its Ref tag's alignment first, the CloneID position as cross-check and fallback.

    Alignment wins a disagreement, since CloneID coordinates may come from
    another assembly version; the disagreement is reported.  A Ref tag with
    several equal-best sites is resolved by the CloneID position when exactly
    one site covers it; otherwise all its sites count as the marker's own.
    Without ``use_clone_position`` (a genome other than the one the CloneIDs
    refer to) the position is never read: an unaligned Ref tag leaves the
    marker unplaced, and equal-best sites stay multi-copy.
    """
    position, status = clone_position(clone_id, contigs) if use_clone_position else (None, "ignored")
    copies = ref_placement.best if ref_placement is not None and ref_placement.sites else ()
    if not copies:
        if position is None:
            return MarkerPlacement((), "unplaced", 0, status)
        contig, base = position
        window = GenomeHit(contig, max(0, base - 1 - tag_length), base + tag_length, 0, 0)
        return MarkerPlacement((window,), "clone_id", 0, "fallback")
    agreeing = [site for site in copies if position is not None and _site_covers(site, *position)]
    check = status if position is None else ("agree" if agreeing else "disagree")
    if len(copies) == 1:
        return MarkerPlacement(copies, "ref_tag", 1, check)
    if len(agreeing) == 1:
        return MarkerPlacement(tuple(agreeing), "ref_tag+clone_id", len(copies), check)
    return MarkerPlacement(copies, "ref_tag_multicopy", len(copies), check)


def _sites_overlap(hit: GenomeHit, site: GenomeHit) -> bool:
    """Same locus: the intervals share at least half of the shorter one."""
    if hit.contig != site.contig:
        return False
    overlap = min(hit.end, site.end) - max(hit.start, site.start)
    return overlap > 0 and 2 * overlap >= min(hit.end - hit.start, site.end - site.start)


def genome_label(
    placement: Optional[GenomePlacement], marker: MarkerPlacement, site_score: Optional[int] = None
) -> GenomeLabel:
    """Label a haplotype by its best score at its marker's own site against its best score anywhere else.

    ``on_target``: the marker's site scores strictly best.  ``paralog``:
    somewhere else scores strictly better (or the haplotype aligns only
    elsewhere).  ``ambiguous``: a tie between the marker's site and another
    site, or a marker whose own Ref tag is multi-copy and unresolved.
    ``unmapped``: no alignment.  ``unplaced``: the marker itself could not be
    located, so the haplotype cannot be judged.  ``site_score`` is the
    sensitive realignment at the marker's site, used when the genome-wide
    search found no hit there.
    """
    sites = placement.sites if placement is not None else ()
    if not marker.sites:
        off = sites
        target_score, target_from = None, "."
    else:
        target = [hit for hit in sites if any(_sites_overlap(hit, site) for site in marker.sites)]
        off = tuple(hit for hit in sites if hit not in target)
        if target:
            target_score, target_from = target[0].score, "genome"
        elif site_score is not None:
            target_score, target_from = site_score, "site"
        else:
            target_score, target_from = None, "."
    off_score = off[0].score if off else None
    off_sites = sum(hit.score == off_score for hit in off)
    detail = (target_score, target_from, off_score, off_sites, off[0] if off else None)
    if not marker.sites:
        label = "unplaced"
    elif target_score is None and off_score is None:
        label = "unmapped"
    elif target_score is None or (off_score is not None and off_score > target_score):
        label = "paralog"
    elif off_score == target_score or len(marker.sites) > 1:
        label = "ambiguous"
    else:
        label = "on_target"
    return GenomeLabel(label, *detail)


def _site_segments(aligner, markers: Mapping[str, MarkerPlacement]) -> Tuple[List[str], Dict[str, List[int]]]:
    """The genome sequence around every marker site, and each marker's segment numbers."""
    segments: List[str] = []
    owned: Dict[str, List[int]] = defaultdict(list)
    for clone, marker in markers.items():
        for site in marker.sites:
            owned[clone].append(len(segments))
            segments.append(aligner.fetch(site.contig, max(0, site.start - SITE_PAD), site.end + SITE_PAD))
    return segments, owned


def _best_by_segment(hits: Iterable[Tuple[int, int]]) -> Dict[int, int]:
    best: Dict[int, int] = {}
    for segment, score in hits:
        best[segment] = max(score, best.get(segment, score))
    return best


def _madc_stream(handle, path: Path):
    """(preamble rows, header row, reader positioned at the first data row) of an open MADC."""
    reader = csv.reader(handle)
    preamble = []
    for row in reader:
        if is_madc_header(row):
            return preamble, row, reader
        preamble.append(row)
    raise MadcError(f"Could not find an AlleleID header row in MADC file: {path}")


def _row_counts(cells: Sequence[str], path: Path) -> List[int]:
    try:
        return list(map(int, cells))
    except ValueError:
        counts = []
        for value in cells:
            try:
                counts.append(int(float(value or 0)))
            except ValueError as exc:
                raise MadcError(f"Non-numeric read count {value!r} in {path}") from exc
        return counts


def _format_site(site: GenomeHit) -> str:
    strand = {1: "+", -1: "-"}.get(site.strand, ".")
    return f"{site.contig}:{site.start + 1}-{site.end}({strand})"


def filter_paralogs(
    madc_path: Path,
    output: Path,
    aligner,
    labels_path: Path,
    markers_path: Path,
    paralog_action: str = "drop",
    ambiguous_action: str = "drop",
    jobs: int = 1,
    provenance: str = "",
    clone_positions: str = "use",
) -> Dict[str, object]:
    """Label every MADC haplotype by where the genome places it; drop or flag paralogous rows.

    Pass 1 reads only row identities.  Every distinct haplotype is aligned to
    the genome once; each marker is located from its Ref row's tag
    (``place_marker``); haplotypes with no genome-wide hit at their marker are
    realigned to the marker sites alone (``SiteIndex``), and each haplotype
    is labelled by ``genome_label``.  Pass 2 streams the table again and
    writes the kept rows unchanged, so a vendor MADC keeps its preamble,
    statistics and column layout, and a HapApp fixed-allele-ID MADC keeps its
    Code_version column.  Canonical Ref and Alt rows are never dropped, since
    downstream dosage tools need both designed alleles; they are reported as
    ``protected`` instead.  In flag mode the PICRef column becomes
    ``GenomeLabel`` on every row.  The labels sidecar has one line per input
    row, dropped rows included.  ``clone_positions="ignore"`` keeps CloneID
    coordinates out of marker placement, for a genome other than the one they
    refer to.
    """
    actions = {"paralog": paralog_action, "ambiguous": ambiguous_action}
    flagging = "flag" in actions.values()

    # Pass 1: row identities, each marker's Ref tag, and the label column.
    haplotypes: Dict[Tuple[str, str], None] = {}
    ref_tags: Dict[str, str] = {}
    tag_lengths: Dict[str, int] = {}
    with open_text(madc_path, "rt") as handle:
        _, header, reader = _madc_stream(handle, madc_path)
        allele_index, clone_index, sequence_index = madc_identity_columns(header, madc_path)
        width = max(allele_index, clone_index, sequence_index) + 1
        label_index = None
        if flagging:
            stripped = [column.strip() for column in header]
            for name in (GENOME_LABEL_HEADER, GENOME_LABEL_COLUMN):
                if name in stripped:
                    label_index = stripped.index(name)
                    break
            else:
                raise MadcError(
                    f"{madc_path} has no {GENOME_LABEL_COLUMN} column to carry the flag; use drop mode, "
                    "or read the labels sidecar"
                )
        for row in reader:
            if len(row) < width or not row[allele_index].strip():
                continue
            clone = row[clone_index].strip()
            sequence = row[sequence_index].strip().upper()
            haplotypes[(clone, sequence)] = None
            tag_lengths.setdefault(clone, len(sequence))
            if clone not in ref_tags and allele_type(row[allele_index].strip()) == "Ref":
                ref_tags[clone] = sequence

    def query(sequence: str) -> str:
        return sequence.replace("-", "")

    sequences = sorted({query(sequence) for _, sequence in haplotypes})
    print(f"Aligning {len(sequences):,} distinct haplotype sequences on {jobs} thread(s)", file=sys.stderr)
    started = time.time()
    placements = map_sequences(aligner, sequences, jobs)
    markers: Dict[str, MarkerPlacement] = {}
    for clone in tag_lengths:
        ref = ref_tags.get(clone)
        markers[clone] = place_marker(
            clone, placements.get(query(ref)) if ref is not None else None, aligner.contigs, tag_lengths[clone],
            use_clone_position=clone_positions == "use",
        )
    labels = {key: genome_label(placements.get(query(key[1])), markers[key[0]]) for key in haplotypes}

    # Realign, at the marker sites alone, every haplotype the genome-wide search found nowhere near its marker.
    missing = [key for key, result in labels.items() if result.target_score is None and markers[key[0]].sites]
    rescued = 0
    if missing:
        segments, owned = _site_segments(aligner, markers)
        site_hits = map_sequences(
            aligner.site_index(segments), sorted({query(sequence) for _, sequence in missing}), jobs, _best_by_segment
        )
        for clone, sequence in missing:
            scores = site_hits[query(sequence)]
            own = [scores[segment] for segment in owned[clone] if segment in scores]
            if own:
                labels[(clone, sequence)] = genome_label(placements.get(query(sequence)), markers[clone], max(own))
                rescued += 1
    alignment_seconds = round(time.time() - started, 1)

    markers_path.parent.mkdir(parents=True, exist_ok=True)
    with markers_path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(GENOME_MARKER_SIDECAR_COLUMNS)
        for clone, marker in markers.items():
            sites = ";".join(_format_site(site) for site in marker.sites) or "."
            writer.writerow([clone, marker.source, marker.ref_tag_sites, marker.clone_check, sites])

    # Pass 2: stream the rows again, writing the kept ones and every row's label.
    def tally():
        return defaultdict(lambda: {"rows": 0, "reads": 0})

    def cell(value) -> object:
        return "." if value is None else value

    by_label, by_outcome = tally(), tally()
    by_category: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(tally)
    overflow = CountOverflowCheck()
    totals = {"input": {"rows": 0, "reads": 0}, "output": {"rows": 0, "reads": 0}}
    output.parent.mkdir(parents=True, exist_ok=True)
    labels_path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(madc_path, "rt") as handle, open_text(output, "wt") as out_handle, labels_path.open(
        "w", newline=""
    ) as label_handle:
        preamble, header, reader = _madc_stream(handle, madc_path)
        writer = csv.writer(out_handle, lineterminator="\n")
        label_writer = csv.writer(label_handle, lineterminator="\n")
        label_writer.writerow(GENOME_LABEL_SIDECAR_COLUMNS)
        if provenance and preamble and preamble[0] and preamble[0][0].startswith("* MADCATS"):
            preamble[0] = [f"{preamble[0][0]} | {provenance}"] + preamble[0][1:]
        writer.writerows(preamble)
        header = list(header)
        if label_index is not None:
            header[label_index] = GENOME_LABEL_HEADER
        writer.writerow(header)
        sample_start = madc_sample_start(header)
        for row in reader:
            if len(row) < width or not row[allele_index].strip():
                continue
            allele_id = row[allele_index].strip()
            clone = row[clone_index].strip()
            sequence = row[sequence_index].strip().upper()
            result = labels[(clone, sequence)]
            category = allele_type(allele_id)
            action = actions.get(result.label, "keep")
            if action == "drop":
                outcome = "protected" if category in ("Ref", "Alt") else "dropped"
            else:
                outcome = "flagged" if action == "flag" else "kept"
            counts = _row_counts(row[sample_start:], madc_path)
            overflow.add(clone, counts)
            reads = sum(counts)
            for bucket in (totals["input"], by_label[result.label], by_outcome[outcome], by_category[category][result.label]):
                bucket["rows"] += 1
                bucket["reads"] += reads
            if outcome != "dropped":
                if label_index is not None:
                    row = list(row) + [""] * max(0, label_index + 1 - len(row))
                    row[label_index] = result.label
                writer.writerow(row)
                totals["output"]["rows"] += 1
                totals["output"]["reads"] += reads
            off_hit = _format_site(result.off_target_hit) if result.off_target_hit is not None else "."
            label_writer.writerow([
                clone, allele_id, sequence, result.label, outcome, cell(result.target_score), result.target_from,
                cell(result.off_target_score), result.off_target_sites, off_hit,
            ])

    warning = overflow.warning(madc_path)
    if warning:
        print(warning, file=sys.stderr)
    sources = Counter(marker.source for marker in markers.values())
    checks = Counter(marker.clone_check for marker in markers.values())
    checked = checks["agree"] + checks["disagree"]
    if checked and checks["disagree"] > 0.1 * checked:
        print(
            f"warning: CloneID positions disagree with the Ref-tag placement for {checks['disagree']:,} of "
            f"{checked:,} markers, so they may refer to another assembly. If they do, rerun with "
            "--clone-positions ignore so they are not used to place markers",
            file=sys.stderr,
        )
    if markers and not sources.keys() - {"unplaced"}:
        print("warning: no marker could be placed on the genome; every row was kept", file=sys.stderr)
    return {
        "markers": {
            "total": len(markers),
            "placement": {source: sources[source] for source in sorted(sources)},
            "clone_id_check": {check: checks[check] for check in sorted(checks)},
            "multicopy_ref_tags": [clone for clone, marker in markers.items() if marker.ref_tag_sites > 1],
            "clone_id_disagreements": [clone for clone, marker in markers.items() if marker.clone_check == "disagree"],
            "unplaced": [clone for clone, marker in markers.items() if marker.source == "unplaced"],
        },
        "haplotypes": {
            "distinct": len(haplotypes),
            "distinct_sequences": len(sequences),
            "labels": dict(Counter(result.label for result in labels.values())),
            "realigned_at_site": {"haplotypes": len(missing), "found_at_site": rescued},
            "alignment_seconds": alignment_seconds,
        },
        "rows": {
            **totals,
            "by_label": {label: dict(by_label[label]) for label in GENOME_LABELS if label in by_label},
            "by_outcome": {outcome: dict(by_outcome[outcome]) for outcome in GENOME_ROW_OUTCOMES if outcome in by_outcome},
            "by_category": {
                category: {label: dict(counts[label]) for label in GENOME_LABELS if label in counts}
                for category, counts in sorted(by_category.items(), key=lambda item: ALLELE_TYPE_ORDER.get(item[0], 99))
            },
        },
        "genome_label_column": GENOME_LABEL_HEADER if label_index is not None else None,
        "count_overflow": overflow.summary(),
    }


def run_filter_paralogs(args: argparse.Namespace) -> int:
    """Clean a finished MADC with a reference genome (see ``filter_paralogs``)."""
    if not args.genome and not args.genome_index:
        raise MadcError("filter-paralogs needs --genome, --genome-index, or both")
    madc = Path(args.madc)
    output = Path(args.output)
    if not madc.is_file():
        raise MadcError(f"MADC not found: {madc}")
    if output.resolve() == madc.resolve():
        raise MadcError("--output must be a new file, not the input MADC")
    genome = Path(args.genome) if args.genome else None
    index = Path(args.genome_index) if args.genome_index else None
    started = time.time()
    aligner = GenomeAligner.open(genome, index, args.jobs)
    labels_path = output.with_suffix(output.suffix + ".genome_labels.csv")
    markers_path = output.with_suffix(output.suffix + ".genome_markers.csv")
    provenance = (
        f"filter-paralogs v{__version__} paralogs={args.paralogs} ambiguous={args.ambiguous} "
        f"clone_positions={args.clone_positions}"
    )
    summary = filter_paralogs(
        madc, output, aligner, labels_path, markers_path,
        paralog_action=args.paralogs, ambiguous_action=args.ambiguous, jobs=args.jobs, provenance=provenance,
        clone_positions=args.clone_positions,
    )
    qc_path = Path(args.qc_output) if args.qc_output else output.with_suffix(output.suffix + ".qc.json")
    qc = {
        "pipeline": "MADCATS",
        "command": "filter-paralogs",
        "version": __version__,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "inputs": {
            "madc": str(madc.resolve()),
            "genome": str(genome.resolve()) if genome else None,
            "genome_index": str(index.resolve()) if index else None,
        },
        "parameters": {
            "paralogs": args.paralogs, "ambiguous": args.ambiguous, "clone_positions": args.clone_positions,
            "jobs": args.jobs,
        },
        "aligner": aligner.description,
        "outputs": {"madc": str(output.resolve()), "labels": str(labels_path), "markers": str(markers_path)},
        "elapsed_seconds": round(time.time() - started, 1),
        **summary,
    }
    _write_json(qc_path, qc)
    rows = summary["rows"]
    dropped = rows["input"]["reads"] - rows["output"]["reads"]
    share = dropped / rows["input"]["reads"] if rows["input"]["reads"] else 0.0
    print(
        f"Kept {rows['output']['rows']:,} of {rows['input']['rows']:,} rows; dropped {dropped:,} reads ({share:.1%})",
        file=sys.stderr,
    )
    print(f"Wrote {output}, its labels {labels_path}, and QC to {qc_path}", file=sys.stderr)
    return 0


def add_filter_paralogs_command(subparsers) -> None:
    filtering = subparsers.add_parser(
        "filter-paralogs",
        help="drop or flag MADC rows that a reference genome places at another locus (needs mappy)",
        description=(
            "Clean a MADC with a reference genome, for panels where paralogous copies of a target co-amplify "
            "with it. Each marker is placed by aligning its Ref tag, with the CloneID position as a "
            "cross-check. Every haplotype is then aligned and labelled on_target, paralog (it aligns better "
            "elsewhere) or ambiguous (it aligns equally well at the marker and elsewhere). Paralog and "
            "ambiguous rows are dropped or flagged; canonical Ref and Alt rows are always kept. Reads MADCs "
            "from generate, vendor DArT MADCs and HapApp fixed-allele-ID MADCs. Needs the mappy aligner: "
            "pip install 'madcats[genome]'."
        ),
    )
    filtering.add_argument(
        "--madc",
        required=True,
        help=(
            "MADC to clean: from generate, a vendor DArT MADC (with or without its starred preamble), or a "
            "HapApp fixed-allele-ID MADC"
        ),
    )
    filtering.add_argument("--genome", help="reference genome FASTA (optionally gzipped)")
    filtering.add_argument(
        "--genome-index",
        help=(
            "minimap2 index (.mmi) built with -x sr; loaded when it exists, otherwise built from --genome and "
            "saved here for the next run"
        ),
    )
    filtering.add_argument(
        "--clone-positions",
        choices=("use", "ignore"),
        default="use",
        help=(
            "CHROM_POS CloneIDs: use them to check placements, choose between equal-best Ref-tag sites and stand "
            "in for Ref tags that do not align, or ignore them when they refer to another assembly than the "
            "genome given (default: use)"
        ),
    )
    filtering.add_argument("--output", required=True, help="cleaned MADC to write")
    filtering.add_argument(
        "--paralogs",
        choices=("drop", "flag"),
        default="drop",
        help="paralog rows: drop them, or keep them labelled in the GenomeLabel column (default: drop)",
    )
    filtering.add_argument(
        "--ambiguous",
        choices=("drop", "flag", "keep"),
        default="drop",
        help="rows with equal-best hits at the marker and elsewhere: drop, flag, or keep unlabelled (default: drop)",
    )
    filtering.add_argument("--qc-output", help="QC JSON (default: OUTPUT.qc.json)")
    filtering.add_argument(
        "--jobs", type=positive_int, default=_default_jobs(), help="alignment threads (default: up to 4)"
    )
    filtering.set_defaults(func=run_filter_paralogs)
