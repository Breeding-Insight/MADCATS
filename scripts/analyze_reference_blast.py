#!/usr/bin/env python3
"""Summarize reference-genome BLAST hits for generated MADC rows.

The input BLAST table is expected to contain these outfmt 6 fields:
qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


MARKER_RE = re.compile(r"^(Chr\d+)_0*(\d+)$", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--blast", required=True, type=Path)
    parser.add_argument("--row-output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    parser.add_argument("--min-identity", type=float, default=80.0)
    parser.add_argument("--min-coverage", type=float, default=90.0)
    parser.add_argument("--max-evalue", type=float, default=1e-5)
    parser.add_argument("--position-tolerance", type=int, default=125)
    return parser.parse_args()


def marker_coordinate(marker: str) -> tuple[str, int]:
    match = MARKER_RE.match(marker)
    if not match:
        raise ValueError(f"Cannot parse marker coordinate: {marker}")
    return match.group(1).lower(), int(match.group(2))


def cluster_locations(hits: list[dict[str, object]], tolerance: int) -> list[tuple[str, float]]:
    clusters: list[tuple[str, float]] = []
    for hit in sorted(hits, key=lambda x: (str(x["subject"]), float(x["midpoint"]))):
        chrom = str(hit["subject"]).lower()
        midpoint = float(hit["midpoint"])
        for index, (old_chrom, old_midpoint) in enumerate(clusters):
            if chrom == old_chrom and abs(midpoint - old_midpoint) <= tolerance:
                clusters[index] = (old_chrom, (old_midpoint + midpoint) / 2)
                break
        else:
            clusters.append((chrom, midpoint))
    return clusters


def summarize(rows: list[dict[str, object]]) -> dict[str, object]:
    categories = Counter(str(row["reference_class"]) for row in rows)
    total_count = sum(int(row["total_count"]) for row in rows)
    weighted = Counter()
    for row in rows:
        weighted[str(row["reference_class"])] += int(row["total_count"])
    multiple = sum(bool(row["multiple_locations"]) for row in rows)
    multiple_count = sum(
        int(row["total_count"]) for row in rows if bool(row["multiple_locations"])
    )
    return {
        "rows": len(rows),
        "total_read_count": total_count,
        "reference_class_rows": dict(sorted(categories.items())),
        "reference_class_read_counts": dict(sorted(weighted.items())),
        "multiple_location_rows": multiple,
        "multiple_location_read_count": multiple_count,
        "multiple_location_row_fraction": multiple / len(rows) if rows else None,
        "multiple_location_read_fraction": multiple_count / total_count if total_count else None,
    }


def main() -> None:
    args = parse_args()
    with args.metadata.open(newline="") as handle:
        metadata = {row["query_id"]: row for row in csv.DictReader(handle)}

    hits: dict[str, list[dict[str, object]]] = defaultdict(list)
    with args.blast.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 12:
                raise ValueError(f"Expected 12 BLAST columns, found {len(fields)}")
            query_id = fields[0]
            hit = {
                "query_length": int(fields[1]),
                "query_start": int(fields[2]),
                "query_end": int(fields[3]),
                "subject": fields[4],
                "subject_length": int(fields[5]),
                "subject_start": int(fields[6]),
                "subject_end": int(fields[7]),
                "alignment_length": int(fields[8]),
                "coverage": float(fields[9]),
                "identity": float(fields[10]),
                "evalue": float(fields[11]),
            }
            hit["midpoint"] = (int(hit["subject_start"]) + int(hit["subject_end"])) / 2
            if (
                float(hit["identity"]) >= args.min_identity
                and float(hit["evalue"]) <= args.max_evalue
            ):
                hits[query_id].append(hit)

    output: list[dict[str, object]] = []
    for query_id, meta in metadata.items():
        marker_chrom, marker_position = marker_coordinate(meta["marker"])
        candidate_hits = hits.get(query_id, [])
        high_quality_hits = [
            hit for hit in candidate_hits if float(hit["coverage"]) >= args.min_coverage
        ]
        design_hits = [
            hit
            for hit in high_quality_hits
            if str(hit["subject"]).lower() == marker_chrom
            and abs(float(hit["midpoint"]) - marker_position) <= args.position_tolerance
        ]
        locations = cluster_locations(high_quality_hits, args.position_tolerance)

        if candidate_hits:
            best = max(
                candidate_hits,
                key=lambda hit: (
                    float(hit["coverage"]) >= args.min_coverage,
                    float(hit["coverage"]),
                    float(hit["identity"]),
                    -float(hit["evalue"]),
                    int(hit["alignment_length"]),
                ),
            )
            best_is_design = (
                float(best["coverage"]) >= args.min_coverage
                and str(best["subject"]).lower() == marker_chrom
                and abs(float(best["midpoint"]) - marker_position) <= args.position_tolerance
            )
            if best_is_design:
                reference_class = "designed_locus_best"
            elif float(best["coverage"]) >= args.min_coverage:
                reference_class = "off_target_best"
            else:
                reference_class = "partial_hit_only"
        else:
            best = None
            reference_class = "no_qualifying_hit"

        output.append(
            {
                **meta,
                "reference_class": reference_class,
                "best_subject": best["subject"] if best else "",
                "best_midpoint": best["midpoint"] if best else "",
                "best_coverage": best["coverage"] if best else "",
                "best_identity": best["identity"] if best else "",
                "best_evalue": best["evalue"] if best else "",
                "has_designed_locus_hit": bool(design_hits),
                "high_quality_locations": len(locations),
                "multiple_locations": len(locations) > 1,
                "qualifying_hits": len(candidate_hits),
            }
        )

    output.sort(key=lambda row: (str(row["population"]), str(row["marker"]), str(row["sequence"])))
    args.row_output.parent.mkdir(parents=True, exist_ok=True)
    with args.row_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)

    populations = sorted({str(row["population"]) for row in output})
    summary: dict[str, object] = {
        "parameters": {
            "min_identity": args.min_identity,
            "min_coverage": args.min_coverage,
            "max_evalue": args.max_evalue,
            "position_tolerance": args.position_tolerance,
        },
        "by_population": {
            population: summarize(
                [row for row in output if row["population"] == population]
            )
            for population in populations
        },
    }
    extras = [row for row in output if row["population"] == "generated_only"]
    extras_by_count = sorted(extras, key=lambda row: int(row["total_count"]), reverse=True)
    summary["generated_only_top_477_by_count"] = summarize(extras_by_count[:477])
    summary["generated_only_after_top_477"] = summarize(extras_by_count[477:])
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
