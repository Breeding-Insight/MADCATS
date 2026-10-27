#!/usr/bin/env python3
"""Characterize generated-only MADC rows and their raw FASTQ quality support."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from madcats import QUALITY_THRESHOLDS, allele_type, hamming_distance, read_madc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generated", required=True, type=Path)
    parser.add_argument("--dart", required=True, type=Path)
    parser.add_argument("--reference-tags", required=True, type=Path)
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument(
        "--cache-fingerprint",
        help="16-hex cache fingerprint to select when the cache directory holds several",
    )
    parser.add_argument("--reference-audit", type=Path)
    parser.add_argument("--row-output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    return parser.parse_args()


def resolve_cache_files(cache_dir: Path, fingerprint: str | None) -> list[Path]:
    """One cache file per sample; mixing fingerprints would double counts or
    map sequences with the wrong panel order."""
    pattern = f"*.{fingerprint}.json.gz" if fingerprint else "*.json.gz"
    files = sorted(cache_dir.glob(pattern))
    if not files:
        sys.exit(f"No cache files matching {pattern!r} in {cache_dir}")
    by_stem: dict[str, list[Path]] = defaultdict(list)
    for path in files:
        by_stem[path.name.split(".")[0]].append(path)
    duplicates = {stem: paths for stem, paths in by_stem.items() if len(paths) > 1}
    if duplicates:
        stem, paths = sorted(duplicates.items())[0]
        sys.exit(
            f"{cache_dir} holds multiple cache fingerprints for {len(duplicates)} sample(s), "
            f"e.g. {stem}: {', '.join(path.name for path in paths)}. "
            "Pass --cache-fingerprint to select one."
        )
    return files


def row_key(row: dict[str, object]) -> tuple[str, str]:
    return str(row["clone_id"]), str(row["sequence"])


def aggregate(rows: list[dict[str, object]]) -> dict[str, object]:
    final_reads = sum(int(row["total_count"]) for row in rows)
    raw_reads = sum(int(row["raw_exact_count"]) for row in rows)
    result: dict[str, object] = {
        "rows": len(rows),
        "final_read_count": final_reads,
        "raw_exact_read_count": raw_reads,
        "supported_10_reads_in_2_samples_rows": sum(
            bool(row["supported_10_reads_in_2_samples"]) for row in rows
        ),
    }
    for index, threshold in enumerate(QUALITY_THRESHOLDS):
        count = sum(int(row[f"raw_q{threshold}_count"]) for row in rows)
        result[f"raw_q{threshold}_count"] = count
        result[f"raw_q{threshold}_fraction"] = count / raw_reads if raw_reads else None
    return result


def main() -> None:
    args = parse_args()
    generated = read_madc(args.generated)
    dart = read_madc(args.dart)
    generated_by_key = {row_key(row): row for row in generated["records"]}
    dart_by_key = {row_key(row): row for row in dart["records"]}
    dart_by_marker: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in dart["records"]:
        dart_by_marker[str(row["clone_id"])].append(row)

    reference_audit: dict[tuple[str, str], dict[str, str]] = {}
    if args.reference_audit:
        with args.reference_audit.open(newline="") as handle:
            reference_audit = {
                (row["marker"], row["sequence"]): row for row in csv.DictReader(handle)
            }

    # Cache payloads key sequences by locus index, which follows the panel
    # build order: the primer-CSV row order, which the design command writes
    # unchanged into the reference-tag CSV.  A reordered or filtered
    # reference-tag file would silently misattribute markers here.
    with args.reference_tags.open(newline="") as handle:
        marker_order = [row["MarkerName"] for row in csv.DictReader(handle)]
    marker_for_index = {index: marker for index, marker in enumerate(marker_order)}

    raw_counts: Counter[tuple[str, str]] = Counter()
    raw_quality: dict[tuple[str, str], list[int]] = defaultdict(
        lambda: [0] * len(QUALITY_THRESHOLDS)
    )
    for cache_path in resolve_cache_files(args.cache_dir, args.cache_fingerprint):
        with gzip.open(cache_path, "rt") as handle:
            payload = json.load(handle)
        for marker_index, sequence, count in payload["counts"]:
            raw_counts[(marker_for_index[int(marker_index)], sequence)] += int(count)
        for marker_index, sequence, *support in payload["quality_support"]:
            key = (marker_for_index[int(marker_index)], sequence)
            for index, count in enumerate(support):
                raw_quality[key][index] += int(count)

    output: list[dict[str, object]] = []
    for key, generated_row in generated_by_key.items():
        if key in dart_by_key:
            population = "common"
        else:
            population = "generated_only"
        marker, sequence = key
        nearest_rows = sorted(
            dart_by_marker[marker],
            key=lambda row: (hamming_distance(sequence, str(row["sequence"])), str(row["sequence"])),
        )
        nearest = nearest_rows[0] if nearest_rows else None
        nearest_sequence = str(nearest["sequence"]) if nearest else ""
        distance = hamming_distance(sequence, nearest_sequence) if nearest else ""
        mismatch_positions = [
            index + 1
            for index, (left, right) in enumerate(zip(sequence, nearest_sequence))
            if left != right
        ]
        substitutions = [
            f"{right}>{left}"
            for left, right in zip(sequence, nearest_sequence)
            if left != right
        ]
        sample_counts = [int(value) for value in generated_row["counts"]]
        record: dict[str, object] = {
            "marker": marker,
            "sequence": sequence,
            "population": population,
            "allele_type": allele_type(str(generated_row["allele_id"])),
            "total_count": sum(sample_counts),
            "nonzero_samples": sum(value > 0 for value in sample_counts),
            "samples_with_at_least_10_reads": sum(value >= 10 for value in sample_counts),
            "supported_10_reads_in_2_samples": sum(value >= 10 for value in sample_counts) >= 2,
            "nearest_dart_sequence": nearest_sequence,
            "nearest_dart_allele_type": allele_type(str(nearest["allele_id"])) if nearest else "",
            "nearest_dart_total_count": sum(int(value) for value in nearest["counts"]) if nearest else 0,
            "nearest_dart_contains_ambiguity": bool(set(nearest_sequence) - set("ACGT")),
            "nearest_dart_distance": distance,
            "mismatch_positions": ";".join(map(str, mismatch_positions)),
            "substitutions_dart_to_generated": ";".join(substitutions),
            "raw_exact_count": raw_counts[key],
            "reference_class": reference_audit.get(key, {}).get("reference_class", ""),
            "multiple_reference_locations": reference_audit.get(key, {}).get(
                "multiple_locations", ""
            ),
        }
        for index, threshold in enumerate(QUALITY_THRESHOLDS):
            record[f"raw_q{threshold}_count"] = raw_quality[key][index]
        output.append(record)

    output.sort(key=lambda row: (str(row["population"]), -int(row["total_count"]), str(row["marker"])))
    args.row_output.parent.mkdir(parents=True, exist_ok=True)
    with args.row_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)

    extras = [row for row in output if row["population"] == "generated_only"]
    common = [row for row in output if row["population"] == "common"]
    extras_by_count = sorted(extras, key=lambda row: int(row["total_count"]), reverse=True)
    distance_one = [row for row in extras if row["nearest_dart_distance"] == 1]
    position_rows: Counter[int] = Counter()
    position_reads: Counter[int] = Counter()
    substitutions_rows: Counter[str] = Counter()
    substitutions_reads: Counter[str] = Counter()
    for row in distance_one:
        position = int(str(row["mismatch_positions"]))
        substitution = str(row["substitutions_dart_to_generated"])
        position_rows[position] += 1
        position_reads[position] += int(row["total_count"])
        substitutions_rows[substitution] += 1
        substitutions_reads[substitution] += int(row["total_count"])

    reference_support_cross_tab: dict[str, dict[str, int]] = defaultdict(
        lambda: {"rows": 0, "final_read_count": 0}
    )
    for row in extras:
        support_group = (
            "supported_10_reads_in_2_samples"
            if row["supported_10_reads_in_2_samples"]
            else "not_supported_10_reads_in_2_samples"
        )
        reference_class = str(row["reference_class"] or "not_tested")
        key = f"{support_group}|{reference_class}"
        reference_support_cross_tab[key]["rows"] += 1
        reference_support_cross_tab[key]["final_read_count"] += int(row["total_count"])

    category_summary: dict[str, dict[str, int]] = defaultdict(
        lambda: {"rows": 0, "final_read_count": 0, "supported_rows": 0}
    )
    for row in extras:
        category = str(row["allele_type"])
        category_summary[category]["rows"] += 1
        category_summary[category]["final_read_count"] += int(row["total_count"])
        category_summary[category]["supported_rows"] += bool(
            row["supported_10_reads_in_2_samples"]
        )

    summary = {
        "all_common": aggregate(common),
        "all_generated_only": aggregate(extras),
        "generated_only_supported_10_reads_in_2_samples": aggregate(
            [row for row in extras if row["supported_10_reads_in_2_samples"]]
        ),
        "generated_only_not_supported_10_reads_in_2_samples": aggregate(
            [row for row in extras if not row["supported_10_reads_in_2_samples"]]
        ),
        "generated_only_top_477_by_count": aggregate(extras_by_count[:477]),
        "generated_only_distance_one": aggregate(distance_one),
        "generated_only_nearest_dart_contains_ambiguity": aggregate(
            [row for row in extras if row["nearest_dart_contains_ambiguity"]]
        ),
        "generated_only_canonical_nearest_dart_contains_ambiguity": aggregate(
            [
                row
                for row in extras
                if row["allele_type"] in ("Ref", "Alt")
                and row["nearest_dart_contains_ambiguity"]
            ]
        ),
        "generated_only_by_allele_type": dict(sorted(category_summary.items())),
        "generated_only_support_by_reference_class": dict(
            sorted(reference_support_cross_tab.items())
        ),
        "distance_one_mismatch_position_rows": dict(sorted(position_rows.items())),
        "distance_one_mismatch_position_reads": dict(sorted(position_reads.items())),
        "distance_one_substitution_rows": dict(sorted(substitutions_rows.items())),
        "distance_one_substitution_reads": dict(sorted(substitutions_reads.items())),
    }
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
