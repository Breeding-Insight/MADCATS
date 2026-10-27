#!/usr/bin/env python3
"""Characterize DArT MADC rows absent from a generated MADC."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from madc_pipeline import QUALITY_THRESHOLDS, allele_type, compatible_mismatches, read_madc


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
    parser.add_argument("--lower-support-generated", type=Path)
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


def quantile(values: list[int], fraction: float) -> Optional[float]:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def summarize(rows: list[dict[str, object]]) -> dict[str, object]:
    dart_counts = [int(row["dart_total_count"]) for row in rows]
    raw_counts = [int(row["raw_assigned_exact_count"]) for row in rows]
    return {
        "rows": len(rows),
        "dart_read_count": sum(dart_counts),
        "raw_assigned_exact_read_count": sum(raw_counts),
        "zero_dart_count_rows": sum(value == 0 for value in dart_counts),
        "positive_dart_count_rows": sum(value > 0 for value in dart_counts),
        "observed_in_assigned_fastq_rows": sum(value > 0 for value in raw_counts),
        "dart_count_median": quantile(dart_counts, 0.5),
        "dart_count_q75": quantile(dart_counts, 0.75),
        "dart_count_q90": quantile(dart_counts, 0.9),
        "dart_count_q95": quantile(dart_counts, 0.95),
        "dart_count_max": max(dart_counts) if dart_counts else None,
    }


def main() -> None:
    args = parse_args()
    generated = read_madc(args.generated)
    dart = read_madc(args.dart)
    generated_by_key = {
        (row["clone_id"], row["sequence"]): row for row in generated["records"]
    }
    generated_canonical = {
        (row["clone_id"], allele_type(row["allele_id"])): row
        for row in generated["records"]
        if allele_type(row["allele_id"]) in ("Ref", "Alt")
    }
    lower_support_keys: set[tuple[str, str]] = set()
    if args.lower_support_generated:
        lower = read_madc(args.lower_support_generated)
        lower_support_keys = {
            (str(row["clone_id"]), str(row["sequence"])) for row in lower["records"]
        }

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

    q20_index = QUALITY_THRESHOLDS.index(20)
    output: list[dict[str, object]] = []
    for row in dart["records"]:
        key = (str(row["clone_id"]), str(row["sequence"]))
        if key in generated_by_key:
            continue
        category = allele_type(str(row["allele_id"]))
        sequence = str(row["sequence"])
        non_acgt = bool(set(sequence) - set("ACGT"))
        counterpart = generated_canonical.get((key[0], category))
        canonical_iupac_compatible = bool(
            counterpart
            and compatible_mismatches(sequence, str(counterpart["sequence"])) == 0
        )
        raw_count = raw_counts[key]
        raw_q20 = raw_quality[key][q20_index]
        dart_counts = [int(value) for value in row["counts"]]
        if canonical_iupac_compatible:
            explanation = "canonical_iupac_representation"
        elif non_acgt:
            explanation = "noncanonical_ambiguity_consensus"
        elif raw_count > 0 and raw_q20 < 30:
            explanation = "observed_below_q20_discovery_support"
        elif raw_count > 0:
            explanation = "observed_but_not_emitted_other"
        elif sum(dart_counts) == 0:
            explanation = "catalogue_zero_not_observed"
        else:
            explanation = "positive_in_dart_not_in_assigned_exact_tags"
        output.append(
            {
                "marker": key[0],
                "allele_type": category,
                "sequence": sequence,
                "contains_non_acgt": non_acgt,
                "canonical_iupac_compatible": canonical_iupac_compatible,
                "dart_total_count": sum(dart_counts),
                "dart_nonzero_samples": sum(value > 0 for value in dart_counts),
                "dart_samples_with_at_least_10_reads": sum(value >= 10 for value in dart_counts),
                "dart_supported_10_reads_in_2_samples": sum(value >= 10 for value in dart_counts) >= 2,
                "raw_assigned_exact_count": raw_count,
                "raw_assigned_q20_count": raw_q20,
                "recovered_at_lower_q20_support": key in lower_support_keys,
                "explanation": explanation,
            }
        )

    output.sort(key=lambda row: (-int(row["dart_total_count"]), str(row["marker"])))
    args.row_output.parent.mkdir(parents=True, exist_ok=True)
    with args.row_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)

    by_explanation: dict[str, dict[str, object]] = {}
    for explanation in sorted({str(row["explanation"]) for row in output}):
        by_explanation[explanation] = summarize(
            [row for row in output if row["explanation"] == explanation]
        )
    by_category: dict[str, dict[str, object]] = {}
    for category in ("Ref", "Alt", "RefMatch", "AltMatch", "Other"):
        by_category[category] = summarize(
            [row for row in output if row["allele_type"] == category]
        )
    supported = [row for row in output if row["dart_supported_10_reads_in_2_samples"]]
    top = output[:100]
    summary = {
        "all_dart_only": summarize(output),
        "by_explanation": by_explanation,
        "by_category": by_category,
        "dart_supported_10_reads_in_2_samples": summarize(supported),
        "recovered_at_lower_q20_support": summarize(
            [row for row in output if row["recovered_at_lower_q20_support"]]
        ),
        "top_100_by_dart_count": summarize(top),
    }
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
