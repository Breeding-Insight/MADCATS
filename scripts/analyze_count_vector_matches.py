#!/usr/bin/env python3
"""Test whether generated-only rows resemble relabeled DArT count vectors."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from madc_pipeline import hamming_distance, natural_key, read_madc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generated", required=True, type=Path)
    parser.add_argument("--dart", required=True, type=Path)
    parser.add_argument("--row-output", required=True, type=Path)
    parser.add_argument("--summary-output", required=True, type=Path)
    return parser.parse_args()


def pearson(left: list[int], right: list[int]) -> Optional[float]:
    n = len(left)
    sx = sum(left)
    sy = sum(right)
    sxx = sum(value * value for value in left)
    syy = sum(value * value for value in right)
    sxy = sum(x * y for x, y in zip(left, right))
    numerator = n * sxy - sx * sy
    denominator = math.sqrt((n * sxx - sx * sx) * (n * syy - sy * sy))
    return numerator / denominator if denominator else None


def main() -> None:
    args = parse_args()
    generated = read_madc(args.generated)
    dart = read_madc(args.dart)
    generated_sample = {name: index for index, name in enumerate(generated["samples"])}
    dart_sample = {name: index for index, name in enumerate(dart["samples"])}
    samples = sorted(set(generated_sample) & set(dart_sample), key=natural_key)

    dart_keys = {(row["clone_id"], row["sequence"]) for row in dart["records"]}
    dart_by_marker: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in dart["records"]:
        dart_by_marker[str(row["clone_id"])].append(row)

    output: list[dict[str, object]] = []
    for generated_row in generated["records"]:
        key = (generated_row["clone_id"], generated_row["sequence"])
        if key in dart_keys:
            continue
        marker, sequence = str(key[0]), str(key[1])
        counts = [generated_row["counts"][generated_sample[sample]] for sample in samples]
        candidates = []
        for dart_row in dart_by_marker[marker]:
            dart_counts = [dart_row["counts"][dart_sample[sample]] for sample in samples]
            correlation = pearson(counts, dart_counts)
            generated_total = sum(counts)
            dart_total = sum(dart_counts)
            total_fold = (
                max(generated_total, dart_total) / min(generated_total, dart_total)
                if generated_total and dart_total
                else math.inf
            )
            candidates.append(
                {
                    "sequence": str(dart_row["sequence"]),
                    "distance": hamming_distance(sequence, str(dart_row["sequence"])),
                    "correlation": correlation,
                    "total": dart_total,
                    "total_fold": total_fold,
                }
            )
        by_correlation = sorted(
            candidates,
            key=lambda candidate: (
                -(candidate["correlation"] if candidate["correlation"] is not None else -2),
                candidate["total_fold"],
                candidate["distance"],
            ),
        )
        best = by_correlation[0] if by_correlation else None
        nearest = min(candidates, key=lambda candidate: (candidate["distance"], candidate["sequence"]))
        output.append(
            {
                "marker": marker,
                "sequence": sequence,
                "generated_total_count": sum(counts),
                "nonzero_samples": sum(value > 0 for value in counts),
                "best_vector_dart_sequence": best["sequence"] if best else "",
                "best_vector_sequence_distance": best["distance"] if best else "",
                "best_vector_pearson": best["correlation"] if best else "",
                "best_vector_dart_total_count": best["total"] if best else "",
                "best_vector_total_fold_difference": best["total_fold"] if best else "",
                "nearest_sequence_distance": nearest["distance"],
                "nearest_sequence_pearson": nearest["correlation"] or "",
                "nearest_sequence_dart_total_count": nearest["total"],
            }
        )

    output.sort(key=lambda row: (-int(row["generated_total_count"]), str(row["marker"])))
    args.row_output.parent.mkdir(parents=True, exist_ok=True)
    with args.row_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)

    groups = {"all": output, "top_477_by_count": output[:477], "after_top_477": output[477:]}
    summary: dict[str, object] = {}
    for name, rows in groups.items():
        thresholds = Counter()
        reads_by_threshold = Counter()
        for row in rows:
            correlation = row["best_vector_pearson"]
            fold = row["best_vector_total_fold_difference"]
            for label, minimum_correlation, maximum_fold in (
                ("pearson_0.99_fold_1.05", 0.99, 1.05),
                ("pearson_0.95_fold_1.20", 0.95, 1.20),
                ("pearson_0.90_fold_2.00", 0.90, 2.00),
            ):
                if (
                    isinstance(correlation, float)
                    and correlation >= minimum_correlation
                    and isinstance(fold, float)
                    and fold <= maximum_fold
                ):
                    thresholds[label] += 1
                    reads_by_threshold[label] += int(row["generated_total_count"])
        summary[name] = {
            "rows": len(rows),
            "generated_read_count": sum(int(row["generated_total_count"]) for row in rows),
            "plausible_relabel_rows": dict(thresholds),
            "plausible_relabel_generated_read_counts": dict(reads_by_threshold),
        }
    args.summary_output.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
