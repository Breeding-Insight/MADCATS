#!/usr/bin/env python3
"""Validate an infer-design panel against a ground-truth reference-tags CSV.

Compares per-locus Ref/Alt tags, variant positions, and Ref/Alt orientation,
and — when the inference QC report is supplied — breaks accuracy down by the
route (exact / fuzzy) and flags each locus was anchored with, so the risky
tiers can be judged separately.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from madcats import load_reference_tags


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inferred", required=True, type=Path, help="reference-tags CSV from infer-design")
    parser.add_argument("--truth", required=True, type=Path, help="reference-tags CSV from the design command")
    parser.add_argument("--inference-report", type=Path, help="infer-design QC JSON for per-route breakdown")
    parser.add_argument("--output", type=Path, help="optional JSON report to write")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inferred = {l.marker: l for l in load_reference_tags(args.inferred)}
    truth = {l.marker: l for l in load_reference_tags(args.truth)}
    routes: dict[str, dict] = {}
    if args.inference_report:
        with args.inference_report.open() as handle:
            report = json.load(handle)
        routes = {entry["marker"]: entry for entry in report.get("per_locus", [])}

    shared = sorted(set(inferred) & set(truth))
    verdicts: dict[str, str] = {}
    for marker in shared:
        got, want = inferred[marker], truth[marker]
        if got.ref_tag == want.ref_tag and got.alt_tag == want.alt_tag:
            verdicts[marker] = (
                "correct" if got.variant_pos == want.variant_pos else "correct_tags_wrong_position"
            )
        elif got.ref_tag == want.alt_tag and got.alt_tag == want.ref_tag:
            verdicts[marker] = "ref_alt_swapped"
        else:
            verdicts[marker] = "wrong_window"
    anchors_ok = sum(
        1 for marker in shared if truth[marker].ref_tag.startswith(inferred[marker].left_primer)
    )

    by_route: dict[str, Counter] = defaultdict(Counter)
    for marker, verdict in verdicts.items():
        route = routes.get(marker, {}).get("route", "unknown")
        if "multi_window" in routes.get(marker, {}).get("flags", []):
            route += "+multi_window_flag"
        by_route[route][verdict] += 1

    summary = {
        "loci": {"truth": len(truth), "inferred": len(inferred), "shared": len(shared)},
        "missing_from_inferred": sorted(set(truth) - set(inferred))[:20],
        "missing_count": len(set(truth) - set(inferred)),
        "extra_in_inferred": sorted(set(inferred) - set(truth)),
        "verdicts": dict(Counter(verdicts.values())),
        "anchor_is_prefix_of_true_ref_tag": anchors_ok,
        "by_route": {route: dict(counts) for route, counts in sorted(by_route.items())},
        "wrong_loci": sorted(
            marker for marker, verdict in verdicts.items() if verdict != "correct"
        ),
    }
    print(json.dumps(summary, indent=2))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(summary, indent=2) + "\n")

    correct = summary["verdicts"].get("correct", 0)
    print(
        f"\n{correct}/{len(shared)} shared loci fully correct "
        f"({correct/len(shared):.2%}); {summary['missing_count']} truth loci missing; "
        f"anchors consistent at {anchors_ok}/{len(shared)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
