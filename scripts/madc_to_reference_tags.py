#!/usr/bin/env python3
"""Convert a delivered MADC's canonical rows into a reference-tags CSV.

A single vendor report reconstructs the panel definition without any design
file: the Ref/Alt rows carry the exact tags, their first difference is the
variant position, and a tag prefix serves as the assignment anchor (the true
primer boundary is unobservable and unnecessary).  Ambiguity-coded canonical
tags (vendor-curated IUPAC) convert as-is; the pipeline classifies
compatible reads as Match categories but such canonical rows accumulate no
exact counts, so they are counted and reported here.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from madc_pipeline import Locus, allele_type, read_madc, write_reference_tags


# Why a marker cannot become a locus, in the words a user needs to act on it.
SKIP_REASONS = {
    "identical_canonical_tags": (
        "Ref and Alt rows carry the same sequence: the target lies outside the window "
        "(a SNP one base past the tag) or the design is presence/absence; the SNP is "
        "unreadable in this window and the marker is dropped"
    ),
    "missing_canonical_row": "no Ref or no Alt row in the MADC; the marker is dropped",
    "length_mismatch": "Ref and Alt rows differ in length; the marker is dropped",
    "variant_at_tag_start": "Ref and Alt differ at the first base, where no anchor can sit; the marker is dropped",
    "tags_diverge_inside_anchor": (
        "Ref and Alt differ inside the assignment anchor; shorten --anchor-length or the marker is dropped"
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--madc", required=True, type=Path, help="delivered MADC CSV")
    parser.add_argument(
        "--anchor-length",
        type=int,
        default=24,
        help="assignment-anchor prefix length, capped below each variant (default: 24)",
    )
    parser.add_argument("--output", required=True, type=Path, help="reference-tags CSV to write")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    parsed = read_madc(args.madc)
    canonical: dict[str, dict[str, str]] = defaultdict(dict)
    for record in parsed["records"]:
        category = allele_type(record["allele_id"])
        if category in ("Ref", "Alt"):
            canonical[record["clone_id"]][category] = record["sequence"]

    loci: list[Locus] = []
    skipped: dict[str, list[str]] = defaultdict(list)
    iupac_loci = 0
    for marker in sorted(canonical):
        pair = canonical[marker]
        if "Ref" not in pair or "Alt" not in pair:
            skipped["missing_canonical_row"].append(marker)
            continue
        ref_tag, alt_tag = pair["Ref"], pair["Alt"]
        if len(ref_tag) != len(alt_tag):
            skipped["length_mismatch"].append(marker)
            continue
        differences = [i for i, (a, b) in enumerate(zip(ref_tag, alt_tag)) if a != b]
        if not differences:
            skipped["identical_canonical_tags"].append(marker)
            continue
        variant_pos = differences[0]
        if variant_pos == 0:
            skipped["variant_at_tag_start"].append(marker)
            continue
        if set(ref_tag + alt_tag) - set("ACGT"):
            iupac_loci += 1
        anchor = ref_tag[: min(args.anchor_length, variant_pos)]
        if not alt_tag.startswith(anchor):
            skipped["tags_diverge_inside_anchor"].append(marker)
            continue
        loci.append(
            Locus(
                index=len(loci),
                marker=marker,
                left_primer=anchor,
                right_primer="",
                ref_base=ref_tag[variant_pos],
                alt_base=alt_tag[variant_pos],
                variant_pos=variant_pos,
                ref_tag=ref_tag,
                alt_tag=alt_tag,
            )
        )

    write_reference_tags(args.output, loci)
    summary = {
        "clones_in_madc": len(canonical),
        "loci_written": len(loci),
        "iupac_canonical_loci": iupac_loci,
        "skipped": {reason: len(markers) for reason, markers in skipped.items()},
    }
    print(json.dumps(summary, indent=2))
    if skipped:
        # Every dropped marker is written beside the tag table so the loss is
        # auditable, and each reason is announced in plain words: a count in a
        # summary is not a warning anyone acts on.
        sidecar = args.output.with_name(args.output.name + ".skipped.tsv")
        with sidecar.open("w") as handle:
            handle.write("marker\treason\texplanation\n")
            for reason, markers in skipped.items():
                for marker in markers:
                    handle.write(f"{marker}\t{reason}\t{SKIP_REASONS.get(reason, reason)}\n")
        for reason, markers in skipped.items():
            print(
                f"WARNING: {len(markers)} marker(s) skipped ({reason}): "
                f"{SKIP_REASONS.get(reason, reason)}; e.g. {', '.join(markers[:3])}",
                file=sys.stderr,
            )
        print(f"WARNING: skipped markers listed in {sidecar}", file=sys.stderr)


if __name__ == "__main__":
    main()
