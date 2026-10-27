#!/usr/bin/env python3
"""Select the highest-count query rows from a BLAST FASTA/metadata pair."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--population", required=True)
    parser.add_argument("--top", required=True, type=int)
    parser.add_argument("--metadata-output", required=True, type=Path)
    parser.add_argument("--fasta-output", required=True, type=Path)
    args = parser.parse_args()

    with args.metadata.open(newline="") as handle:
        rows = [
            row for row in csv.DictReader(handle) if row["population"] == args.population
        ]
    rows.sort(key=lambda row: (-int(row["total_count"]), row["query_id"]))
    selected = rows[: args.top]
    query_ids = {row["query_id"] for row in selected}

    sequences: dict[str, str] = {}
    current = ""
    with args.fasta.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in query_ids:
                    sequences[current] = ""
            elif current in sequences:
                sequences[current] += line

    args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
    with args.metadata_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(selected[0]))
        writer.writeheader()
        writer.writerows(selected)
    with args.fasta_output.open("w") as handle:
        for row in selected:
            query_id = row["query_id"]
            handle.write(f">{query_id}\n{sequences[query_id]}\n")


if __name__ == "__main__":
    main()
