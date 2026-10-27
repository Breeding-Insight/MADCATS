"""Smoke-test a built madcats executable against the source module.

Writes a small synthetic panel and FASTQ set, runs ``generate`` through the
executable with several workers (the frozen-multiprocessing path) and through
``madcats.py`` with one, and requires byte-identical MADC tables.

    python packaging/smoke_test.py dist/madcats
"""
import csv
import gzip
import random
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def write_fixture(root: Path) -> None:
    rng = random.Random(7)
    fastq_dir = root / "fastq"
    fastq_dir.mkdir()
    markers = []
    for index in range(6):
        primer = "".join(rng.choice("ACGT") for _ in range(20))
        flank = "".join(rng.choice("ACGT") for _ in range(60))
        ref = primer + flank[:20] + "C" + flank[20:60]
        alt = primer + flank[:20] + "A" + flank[20:60]
        markers.append((f"m{index}", primer, ref[:81], alt[:81]))
    with (root / "tags.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["MarkerName", "LeftPrimerSequence", "RefSequence", "AltSequence"])
        writer.writerows(markers)
    for sample in range(8):
        with gzip.open(fastq_dir / f"s{sample}.fastq.gz", "wt") as handle:
            read = 0
            for _, _, ref, alt in markers:
                for _ in range(rng.randint(20, 60)):
                    sequence = rng.choice([ref, alt])
                    if rng.random() < 0.3:
                        position = rng.randrange(25, 81)
                        sequence = sequence[:position] + rng.choice("ACGT") + sequence[position + 1:]
                    handle.write(f"@r{read}\n{sequence}\n+\n{'I' * len(sequence)}\n")
                    read += 1


def main() -> int:
    executable = Path(sys.argv[1]).resolve()
    subprocess.run([str(executable), "--version"], check=True)
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        write_fixture(root)
        common = ["generate", "--fastq-dir", str(root / "fastq"), "--reference-tags", str(root / "tags.csv")]
        subprocess.run([str(executable), *common, "--jobs", "4", "--output", str(root / "frozen.csv")],
                       check=True, timeout=300)
        subprocess.run([sys.executable, str(REPO / "madcats.py"), *common, "--jobs", "1",
                        "--output", str(root / "source.csv")], check=True, timeout=300)
        frozen = (root / "frozen.csv").read_bytes()
        if frozen != (root / "source.csv").read_bytes():
            print("FAIL: executable and source module produced different MADC tables", file=sys.stderr)
            return 1
        rows = len(frozen.splitlines()) - 1
        print(f"OK: identical MADC tables ({rows} allele rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
