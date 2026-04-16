#!/usr/bin/env python3

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


DEFAULT_DB_PATH = "PATH_TO_DB"


def run_blastn_for_file(fasta_path: Path, db_path: str) -> None:
    output_path = fasta_path.with_name(f"{fasta_path.stem}_blasted.txt")
    command = ["blastn", "-query", str(fasta_path), "-db", db_path]

    with output_path.open("w", encoding="utf-8") as output_handle:
        subprocess.run(command, stdout=output_handle, check=True)


def find_fasta_files(input_dir: Path) -> list[Path]:
    return sorted(path for path in input_dir.iterdir() if path.is_file() and path.suffix == ".fasta")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run blastn on every .fasta file in a directory and save stdout to *_blasted.txt files."
    )
    parser.add_argument("directory", help="Directory containing .fasta files")
    parser.add_argument(
        "--db",
        default=DEFAULT_DB_PATH,
        help="Path to the BLAST database. Replace the default or pass --db explicitly.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_dir = Path(args.directory).expanduser().resolve()

    if not input_dir.is_dir():
        raise SystemExit(f"Not a directory: {input_dir}")

    fasta_files = find_fasta_files(input_dir)
    if not fasta_files:
        print(f"No .fasta files found in {input_dir}")
        return 0

    for fasta_path in fasta_files:
        print(f"Running blastn for {fasta_path.name}")
        run_blastn_for_file(fasta_path, args.db)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
