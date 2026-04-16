#!/usr/bin/env python3

import argparse
import re
import sys
from pathlib import Path

import pysam


def parse_region(region_str):
    """
    Parse a region like chr1:1000-2000 into (chrom, start0, end0).
    Input is treated as 1-based inclusive.
    Output is 0-based half-open for pysam.fetch().
    """
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", region_str.replace(",", ""))
    if not match:
        raise ValueError(
            f"Invalid region format: {region_str}. Expected format: chrom:start-end"
        )

    chrom, start_str, end_str = match.groups()
    start = int(start_str)
    end = int(end_str)

    if start < 1 or end < start:
        raise ValueError(
            f"Invalid coordinates in region: {region_str}. Require 1 <= start <= end."
        )

    return chrom, start - 1, end


def sanitize_filename(name):
    """
    Replace characters that are unsafe in filenames.
    """
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def wrap_fasta_sequence(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract primary-aligning reads overlapping a genomic region "
            "and save each read as an individual FASTA file."
        )
    )
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument(
        "--bai",
        default=None,
        help=(
            "Optional BAM index (.bai). If omitted, pysam will try to find it "
            "automatically."
        ),
    )
    parser.add_argument(
        "--region",
        required=True,
        help=(
            "Genomic region in format chrom:start-end (1-based inclusive), "
            "e.g. chr1:1000-2000"
        ),
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for per-read FASTA files",
    )
    parser.add_argument(
        "--mapq-min",
        type=int,
        default=0,
        help="Minimum MAPQ to keep a read (default: 0)",
    )

    args = parser.parse_args()

    bam_path = Path(args.bam)
    bai_path = Path(args.bai) if args.bai else None
    outdir = Path(args.outdir)

    if not bam_path.exists():
        sys.exit(f"ERROR: BAM file not found: {bam_path}")

    if bai_path is not None and not bai_path.exists():
        sys.exit(f"ERROR: BAI file not found: {bai_path}")

    outdir.mkdir(parents=True, exist_ok=True)

    try:
        chrom, start0, end0 = parse_region(args.region)
    except ValueError as e:
        sys.exit(f"ERROR: {e}")

    seen_reads = set()
    written = 0
    skipped_no_seq = 0

    try:
        if bai_path is not None:
            bam = pysam.AlignmentFile(str(bam_path), "rb", index_filename=str(bai_path))
        else:
            bam = pysam.AlignmentFile(str(bam_path), "rb")
    except Exception as e:
        sys.exit(f"ERROR: Failed to open BAM/index: {e}")

    try:
        for read in bam.fetch(chrom, start0, end0):
            if read.is_unmapped:
                continue
            if read.is_secondary:
                continue
            if read.is_supplementary:
                continue
            if read.mapping_quality < args.mapq_min:
                continue

            read_name = read.query_name

            if read_name in seen_reads:
                continue
            seen_reads.add(read_name)

            seq = read.query_sequence
            if not seq:
                skipped_no_seq += 1
                continue

            safe_name = sanitize_filename(read_name)
            fasta_path = outdir / f"{safe_name}.fasta"

            # Handle rare filename collisions after sanitization.
            if fasta_path.exists():
                counter = 2
                while True:
                    candidate = outdir / f"{safe_name}_{counter}.fasta"
                    if not candidate.exists():
                        fasta_path = candidate
                        break
                    counter += 1

            header = (
                f">{read_name} {chrom}:{read.reference_start + 1}-{read.reference_end} "
                f"MAPQ={read.mapping_quality}"
            )
            fasta_text = header + "\n" + wrap_fasta_sequence(seq) + "\n"

            with open(fasta_path, "w") as fh:
                fh.write(fasta_text)

            written += 1

    except ValueError as e:
        sys.exit(f"ERROR while fetching region {args.region}: {e}")
    finally:
        bam.close()

    print(f"Reads written: {written}")
    print(f"Reads skipped due to missing sequence: {skipped_no_seq}")
    print(f"Output directory: {outdir}")


if __name__ == "__main__":
    main()
