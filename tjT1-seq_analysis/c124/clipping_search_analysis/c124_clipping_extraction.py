#!/usr/bin/env python3

import argparse
import configparser
import csv
import hashlib
import re
import secrets
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pysam


SCRIPT_LABEL = "clipping_extraction"
SCRIPT_VERSION = "0.1.0"
RUN_ID_LENGTH = 20


def build_parser():
    """Build the command-line parser for clipping extraction."""
    parser = argparse.ArgumentParser(
        description=(
            "Extract read IDs and clipping-related information for reads "
            "overlapping a target genomic interval."
        )
    )
    parser.add_argument(
        "-f",
        "--bam-file",
        help="Input BAM file.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        help="Parent output directory. Must already exist.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        help="Required run prefix/ID, for example: tjt1-c-end",
    )
    parser.add_argument(
        "-r",
        "--range",
        nargs=2,
        type=int,
        metavar=("START", "END"),
        help="Target interval coordinates, for example: -r 8644000 8644500",
    )
    parser.add_argument(
        "-c",
        "--chromosome",
        help="Chromosome name for the search interval.",
    )
    parser.add_argument(
        "--clipping-side",
        choices=["left", "right", "both"],
        default="both",
        help="Which clipping side(s) to extract. Default: both",
    )
    parser.add_argument(
        "--alignment-scope",
        choices=["primary", "primary_and_supplementary"],
        default="primary_and_supplementary",
        help=(
            "Which alignments to consider in the region scan. "
            "Default: primary_and_supplementary"
        ),
    )
    parser.add_argument(
        "--min-clip-length",
        type=int,
        default=100,
        help="Minimum clipping length to count as clipped. Default: 100",
    )
    parser.add_argument(
        "--write-whole-fastas",
        action="store_true",
        help="Also write FASTA records for full reads.",
    )
    parser.add_argument(
        "--make-config",
        help="Write a sample INI config file to this path and exit.",
    )
    parser.add_argument(
        "--use-config",
        help=(
            "Read arguments from an INI config file and run with those settings. "
            "This flag cannot be combined with any other run-setting flag."
        ),
    )
    return parser


def create_sample_config(config_path):
    """Write a sample INI config file for this script."""
    config_file = Path(config_path)
    if config_file.exists():
        raise ValueError(f"Refusing to overwrite existing config file: {config_file}")

    config_text = """[clipping_extraction]
bam_file = /path/to/input.bam
output_dir = /path/to/existing_output_parent
prefix = tjt1-c-end
chromosome = V
range_start = 8644000
range_end = 8644500
clipping_side = both
alignment_scope = primary_and_supplementary
write_whole_fastas = false
min_clip_length = 100
"""
    config_file.parent.mkdir(parents=True, exist_ok=True)
    config_file.write_text(config_text)


def validate_make_config_prefix(prefix_value):
    """Validate the prefix used to build a sample config filename."""
    if "." in prefix_value:
        raise ValueError(
            "--make-config expects only a prefix, without periods or a filename extension."
        )
    if not prefix_value.strip():
        raise ValueError("--make-config prefix cannot be empty.")
    return prefix_value


def resolve_make_config_path(prefix_value):
    prefix = validate_make_config_prefix(prefix_value)
    return Path(f"{prefix}_config.ini").resolve()


def parse_config(config_path):
    parser = configparser.ConfigParser()
    read_ok = parser.read(config_path)
    if not read_ok:
        raise ValueError(f"Could not read config file: {config_path}")
    if "clipping_extraction" not in parser:
        raise ValueError(
            "Config file must contain a [clipping_extraction] section."
        )

    section = parser["clipping_extraction"]
    bam_file = section.get("bam_file", "").strip()
    output_dir = section.get("output_dir", "").strip()
    prefix = section.get("prefix", "").strip()
    chromosome = section.get("chromosome", "").strip()
    range_start = section.get("range_start", "").strip()
    range_end = section.get("range_end", "").strip()
    clipping_side = section.get("clipping_side", "both").strip()
    alignment_scope = section.get("alignment_scope", "primary_and_supplementary").strip()
    write_whole_fastas = section.getboolean("write_whole_fastas", fallback=False)
    min_clip_length = section.getint("min_clip_length", fallback=100)

    if not bam_file:
        raise ValueError("Config file must define clipping_extraction.bam_file")
    if not output_dir:
        raise ValueError("Config file must define clipping_extraction.output_dir")
    if not prefix:
        raise ValueError("Config file must define clipping_extraction.prefix")
    if not chromosome:
        raise ValueError("Config file must define clipping_extraction.chromosome")
    if not range_start or not range_end:
        raise ValueError(
            "Config file must define clipping_extraction.range_start and range_end"
        )
    if clipping_side not in {"left", "right", "both"}:
        raise ValueError(
            "Config value clipping_side must be one of: left, right, both"
        )
    if alignment_scope not in {"primary", "primary_and_supplementary"}:
        raise ValueError(
            "Config value alignment_scope must be one of: "
            "primary, primary_and_supplementary"
        )
    if min_clip_length < 0:
        raise ValueError("Config value min_clip_length must be >= 0")

    return {
        "bam_file": bam_file,
        "output_dir": output_dir,
        "prefix": prefix,
        "chromosome": chromosome,
        "range": [int(range_start), int(range_end)],
        "clipping_side": clipping_side,
        "alignment_scope": alignment_scope,
        "write_whole_fastas": write_whole_fastas,
        "min_clip_length": min_clip_length,
        "config_source_path": str(Path(config_path).resolve()),
    }


def validate_args(args):
    if args.make_config:
        other_args_used = any(
            value is not None
            for value in [
                args.bam_file,
                args.output_dir,
                args.prefix,
                args.range,
                args.use_config,
            ]
        )
        if other_args_used or args.chromosome is not None:
            raise ValueError(
                "--make-config cannot be combined with any other run-setting flag."
            )
        if (
            args.clipping_side != "both"
            or args.alignment_scope != "primary_and_supplementary"
            or args.write_whole_fastas
            or args.min_clip_length != 100
        ):
            raise ValueError(
                "--make-config cannot be combined with any other run-setting flag."
            )
        return {"make_config": str(resolve_make_config_path(args.make_config))}

    if args.use_config:
        other_args_used = any(
            value is not None
            for value in [
                args.bam_file,
                args.output_dir,
                args.prefix,
                args.range,
                args.make_config,
            ]
        )
        if other_args_used or args.chromosome is not None:
            raise ValueError(
                "--use-config cannot be combined with any other run-setting flag."
            )
        if (
            args.clipping_side != "both"
            or args.alignment_scope != "primary_and_supplementary"
            or args.write_whole_fastas
            or args.min_clip_length != 100
        ):
            raise ValueError(
                "--use-config cannot be combined with any other run-setting flag."
            )
        return parse_config(args.use_config)

    if (
        not args.bam_file
        or not args.output_dir
        or not args.prefix
        or not args.chromosome
        or not args.range
    ):
        raise ValueError(
            "You must provide -f/--bam-file, -o/--output-dir, "
            "-p/--prefix, -c/--chromosome, and -r/--range unless using --use-config."
        )

    return {
        "bam_file": args.bam_file,
        "output_dir": args.output_dir,
        "prefix": args.prefix,
        "chromosome": args.chromosome,
        "range": args.range,
        "clipping_side": args.clipping_side,
        "alignment_scope": args.alignment_scope,
        "write_whole_fastas": args.write_whole_fastas,
        "min_clip_length": args.min_clip_length,
        "config_source_path": None,
    }


def clip_passes_threshold(clip_length, settings):
    return clip_length >= settings["min_clip_length"]


def should_keep_alignment(read, settings):
    if read.is_unmapped:
        return False
    if read.is_secondary:
        return False
    if settings["alignment_scope"] == "primary":
        return not read.is_supplementary
    if settings["alignment_scope"] == "primary_and_supplementary":
        return True
    raise ValueError(f"Unsupported alignment_scope: {settings['alignment_scope']}")


def validate_range(start, end):
    if start < 1 or end < start:
        raise ValueError("Range must satisfy 1 <= start <= end.")


def build_run_directory_name(prefix, chromosome, start, end):
    return f"{prefix}_{SCRIPT_LABEL}_{chromosome.lower()}{start}-{end}"


def sanitize_filename(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def generate_run_uid():
    alphabet = "abcdefghijklmnopqrstuvwxyz0123456789"
    return "".join(secrets.choice(alphabet) for _ in range(RUN_ID_LENGTH))


def get_script_sha256(script_path):
    digest = hashlib.sha256()
    with open(script_path, "rb") as handle:
        digest.update(handle.read())
    return digest.hexdigest()


def get_git_commit(script_path):
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=script_path.parent,
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return None
    return completed.stdout.strip()


def prepare_run_directory(settings, script_path):
    output_dir = Path(settings["output_dir"]).resolve()
    if not output_dir.is_dir():
        raise ValueError(
            f"Output directory must already exist and be a directory: {output_dir}"
        )

    start, end = settings["range"]
    run_dir_name = build_run_directory_name(
        settings["prefix"],
        settings["chromosome"],
        start,
        end,
    )
    run_dir = output_dir / run_dir_name
    if run_dir.exists():
        raise ValueError(f"Refusing to run because output directory exists: {run_dir}")

    run_metadata = {
        "run_uid": generate_run_uid(),
        "script_label": SCRIPT_LABEL,
        "script_version": SCRIPT_VERSION,
        "script_sha256": get_script_sha256(script_path),
        "git_commit": get_git_commit(script_path),
        "run_dir": run_dir,
    }
    return run_metadata


def get_clip_lengths(read):
    if not read.cigartuples:
        return 0, 0, False, False

    left_op, left_len = read.cigartuples[0]
    right_op, right_len = read.cigartuples[-1]

    left_clip = left_len if left_op in {4, 5} else 0
    right_clip = right_len if right_op in {4, 5} else 0
    left_soft = left_op == 4
    right_soft = right_op == 4

    return left_clip, right_clip, left_soft, right_soft


def parse_cigar_end_clipping(cigar_string):
    cigar_ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar_string)
    if not cigar_ops:
        return {
            "sa_left_clip_length": 0,
            "sa_right_clip_length": 0,
            "sa_left_clip_is_soft": False,
            "sa_right_clip_is_soft": False,
        }

    left_len_str, left_op = cigar_ops[0]
    right_len_str, right_op = cigar_ops[-1]

    left_clip_length = int(left_len_str) if left_op in {"S", "H"} else 0
    right_clip_length = int(right_len_str) if right_op in {"S", "H"} else 0

    return {
        "sa_left_clip_length": left_clip_length,
        "sa_right_clip_length": right_clip_length,
        "sa_left_clip_is_soft": left_op == "S",
        "sa_right_clip_is_soft": right_op == "S",
    }


def get_reference_consumed_bases(cigar_string):
    cigar_ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar_string)
    reference_bases = 0
    for length_str, op in cigar_ops:
        if op in {"M", "D", "N", "=", "X"}:
            reference_bases += int(length_str)
    return reference_bases


def wrap_fasta_sequence(sequence, width=80):
    return "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))


def get_clipped_sequences(read):
    sequence = read.query_sequence or ""
    if not sequence or not read.cigartuples:
        return "", ""

    left_op, left_len = read.cigartuples[0]
    right_op, right_len = read.cigartuples[-1]

    left_sequence = sequence[:left_len] if left_op in {4, 5} and left_op == 4 else ""
    right_sequence = sequence[-right_len:] if right_op in {4, 5} and right_op == 4 else ""

    return left_sequence, right_sequence


def collect_primary_region_rows(settings):
    bam_path = Path(settings["bam_file"]).resolve()
    chromosome = settings["chromosome"]
    start, end = settings["range"]
    start0 = start - 1
    rows = []

    with pysam.AlignmentFile(str(bam_path), "rb") as bam_file:
        for read in bam_file.fetch(chromosome, start0, end):
            if not should_keep_alignment(read, settings):
                continue

            left_clip, right_clip, left_soft, right_soft = get_clip_lengths(read)
            left_sequence, right_sequence = get_clipped_sequences(read)
            has_sa = read.has_tag("SA")
            sa_entries = parse_sa_tag(read.get_tag("SA")) if has_sa else []

            rows.append(
                {
                    "read_id": read.query_name,
                    "chromosome": chromosome,
                    "region_start": start,
                    "region_end": end,
                    "reference_start": read.reference_start + 1,
                    "reference_end": read.reference_end,
                    "mapping_quality": read.mapping_quality,
                    "is_reverse": read.is_reverse,
                    "alignment_type": (
                        "supplementary" if read.is_supplementary else "primary"
                    ),
                    "query_length": read.query_length or 0,
                    "query_sequence": read.query_sequence or "",
                    "left_clip_length": left_clip,
                    "right_clip_length": right_clip,
                    "left_clip_passes_threshold": clip_passes_threshold(left_clip, settings),
                    "right_clip_passes_threshold": clip_passes_threshold(right_clip, settings),
                    "left_soft_clip": left_soft,
                    "right_soft_clip": right_soft,
                    "left_clip_sequence": left_sequence,
                    "right_clip_sequence": right_sequence,
                    "has_sa_tag": has_sa,
                    "sa_entries": sa_entries,
                }
            )

    return rows


def write_primary_region_tsv(rows, run_metadata):
    run_dir = run_metadata["run_dir"]
    output_path = run_dir / f"{run_dir.name}_primary_region_reads.tsv"
    fieldnames = [
        "read_id",
        "chromosome",
        "region_start",
        "region_end",
        "reference_start",
        "reference_end",
        "mapping_quality",
        "is_reverse",
        "alignment_type",
        "query_length",
        "left_clip_length",
        "right_clip_length",
        "left_clip_passes_threshold",
        "right_clip_passes_threshold",
        "left_soft_clip",
        "right_soft_clip",
        "has_sa_tag",
    ]

    with open(output_path, "w", newline="") as handle:
        handle.write(f"# run_uid={run_metadata['run_uid']}\n")
        handle.write(f"# script_version={run_metadata['script_version']}\n")
        handle.write(f"# script_sha256={run_metadata['script_sha256']}\n")
        if run_metadata["git_commit"]:
            handle.write(f"# git_commit={run_metadata['git_commit']}\n")

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    return output_path


def should_write_left_clip(settings):
    return settings["clipping_side"] in {"left", "both"}


def should_write_right_clip(settings):
    return settings["clipping_side"] in {"right", "both"}


def build_clip_fasta_header(row, side):
    if side == "left":
        clip_length = row["left_clip_length"]
        is_soft = row["left_soft_clip"]
    else:
        clip_length = row["right_clip_length"]
        is_soft = row["right_soft_clip"]

    soft_label = "soft" if is_soft else "hard_or_none"
    return (
        f">{row['read_id']} "
        f"side={side} "
        f"clip_length={clip_length} "
        f"clip_type={soft_label} "
        f"ref={row['chromosome']}:{row['reference_start']}-{row['reference_end']} "
        f"mapq={row['mapping_quality']} "
        f"sa_tag={str(row['has_sa_tag']).lower()}"
    )


def parse_sa_tag(sa_tag_value):
    entries = []
    for raw_entry in sa_tag_value.split(";"):
        raw_entry = raw_entry.strip()
        if not raw_entry:
            continue

        fields = raw_entry.split(",")
        if len(fields) < 6:
            continue

        chrom, pos, strand, cigar, mapq, nm = fields[:6]
        try:
            pos = int(pos)
        except ValueError:
            continue

        reference_span = get_reference_consumed_bases(cigar)
        reference_end = pos + reference_span - 1 if reference_span > 0 else pos

        entries.append(
            {
                "sa_chromosome": chrom,
                "sa_position": pos,
                "sa_reference_start": pos,
                "sa_reference_end": reference_end,
                "sa_reference_span": reference_span,
                "sa_strand": strand,
                "sa_cigar": cigar,
                "sa_mapq": mapq,
                "sa_nm": nm,
                "sa_raw_entry": raw_entry,
                **parse_cigar_end_clipping(cigar),
            }
        )

    return entries


def sa_entry_supports_left_clip(sa_entry):
    return sa_entry["sa_right_clip_length"] > 0


def sa_entry_supports_right_clip(sa_entry):
    return sa_entry["sa_left_clip_length"] > 0


def collect_sa_side_rows(rows, settings):
    sa_rows = []

    for row in rows:
        if not row["has_sa_tag"]:
            continue

        if should_write_left_clip(settings) and row["left_clip_passes_threshold"]:
            for sa_entry in row["sa_entries"]:
                if not sa_entry_supports_left_clip(sa_entry):
                    continue
                sa_rows.append(
                    {
                        "read_id": row["read_id"],
                        "clip_side": "left",
                        "clip_length": row["left_clip_length"],
                        "clip_is_soft": row["left_soft_clip"],
                        "read_chromosome": row["chromosome"],
                        "read_reference_start": row["reference_start"],
                        "read_reference_end": row["reference_end"],
                        **sa_entry,
                    }
                )

        if should_write_right_clip(settings) and row["right_clip_passes_threshold"]:
            for sa_entry in row["sa_entries"]:
                if not sa_entry_supports_right_clip(sa_entry):
                    continue
                sa_rows.append(
                    {
                        "read_id": row["read_id"],
                        "clip_side": "right",
                        "clip_length": row["right_clip_length"],
                        "clip_is_soft": row["right_soft_clip"],
                        "read_chromosome": row["chromosome"],
                        "read_reference_start": row["reference_start"],
                        "read_reference_end": row["reference_end"],
                        **sa_entry,
                    }
                )

    return sa_rows


def write_sa_tag_tsvs(rows, settings, run_metadata):
    sa_rows = collect_sa_side_rows(rows, settings)
    if not sa_rows:
        return []

    fieldnames = [
        "read_id",
        "clip_side",
        "clip_length",
        "clip_is_soft",
        "read_chromosome",
        "read_reference_start",
        "read_reference_end",
        "sa_chromosome",
        "sa_position",
        "sa_reference_start",
        "sa_reference_end",
        "sa_reference_span",
        "sa_strand",
        "sa_cigar",
        "sa_left_clip_length",
        "sa_right_clip_length",
        "sa_left_clip_is_soft",
        "sa_right_clip_is_soft",
        "sa_mapq",
        "sa_nm",
        "sa_raw_entry",
    ]

    output_paths = []
    sa_rows_by_chromosome = {}
    for row in sa_rows:
        sa_rows_by_chromosome.setdefault(row["sa_chromosome"], []).append(row)

    for sa_chromosome, chromosome_rows in sorted(sa_rows_by_chromosome.items()):
        safe_chromosome = sanitize_filename(sa_chromosome)
        output_path = (
            run_metadata["run_dir"]
            / f"{run_metadata['run_dir'].name}_sa_tags_{safe_chromosome}.tsv"
        )

        with open(output_path, "w", newline="") as handle:
            handle.write(f"# run_uid={run_metadata['run_uid']}\n")
            handle.write(f"# script_version={run_metadata['script_version']}\n")
            handle.write(f"# script_sha256={run_metadata['script_sha256']}\n")
            if run_metadata["git_commit"]:
                handle.write(f"# git_commit={run_metadata['git_commit']}\n")

            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(chromosome_rows)

        output_paths.append(output_path)

    return output_paths


def write_clipped_fastas(rows, settings, run_metadata):
    run_dir = run_metadata["run_dir"]
    left_path = run_dir / f"{run_dir.name}_left_clipping_sequences.fasta"
    right_path = run_dir / f"{run_dir.name}_right_clipping_sequences.fasta"

    left_written = 0
    right_written = 0
    left_output_path = None
    right_output_path = None

    if should_write_left_clip(settings):
        with open(left_path, "w") as handle:
            handle.write(f"# run_uid={run_metadata['run_uid']}\n")
            handle.write(f"# script_version={run_metadata['script_version']}\n")
            handle.write(f"# script_sha256={run_metadata['script_sha256']}\n")
            if run_metadata["git_commit"]:
                handle.write(f"# git_commit={run_metadata['git_commit']}\n")

            for row in rows:
                sequence = row["left_clip_sequence"]
                if not sequence or not row["left_clip_passes_threshold"]:
                    continue
                handle.write(build_clip_fasta_header(row, "left") + "\n")
                handle.write(wrap_fasta_sequence(sequence) + "\n")
                left_written += 1

        left_output_path = left_path

    if should_write_right_clip(settings):
        with open(right_path, "w") as handle:
            handle.write(f"# run_uid={run_metadata['run_uid']}\n")
            handle.write(f"# script_version={run_metadata['script_version']}\n")
            handle.write(f"# script_sha256={run_metadata['script_sha256']}\n")
            if run_metadata["git_commit"]:
                handle.write(f"# git_commit={run_metadata['git_commit']}\n")

            for row in rows:
                sequence = row["right_clip_sequence"]
                if not sequence or not row["right_clip_passes_threshold"]:
                    continue
                handle.write(build_clip_fasta_header(row, "right") + "\n")
                handle.write(wrap_fasta_sequence(sequence) + "\n")
                right_written += 1

        right_output_path = right_path

    return {
        "left_output_path": left_output_path,
        "right_output_path": right_output_path,
        "left_written": left_written,
        "right_written": right_written,
    }


def get_clipped_sides_label(row):
    has_left = row["left_clip_passes_threshold"]
    has_right = row["right_clip_passes_threshold"]

    if has_left and has_right:
        return "both"
    if has_left:
        return "left"
    if has_right:
        return "right"
    return "none"


def row_matches_read_id_group(row, clipped_side):
    has_left = row["left_clip_passes_threshold"]
    has_right = row["right_clip_passes_threshold"]

    if clipped_side == "left":
        return has_left
    if clipped_side == "right":
        return has_right
    if clipped_side == "both":
        return has_left and has_right
    if clipped_side == "none":
        return not has_left and not has_right

    raise ValueError(f"Unsupported clipped-side group: {clipped_side}")


def unique_read_ids_for_clipped_side(rows, clipped_side):
    read_ids = set()
    for row in rows:
        if row_matches_read_id_group(row, clipped_side):
            read_ids.add(row["read_id"])
    return sorted(read_ids)


def write_read_id_list(output_path, read_ids, run_metadata):
    with open(output_path, "w") as handle:
        handle.write(f"# run_uid={run_metadata['run_uid']}\n")
        handle.write(f"# script_version={run_metadata['script_version']}\n")
        handle.write(f"# script_sha256={run_metadata['script_sha256']}\n")
        if run_metadata["git_commit"]:
            handle.write(f"# git_commit={run_metadata['git_commit']}\n")
        for read_id in read_ids:
            handle.write(f"{read_id}\n")


def write_unique_read_id_files(rows, run_metadata):
    run_dir = run_metadata["run_dir"]
    output_paths = {}

    for clipped_side in ["left", "right", "both", "none"]:
        read_ids = unique_read_ids_for_clipped_side(rows, clipped_side)
        side_label = "unclipped" if clipped_side == "none" else f"{clipped_side}_clipped"
        output_path = run_dir / f"{run_dir.name}_{side_label}_read_ids.txt"
        write_read_id_list(output_path, read_ids, run_metadata)
        output_paths[clipped_side] = {
            "output_path": output_path,
            "count": len(read_ids),
        }

    return output_paths


def write_unclipped_summary(rows, run_metadata):
    unclipped_rows = []
    for row in rows:
        if get_clipped_sides_label(row) != "none":
            continue
        unclipped_rows.append(
            {
                "read_id": row["read_id"],
                "chromosome": row["chromosome"],
                "reference_start": row["reference_start"],
                "reference_end": row["reference_end"],
                "mapping_quality": row["mapping_quality"],
                "is_reverse": row["is_reverse"],
                "query_length": row["query_length"],
                "has_sa_tag": row["has_sa_tag"],
            }
        )

    output_path = run_metadata["run_dir"] / f"{run_metadata['run_dir'].name}_unclipped_summary.tsv"
    fieldnames = [
        "read_id",
        "chromosome",
        "reference_start",
        "reference_end",
        "mapping_quality",
        "is_reverse",
        "query_length",
        "has_sa_tag",
    ]

    with open(output_path, "w", newline="") as handle:
        handle.write(f"# run_uid={run_metadata['run_uid']}\n")
        handle.write(f"# script_version={run_metadata['script_version']}\n")
        handle.write(f"# script_sha256={run_metadata['script_sha256']}\n")
        if run_metadata["git_commit"]:
            handle.write(f"# git_commit={run_metadata['git_commit']}\n")

        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(unclipped_rows)

    return {
        "output_path": output_path,
        "count": len(unclipped_rows),
    }


def count_unique_read_ids(rows):
    return len({row["read_id"] for row in rows})


def build_summary_rows(rows):
    summary_rows = []

    total_overlapping_rows = len(rows)
    primary_rows = [row for row in rows if row["alignment_type"] == "primary"]
    supplementary_rows = [
        row for row in rows if row["alignment_type"] == "supplementary"
    ]
    unique_overlapping_reads = count_unique_read_ids(rows)
    unique_reads_with_primary_alignments = len({row["read_id"] for row in primary_rows})
    unique_reads_with_supplementary_alignments = len(
        {row["read_id"] for row in supplementary_rows}
    )
    left_reads = unique_read_ids_for_clipped_side(rows, "left")
    right_reads = unique_read_ids_for_clipped_side(rows, "right")
    both_reads = unique_read_ids_for_clipped_side(rows, "both")
    unclipped_reads = unique_read_ids_for_clipped_side(rows, "none")
    reads_with_sa = [row for row in rows if row["has_sa_tag"]]
    unique_reads_with_sa = len({row["read_id"] for row in reads_with_sa})

    summary_rows.extend(
        [
            {"metric": "total_overlapping_rows", "value": total_overlapping_rows},
            {"metric": "primary_overlapping_rows", "value": len(primary_rows)},
            {
                "metric": "supplementary_overlapping_rows",
                "value": len(supplementary_rows),
            },
            {"metric": "unique_overlapping_reads", "value": unique_overlapping_reads},
            {
                "metric": "unique_reads_with_primary_alignments",
                "value": unique_reads_with_primary_alignments,
            },
            {
                "metric": "unique_reads_with_supplementary_alignments",
                "value": unique_reads_with_supplementary_alignments,
            },
            {"metric": "unique_reads_with_left_clipping", "value": len(left_reads)},
            {"metric": "unique_reads_with_right_clipping", "value": len(right_reads)},
            {"metric": "unique_reads_with_both_side_clipping", "value": len(both_reads)},
            {"metric": "unique_unclipped_reads", "value": len(unclipped_reads)},
            {"metric": "overlapping_rows_with_sa_tag", "value": len(reads_with_sa)},
            {"metric": "unique_reads_with_sa_tag", "value": unique_reads_with_sa},
        ]
    )

    sa_entry_counts_by_chromosome = {}
    unique_sa_reads_by_chromosome = {}

    for row in rows:
        if not row["has_sa_tag"]:
            continue
        for sa_entry in row["sa_entries"]:
            chromosome = sa_entry["sa_chromosome"]
            sa_entry_counts_by_chromosome[chromosome] = (
                sa_entry_counts_by_chromosome.get(chromosome, 0) + 1
            )
            unique_sa_reads_by_chromosome.setdefault(chromosome, set()).add(row["read_id"])

    for chromosome in sorted(sa_entry_counts_by_chromosome):
        safe_chromosome = sanitize_filename(chromosome)
        summary_rows.append(
            {
                "metric": f"sa_entries_to_{safe_chromosome}",
                "value": sa_entry_counts_by_chromosome[chromosome],
            }
        )
        summary_rows.append(
            {
                "metric": f"unique_reads_with_sa_to_{safe_chromosome}",
                "value": len(unique_sa_reads_by_chromosome[chromosome]),
            }
        )

    return summary_rows


def write_summary_tsv(rows, run_metadata):
    output_path = run_metadata["run_dir"] / f"{run_metadata['run_dir'].name}_summary.tsv"
    summary_rows = build_summary_rows(rows)

    with open(output_path, "w", newline="") as handle:
        handle.write(f"# run_uid={run_metadata['run_uid']}\n")
        handle.write(f"# script_version={run_metadata['script_version']}\n")
        handle.write(f"# script_sha256={run_metadata['script_sha256']}\n")
        if run_metadata["git_commit"]:
            handle.write(f"# git_commit={run_metadata['git_commit']}\n")
        handle.write(
            "# summary_scope=top-level counts reflect the collected alignment rows "
            "(primary or primary_and_supplementary depending on alignment_scope)\n"
        )
        handle.write(
            "# sa_counts_scope=per-chromosome SA metrics are counted from parsed SA entries "
            "present on those collected rows\n"
        )

        writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    return {
        "output_path": output_path,
        "row_count": len(summary_rows),
    }


def get_log_timestamp():
    return datetime.now(timezone.utc).isoformat()


def initialize_run_log(settings, run_metadata):
    output_path = run_metadata["run_dir"] / f"{run_metadata['run_dir'].name}_run.log"
    with open(output_path, "w") as handle:
        handle.write(f"{get_log_timestamp()}\tSTART\trun_uid\t{run_metadata['run_uid']}\n")
        handle.write(
            f"{get_log_timestamp()}\tINFO\tscript_label\t{run_metadata['script_label']}\n"
        )
        handle.write(
            f"{get_log_timestamp()}\tINFO\tscript_version\t{run_metadata['script_version']}\n"
        )
        handle.write(
            f"{get_log_timestamp()}\tINFO\tscript_sha256\t{run_metadata['script_sha256']}\n"
        )
        handle.write(
            f"{get_log_timestamp()}\tINFO\tgit_commit\t"
            f"{run_metadata['git_commit'] if run_metadata['git_commit'] else 'unavailable'}\n"
        )
        handle.write(
            f"{get_log_timestamp()}\tINFO\tbam_file\t{Path(settings['bam_file']).resolve()}\n"
        )
        handle.write(
            f"{get_log_timestamp()}\tINFO\toutput_dir\t{Path(settings['output_dir']).resolve()}\n"
        )
        handle.write(f"{get_log_timestamp()}\tINFO\trun_dir\t{run_metadata['run_dir']}\n")
        handle.write(f"{get_log_timestamp()}\tINFO\tprefix\t{settings['prefix']}\n")
        handle.write(f"{get_log_timestamp()}\tINFO\tchromosome\t{settings['chromosome']}\n")
        handle.write(f"{get_log_timestamp()}\tINFO\trange_start\t{settings['range'][0]}\n")
        handle.write(f"{get_log_timestamp()}\tINFO\trange_end\t{settings['range'][1]}\n")
        handle.write(
            f"{get_log_timestamp()}\tINFO\tclipping_side\t{settings['clipping_side']}\n"
        )
        handle.write(
            f"{get_log_timestamp()}\tINFO\talignment_scope\t{settings['alignment_scope']}\n"
        )
        handle.write(
            f"{get_log_timestamp()}\tINFO\twrite_whole_fastas\t"
            f"{str(settings['write_whole_fastas']).lower()}\n"
        )
        handle.write(
            f"{get_log_timestamp()}\tINFO\tmin_clip_length\t{settings['min_clip_length']}\n"
        )
    return output_path


def append_run_log(log_path, status, key, value):
    with open(log_path, "a") as handle:
        handle.write(f"{get_log_timestamp()}\t{status}\t{key}\t{value}\n")


def write_whole_read_fastas(rows, settings, run_metadata):
    if not settings["write_whole_fastas"]:
        return {
            "output_path": None,
            "written": 0,
        }

    output_path = run_metadata["run_dir"] / f"{run_metadata['run_dir'].name}_whole_reads.fasta"
    written = 0

    with open(output_path, "w") as handle:
        handle.write(f"# run_uid={run_metadata['run_uid']}\n")
        handle.write(f"# script_version={run_metadata['script_version']}\n")
        handle.write(f"# script_sha256={run_metadata['script_sha256']}\n")
        if run_metadata["git_commit"]:
            handle.write(f"# git_commit={run_metadata['git_commit']}\n")

        seen_read_ids = set()
        for row in rows:
            if row["read_id"] in seen_read_ids:
                continue

            sequence = row["query_sequence"]
            if not sequence:
                continue

            seen_read_ids.add(row["read_id"])
            clipped_sides = get_clipped_sides_label(row)
            is_clipped = "true" if clipped_sides != "none" else "false"

            handle.write(
                f">{row['read_id']} "
                f"ref={row['chromosome']}:{row['reference_start']}-{row['reference_end']} "
                f"mapq={row['mapping_quality']} "
                f"is_clipped={is_clipped} "
                f"clipped_sides={clipped_sides} "
                f"left_clip={row['left_clip_length']} "
                f"right_clip={row['right_clip_length']} "
                f"sa_tag={str(row['has_sa_tag']).lower()}\n"
            )
            handle.write(wrap_fasta_sequence(sequence) + "\n")
            written += 1

    return {
        "output_path": output_path,
        "written": written,
    }


def build_run_config_text(settings):
    start, end = settings["range"]
    return """[clipping_extraction]
bam_file = {bam_file}
output_dir = {output_dir}
prefix = {prefix}
chromosome = {chromosome}
range_start = {range_start}
range_end = {range_end}
clipping_side = {clipping_side}
alignment_scope = {alignment_scope}
write_whole_fastas = {write_whole_fastas}
min_clip_length = {min_clip_length}
""".format(
        bam_file=Path(settings["bam_file"]).resolve(),
        output_dir=Path(settings["output_dir"]).resolve(),
        prefix=settings["prefix"],
        chromosome=settings["chromosome"],
        range_start=start,
        range_end=end,
        clipping_side=settings["clipping_side"],
        alignment_scope=settings["alignment_scope"],
        write_whole_fastas=str(settings["write_whole_fastas"]).lower(),
        min_clip_length=settings["min_clip_length"],
    )


def write_run_config(settings, run_metadata):
    run_dir = run_metadata["run_dir"]
    if settings["config_source_path"] is not None:
        source_path = Path(settings["config_source_path"]).resolve()
        config_path = run_dir / source_path.name
        shutil.copy2(source_path, config_path)
    else:
        config_path = run_dir / f"{settings['prefix']}_config.ini"
        config_path.write_text(build_run_config_text(settings))

    return config_path


def finalize_original_config_removal(settings, copied_config_path):
    return


def main():
    args = build_parser().parse_args()
    log_path = None

    try:
        settings = validate_args(args)
        if "make_config" in settings:
            create_sample_config(settings["make_config"])
            print(f"Wrote config template: {settings['make_config']}")
            return

        bam_path = Path(settings["bam_file"]).resolve()
        if not bam_path.is_file():
            raise ValueError(f"BAM file does not exist: {bam_path}")

        start, end = settings["range"]
        validate_range(start, end)

        script_path = Path(__file__).resolve()
        run_metadata = prepare_run_directory(settings, script_path)
        run_metadata["run_dir"].mkdir(parents=True, exist_ok=False)
        log_path = initialize_run_log(settings, run_metadata)
        append_run_log(log_path, "DONE", "run_directory_created", run_metadata["run_dir"])

        append_run_log(log_path, "START", "write_run_config", "begin")
        run_config_path = write_run_config(settings, run_metadata)
        append_run_log(log_path, "DONE", "run_config_path", run_config_path)

        append_run_log(log_path, "START", "collect_primary_region_rows", "begin")
        primary_rows = collect_primary_region_rows(settings)
        append_run_log(log_path, "DONE", "primary_rows_found", len(primary_rows))

        append_run_log(log_path, "START", "write_primary_region_tsv", "begin")
        primary_tsv_path = write_primary_region_tsv(primary_rows, run_metadata)
        append_run_log(log_path, "DONE", "primary_region_tsv", primary_tsv_path)

        append_run_log(log_path, "START", "write_clipped_fastas", "begin")
        clip_fasta_outputs = write_clipped_fastas(primary_rows, settings, run_metadata)
        append_run_log(
            log_path,
            "DONE",
            "left_clipping_fasta",
            clip_fasta_outputs["left_output_path"]
            if clip_fasta_outputs["left_output_path"] is not None
            else "not_written",
        )
        append_run_log(
            log_path,
            "DONE",
            "right_clipping_fasta",
            clip_fasta_outputs["right_output_path"]
            if clip_fasta_outputs["right_output_path"] is not None
            else "not_written",
        )

        append_run_log(log_path, "START", "write_sa_tag_tsvs", "begin")
        sa_tsv_paths = write_sa_tag_tsvs(primary_rows, settings, run_metadata)
        append_run_log(log_path, "DONE", "sa_tag_tsv_count", len(sa_tsv_paths))

        append_run_log(log_path, "START", "write_unique_read_id_files", "begin")
        unique_read_id_outputs = write_unique_read_id_files(primary_rows, run_metadata)
        append_run_log(
            log_path,
            "DONE",
            "unique_reads_with_left_clipping",
            unique_read_id_outputs["left"]["count"],
        )
        append_run_log(
            log_path,
            "DONE",
            "unique_reads_with_right_clipping",
            unique_read_id_outputs["right"]["count"],
        )
        append_run_log(
            log_path,
            "DONE",
            "unique_reads_with_both_side_clipping",
            unique_read_id_outputs["both"]["count"],
        )
        append_run_log(
            log_path,
            "DONE",
            "unique_unclipped_reads",
            unique_read_id_outputs["none"]["count"],
        )

        append_run_log(log_path, "START", "write_unclipped_summary", "begin")
        unclipped_summary_output = write_unclipped_summary(primary_rows, run_metadata)
        append_run_log(
            log_path,
            "DONE",
            "unclipped_summary_path",
            unclipped_summary_output["output_path"],
        )

        append_run_log(log_path, "START", "write_summary_tsv", "begin")
        summary_output = write_summary_tsv(primary_rows, run_metadata)
        append_run_log(log_path, "DONE", "run_summary_path", summary_output["output_path"])

        append_run_log(log_path, "START", "write_whole_read_fastas", "begin")
        whole_read_fasta_output = write_whole_read_fastas(
            primary_rows, settings, run_metadata
        )
        append_run_log(
            log_path,
            "DONE",
            "whole_read_fasta",
            whole_read_fasta_output["output_path"]
            if whole_read_fasta_output["output_path"] is not None
            else "not_written",
        )
        append_run_log(
            log_path,
            "DONE",
            "whole_reads_written",
            whole_read_fasta_output["written"],
        )
        append_run_log(log_path, "SUCCESS", "status", "completed")

        print("Argument validation complete.")
        print(f"BAM file: {bam_path}")
        print(
            f"Region: {settings['chromosome']}:{start}-{end} "
            f"(clipping side: {settings['clipping_side']})"
        )
        print(f"Run directory created: {run_metadata['run_dir']}")
        print(f"Run config saved: {run_config_path}")
        print(f"Primary-read TSV saved: {primary_tsv_path}")
        print(f"Primary reads found: {len(primary_rows)}")
        if clip_fasta_outputs["left_output_path"] is not None:
            print(f"Left-clipping FASTA saved: {clip_fasta_outputs['left_output_path']}")
            print(f"Left clipped sequences written: {clip_fasta_outputs['left_written']}")
        if clip_fasta_outputs["right_output_path"] is not None:
            print(f"Right-clipping FASTA saved: {clip_fasta_outputs['right_output_path']}")
            print(f"Right clipped sequences written: {clip_fasta_outputs['right_written']}")
        if sa_tsv_paths:
            print(f"SA-tag chromosome TSVs written: {len(sa_tsv_paths)}")
            for sa_tsv_path in sa_tsv_paths:
                print(f"SA-tag TSV saved: {sa_tsv_path}")
        print(
            f"Unique reads with left clipping: {unique_read_id_outputs['left']['count']} "
            f"({unique_read_id_outputs['left']['output_path']})"
        )
        print(
            f"Unique reads with right clipping: {unique_read_id_outputs['right']['count']} "
            f"({unique_read_id_outputs['right']['output_path']})"
        )
        print(
            f"Unique reads with both-side clipping: {unique_read_id_outputs['both']['count']} "
            f"({unique_read_id_outputs['both']['output_path']})"
        )
        print(
            f"Unique unclipped reads: {unique_read_id_outputs['none']['count']} "
            f"({unique_read_id_outputs['none']['output_path']})"
        )
        print(
            f"Unclipped summary saved: {unclipped_summary_output['output_path']} "
            f"({unclipped_summary_output['count']} rows)"
        )
        print(
            f"Run summary saved: {summary_output['output_path']} "
            f"({summary_output['row_count']} metrics)"
        )
        if whole_read_fasta_output["output_path"] is not None:
            print(f"Whole-read FASTA saved: {whole_read_fasta_output['output_path']}")
            print(f"Whole reads written: {whole_read_fasta_output['written']}")
        print(f"Run log saved: {log_path}")
        print(f"Run UID: {run_metadata['run_uid']}")
        print(f"Script version: {run_metadata['script_version']}")
        print(f"Script SHA256: {run_metadata['script_sha256']}")
        if run_metadata["git_commit"]:
            print(f"Git commit: {run_metadata['git_commit']}")
        else:
            print("Git commit: unavailable")

        finalize_original_config_removal(settings, run_config_path)

    except ValueError as exc:
        if log_path is not None:
            append_run_log(log_path, "ERROR", "message", str(exc))
        sys.exit(f"ERROR: {exc}")


if __name__ == "__main__":
    main()
