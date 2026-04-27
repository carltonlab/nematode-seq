#!/usr/bin/env python3

import argparse
import configparser
import hashlib
import secrets
import subprocess
import sys
from pathlib import Path


SCRIPT_LABEL = "clipping_extraction"
SCRIPT_VERSION = "0.1.0"
DEFAULT_CHROMOSOME = "V"
RUN_ID_LENGTH = 20


def build_parser():
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
        default=DEFAULT_CHROMOSOME,
        help=f"Chromosome name for the search interval. Default: {DEFAULT_CHROMOSOME}",
    )
    parser.add_argument(
        "--clipping-side",
        choices=["left", "right", "both"],
        default="both",
        help="Which clipping side(s) to extract. Default: both",
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
write_whole_fastas = false
"""
    config_file.parent.mkdir(parents=True, exist_ok=True)
    config_file.write_text(config_text)


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
    chromosome = section.get("chromosome", DEFAULT_CHROMOSOME).strip()
    range_start = section.get("range_start", "").strip()
    range_end = section.get("range_end", "").strip()
    clipping_side = section.get("clipping_side", "both").strip()
    write_whole_fastas = section.getboolean("write_whole_fastas", fallback=False)

    if not bam_file:
        raise ValueError("Config file must define clipping_extraction.bam_file")
    if not output_dir:
        raise ValueError("Config file must define clipping_extraction.output_dir")
    if not prefix:
        raise ValueError("Config file must define clipping_extraction.prefix")
    if not range_start or not range_end:
        raise ValueError(
            "Config file must define clipping_extraction.range_start and range_end"
        )
    if clipping_side not in {"left", "right", "both"}:
        raise ValueError(
            "Config value clipping_side must be one of: left, right, both"
        )

    return {
        "bam_file": bam_file,
        "output_dir": output_dir,
        "prefix": prefix,
        "chromosome": chromosome,
        "range": [int(range_start), int(range_end)],
        "clipping_side": clipping_side,
        "write_whole_fastas": write_whole_fastas,
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
        if other_args_used or args.chromosome != DEFAULT_CHROMOSOME:
            raise ValueError(
                "--make-config cannot be combined with any other run-setting flag."
            )
        if args.clipping_side != "both" or args.write_whole_fastas:
            raise ValueError(
                "--make-config cannot be combined with any other run-setting flag."
            )
        return {"make_config": args.make_config}

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
        if other_args_used or args.chromosome != DEFAULT_CHROMOSOME:
            raise ValueError(
                "--use-config cannot be combined with any other run-setting flag."
            )
        if args.clipping_side != "both" or args.write_whole_fastas:
            raise ValueError(
                "--use-config cannot be combined with any other run-setting flag."
            )
        return parse_config(args.use_config)

    if not args.bam_file or not args.output_dir or not args.prefix or not args.range:
        raise ValueError(
            "You must provide -f/--bam-file, -o/--output-dir, "
            "-p/--prefix, and -r/--range unless using --use-config."
        )

    return {
        "bam_file": args.bam_file,
        "output_dir": args.output_dir,
        "prefix": args.prefix,
        "chromosome": args.chromosome,
        "range": args.range,
        "clipping_side": args.clipping_side,
        "write_whole_fastas": args.write_whole_fastas,
    }


def validate_range(start, end):
    if start < 1 or end < start:
        raise ValueError("Range must satisfy 1 <= start <= end.")


def build_run_directory_name(prefix, chromosome, start, end):
    return f"{prefix}_{SCRIPT_LABEL}_{chromosome.lower()}{start}-{end}"


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


def main():
    args = build_parser().parse_args()

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

        print("Argument validation complete.")
        print(f"BAM file: {bam_path}")
        print(
            f"Region: {settings['chromosome']}:{start}-{end} "
            f"(clipping side: {settings['clipping_side']})"
        )
        print(f"Run directory to create: {run_metadata['run_dir']}")
        print(f"Run UID: {run_metadata['run_uid']}")
        print(f"Script version: {run_metadata['script_version']}")
        print(f"Script SHA256: {run_metadata['script_sha256']}")
        if run_metadata["git_commit"]:
            print(f"Git commit: {run_metadata['git_commit']}")
        else:
            print("Git commit: unavailable")

    except ValueError as exc:
        sys.exit(f"ERROR: {exc}")


if __name__ == "__main__":
    main()
