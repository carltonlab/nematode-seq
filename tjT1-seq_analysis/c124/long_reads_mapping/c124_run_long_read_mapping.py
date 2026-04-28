#!/usr/bin/env python3
import argparse
import configparser
import os
import subprocess
import sys
from pathlib import Path


def build_parser():
    parser = argparse.ArgumentParser(
        description="Launch the long-read minimap2 mapping workflow."
    )
    parser.add_argument(
        "--snakefile",
        help=(
            "Path to the Snakemake workflow file. "
            "If omitted, the script looks in its own directory."
        ),
    )
    parser.add_argument(
        "-i",
        "--inputs",
        nargs="+",
        help="One or more long-read FASTQ files to map.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        help="Reference FASTA file.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        help="Output directory for long-read mapping results.",
    )
    parser.add_argument(
        "--config-file",
        help="INI configuration file describing inputs, reference, and output directory.",
    )
    parser.add_argument(
        "--create-config",
        help="Write a sample INI config file to the given path and exit.",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="Number of Snakemake cores and minimap2 threads to use.",
    )
    parser.add_argument(
        "--preset",
        default="map-ont",
        help="Minimap2 preset to use. Default: map-ont",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run Snakemake in dry-run mode.",
    )
    return parser


def parse_config_file(config_path):
    parser = configparser.ConfigParser()
    read_files = parser.read(config_path)
    if not read_files:
        raise ValueError(f"Could not read config file: {config_path}")
    if "workflow" not in parser:
        raise ValueError("Config file must contain a [workflow] section")

    workflow = parser["workflow"]
    inputs_value = workflow.get("inputs", "").strip()
    reference = workflow.get("reference", "").strip()
    output_dir = workflow.get("output_dir", "").strip()
    cores = workflow.get("cores", "").strip()
    preset = workflow.get("preset", "map-ont").strip()
    dry_run = workflow.get("dry_run", "").strip().lower()

    inputs = [item.strip() for item in inputs_value.split(",") if item.strip()]
    if not inputs:
        raise ValueError("Config file must define workflow.inputs")
    if not reference:
        raise ValueError("Config file must define workflow.reference")
    if not output_dir:
        raise ValueError("Config file must define workflow.output_dir")

    return {
        "inputs": inputs,
        "reference": reference,
        "output_dir": output_dir,
        "cores": int(cores) if cores else 1,
        "preset": preset,
        "dry_run": dry_run in {"1", "true", "yes", "on"},
    }


def create_sample_config(config_path):
    path = Path(config_path)
    if path.exists():
        raise ValueError(f"Refusing to overwrite existing config file: {config_path}")

    content = """[workflow]
inputs = /path/to/M5_merge.fastq.gz
reference = /path/to/reference.fa
output_dir = /path/to/long_read_mapping_results
cores = 16
preset = map-ont
dry_run = true
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


def validate_args(args):
    if args.create_config and (
        args.config_file is not None
        or args.inputs is not None
        or args.reference is not None
        or args.output_dir is not None
        or args.cores != 1
        or args.preset != "map-ont"
        or args.dry_run
    ):
        raise ValueError(
            "You cannot use any other flags together with --create-config."
        )
    if args.create_config:
        return {"create_config": args.create_config}

    if args.config_file and (
        args.inputs is not None
        or args.reference is not None
        or args.output_dir is not None
        or args.cores != 1
        or args.preset != "map-ont"
        or args.create_config is not None
        or args.dry_run
    ):
        raise ValueError(
            "You cannot use any other flags together with --config-file."
        )

    if args.config_file:
        settings = parse_config_file(args.config_file)
        settings["snakefile"] = args.snakefile
        return settings

    if not args.inputs or not args.reference or not args.output_dir:
        raise ValueError(
            "Without --config-file, you must provide -i/--inputs, "
            "-r/--reference, and -o/--output-dir."
        )

    return {
        "snakefile": args.snakefile,
        "inputs": args.inputs,
        "reference": args.reference,
        "output_dir": args.output_dir,
        "cores": args.cores,
        "preset": args.preset,
        "dry_run": args.dry_run,
    }


def resolve_snakefile(user_path=None):
    if user_path:
        snakefile = Path(user_path).resolve()
        if not snakefile.is_file():
            raise ValueError(f"Snakemake workflow file does not exist: {snakefile}")
        return snakefile

    script_dir = Path(__file__).resolve().parent
    candidates = [
        script_dir / "Snakefile",
        script_dir / "snakefile",
        script_dir / "c124_long_read_mapping_Snakemake",
    ]

    for candidate in candidates:
        if candidate.is_file():
            return candidate

    raise ValueError(
        "Could not find a Snakemake workflow file next to the runner script. "
        "Use --snakefile to provide one explicitly."
    )


def build_snakemake_command(settings, snakefile):
    command = [
        "snakemake",
        "-s",
        str(snakefile),
        "--cores",
        str(settings["cores"]),
        "--config",
        f"input_fastqs={','.join(settings['inputs'])}",
        f"reference={settings['reference']}",
        f"output_dir={settings['output_dir']}",
        f"cores={settings['cores']}",
        f"preset={settings['preset']}",
    ]
    if settings["dry_run"]:
        command.append("-n")
    return command


def main():
    args = build_parser().parse_args()
    settings = validate_args(args)
    if "create_config" in settings:
        create_sample_config(settings["create_config"])
        return

    resolved_inputs = [Path(path).resolve() for path in settings["inputs"]]
    for path in resolved_inputs:
        if not path.is_file():
            raise ValueError(f"Input FASTQ does not exist: {path}")
    reference = Path(settings["reference"]).resolve()
    if not reference.is_file():
        raise ValueError(f"Reference FASTA does not exist: {reference}")

    settings["inputs"] = [str(path) for path in resolved_inputs]
    settings["reference"] = str(reference)
    settings["output_dir"] = str(Path(settings["output_dir"]).resolve())

    Path(settings["output_dir"]).mkdir(parents=True, exist_ok=True)
    snakefile = resolve_snakefile(settings.get("snakefile"))

    snakemake_env = {
        "XDG_CACHE_HOME": "/tmp/.cache",
        "SNAKEMAKE_OUTPUT_CACHE": "/tmp/.snakemake",
    }

    command = build_snakemake_command(settings, snakefile)
    subprocess.run(command, check=True, env={**os.environ, **snakemake_env})


if __name__ == "__main__":
    try:
        main()
    except ValueError as exc:
        print(f"ValueError: {exc}", file=sys.stderr)
        sys.exit(1)
