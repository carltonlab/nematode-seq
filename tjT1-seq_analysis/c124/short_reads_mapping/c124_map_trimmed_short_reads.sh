#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  c124_map_trimmed_short_reads.sh --init-config <config.sh>
  c124_map_trimmed_short_reads.sh --config <config.sh>

Modes:
  --init-config   Create a template config file for editing.
  --config        Run short-read mapping using values from the config file.

Config variables:
  SAMPLE_NAME     Sample label used for output naming.
  R1              Path to trimmed read 1 FASTQ.
  R2              Path to trimmed read 2 FASTQ.
  REF             Path to reference FASTA.
  OUTDIR          Directory for BAM and QC outputs.
  THREADS         Number of threads to use.
  ALIGNER         bwa or bwa-mem2.
  READ_GROUP_ID   Read group ID.
  READ_GROUP_SM   Read group sample name.
  READ_GROUP_LB   Read group library name.
  READ_GROUP_PL   Read group platform, e.g. ILLUMINA.
  SORT_MEMORY     samtools sort memory per thread, e.g. 1G.
EOF
}

init_config() {
  local config_path="$1"

  if [[ -e "$config_path" ]]; then
    echo "Refusing to overwrite existing file: $config_path" >&2
    exit 1
  fi

  cat >"$config_path" <<'EOF'
#!/usr/bin/env bash

# Fill in all required values before running the mapping script.

SAMPLE_NAME=""
R1=""
R2=""
REF=""
OUTDIR=""

# Tooling
THREADS=8
ALIGNER="bwa-mem2"
SORT_MEMORY="1G"

# Read group
READ_GROUP_ID=""
READ_GROUP_SM=""
READ_GROUP_LB=""
READ_GROUP_PL="ILLUMINA"
EOF

  echo "Created config template: $config_path"
  echo "Edit the file, then run:"
  echo "  $(basename "$0") --config $config_path"
}

require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Required command not found: $cmd" >&2
    exit 1
  fi
}

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "$path" ]]; then
    echo "$label does not exist: $path" >&2
    exit 1
  fi
}

require_var() {
  local name="$1"
  if [[ -z "${!name:-}" ]]; then
    echo "Required config variable is empty: $name" >&2
    exit 1
  fi
}

ensure_bwa_index() {
  local ref="$1"
  if [[ -f "${ref}.bwt" && -f "${ref}.pac" && -f "${ref}.ann" && -f "${ref}.amb" && -f "${ref}.sa" ]]; then
    echo "BWA index already present for $ref"
    return
  fi

  echo "Creating BWA index for $ref"
  bwa index "$ref"
}

ensure_bwa_mem2_index() {
  local ref="$1"
  if [[ -f "${ref}.bwt.2bit.64" ]]; then
    echo "bwa-mem2 index already present for $ref"
    return
  fi

  echo "Creating bwa-mem2 index for $ref"
  bwa-mem2 index "$ref"
}

run_mapping() {
  local config_path="$1"

  if [[ ! -f "$config_path" ]]; then
    echo "Config file not found: $config_path" >&2
    exit 1
  fi

  # shellcheck disable=SC1090
  source "$config_path"

  require_var SAMPLE_NAME
  require_var R1
  require_var R2
  require_var REF
  require_var OUTDIR
  require_var THREADS
  require_var ALIGNER
  require_var READ_GROUP_ID
  require_var READ_GROUP_SM
  require_var READ_GROUP_LB
  require_var READ_GROUP_PL
  require_var SORT_MEMORY

  require_file "$R1" "R1 FASTQ"
  require_file "$R2" "R2 FASTQ"
  require_file "$REF" "Reference FASTA"
  mkdir -p "$OUTDIR"

  require_command samtools

  local bam_prefix="${OUTDIR}/${SAMPLE_NAME}_vs_$(basename "${REF%.*}")"
  local sorted_bam="${bam_prefix}.sorted.bam"
  local flagstat_txt="${bam_prefix}.flagstat.txt"
  local idxstats_txt="${bam_prefix}.idxstats.txt"
  local read_group

  read_group="@RG\tID:${READ_GROUP_ID}\tSM:${READ_GROUP_SM}\tLB:${READ_GROUP_LB}\tPL:${READ_GROUP_PL}"

  case "$ALIGNER" in
    bwa)
      require_command bwa
      ensure_bwa_index "$REF"
      echo "Running alignment with bwa mem"
      bwa mem -t "$THREADS" -R "$read_group" "$REF" "$R1" "$R2" \
        | samtools sort -@ "$THREADS" -m "$SORT_MEMORY" -o "$sorted_bam" -
      ;;
    bwa-mem2)
      require_command bwa-mem2
      ensure_bwa_mem2_index "$REF"
      echo "Running alignment with bwa-mem2 mem"
      bwa-mem2 mem -t "$THREADS" -R "$read_group" "$REF" "$R1" "$R2" \
        | samtools sort -@ "$THREADS" -m "$SORT_MEMORY" -o "$sorted_bam" -
      ;;
    *)
      echo "Unsupported ALIGNER: $ALIGNER" >&2
      echo "Supported values: bwa, bwa-mem2" >&2
      exit 1
      ;;
  esac

  samtools index "$sorted_bam"
  samtools quickcheck "$sorted_bam"
  samtools flagstat "$sorted_bam" >"$flagstat_txt"
  samtools idxstats "$sorted_bam" >"$idxstats_txt"

  echo "Finished."
  echo "Sorted BAM: $sorted_bam"
  echo "BAM index: ${sorted_bam}.bai"
  echo "Flagstat: $flagstat_txt"
  echo "Idxstats: $idxstats_txt"
}

if [[ $# -ne 2 ]]; then
  usage
  exit 1
fi

case "$1" in
  --init-config)
    init_config "$2"
    ;;
  --config)
    run_mapping "$2"
    ;;
  -h|--help)
    usage
    ;;
  *)
    usage
    exit 1
    ;;
esac
