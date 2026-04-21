#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  c124_map_clipped_fasta_to_reference.sh <softclip_fasta> <reference_fasta> <partner_chr> <output_prefix> [aligner_preset] [min_mapq]
EOF
}

require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Required command not found: $cmd" >&2
    exit 1
  fi
}

if [[ $# -lt 4 || $# -gt 6 ]]; then
  usage
  exit 1
fi

require_command minimap2
require_command samtools
require_command awk
require_command sort

SOFTCLIP_FASTA="$1"
REFERENCE_FASTA="$2"
PARTNER_CHR="$3"
OUTPUT_PREFIX="$4"
ALIGNER_PRESET="${5:-map-ont}"
MIN_MAPQ="${6:-20}"

if [[ ! -f "$SOFTCLIP_FASTA" ]]; then
  echo "Soft-clipped FASTA not found: $SOFTCLIP_FASTA" >&2
  exit 1
fi

if [[ ! -f "$REFERENCE_FASTA" ]]; then
  echo "Reference FASTA not found: $REFERENCE_FASTA" >&2
  exit 1
fi

if [[ ! "$MIN_MAPQ" =~ ^[0-9]+$ ]]; then
  echo "min_mapq must be a non-negative integer: $MIN_MAPQ" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

SORTED_BAM="${OUTPUT_PREFIX}.sorted.bam"
IDXSTATS_OUT="${OUTPUT_PREFIX}.idxstats.txt"
FLAGSTAT_OUT="${OUTPUT_PREFIX}.flagstat.txt"
CHROM_COUNTS_OUT="${OUTPUT_PREFIX}.chromosome_counts.tsv"
PARTNER_HIGH_MAPQ_OUT="${OUTPUT_PREFIX}.${PARTNER_CHR}_high_mapq.tsv"
PARTNER_POSITIONS_OUT="${OUTPUT_PREFIX}.${PARTNER_CHR}_positions.tsv"

minimap2 -a -x "$ALIGNER_PRESET" "$REFERENCE_FASTA" "$SOFTCLIP_FASTA" \
  | samtools sort -o "$SORTED_BAM" -

samtools index "$SORTED_BAM"
samtools flagstat "$SORTED_BAM" >"$FLAGSTAT_OUT"
samtools idxstats "$SORTED_BAM" >"$IDXSTATS_OUT"

samtools view "$SORTED_BAM" | awk '
BEGIN {
  OFS = "\t"
  print "chromosome", "mapped_segments"
}
$3 != "*" {
  chr_count[$3]++
}
END {
  for (chr in chr_count) {
    print chr, chr_count[chr]
  }
}
' | sort -k2,2nr -k1,1 >"$CHROM_COUNTS_OUT"

samtools view "$SORTED_BAM" "$PARTNER_CHR" | awk -v min_mapq="$MIN_MAPQ" '
BEGIN {
  OFS = "\t"
  print "read_name", "flag", "chr", "pos", "mapq", "cigar"
}
$5 >= min_mapq {
  print $1, $2, $3, $4, $5, $6
}
' >"$PARTNER_HIGH_MAPQ_OUT"

samtools view "$SORTED_BAM" "$PARTNER_CHR" | awk -v min_mapq="$MIN_MAPQ" -v partner_chr="$PARTNER_CHR" '
$5 >= min_mapq {
  pos_count[$4]++
  read_seen[$4, $1] = 1
}
END {
  OFS = "\t"
  print "chr", "pos", "supporting_alignments", "supporting_unique_reads"
  for (pos in pos_count) {
    unique_reads = 0
    for (k in read_seen) {
      split(k, parts, SUBSEP)
      if (parts[1] == pos) {
        unique_reads++
      }
    }
    print partner_chr, pos, pos_count[pos], unique_reads
  }
}
' | sort -k4,4nr -k3,3nr -k2,2n >"$PARTNER_POSITIONS_OUT"

echo "Sorted BAM: $SORTED_BAM"
echo "BAM index: ${SORTED_BAM}.bai"
echo "Flagstat: $FLAGSTAT_OUT"
echo "Idxstats: $IDXSTATS_OUT"
echo "Chromosome counts: $CHROM_COUNTS_OUT"
echo "Partner chromosome high-MAPQ hits: $PARTNER_HIGH_MAPQ_OUT"
echo "Partner chromosome positions: $PARTNER_POSITIONS_OUT"
