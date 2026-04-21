#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  c124_collect_breakpoint_evidence.sh <bam> <region> <partner_chr> <output_prefix> [min_clip_len]

Outputs:
  <output_prefix>.per_read.tsv
  <output_prefix>.summary.tsv
  <output_prefix>.sa_partner_positions.tsv
  <output_prefix>.softclipped_segments.tsv
  <output_prefix>.softclipped_segments.fasta
EOF
}

require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Required command not found: $cmd" >&2
    exit 1
  fi
}

if [[ $# -lt 4 || $# -gt 5 ]]; then
  usage
  exit 1
fi

require_command samtools
require_command awk

BAM="$1"
REGION="$2"
PARTNER_CHR="$3"
OUTPUT_PREFIX="$4"
MIN_CLIP_LEN="${5:-100}"

if [[ ! -f "$BAM" ]]; then
  echo "BAM not found: $BAM" >&2
  exit 1
fi

if [[ ! -f "${BAM}.bai" ]]; then
  echo "BAM index not found: ${BAM}.bai" >&2
  exit 1
fi

if [[ ! "$MIN_CLIP_LEN" =~ ^[0-9]+$ ]] || [[ "$MIN_CLIP_LEN" -lt 1 ]]; then
  echo "min_clip_len must be a positive integer: $MIN_CLIP_LEN" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

PER_READ_OUT="${OUTPUT_PREFIX}.per_read.tsv"
SUMMARY_OUT="${OUTPUT_PREFIX}.summary.tsv"
SA_POS_OUT="${OUTPUT_PREFIX}.sa_partner_positions.tsv"
SOFTCLIP_TSV_OUT="${OUTPUT_PREFIX}.softclipped_segments.tsv"
SOFTCLIP_FASTA_OUT="${OUTPUT_PREFIX}.softclipped_segments.fasta"

samtools view "$BAM" "$REGION" | awk \
  -v region="$REGION" \
  -v partner_chr="$PARTNER_CHR" \
  -v min_clip_len="$MIN_CLIP_LEN" \
  -v per_read_out="$PER_READ_OUT" \
  -v summary_out="$SUMMARY_OUT" \
  -v sa_pos_out="$SA_POS_OUT" \
  -v softclip_tsv_out="$SOFTCLIP_TSV_OUT" \
  -v softclip_fasta_out="$SOFTCLIP_FASTA_OUT" '
BEGIN {
  OFS = "\t"
  partner_label = "SA_to_" partner_chr

  print "read_name", "flag", "primary_chr", "primary_pos", "mapq", "cigar", \
        "has_softclip", "left_clip_len", "right_clip_len", \
        "has_sa", "class", "sa_chromosomes", "sa_positions", "sa_tag" > per_read_out

  print "read_name", "flag", "primary_chr", "primary_pos", "mapq", "cigar", \
        "clip_side", "clip_len", "clip_seq_len" > softclip_tsv_out
}

function append_unique(list, value, key) {
  if (value == "") {
    return list
  }
  if (!(key in seen_values)) {
    seen_values[key] = 1
    if (list == "") {
      return value
    }
    return list "," value
  }
  return list
}

function left_clip_len(cigar,    m) {
  if (match(cigar, /^[0-9]+S/)) {
    m = substr(cigar, RSTART, RLENGTH)
    sub(/S$/, "", m)
    return m + 0
  }
  return 0
}

function right_clip_len(cigar,    m) {
  if (match(cigar, /[0-9]+S$/)) {
    m = substr(cigar, RSTART, RLENGTH)
    sub(/S$/, "", m)
    return m + 0
  }
  return 0
}

{
  total_alignments++

  qname = $1
  flag = $2
  rname = $3
  pos = $4
  mapq = $5
  cigar = $6
  seq = $10

  has_sa = 0
  has_sa_to_partner = 0
  sa_value = ""
  sa_chromosomes = ""
  sa_positions = ""

  delete seen_values

  lclip = left_clip_len(cigar)
  rclip = right_clip_len(cigar)
  has_softclip = (lclip > 0 || rclip > 0)

  for (i = 12; i <= NF; i++) {
    if ($i ~ /^SA:Z:/) {
      has_sa = 1
      sa_value = substr($i, 6)
      n = split(sa_value, sa_entries, ";")
      for (j = 1; j <= n; j++) {
        if (sa_entries[j] == "") {
          continue
        }

        split(sa_entries[j], sa_fields, ",")
        sa_chr = sa_fields[1]
        sa_pos = sa_fields[2]

        sa_chromosomes = append_unique(sa_chromosomes, sa_chr, "chr" SUBSEP sa_chr)
        sa_positions = append_unique(sa_positions, sa_chr ":" sa_pos, "pos" SUBSEP sa_chr ":" sa_pos)

        chr_count[sa_chr]++
        pos_count[sa_chr SUBSEP sa_pos]++
        pos_reads[sa_chr SUBSEP sa_pos SUBSEP qname] = 1

        if (sa_chr == partner_chr) {
          has_sa_to_partner = 1
        }
      }
      break
    }
  }

  if (!(qname in seen_read)) {
    seen_read[qname] = 1
    unique_reads++
  }

  if (has_sa_to_partner) {
    class_name = partner_label
  } else if (has_sa) {
    class_name = "SA_to_other"
  } else if (has_softclip) {
    class_name = "softclipped_no_SA"
  } else {
    class_name = "other"
  }

  class_count[class_name]++

  unique_class_key = class_name SUBSEP qname
  if (!(unique_class_key in seen_unique_class)) {
    seen_unique_class[unique_class_key] = 1
    unique_class_count[class_name]++
  }

  print qname, flag, rname, pos, mapq, cigar, has_softclip, lclip, rclip, has_sa, \
        class_name, sa_chromosomes, sa_positions, sa_value > per_read_out

  if (lclip >= min_clip_len && seq != "*" && length(seq) >= lclip) {
    left_seq = substr(seq, 1, lclip)
    print qname, flag, rname, pos, mapq, cigar, "left", lclip, length(left_seq) >> softclip_tsv_out
    print ">" qname "|left|" rname ":" pos "|clip=" lclip >> softclip_fasta_out
    print left_seq >> softclip_fasta_out
  }

  if (rclip >= min_clip_len && seq != "*" && length(seq) >= rclip) {
    right_seq = substr(seq, length(seq) - rclip + 1, rclip)
    print qname, flag, rname, pos, mapq, cigar, "right", rclip, length(right_seq) >> softclip_tsv_out
    print ">" qname "|right|" rname ":" pos "|clip=" rclip >> softclip_fasta_out
    print right_seq >> softclip_fasta_out
  }
}

END {
  print "query_region", "metric", "count" > summary_out
  print region, "total_alignments", total_alignments + 0 >> summary_out
  print region, "unique_reads", unique_reads + 0 >> summary_out
  print region, "unique_reads_" partner_label, unique_class_count[partner_label] + 0 >> summary_out
  print region, "unique_reads_SA_to_other", unique_class_count["SA_to_other"] + 0 >> summary_out
  print region, "unique_reads_softclipped_no_SA", unique_class_count["softclipped_no_SA"] + 0 >> summary_out
  print region, "unique_reads_other", unique_class_count["other"] + 0 >> summary_out
  print region, "alignments_" partner_label, class_count[partner_label] + 0 >> summary_out
  print region, "alignments_SA_to_other", class_count["SA_to_other"] + 0 >> summary_out
  print region, "alignments_softclipped_no_SA", class_count["softclipped_no_SA"] + 0 >> summary_out
  print region, "alignments_other", class_count["other"] + 0 >> summary_out

  print "query_region", "partner_chr", "partner_pos", "supporting_sa_entries", "supporting_unique_reads" > sa_pos_out
  for (key in pos_count) {
    split(key, parts, SUBSEP)
    partner_chr_name = parts[1]
    partner_pos = parts[2]

    unique_support = 0
    for (read_key in pos_reads) {
      split(read_key, read_parts, SUBSEP)
      if (read_parts[1] == partner_chr_name && read_parts[2] == partner_pos) {
        unique_support++
      }
    }

    print region, partner_chr_name, partner_pos, pos_count[key], unique_support >> sa_pos_out
  }
}
'

echo "Per-read table: $PER_READ_OUT"
echo "Summary: $SUMMARY_OUT"
echo "SA partner positions: $SA_POS_OUT"
echo "Soft-clipped segments table: $SOFTCLIP_TSV_OUT"
echo "Soft-clipped segments FASTA: $SOFTCLIP_FASTA_OUT"
