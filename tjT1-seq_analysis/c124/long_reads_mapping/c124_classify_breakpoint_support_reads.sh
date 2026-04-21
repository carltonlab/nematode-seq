#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  c124_classify_breakpoint_support_reads.sh <bam> <region> <output_prefix>

Example:
  c124_classify_breakpoint_support_reads.sh \
    /path/to/long_reads.sorted.bam \
    V:8644820-8644860 \
    /path/to/output/v_breakpoint_support

Outputs:
  <output_prefix>.per_read.tsv
  <output_prefix>.summary.tsv
  <output_prefix>.partner_positions.tsv
EOF
}

require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Required command not found: $cmd" >&2
    exit 1
  fi
}

if [[ $# -ne 3 ]]; then
  usage
  exit 1
fi

require_command samtools
require_command awk

BAM="$1"
REGION="$2"
OUTPUT_PREFIX="$3"

if [[ ! -f "$BAM" ]]; then
  echo "BAM not found: $BAM" >&2
  exit 1
fi

if [[ ! -f "${BAM}.bai" ]]; then
  echo "BAM index not found: ${BAM}.bai" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

PER_READ_OUT="${OUTPUT_PREFIX}.per_read.tsv"
SUMMARY_OUT="${OUTPUT_PREFIX}.summary.tsv"
POSITIONS_OUT="${OUTPUT_PREFIX}.partner_positions.tsv"

samtools view "$BAM" "$REGION" | awk -v region="$REGION" -v per_read_out="$PER_READ_OUT" -v summary_out="$SUMMARY_OUT" -v positions_out="$POSITIONS_OUT" '
BEGIN {
  OFS = "\t"
  print "read_name", "flag", "primary_chr", "primary_pos", "mapq", "cigar", "has_softclip", "has_sa", "class", "sa_chromosomes", "sa_positions", "sa_tag" > per_read_out
}

function append_unique(list, value, seen_map_key) {
  if (value == "") {
    return list
  }
  if (!(seen_map_key in seen_values)) {
    seen_values[seen_map_key] = 1
    if (list == "") {
      return value
    }
    return list "," value
  }
  return list
}

{
  total_alignments++

  qname = $1
  flag = $2
  rname = $3
  pos = $4
  mapq = $5
  cigar = $6

  has_softclip = (cigar ~ /[0-9]+S/)
  has_sa = 0
  sa_value = ""
  sa_chromosomes = ""
  sa_positions = ""

  delete seen_values

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
      }
      break
    }
  }

  if (has_sa) {
    class_name = "has_SA"
  } else if (has_softclip) {
    class_name = "softclipped_only"
  } else {
    class_name = "local_only"
  }

  class_count[class_name]++

  if (!(qname in seen_read)) {
    seen_read[qname] = 1
    unique_reads++
  }

  unique_class_key = class_name SUBSEP qname
  if (!(unique_class_key in seen_unique_class)) {
    seen_unique_class[unique_class_key] = 1
    unique_class_count[class_name]++
  }

  print qname, flag, rname, pos, mapq, cigar, has_softclip, has_sa, class_name, sa_chromosomes, sa_positions, sa_value > per_read_out
}

END {
  print "query_region", "metric", "count" > summary_out
  print region, "total_alignments", total_alignments + 0 >> summary_out
  print region, "unique_reads", unique_reads + 0 >> summary_out
  print region, "unique_reads_has_SA", unique_class_count["has_SA"] + 0 >> summary_out
  print region, "unique_reads_softclipped_only", unique_class_count["softclipped_only"] + 0 >> summary_out
  print region, "unique_reads_local_only", unique_class_count["local_only"] + 0 >> summary_out
  print region, "alignments_has_SA", class_count["has_SA"] + 0 >> summary_out
  print region, "alignments_softclipped_only", class_count["softclipped_only"] + 0 >> summary_out
  print region, "alignments_local_only", class_count["local_only"] + 0 >> summary_out

  print "query_region", "partner_chr", "partner_pos", "supporting_sa_entries", "supporting_unique_reads" > positions_out
  for (key in pos_count) {
    split(key, parts, SUBSEP)
    partner_chr = parts[1]
    partner_pos = parts[2]

    unique_support = 0
    for (read_key in pos_reads) {
      split(read_key, read_parts, SUBSEP)
      if (read_parts[1] == partner_chr && read_parts[2] == partner_pos) {
        unique_support++
      }
    }

    print region, partner_chr, partner_pos, pos_count[key], unique_support >> positions_out
  }
}
'

echo "Per-read classification: $PER_READ_OUT"
echo "Summary: $SUMMARY_OUT"
echo "Partner positions: $POSITIONS_OUT"
