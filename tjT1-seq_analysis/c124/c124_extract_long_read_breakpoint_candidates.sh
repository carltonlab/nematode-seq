#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  c124_extract_long_read_breakpoint_candidates.sh <bam> <region> <partner_chr> <output_prefix>

Example:
  c124_extract_long_read_breakpoint_candidates.sh \
    /path/to/long_reads.sorted.bam \
    V:9500000-9600000 \
    IV \
    /path/to/output/c124_V_breakpoint_to_IV

Outputs:
  <output_prefix>.candidate_read_names.txt
  <output_prefix>.candidate_summary.tsv
  <output_prefix>.candidate_alignments.sam
  <output_prefix>.candidate_reads.fasta

Candidate reads are those with an alignment overlapping <region> and an SA tag
that points to <partner_chr>. The partner chromosome position is discovered from
the SA tag and reported in the summary table.
EOF
}

require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Required command not found: $cmd" >&2
    exit 1
  fi
}

if [[ $# -ne 4 ]]; then
  usage
  exit 1
fi

require_command samtools
require_command awk
require_command sort

BAM="$1"
REGION="$2"
PARTNER_CHR="$3"
OUTPUT_PREFIX="$4"

if [[ ! -f "$BAM" ]]; then
  echo "BAM not found: $BAM" >&2
  exit 1
fi

if [[ ! -f "${BAM}.bai" ]]; then
  echo "BAM index not found: ${BAM}.bai" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

NAMES_OUT="${OUTPUT_PREFIX}.candidate_read_names.txt"
SUMMARY_OUT="${OUTPUT_PREFIX}.candidate_summary.tsv"
ALIGNMENTS_OUT="${OUTPUT_PREFIX}.candidate_alignments.sam"
FASTA_OUT="${OUTPUT_PREFIX}.candidate_reads.fasta"
TMP_REGION_SAM="${OUTPUT_PREFIX}.region.tmp.sam"

samtools view "$BAM" "$REGION" >"$TMP_REGION_SAM"

awk -v partner="$PARTNER_CHR" '
BEGIN {
  OFS = "\t"
  print "read_name", "flag", "region_chr", "region_pos", "mapq", "cigar", "partner_chr", "partner_pos", "partner_strand", "partner_cigar", "partner_mapq", "sa_tag"
}
{
  qname = $1
  flag = $2
  rname = $3
  pos = $4
  mapq = $5
  cigar = $6
  sa_value = ""

  for (i = 12; i <= NF; i++) {
    if ($i ~ /^SA:Z:/) {
      sa_value = substr($i, 6)
      n = split(sa_value, sa_entries, ";")
      for (j = 1; j <= n; j++) {
        if (sa_entries[j] == "") {
          continue
        }
        split(sa_entries[j], sa_fields, ",")
        if (sa_fields[1] == partner) {
          print qname, flag, rname, pos, mapq, cigar, sa_fields[1], sa_fields[2], sa_fields[3], sa_fields[4], sa_fields[5], sa_value
        }
      }
    }
  }
}
' "$TMP_REGION_SAM" >"$SUMMARY_OUT"

tail -n +2 "$SUMMARY_OUT" | cut -f1 | sort -u >"$NAMES_OUT"

if [[ ! -s "$NAMES_OUT" ]]; then
  samtools view -H "$BAM" >"$ALIGNMENTS_OUT"
  : >"$FASTA_OUT"
  rm -f "$TMP_REGION_SAM"
  echo "No candidate reads found for $REGION with SA tags to $PARTNER_CHR"
  echo "Created outputs:"
  echo "  $NAMES_OUT"
  echo "  $SUMMARY_OUT"
  echo "  $ALIGNMENTS_OUT"
  echo "  $FASTA_OUT"
  exit 0
fi

samtools view -h "$BAM" | awk '
NR == FNR {
  keep[$1] = 1
  next
}
/^@/ {
  print
  next
}
keep[$1]
' "$NAMES_OUT" - >"$ALIGNMENTS_OUT"

samtools view "$BAM" "$REGION" | awk '
NR == FNR {
  keep[$1] = 1
  next
}
keep[$1] && $10 != "*" && !seen[$1]++ {
  print ">" $1
  print $10
}
' "$NAMES_OUT" - >"$FASTA_OUT"

rm -f "$TMP_REGION_SAM"

echo "Candidate read names: $NAMES_OUT"
echo "Candidate summary: $SUMMARY_OUT"
echo "Candidate alignments: $ALIGNMENTS_OUT"
echo "Candidate FASTA: $FASTA_OUT"
