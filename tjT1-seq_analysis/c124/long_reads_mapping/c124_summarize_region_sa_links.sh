#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  c124_summarize_region_sa_links.sh <bam> <region> <output_prefix> [chromosomes]

Example:
  c124_summarize_region_sa_links.sh \
    /path/to/long_reads.sorted.bam \
    V:8644000-8646000 \
    /path/to/output/c124_V_8644k_8646k \
    I,II,III,IV,V,X

Outputs:
  <output_prefix>.reads_in_region.tsv
  <output_prefix>.chromosome_summary.tsv
  <output_prefix>.partner_positions.tsv

Definitions:
  reads_in_region.tsv
    One row per alignment overlapping <region>, with one 0/1 column per chromosome
    indicating whether that chromosome appears in the read's SA tag.

  chromosome_summary.tsv
    Aggregate counts for the queried region, including total alignments, unique reads,
    reads with any SA tag, and counts of reads/SA entries per chromosome.

  partner_positions.tsv
    Counts of SA-supported partner positions by chromosome and coordinate.

If [chromosomes] is omitted, the script uses:
  I,II,III,IV,V,X
EOF
}

require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Required command not found: $cmd" >&2
    exit 1
  fi
}

if [[ $# -lt 3 || $# -gt 4 ]]; then
  usage
  exit 1
fi

require_command samtools
require_command awk
require_command sort

BAM="$1"
REGION="$2"
OUTPUT_PREFIX="$3"
CHROMS="${4:-I,II,III,IV,V,X}"

if [[ ! -f "$BAM" ]]; then
  echo "BAM not found: $BAM" >&2
  exit 1
fi

if [[ ! -f "${BAM}.bai" ]]; then
  echo "BAM index not found: ${BAM}.bai" >&2
  exit 1
fi

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

READS_OUT="${OUTPUT_PREFIX}.reads_in_region.tsv"
SUMMARY_OUT="${OUTPUT_PREFIX}.chromosome_summary.tsv"
POSITIONS_OUT="${OUTPUT_PREFIX}.partner_positions.tsv"
TMP_SAM="${OUTPUT_PREFIX}.region.tmp.sam"

samtools view "$BAM" "$REGION" >"$TMP_SAM"

awk -v chroms="$CHROMS" -v region="$REGION" -v reads_out="$READS_OUT" -v summary_out="$SUMMARY_OUT" -v positions_out="$POSITIONS_OUT" '
BEGIN {
  OFS = "\t"
  nchrom = split(chroms, chrom, ",")

  print "read_name", "flag", "primary_chr", "primary_pos", "mapq", "cigar", "has_sa", "sa_tag", build_chr_headers() > reads_out
  print "query_region", "chromosome", "unique_reads_with_sa_to_chr", "sa_entries_to_chr" > summary_out
  print "query_region", "partner_chr", "partner_pos", "supporting_sa_entries", "supporting_unique_reads" > positions_out
}

function build_chr_headers(    out, i) {
  out = ""
  for (i = 1; i <= nchrom; i++) {
    out = out OFS "sa_has_" chrom[i]
  }
  return out
}

{
  total_alignments++
  qname = $1
  flag = $2
  rname = $3
  pos = $4
  mapq = $5
  cigar = $6
  sa_value = ""
  has_sa = 0

  if (!(qname in seen_read)) {
    seen_read[qname] = 1
    total_unique_reads++
  }

  delete sa_has
  for (i = 1; i <= nchrom; i++) {
    sa_has[chrom[i]] = 0
  }

  for (i = 12; i <= NF; i++) {
    if ($i ~ /^SA:Z:/) {
      has_sa = 1
      sa_value = substr($i, 6)
      break
    }
  }

  if (has_sa) {
    reads_with_sa_alignments++
    if (!(qname in seen_sa_read)) {
      seen_sa_read[qname] = 1
      unique_reads_with_sa++
    }

    n = split(sa_value, sa_entries, ";")
    for (j = 1; j <= n; j++) {
      if (sa_entries[j] == "") {
        continue
      }
      split(sa_entries[j], sa_fields, ",")
      partner_chr = sa_fields[1]
      partner_pos = sa_fields[2]

      sa_entries_total++
      pos_key = partner_chr OFS partner_pos
      pos_count[pos_key]++
      pos_reads[pos_key, qname] = 1

      for (k = 1; k <= nchrom; k++) {
        chr_name = chrom[k]
        if (partner_chr == chr_name) {
          sa_has[chr_name] = 1
          sa_entries_by_chr[chr_name]++
          read_chr_key = chr_name SUBSEP qname
          if (!(read_chr_key in seen_chr_read)) {
            seen_chr_read[read_chr_key] = 1
            unique_reads_by_chr[chr_name]++
          }
        }
      }
    }
  }

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", qname, flag, rname, pos, mapq, cigar, has_sa, sa_value > reads_out
  for (i = 1; i <= nchrom; i++) {
    printf "\t%d", sa_has[chrom[i]] > reads_out
  }
  printf "\n" > reads_out
}

END {
  print region, "ALL", total_unique_reads, total_alignments >> summary_out
  print region, "ANY_SA", unique_reads_with_sa + 0, reads_with_sa_alignments + 0 >> summary_out
  for (i = 1; i <= nchrom; i++) {
    chr_name = chrom[i]
    print region, chr_name, unique_reads_by_chr[chr_name] + 0, sa_entries_by_chr[chr_name] + 0 >> summary_out
  }

  for (key in pos_count) {
    split(key, fields, OFS)
    partner_chr = fields[1]
    partner_pos = fields[2]
    unique_support = 0
    for (read_key in pos_reads) {
      split(read_key, pair, SUBSEP)
      if (pair[1] == key) {
        unique_support++
      }
    }
    print region, partner_chr, partner_pos, pos_count[key], unique_support >> positions_out
  }
}
' "$TMP_SAM"

rm -f "$TMP_SAM"

echo "Reads in region: $READS_OUT"
echo "Chromosome summary: $SUMMARY_OUT"
echo "Partner positions: $POSITIONS_OUT"
