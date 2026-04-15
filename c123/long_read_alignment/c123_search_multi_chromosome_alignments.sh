#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <input.bam> <output.tsv>" >&2
  exit 1
fi

BAM="$1"
OUTPUT="$2"

samtools view "$BAM" | awk '
BEGIN{OFS="\t"}
{
  read=$1
  chr=$3
  if (chr=="*") next
  print read, chr
}' | sort -u | awk '
BEGIN{OFS="\t"}
{
  if ($1 != prev_read && NR > 1) {
    if (nchr > 1) print prev_read, chrs
    chrs=$2
    nchr=1
  } else if ($1 == prev_read) {
    if (index("," chrs ",", "," $2 ",") == 0) {
      chrs=chrs "," $2
      nchr++
    }
  } else {
    chrs=$2
    nchr=1
  }
  prev_read=$1
}
END{
  if (NR > 0 && nchr > 1) print prev_read, chrs
}' > "$OUTPUT"

echo "Wrote $OUTPUT"
