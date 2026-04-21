#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  c124_blast_clipped_fasta_to_reference.sh <query_fasta> <blast_db> <output_dir> [max_target_seqs] [num_threads]

Outputs in <output_dir>:
  blast_hits.tsv
  blast_hits_summary.tsv

Notes:
  - <blast_db> should be the value you would pass directly to blastn -db
  - Output format is tabular with standard alignment fields
EOF
}

require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Required command not found: $cmd" >&2
    exit 1
  fi
}

if [[ $# -lt 3 || $# -gt 5 ]]; then
  usage
  exit 1
fi

require_command blastn
require_command awk
require_command sort

QUERY_FASTA="$1"
BLAST_DB="$2"
OUTPUT_DIR="$3"
MAX_TARGET_SEQS="${4:-20}"
NUM_THREADS="${5:-1}"

if [[ ! -f "$QUERY_FASTA" ]]; then
  echo "Query FASTA not found: $QUERY_FASTA" >&2
  exit 1
fi

if [[ ! "$MAX_TARGET_SEQS" =~ ^[0-9]+$ ]] || [[ "$MAX_TARGET_SEQS" -lt 1 ]]; then
  echo "max_target_seqs must be a positive integer: $MAX_TARGET_SEQS" >&2
  exit 1
fi

if [[ ! "$NUM_THREADS" =~ ^[0-9]+$ ]] || [[ "$NUM_THREADS" -lt 1 ]]; then
  echo "num_threads must be a positive integer: $NUM_THREADS" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

HITS_OUT="${OUTPUT_DIR}/blast_hits.tsv"
SUMMARY_OUT="${OUTPUT_DIR}/blast_hits_summary.tsv"

blastn \
  -query "$QUERY_FASTA" \
  -db "$BLAST_DB" \
  -task blastn \
  -max_target_seqs "$MAX_TARGET_SEQS" \
  -num_threads "$NUM_THREADS" \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
  -out "$HITS_OUT"

{
  printf "subject_id\thit_count\tunique_queries\tbest_bitscore\tbest_evalue\n"
  awk '
  {
    subject = $2
    query = $1
    bitscore = $12 + 0
    evalue = $11

    hit_count[subject]++
    query_seen[subject, query] = 1

    if (!(subject in best_bitscore) || bitscore > best_bitscore[subject]) {
      best_bitscore[subject] = bitscore
      best_evalue[subject] = evalue
    }
  }
  END {
    OFS = "\t"
    for (subject in hit_count) {
      unique_queries = 0
      for (k in query_seen) {
        split(k, parts, SUBSEP)
        if (parts[1] == subject) {
          unique_queries++
        }
      }
      print subject, hit_count[subject], unique_queries, best_bitscore[subject], best_evalue[subject]
    }
  }
  ' "$HITS_OUT" | sort -k3,3nr -k2,2nr -k4,4nr -k1,1
} >"$SUMMARY_OUT"

echo "BLAST hits: $HITS_OUT"
echo "BLAST summary: $SUMMARY_OUT"
