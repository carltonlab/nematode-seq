#!/bin/bash
# Find truly unique mutations with mutation type filtering
# Rules:
# 1. Only count mutations with count > 1
# 2. All insertions (any length) = one mutation type '+'
# 3. All deletions (any length) = one mutation type '-'
# 4. Mutation is unique if that TYPE is absent in ancestor
# 5. Filter out positions with 3+ total mutation types

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <evolved_file.txt> <ancestor_file.txt> <output_prefix>"
    echo "Example: $0 generation_50.txt generation_1.txt unique_mutations"
    exit 1
fi

EVOLVED_FILE="$1"
ANCESTOR_FILE="$2"
OUTPUT_PREFIX="$3"

if [ ! -f "$EVOLVED_FILE" ]; then
    echo "Error: Evolved file '$EVOLVED_FILE' not found!"
    exit 1
fi

if [ ! -f "$ANCESTOR_FILE" ]; then
    echo "Error: Ancestor file '$ANCESTOR_FILE' not found!"
    exit 1
fi

echo "Finding truly unique mutations..."
echo "Rules:"
echo "  - Only count mutations with count > 1"
echo "  - All insertions (+1, +2, +4...) = one type '+'"
echo "  - All deletions (-1, -2, -3...) = one type '-'"
echo "  - Filter positions with 3+ mutation types"
echo "================================"
echo "Evolved sample: $EVOLVED_FILE"
echo "Ancestor sample: $ANCESTOR_FILE"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""

OUTPUT_UNIQUE="${OUTPUT_PREFIX}_truly_unique_mutations.txt"
OUTPUT_SUMMARY="${OUTPUT_PREFIX}_summary.txt"
TEMP_FILE=$(mktemp)

awk '
BEGIN {
    FS = "\t"
    OFS = "\t"
    print "Processing files..." > "/dev/stderr"
}

# Function to count mutation types at a position
function count_mutation_types(bases_found, insertions, deletions) {
    mutation_types = 0
    
    # Count number of different bases
    n_bases = split(bases_found, bases_array, ",")
    mutation_types += n_bases
    
    # Check for ANY insertions (regardless of count)
    if (insertions != "" && insertions != "-" && insertions !~ /^-IN/ && insertions !~ /^-NN/ && insertions !~ /COMPLEX/) {
        mutation_types += 1
    }
    
    # Check for ANY deletions (regardless of count)
    if (deletions != "" && deletions != "-" && deletions !~ /^-IN/ && deletions !~ /^-NN/ && deletions !~ /COMPLEX/) {
        mutation_types += 1
    }
    
    return mutation_types
}

# Function to check if there are ANY insertions with count > 1
function has_insertions_gt1(indel_string) {
    if (indel_string == "" || indel_string == "-" || indel_string ~ /^-IN/ || indel_string ~ /^-NN/ || indel_string ~ /COMPLEX/) {
        return 0
    }
    temp = indel_string
    while (match(temp, /\+[0-9]+[ACGTN]+\(([0-9]+)\)/, arr)) {
        if (arr[1] > 1) return 1
        sub(/\+[0-9]+[ACGTN]+\([0-9]+\)/, "", temp)
    }
    return 0
}

# Function to check if there are ANY deletions with count > 1
function has_deletions_gt1(indel_string) {
    if (indel_string == "" || indel_string == "-" || indel_string ~ /^-IN/ || indel_string ~ /^-NN/ || indel_string ~ /COMPLEX/) {
        return 0
    }
    temp = indel_string
    while (match(temp, /-[0-9]+[ACGTN]+\(([0-9]+)\)/, arr)) {
        if (arr[1] > 1) return 1
        sub(/-[0-9]+[ACGTN]+\([0-9]+\)/, "", temp)
    }
    return 0
}

# Function to extract bases and their counts
function parse_bases(bases_found, base_counts, bases_array, counts_array) {
    delete bases_array
    delete counts_array
    
    n_bases = split(bases_found, bases_list, ",")
    
    temp = base_counts
    for (i = 1; i <= n_bases; i++) {
        base = bases_list[i]
        if (match(temp, base "\\(([0-9]+)\\)", arr)) {
            bases_array[base] = 1
            counts_array[base] = arr[1]
        }
    }
}

# Read ancestor file
NR == FNR {
    if ($1 == "CHROM") next
    
    chrom = $1
    pos = $2
    bases_found = $3
    base_counts = $4
    insertions = $5
    deletions = $6
    
    key = chrom":"pos
    
    # Parse bases in ancestor (only count if > 1)
    parse_bases(bases_found, base_counts, ancestor_bases, ancestor_counts)
    for (base in ancestor_bases) {
        if (ancestor_counts[base] > 1) {
            ancestor_has[key":"base] = 1
            ancestor_count[key":"base] = ancestor_counts[base]
        }
    }
    
    # Check if ancestor has ANY insertions with count > 1
    if (has_insertions_gt1(insertions)) {
        ancestor_has_ins[key] = 1
    }
    
    # Check if ancestor has ANY deletions with count > 1
    if (has_deletions_gt1(deletions)) {
        ancestor_has_del[key] = 1
    }
    
    ancestor_data[key] = $0
    next
}

# Read evolved file
{
    if ($1 == "CHROM") next
    
    chrom = $1
    pos = $2
    bases_found = $3
    base_counts = $4
    insertions = $5
    deletions = $6
    zygosity = $7
    
    key = chrom":"pos
    
    # FIRST: Check if this position has ≤2 mutation types
    total_mutation_types = count_mutation_types(bases_found, insertions, deletions)
    
    if (total_mutation_types > 2) {
        filtered_complex++
        next
    }
    
    # Parse bases in evolved sample
    parse_bases(bases_found, base_counts, evolved_bases, evolved_counts)
    
    # Check for truly new mutations
    found_unique = 0
    unique_snps = ""
    unique_insertions = ""
    unique_deletions = ""
    
    snp_count = 0
    ins_count = 0
    del_count = 0
    
    # Check each base in evolved sample (only if count > 1)
    for (base in evolved_bases) {
        count = evolved_counts[base]
        
        if (count > 1) {
            # Check if this base is ABSENT in ancestor (or has count = 1 in ancestor)
            if (!(key":"base in ancestor_has)) {
                found_unique = 1
                if (unique_snps != "") unique_snps = unique_snps ","
                unique_snps = unique_snps base "(" count ")"
                snp_count++
            }
        }
    }
    
    # Check if evolved has insertions (count > 1) AND ancestor has NO insertions (count > 1)
    if (has_insertions_gt1(insertions) && !(key in ancestor_has_ins)) {
        found_unique = 1
        unique_insertions = insertions
        ins_count++
    }
    
    # Check if evolved has deletions (count > 1) AND ancestor has NO deletions (count > 1)
    if (has_deletions_gt1(deletions) && !(key in ancestor_has_del)) {
        found_unique = 1
        unique_deletions = deletions
        del_count++
    }
    
    # If we found unique mutations at this position, output the line
    if (found_unique) {
        mutation_type = ""
        if (snp_count > 0) mutation_type = mutation_type "SNP"
        if (ins_count > 0) mutation_type = mutation_type (mutation_type != "" ? "+" : "") "INS"
        if (del_count > 0) mutation_type = mutation_type (mutation_type != "" ? "+" : "") "DEL"
        
        new_info = "[" mutation_type ": "
        if (unique_snps != "") new_info = new_info unique_snps " "
        if (unique_insertions != "") new_info = new_info "+" unique_insertions " "
        if (unique_deletions != "") new_info = new_info "-" unique_deletions " "
        new_info = new_info "]"
        
        print $0 "\t" new_info >> "'$TEMP_FILE'"
        
        total_unique++
        if (snp_count > 0) total_snps++
        if (ins_count > 0) total_ins++
        if (del_count > 0) total_dels++
    }
}

END {
    print "" > "/dev/stderr"
    print "================================" > "/dev/stderr"
    print "ANALYSIS COMPLETE" > "/dev/stderr"
    print "================================" > "/dev/stderr"
    print "Filtered out (3+ mutation types): " filtered_complex > "/dev/stderr"
    print "Total unique mutation sites: " total_unique > "/dev/stderr"
    print "  Sites with new SNPs: " total_snps > "/dev/stderr"
    print "  Sites with new insertions: " total_ins > "/dev/stderr"
    print "  Sites with new deletions: " total_dels > "/dev/stderr"
    print "" > "/dev/stderr"
}
' "$ANCESTOR_FILE" "$EVOLVED_FILE"

# Sort and add header
echo "Sorting results..."
echo -e "CHROM\tPOS\tBASES_FOUND\tBASE_COUNTS\tINSERTIONS\tDELETIONS\tZYGOSITY\tNEW_MUTATIONS" > "$OUTPUT_UNIQUE"
sort -k1,1V -k2,2n "$TEMP_FILE" >> "$OUTPUT_UNIQUE"

# Create summary
{
    echo "================================"
    echo "UNIQUE MUTATIONS SUMMARY"
    echo "================================"
    echo "Evolved sample: $EVOLVED_FILE"
    echo "Ancestor sample: $ANCESTOR_FILE"
    echo ""
    echo "Rules applied:"
    echo "  - Only mutations with count > 1 considered"
    echo "  - Positions with 3+ mutation types filtered out"
    echo "  - All insertions treated as '+' (any length)"
    echo "  - All deletions treated as '-' (any length)"
    echo ""
    echo "Total unique mutation sites: $(tail -n +2 "$OUTPUT_UNIQUE" | wc -l)"
    echo ""
    echo "Breakdown by chromosome:"
    tail -n +2 "$OUTPUT_UNIQUE" | awk '{print $1}' | sort | uniq -c | awk '{printf "  %-5s: %5d sites\n", $2, $1}'
    echo ""
    echo "Mutation types:"
    tail -n +2 "$OUTPUT_UNIQUE" | grep -c 'SNP:' | awk '{print "  New SNPs only: " $1}'
    tail -n +2 "$OUTPUT_UNIQUE" | grep -c '\[INS:' | awk '{print "  New insertions only: " $1}'
    tail -n +2 "$OUTPUT_UNIQUE" | grep -c '\[DEL:' | awk '{print "  New deletions only: " $1}'
    tail -n +2 "$OUTPUT_UNIQUE" | grep -c 'SNP+INS\|INS+SNP' | awk '{print "  SNP + Insertion: " $1}'
    tail -n +2 "$OUTPUT_UNIQUE" | grep -c 'SNP+DEL\|DEL+SNP' | awk '{print "  SNP + Deletion: " $1}'
} > "$OUTPUT_SUMMARY"

rm "$TEMP_FILE"

echo ""
echo "Output files:"
echo "  Mutations: $OUTPUT_UNIQUE"
echo "  Summary: $OUTPUT_SUMMARY"
echo ""
cat "$OUTPUT_SUMMARY"
echo ""
echo "First 20 unique mutations:"
head -21 "$OUTPUT_UNIQUE"
