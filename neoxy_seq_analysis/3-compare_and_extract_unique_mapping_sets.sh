#!/bin/bash
# Compare two mapping files to find truly unique mutations (SNPs + Indels)
# A mutation is unique only if:
# 1. The alternative base/indel has count > 1 in one sample
# 2. That alternative base/indel is completely ABSENT in the other sample

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <evolved_file.txt> <ancestor_file.txt> <output_prefix>"
    echo "Example: $0 generation_50.txt generation_1.txt unique_mutations"
    echo ""
    echo "Input format expected:"
    echo "CHROM  POS  BASES_FOUND  BASE_COUNTS  INSERTIONS  DELETIONS  ZYGOSITY"
    exit 1
fi

EVOLVED_FILE="$1"
ANCESTOR_FILE="$2"
OUTPUT_PREFIX="$3"

# Check if input files exist
if [ ! -f "$EVOLVED_FILE" ]; then
    echo "Error: File '$EVOLVED_FILE' not found!"
    exit 1
fi

if [ ! -f "$ANCESTOR_FILE" ]; then
    echo "Error: File '$ANCESTOR_FILE' not found!"
    exit 1
fi

echo "Finding truly unique mutations (SNPs + Indels)..."
echo "================================"
echo "Evolved sample: $EVOLVED_FILE"
echo "Ancestor sample: $ANCESTOR_FILE"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""

# Create output files
OUTPUT_UNIQUE="${OUTPUT_PREFIX}_truly_unique_mutations.txt"
OUTPUT_SUMMARY="${OUTPUT_PREFIX}_summary.txt"
TEMP_FILE=$(mktemp)

# Run the comparison using awk
awk '
BEGIN {
    FS = "\t"
    OFS = "\t"
    print "Processing files..." > "/dev/stderr"
}

# Function to parse insertions and extract individual insertions with counts
function parse_indels(indel_string, indels_array, counts_array) {
    delete indels_array
    delete counts_array
    
    if (indel_string == "" || indel_string == "-" || indel_string ~ /^-IN/ || indel_string ~ /^-NN/) {
        return 0
    }
    
    # Handle multiple indels: +1T(2),+2NN(1) or -1N(5),-2NN(3)
    n = split(indel_string, parts, ",")
    for (i = 1; i <= n; i++) {
        # Extract indel and count: +1T(2) -> indel=+1T, count=2
        if (match(parts[i], /([+-][0-9]+[ACGTN]+)\(([0-9]+)\)/, arr)) {
            indel = arr[1]
            count = arr[2]
            indels_array[indel] = 1
            counts_array[indel] = count
        }
    }
    return length(indels_array)
}

# Function to extract bases and their counts
function parse_bases(bases_found, base_counts, bases_array, counts_array) {
    delete bases_array
    delete counts_array
    
    # Split BASES_FOUND: e.g., "A,G,T" or "A,C"
    n_bases = split(bases_found, bases_list, ",")
    
    # Extract counts from BASE_COUNTS: e.g., "A(17),G(1)" -> A=17, G=1
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
    if ($1 == "CHROM") next  # Skip header
    
    chrom = $1
    pos = $2
    bases_found = $3
    base_counts = $4
    insertions = $5
    deletions = $6
    
    key = chrom":"pos
    
    # Parse bases present in ancestor
    parse_bases(bases_found, base_counts, ancestor_bases, ancestor_counts)
    for (base in ancestor_bases) {
        ancestor_has[key":"base] = 1
        ancestor_count[key":"base] = ancestor_counts[base]
    }
    
    # Parse insertions in ancestor
    parse_indels(insertions, ancestor_ins, ancestor_ins_counts)
    for (indel in ancestor_ins) {
        ancestor_has_ins[key":"indel] = 1
        ancestor_ins_count[key":"indel] = ancestor_ins_counts[indel]
    }
    
    # Parse deletions in ancestor
    parse_indels(deletions, ancestor_del, ancestor_del_counts)
    for (indel in ancestor_del) {
        ancestor_has_del[key":"indel] = 1
        ancestor_del_count[key":"indel] = ancestor_del_counts[indel]
    }
    
    ancestor_data[key] = $0
    next
}

# Read evolved file
{
    if ($1 == "CHROM") next  # Skip header
    
    chrom = $1
    pos = $2
    bases_found = $3
    base_counts = $4
    insertions = $5
    deletions = $6
    zygosity = $7
    
    key = chrom":"pos
    
    # Parse bases in evolved sample
    parse_bases(bases_found, base_counts, evolved_bases, evolved_counts)
    
    # Parse indels in evolved sample
    parse_indels(insertions, evolved_ins, evolved_ins_counts)
    parse_indels(deletions, evolved_del, evolved_del_counts)
    
    # Check for truly new mutations
    found_unique = 0
    unique_snps = ""
    unique_insertions = ""
    unique_deletions = ""
    
    snp_count = 0
    ins_count = 0
    del_count = 0
    
    # Check each base in evolved sample
    for (base in evolved_bases) {
        count = evolved_counts[base]
        
        # Only consider if count > 1 (avoid sequencing errors)
        if (count > 1) {
            # Check if this base is ABSENT in ancestor
            if (!(key":"base in ancestor_has)) {
                found_unique = 1
                if (unique_snps != "") unique_snps = unique_snps ","
                unique_snps = unique_snps base "(" count ")"
                snp_count++
            }
        }
    }
    
    # Check each insertion in evolved sample
    for (indel in evolved_ins) {
        count = evolved_ins_counts[indel]
        
        # Only consider if count > 1
        if (count > 1) {
            # Check if this insertion is ABSENT in ancestor
            if (!(key":"indel in ancestor_has_ins)) {
                found_unique = 1
                if (unique_insertions != "") unique_insertions = unique_insertions ","
                unique_insertions = unique_insertions indel "(" count ")"
                ins_count++
            }
        }
    }
    
    # Check each deletion in evolved sample
    for (indel in evolved_del) {
        count = evolved_del_counts[indel]
        
        # Only consider if count > 1
        if (count > 1) {
            # Check if this deletion is ABSENT in ancestor
            if (!(key":"indel in ancestor_has_del)) {
                found_unique = 1
                if (unique_deletions != "") unique_deletions = unique_deletions ","
                unique_deletions = unique_deletions indel "(" count ")"
                del_count++
            }
        }
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
echo "Creating summary..."
{
    echo "================================"
    echo "UNIQUE MUTATIONS SUMMARY"
    echo "================================"
    echo "Evolved sample: $EVOLVED_FILE"
    echo "Ancestor sample: $ANCESTOR_FILE"
    echo ""
    echo "Total unique mutation sites: $(tail -n +2 "$OUTPUT_UNIQUE" | wc -l)"
    echo ""
    echo "Breakdown by chromosome:"
    tail -n +2 "$OUTPUT_UNIQUE" | awk '{print $1}' | sort | uniq -c | awk '{printf "  %-5s: %5d sites\n", $2, $1}'
    echo ""
    echo "Mutation types:"
    tail -n +2 "$OUTPUT_UNIQUE" | grep -o '\[SNP[^]]*\]' | wc -l | awk '{print "  New SNPs: " $1}'
    tail -n +2 "$OUTPUT_UNIQUE" | grep -o 'INS' | wc -l | awk '{print "  New insertions: " $1}'
    tail -n +2 "$OUTPUT_UNIQUE" | grep -o 'DEL' | wc -l | awk '{print "  New deletions: " $1}'
} > "$OUTPUT_SUMMARY"

# Cleanup
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
