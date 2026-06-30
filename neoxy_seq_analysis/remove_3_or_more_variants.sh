#!/bin/bash
# Filter positions with 3+ mutation types
# Count ALL insertions as 1 type and ALL deletions as 1 type (regardless of count)

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.txt> <output.txt>"
    exit 1
fi

INPUT="$1"
OUTPUT="$2"

awk '
BEGIN {
    FS = "\t"
    OFS = "\t"
    total = 0
    kept = 0
    filtered = 0
}

NR == 1 {
    print
    next
}

{
    total++
    
    bases_found = $3
    insertions = $5
    deletions = $6
    
    # Count mutation types
    mutation_types = 0
    
    # 1. Count number of DIFFERENT bases
    n_bases = split(bases_found, bases_array, ",")
    mutation_types += n_bases
    
    # 2. Check if there are ANY insertions (regardless of count)
    if (insertions != "" && insertions != "-" && insertions !~ /^-IN/ && insertions !~ /^-NN/ && insertions !~ /COMPLEX/) {
        mutation_types += 1
    }
    
    # 3. Check if there are ANY deletions (regardless of count)
    if (deletions != "" && deletions != "-" && deletions !~ /^-IN/ && deletions !~ /^-NN/ && deletions !~ /COMPLEX/) {
        mutation_types += 1
    }
    
    # Keep if ≤2 mutation types
    if (mutation_types <= 2) {
        print
        kept++
    } else {
        filtered++
    }
}

END {
    print "" > "/dev/stderr"
    print "================================" > "/dev/stderr"
    print "FILTERING SUMMARY" > "/dev/stderr"
    print "================================" > "/dev/stderr"
    print "Total variants: " total > "/dev/stderr"
    print "Kept (≤2 mutation types): " kept > "/dev/stderr"
    print "Filtered (3+ mutation types): " filtered > "/dev/stderr"
    print "Percentage kept: " sprintf("%.2f%%", (kept/total)*100) > "/dev/stderr"
}
' "$INPUT" > "$OUTPUT"

echo ""
echo "Output saved to: $OUTPUT"
echo ""
echo "First 20 filtered variants:"
head -21 "$OUTPUT"
