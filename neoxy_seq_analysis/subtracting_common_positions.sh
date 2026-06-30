#!/bin/bash
# Script to subtract common positions between two variant files
# Keeps only positions unique to the first file

# Usage check
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <1st_gen.txt> <50th_gen.txt> <output_unique.txt>"
    echo "Example: $0 xy_1st_all_variants.txt xy_50th_all_variants.txt xy_1st_unique.txt"
    exit 1
fi

FIRST_GEN="$1"
FIFTIETH_GEN="$2"
OUTPUT="$3"

# Check if input files exist
if [ ! -f "$FIRST_GEN" ]; then
    echo "Error: First generation file '$FIRST_GEN' not found!"
    exit 1
fi

if [ ! -f "$FIFTIETH_GEN" ]; then
    echo "Error: 50th generation file '$FIFTIETH_GEN' not found!"
    exit 1
fi

echo "Processing variant files..."
echo "================================"
echo "1st generation file: $FIRST_GEN"
echo "50th generation file: $FIFTIETH_GEN"
echo "Output file: $OUTPUT"
echo ""

# Create temporary files
TEMP_1ST=$(mktemp)
TEMP_50TH=$(mktemp)

# Extract and normalize positions from 1st gen (CHROM:POS or CHROM POS or CHROM\tPOS)
echo "Extracting positions from 1st generation..."
awk '{
    # Handle different separators (tab, space, colon)
    if ($0 ~ /:/) {
        split($0, a, /[:\t ]/)
        print a[1]":"a[2]
    } else {
        print $1":"$2
    }
}' "$FIRST_GEN" | sort -u > "$TEMP_1ST"

FIRST_COUNT=$(wc -l < "$TEMP_1ST")
echo "  Total unique positions in 1st gen: $FIRST_COUNT"

# Extract and normalize positions from 50th gen
echo "Extracting positions from 50th generation..."
awk '{
    # Handle different separators (tab, space, colon)
    if ($0 ~ /:/) {
        split($0, a, /[:\t ]/)
        print a[1]":"a[2]
    } else {
        print $1":"$2
    }
}' "$FIFTIETH_GEN" | sort -u > "$TEMP_50TH"

FIFTIETH_COUNT=$(wc -l < "$TEMP_50TH")
echo "  Total unique positions in 50th gen: $FIFTIETH_COUNT"

# Find positions unique to 1st generation (subtract 50th from 1st)
echo ""
echo "Finding positions unique to 1st generation..."
comm -23 "$TEMP_1ST" "$TEMP_50TH" > "${OUTPUT}.positions_only"

UNIQUE_COUNT=$(wc -l < "${OUTPUT}.positions_only")
COMMON_COUNT=$((FIRST_COUNT - UNIQUE_COUNT))

echo "  Common positions (in both): $COMMON_COUNT"
echo "  Unique to 1st gen: $UNIQUE_COUNT"

# Now extract the full variant information (CHROM POS REF ALT) for unique positions
echo ""
echo "Extracting full variant information..."

# Create associative array lookup for unique positions
awk 'NR==FNR {
    unique[$1]=1
    next
}
{
    # Extract CHROM and POS
    if ($0 ~ /:/) {
        split($0, a, /[:\t ]/)
        chrom = a[1]
        pos = a[2]
        ref = a[3]
        alt = a[4]
    } else {
        chrom = $1
        pos = $2
        ref = $3
        alt = $4
    }
    
    key = chrom":"pos
    
    # If this position is unique to 1st gen, print full info
    if (key in unique) {
        print chrom"\t"pos"\t"ref"\t"alt
    }
}' "${OUTPUT}.positions_only" "$FIRST_GEN" > "$OUTPUT"

# Clean up temp files
rm "$TEMP_1ST" "$TEMP_50TH" "${OUTPUT}.positions_only"

echo ""
echo "================================"
echo "COMPLETE!"
echo "================================"
echo "Output saved to: $OUTPUT"
echo "Format: CHROM  POS  REF  ALT (tab-separated)"
echo ""
echo "Summary:"
echo "  1st generation total: $FIRST_COUNT positions"
echo "  50th generation total: $FIFTIETH_COUNT positions"
echo "  Common (removed): $COMMON_COUNT positions"
echo "  Unique to 1st (kept): $UNIQUE_COUNT positions"
