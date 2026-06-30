#!/usr/bin/bash

# Script to extract CHROM, POS, REF, ALT from VCF.gz and save to TXT
# Usage: ./vcf_to_txt.sh <input_vcf.gz> <output_txt>

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_vcf.gz> <output_txt>"
    echo "Example: $0 variants.vcf.gz variants.txt"
    exit 1
fi

VCF_FILE="$1"
TXT_FILE="$2"

# Write header to TXT
echo "CHROM:POS:REF:ALT" > "$TXT_FILE"

# Extract variants from VCF.gz (skip header lines starting with #)
zcat "$VCF_FILE" | grep -v "^#" | cut -f1,2,4,5 | tr '\t' ':' >> "$TXT_FILE"

echo "Done! Extracted variants saved to: $TXT_FILE"
echo "Total variants: $(( $(wc -l < "$TXT_FILE") - 0 ))"
