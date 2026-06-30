#!/usr/bin/env python3

"""
Enhanced variant analysis script with indel tracking
Analyzes BAM files at variant positions and extracts:
- Unique bases (A, T, G, C)
- Insertions with counts
- Deletions with counts
"""

import sys
import subprocess
import re
from collections import defaultdict, Counter

def parse_mpileup_bases(bases_string):
    """
    Parse mpileup bases string to extract:
    - Pure bases (A, T, G, C)
    - Insertions (with sequences)
    - Deletions (with sequences)
    
    Returns: dict with 'bases', 'insertions', 'deletions'
    """
    bases_upper = bases_string.upper()
    
    pure_bases = []
    insertions = []
    deletions = []
    
    i = 0
    length = len(bases_upper)
    
    while i < length:
        char = bases_upper[i]
        
        if char in 'ATCG':
            # Direct base call
            pure_bases.append(char)
            i += 1
            
        elif char in '.,':
            # Reference match - skip for now (or could add ref base)
            i += 1
            
        elif char == '^':
            # Start of read segment, skip mapping quality char
            i += 2
            
        elif char == '$':
            # End of read segment
            i += 1
            
        elif char == '+':
            # Insertion
            i += 1
            # Parse the length
            length_str = ''
            while i < length and bases_upper[i].isdigit():
                length_str += bases_upper[i]
                i += 1
            
            ins_len = int(length_str)
            # Extract insertion sequence
            ins_seq = bases_upper[i:i+ins_len]
            insertions.append(f"+{ins_len}{ins_seq}")
            i += ins_len
            
        elif char == '-':
            # Deletion
            i += 1
            # Parse the length
            length_str = ''
            while i < length and bases_upper[i].isdigit():
                length_str += bases_upper[i]
                i += 1
            
            del_len = int(length_str)
            # Extract deletion sequence
            del_seq = bases_upper[i:i+del_len]
            deletions.append(f"-{del_len}{del_seq}")
            i += del_len
            
        elif char in '*N#':
            # Deletion placeholder, unknown, or padding
            i += 1
            
        else:
            # Skip any other character
            i += 1
    
    return {
        'bases': pure_bases,
        'insertions': insertions,
        'deletions': deletions
    }

def get_mpileup_at_position(bam_file, chrom, pos):
    """
    Run samtools mpileup at a specific position
    Returns the bases string or None
    """
    try:
        cmd = [
            'samtools', 'mpileup',
            '-aa',  # Output all positions
            '-q', '0',  # Min mapping quality
            '-Q', '0',  # Min base quality
            '-r', f"{chrom}:{pos}-{pos}",
            bam_file
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.returncode != 0 or not result.stdout.strip():
            return None
        
        # Parse mpileup output (tab-separated)
        # Format: chrom pos ref depth bases qualities
        fields = result.stdout.strip().split('\t')
        
        if len(fields) >= 5:
            return fields[4]  # bases column
        
        return None
        
    except Exception as e:
        print(f"Error running mpileup: {e}", file=sys.stderr)
        return None

def format_indel_counts(indel_list):
    """
    Convert list of indels to count format
    Example: ['+3ATG', '+3ATG', '+2GC'] -> '+3ATG(2),+2GC(1)'
    """
    if not indel_list:
        return ''
    
    counts = Counter(indel_list)
    return ','.join([f"{indel}({count})" for indel, count in counts.most_common()])

def determine_zygosity(base_counts, insertions, deletions, ref_base, alt_base, hetero_min=0.20, hetero_max=0.80):
    """
    Determine zygosity based on allele frequencies
    
    Parameters:
    - base_counts: Counter object with base frequencies
    - insertions: list of insertions
    - deletions: list of deletions
    - ref_base: reference allele
    - alt_base: alternative allele from VCF
    - hetero_min: minimum frequency for heterozygous allele (default 0.20 = 20%)
    - hetero_max: maximum frequency for heterozygous allele (default 0.80 = 80%)
    
    Returns: string ('HOMOZYGOUS_REF', 'HOMOZYGOUS_ALT', 'HETEROZYGOUS', 'COMPLEX', 'NO_COVERAGE')
    """
    total_bases = sum(base_counts.values())
    
    if total_bases == 0:
        return 'NO_COVERAGE'
    
    # Get the most common allele and its frequency
    if base_counts:
        most_common_base, most_common_count = base_counts.most_common(1)[0]
        most_common_freq = most_common_count / total_bases
    else:
        return 'NO_COVERAGE'
    
    # Count number of different alleles present
    num_alleles = len(base_counts)
    has_indels = len(insertions) > 0 or len(deletions) > 0
    
    # HOMOZYGOUS: Only one allele at 100% frequency, no indels
    if num_alleles == 1 and not has_indels and most_common_freq == 1.0:
        if most_common_base.upper() == ref_base.upper():
            return 'HOMOZYGOUS_REF'
        elif most_common_base.upper() == alt_base.upper():
            return 'HOMOZYGOUS_ALT'
        else:
            return 'HOMOZYGOUS_OTHER'
    
    # HETEROZYGOUS: Two alleles (ideally REF and ALT)
    if num_alleles == 2 and not has_indels:
        # Check if ref and alt are the two alleles
        bases_present = set(b.upper() for b in base_counts.keys())
        
        # Get counts for each allele (case-insensitive)
        allele_counts = {}
        for base, count in base_counts.items():
            base_upper = base.upper()
            allele_counts[base_upper] = allele_counts.get(base_upper, 0) + count
        
        if ref_base.upper() in allele_counts and alt_base.upper() in allele_counts:
            # Classic heterozygous (ref/alt)
            ref_count = allele_counts[ref_base.upper()]
            alt_count = allele_counts[alt_base.upper()]
            ref_freq = ref_count / total_bases
            alt_freq = alt_count / total_bases
            
            # Both alleles present at reasonable frequencies
            if hetero_min <= ref_freq <= hetero_max and hetero_min <= alt_freq <= hetero_max:
                return 'HETEROZYGOUS'
            # One allele dominates but not 100%
            elif ref_freq > hetero_max:
                return 'MOSTLY_REF'
            elif alt_freq > hetero_max:
                return 'MOSTLY_ALT'
        
        # Two alleles but not the expected ref/alt
        return 'HETEROZYGOUS_OTHER'
    
    # COMPLEX: Multiple alleles (>2) or presence of indels or mixed frequencies
    if num_alleles > 2 or has_indels or (num_alleles > 1 and most_common_freq < 1.0):
        return 'COMPLEX'
    
    # Edge case: single allele but with indels
    if num_alleles == 1 and has_indels:
        return 'COMPLEX'
    
    # Default for edge cases
    return 'COMPLEX'

def main():
    if len(sys.argv) != 4:
        print("Usage: python3 analyze_variants_with_indels.py <variants_txt> <bam_file> <output_txt>")
        print("Example: python3 analyze_variants_with_indels.py variants.txt sample.bam results.txt")
        sys.exit(1)
    
    variants_file = sys.argv[1]
    bam_file = sys.argv[2]
    output_file = sys.argv[3]
    
    print("Variant Analysis with Indel Tracking")
    print("=" * 50)
    print(f"Variants file: {variants_file}")
    print(f"BAM file: {bam_file}")
    print(f"Output file: {output_file}")
    print("=" * 50)
    
    processed = 0
    
    with open(variants_file, 'r') as vf, open(output_file, 'w') as of:
        # Write header
        of.write("CHROM\tPOS\tREF\tALT\tBASES_FOUND\tBASE_COUNTS\tINSERTIONS\tDELETIONS\tZYGOSITY\n")
        
        # Skip header in input file
        next(vf)
        
        for line in vf:
            line = line.strip()
            if not line:
                continue
            
            # Parse the variant line (format: CHROM:POS:REF:ALT)
            parts = line.split(':')
            if len(parts) < 4:
                print(f"Warning: Skipping malformed line: {line}")
                continue
            
            chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
            
            # Get mpileup at this position
            bases_string = get_mpileup_at_position(bam_file, chrom, pos)
            
            if bases_string is None:
                print(f"Warning: No coverage at {chrom}:{pos}")
                continue
            
            # Parse the bases
            parsed = parse_mpileup_bases(bases_string)
            
            # Get unique bases and counts
            base_counts = Counter(parsed['bases'])
            unique_bases = sorted(base_counts.keys())
            bases_str = ','.join(unique_bases)
            
            # Format base counts
            base_counts_str = ','.join([f"{base}({count})" for base, count in base_counts.most_common()])
            
            # Format indels with counts
            insertions_str = format_indel_counts(parsed['insertions'])
            deletions_str = format_indel_counts(parsed['deletions'])
            
            # Determine zygosity
            zygosity = determine_zygosity(
                base_counts, 
                parsed['insertions'], 
                parsed['deletions'],
                ref,
                alt
            )
            
            # Write output
            of.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{bases_str}\t{base_counts_str}\t{insertions_str}\t{deletions_str}\t{zygosity}\n")
            
            processed += 1
            
            # Display progress
            print(f"[{processed}] {chrom}:{pos} - {zygosity}")
            if bases_str:
                print(f"    Bases: {{{base_counts_str}}}")
            if insertions_str:
                print(f"    Insertions: {insertions_str}")
            if deletions_str:
                print(f"    Deletions: {deletions_str}")
    
    print()
    print("=" * 50)
    print(f"Analysis complete!")
    print(f"Positions processed: {processed}")
    print(f"Results saved to: {output_file}")

if __name__ == '__main__':
    main()
