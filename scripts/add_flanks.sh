#!/bin/bash

# Define the input files
BED_FILE="$1"
CONTIG_SIZES_FILE="$2"
OUTPUT_FILE="elements_with_flanks_NEW.bed"

# Create an associative array to store contig sizes
declare -A contig_sizes

# Read the contig sizes into the array
while read -r contig size; do
    contig_sizes[$contig]=$size
done < "$CONTIG_SIZES_FILE"

# Process the BED file and add flanking regions
awk -v OFS="\t" -v flanksize=15000 -v contig_sizes_file="$CONTIG_SIZES_FILE" '
BEGIN {
    # Read contig sizes into an array
    while ((getline < contig_sizes_file) > 0) {
        contig_sizes[$1] = $2
    }
}
{
    contig = $1
    start = $2
    end = $3
    name = $4

    # Calculate combined flank range
    left_start = (start - flanksize < 1) ? 1 : start - flanksize
    right_end = (end + flanksize > contig_sizes[contig]) ? contig_sizes[contig] : end + flanksize

    # Print combined entry
    print contig, left_start, right_end, name "_combined_flank"
}' "$BED_FILE" > "$OUTPUT_FILE"

echo "Output written to $OUTPUT_FILE"
