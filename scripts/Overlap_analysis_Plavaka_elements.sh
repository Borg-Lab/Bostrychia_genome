#!/bin/sh

FASTA="$1"
ADOMET="total_AdoMetSet_sequences.fa"
OUTFILE="$2"
PREFIX="$3"

CWD=$(pwd)

# Check if sufficient arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fasta_file> <outfile> <prefix>"
    exit 1
fi

# Copy fasta file and repeats out file in respective directory
grep ">" "$ADOMET" > total_AdometSET_sequences.txt
grep "$PREFIX" total_AdometSET_sequences.txt > "${PREFIX}_AdometSET_seqs.txt"
sed 's/^.\{6\}//' "${PREFIX}_AdometSET_seqs.txt" > "${PREFIX}_AdometSET_seqs_trim.txt"

# Edit out file
awk 'NR > 3 {print $5, $6, $7, $10, $11}' "$OUTFILE" > "${OUTFILE}_edit"
awk -v OFS="\t" '$1=$1' "${OUTFILE}_edit" > "${OUTFILE}_edit.txt"
awk -F'\t' -v prefix="$PREFIX" 'BEGIN {OFS="\t"} { $4 = prefix "_" $4 "#" $5; $5 = ""; print $1, $2, $3, $4 }' "${OUTFILE}_edit.txt" > "${OUTFILE}_edit_merged.txt"

tr -d '\r' < "${PREFIX}_AdometSET_seqs_trim.txt" > "${PREFIX}_AdometSET_seqs_trim_clean.txt"
grep -w -f "${PREFIX}_AdometSET_seqs_trim_clean.txt" "${OUTFILE}_edit_merged.txt" > "${PREFIX}_AdometSET.bed"

# Prepare contig length file
awk '/^>/ {if (seqlen) print name, seqlen; name=substr($0, 2); seqlen=0; next} {seqlen += length($0)} END {print name, seqlen}' "$FASTA" > "${PREFIX}_contig_lengths.txt"

# Add flanks on both sides (30kb total)
scripts/add_flanks.sh "$CWD/${PREFIX}_AdometSET.bed" "$CWD/${PREFIX}_contig_lengths.txt"
#output: elements_with_flanks_NEW.bed
awk '$3 >= $2' elements_with_flanks_NEW.bed > "${PREFIX}_flanks_NEW_edit.bed"

# Generating fasta based on bed file
bedtools getfasta -fi "$FASTA" -bed "${PREFIX}_flanks_NEW_edit.bed" -name > "${PREFIX}_flanks_NEW_edit.fasta"
blastx -query "${PREFIX}_flanks_NEW_edit.fasta" -db ../database/total_prots_redAlgae_CDHIT.fa -evalue 0.001 -outfmt "6 qseqid qstart qend sseqid sstart send pident slen qlen length mismatch gapopen evalue" > out_NEW
awk '{print $1, $4}' out_NEW > out_NEW_filtered
sort out_NEW_filtered | uniq > output_NEW.txt
sed -i 's/ /\t/' output_NEW.txt
awk -F'\t' 'BEGIN {OFS=FS} {split($2, a, "_"); $2=a[1]; $3=a[2]; print $1, $2}' output_NEW.txt > output_NEW.tsv
sort output_NEW.tsv | uniq > output_final_NEW.txt

# Process output file
python scripts/process_ids.py > analysis_out
