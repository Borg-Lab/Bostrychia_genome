
#!/bin/bash

# Define directories
PROT_DIR="output_fastas"   # Directory where protein FASTA files are stored
CDS_DIR="output_fastas"    # Directory where CDS FASTA files are stored
OUTPUT_DIR="aligned_outputs"  # Directory for the output files

# Make sure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through all protein FASTA files
for prot_fasta in "$PROT_DIR"/*_prot.fasta; do
    # Extract the base name of the file (without the directory and extension)
    base_name=$(basename "$prot_fasta" _prot.fasta)

    # Corresponding CDS FASTA file
    cds_fasta="$CDS_DIR/${base_name}_cds.fasta"

    # Check if corresponding CDS file exists
    if [ ! -f "$cds_fasta" ]; then
        echo "CDS file for $base_name not found!"
        continue
    fi

    # MAFFT alignment (protein sequences)
    aln_file="$OUTPUT_DIR/${base_name}_aln.fasta"
    echo "Running MAFFT on $prot_fasta..."
    mafft --thread 50 --auto "$prot_fasta" > "$aln_file"

    # PAL2NAL conversion (using aligned protein and CDS sequences)
    pal2nal_output="$OUTPUT_DIR/${base_name}.pal2nal"
    echo "Running PAL2NAL on $aln_file and $cds_fasta..."
    pal2nal.pl "$aln_file" "$cds_fasta" -output paml -nogap > "$pal2nal_output"

    echo "Finished processing $base_name"
done
