from Bio import SeqIO
import sys
import os

# Get input directory from command-line arguments
input_dir = sys.argv[1]

# List all files in the input directory
all_files = os.listdir(input_dir)

# Extract unique species names by removing file extensions
all_species = set()
for filename in all_files:
    species_name = filename.split('.')[0]  # Assuming format: species.fasta
    all_species.add(species_name)

# Dictionary to store alignment records per species
alignment_records = {species: {} for species in all_species}

# Read each FASTA file and store sequences by gene
for filename in all_files:
    species_name = filename.split('.')[0]
    file_path = os.path.join(input_dir, filename)
    
    # Parse the FASTA file
    for record in SeqIO.parse(file_path, "fasta"):
        gene_name = record.id.split('_')[0]  # Extract gene identifier
        
        # Store sequence record under the corresponding species and gene
        if gene_name not in alignment_records[species_name]:
            alignment_records[species_name][gene_name] = []
        alignment_records[species_name][gene_name].append(record)

# Print collected data summary
for species, genes in alignment_records.items():
    print(f"Species: {species}, Total genes: {len(genes)}")
