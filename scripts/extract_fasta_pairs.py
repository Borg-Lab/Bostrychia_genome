from Bio import SeqIO

def load_gene_pairs(gene_pairs_file):
    """Load gene pairs from the txt file"""
    pairs = []
    with open(gene_pairs_file, 'r') as f:
        for line in f:
            male_id, female_id = line.strip().split()
            pairs.append((male_id, female_id))
    return pairs

def extract_sequences(fasta_file, gene_ids):
    """Extract sequences for the given gene IDs from a FASTA file"""
    sequences = {gene_id: None for gene_id in gene_ids}

    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in sequences:
            sequences[record.id] = record
            if all(sequences.values()):
                break  # Stop when we have found all sequences

    return sequences

def write_fasta(output_file, sequences):
    """Write the given sequences to a FASTA file"""
    with open(output_file, 'w') as f:
        SeqIO.write(sequences, f, "fasta")

def process_gene_pairs(gene_pairs_file, prot_fasta, cds_fasta, output_dir):
    """Process each gene pair and write the corresponding FASTA files"""
    gene_pairs = load_gene_pairs(gene_pairs_file)

    for male_id, female_id in gene_pairs:
        # Extract sequences from protein FASTA
        prot_seqs = extract_sequences(prot_fasta, [male_id, female_id])
        # Write protein FASTA for the pair
        output_prot_file = f"{output_dir}/{male_id}_{female_id}_prot.fasta"
        write_fasta(output_prot_file, prot_seqs.values())

        # Extract sequences from CDS FASTA
        cds_seqs = extract_sequences(cds_fasta, [male_id, female_id])
        # Write CDS FASTA for the pair
        output_cds_file = f"{output_dir}/{male_id}_{female_id}_cds.fasta"
        write_fasta(output_cds_file, cds_seqs.values())

# Usage
gene_pairs_file = 'Matching.tsv'  # File with gene pairs (male and female IDs)
prot_fasta = 'MF_prot.fa'       # Protein FASTA file
cds_fasta = 'MF_cds.fa'             # CDS FASTA file
output_dir = 'output_fastas'        # Directory to save the resulting FASTA files

# Make sure output directory exists
import os
os.makedirs(output_dir, exist_ok=True)

process_gene_pairs(gene_pairs_file, prot_fasta, cds_fasta, output_dir)
