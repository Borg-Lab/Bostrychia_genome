from collections import defaultdict

def parse_fasta(fasta_file):
    gene_ids = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                gene_id = line[1:].strip()
                gene_ids.append(gene_id)
    return gene_ids

def parse_blast_results(blast_file):
    best_hits = {}
    with open(blast_file) as f:
        for line in f:
            query, subject, *rest = line.strip().split("\t")
            if query not in best_hits:
                best_hits[query] = subject
    return best_hits

def find_reciprocal_best_hits(blast1, blast2):
    rbh = {}
    for gene1, gene2 in blast1.items():
        if gene2 in blast2 and blast2[gene2] == gene1:
            rbh[gene1] = gene2
        else:
            rbh[gene1] = '.'
    return rbh

# Parse the FASTA files to get all male and female gene IDs
male_gene_ids = parse_fasta("Bm3235_M_v3_annoV5_iso_TErem.codingseq")
female_gene_ids = parse_fasta("Bm3235_F_assembly_annoV5_iso_TErem_v2.codingseq")

# Parse BLAST results
blast1_results = parse_blast_results("Male_Female_blast.out")
blast2_results = parse_blast_results("Female_Male_blast.out")

# Find reciprocal best hits
reciprocal_best_hits_male_female = find_reciprocal_best_hits(blast1_results, blast2_results)
reciprocal_best_hits_female_male = find_reciprocal_best_hits(blast2_results, blast1_results)

# Write output for male-to-female
with open("reciprocal_best_hits_male_female.tsv", "w") as f:
    for gene1 in male_gene_ids:
        gene2 = reciprocal_best_hits_male_female.get(gene1, '.')
        f.write(f"{gene1}\t{gene2}\n")

# Write output for female-to-male
with open("reciprocal_best_hits_female_male.tsv", "w") as f:
    for gene1 in female_gene_ids:
        gene2 = reciprocal_best_hits_female_male.get(gene1, '.')
        f.write(f"{gene1}\t{gene2}\n")
