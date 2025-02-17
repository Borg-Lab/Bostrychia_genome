# Read reciprocal best hits
def read_reciprocal_hits(file):
    hits = {}
    with open(file) as f:
        for line in f:
            gene1, gene2 = line.strip().split("\t")
            hits[gene1] = gene2
    return hits

# Load the best hit files
male_to_female_hits = read_reciprocal_hits("reciprocal_best_hits_male_female.tsv")
female_to_male_hits = read_reciprocal_hits("reciprocal_best_hits_female_male.tsv")

# Summarize results
male_specific = set()
female_specific = set()
matching_genes = set()

for male_gene, female_gene in male_to_female_hits.items():
    if female_gene == ".":
        male_specific.add(male_gene)
    elif female_gene in female_to_male_hits and female_to_male_hits[female_gene] == male_gene:
        matching_genes.add((male_gene, female_gene))

for female_gene, male_gene in female_to_male_hits.items():
    if male_gene == ".":
        female_specific.add(female_gene)

# Write summary to file
with open("reciprocal_best_hits_summary.tsv", "w") as f:
    f.write("Gene Type\tGene 1\tGene 2\n")
    for male_gene, female_gene in matching_genes:
        f.write(f"Matching\t{male_gene}\t{female_gene}\n")
    for male_gene in male_specific:
        f.write(f"Male Specific\t{male_gene}\t.\n")
    for female_gene in female_specific:
        f.write(f"Female Specific\t.\t{female_gene}\n")
