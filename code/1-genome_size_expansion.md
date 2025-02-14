# Generate red algal phylogeny and perform ancestral state reconstruction of genome sizes

To perform the ancestral state reconstruction of genome size estimates from Kapraun &
Freshwater (2012), a concatenated gene tree was generated using conserved plastid, mitochondrial
and nuclear marker genes (atpA, atpB, cox1, psaA, psbA, rbcL and EF2) from red algae and selected
outgroup species. Gene sequences were retrieved from the NCBI Protein
database by searching with the keywords “Rhodophyta” and the respective gene name.

## 1. Generating concatenated gene tree

### 1.1 Prepare input data
First the genes were combined by species and redundant accessions were removed. 
```r
# Load packages
library(dplyr)

# Define gene names
genes <- c("atpA", "atpB", "cox1", "rbcL", "psaA", "psbA", "EF2")

# Read and preprocess all gene tables
data_list <- lapply(genes, function(gene) {
  df <- read.table(paste0("table_sort_", gene, ".txt"), sep=";", header=TRUE, stringsAsFactors=FALSE)[,1:2]
  colnames(df) <- c("Species", gene)  
  df %>% distinct(Species, .keep_all = TRUE)
})

# Merge data frames
Combined_genes <- Reduce(function(x, y) merge(x, y, by="Species", all=TRUE), data_list)

write.csv(Combined_genes, file="Merged_datasets_genes.csv")
write.table(Combined_genes,file="Merged_datasets_genes.txt")
```
Afterwards, the entries were filtered manually to keep just the species with accessions for at least three genes.

Using `seqtk` to get fasta files for each gene. 
```bash
# For each gene:
seqtk subseq "gene".fasta seq_"gene".txt > seq_"gene".fa 
```
Generate individual gene alignments.
```bash
# For each gene:
 mafft --auto seq_"gene"_total.fa > aln_"gene"_total.fa
```
Concatenate alignments using a custom script.
```bash
 ./concatenate_gene_alignments.py aln_EF2_total.fa aln_atpA_total.fa aln_atpB_total.fa aln_cox1_total.fa aln_psaA_total.fa aln_psbA_total.fa aln_rbcL_total.fa 
```



