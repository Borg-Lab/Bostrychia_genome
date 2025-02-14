# Generate red algal phylogeny and perform ancestral state reconstruction of genome sizes

To perform the ancestral state reconstruction of genome size estimates from Kapraun &
Freshwater (2012), a concatenated gene tree was generated using conserved plastid, mitochondrial
and nuclear marker genes (atpA, atpB, cox1, psaA, psbA, rbcL and EF2) from red algae and selected
outgroup species. Gene sequences were retrieved from the NCBI Protein
database by searching with the keywords “Rhodophyta” and the respective gene name.



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
full_all_crypto <- Reduce(function(x, y) merge(x, y, by="Species", all=TRUE), data_list)



