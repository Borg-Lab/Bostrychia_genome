# Generate red algal phylogeny and perform ancestral state reconstruction of genome sizes

To perform the ancestral state reconstruction of genome size estimates from Kapraun &
Freshwater (2012), a concatenated gene tree was generated using conserved plastid, mitochondrial
and nuclear marker genes (atpA, atpB, cox1, psaA, psbA, rbcL and EF2) from red algae and selected
outgroup species. Gene sequences were retrieved from the NCBI Protein
database by searching with the keywords “Rhodophyta” and the respective gene name.

## 1. Generating concatenated gene tree

### 1.1 Prepare input alignment
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
Concatenate alignments using the custom script
[concatenate_gene_alignments.py](https://github.com/Borg-Lab/Bostrychia_genome/edit/main/scripts/concatenate_gene_alignments.py)
```bash
 ./concatenate_gene_alignments.py aln_EF2_total.fa aln_atpA_total.fa aln_atpB_total.fa aln_cox1_total.fa aln_psaA_total.fa aln_psbA_total.fa aln_rbcL_total.fa 
```

### 1.2 Determine the optimal model of evolution per gene partition
```bash
modeltest-ng-static -d aa -i concatenated_gene_alignments.fasta -o output_modelfinder -q partitions -T raxml
```

### 1.3 Generate maximum likelihood phylogeny by running RAxML
```bash
raxmlHPC-PTHREADS-AVX -m PROTGAMMALG -T 45 -x 12342 -# 100 -q partitions \
-k -s concatenated_gene_alignments.fasta -p 27111997 \
-n rhodo -w run1
```
Using bootstraps to build final tree
```bash
raxmlHPC-PTHREADS-AVX -m PROTGAMMALG -T 20 -z RAxML_bootstrap.rhodo \
-q partitions \
-s concatenated_gene_alignments.fasta \
-p 27111997 -n rhodo -w run2
```

## 2. Perform ancestral state reconstruction of genome size estimates
```R
library(ape)
library(phytools)
library(picante)
library(diversitree)
library(corHMM)
library(qgraph)   
library(geiger)
library(coda)
library(plotrix)
library(RColorBrewer)
library(ggsci)
library(scales)
library(ggtree)
library(dplyr)

## 1. Load tree
tree_unrooted<-read.tree("RAxML_bestTree.RHODO")

## Root tree
tree_root<-root(tree_unrooted,outgroup = "Chlamydomonas_reinhardtii", resolve.root = T)
write.tree(tree_root, file="RAxML_bestTree_rooted.tre")

## 2. Load bootstrap replicates
trees_total<-read.tree("RAxML_bootstrap.rhodo")
treesBS_root<-root(trees_total,outgroup = "Chlamydomonas_reinhardtii", resolve.root = T)

## 3. Load table
traits <- read.table("Genome_sizes_red_algae.csv", sep=";", h=T, stringsAsFactors=F)
rownames(traits)<-traits[,1]
traits$V2<-as.double(traits$V2)

## 4. Prune tree to keep species with genome size estimates
tips<-tree_root$tip.label
trait <- traits[,1]
diff<-setdiff(tips,trait)
diff
pruned_tree<-drop.tip(tree_root,diff)
write.tree(pruned_tree, file="RAxML_bestTree_pruned.tre")

## 5. Find the appropriate model of trait evolution
ancestral_states_function=function(trait_number){
  
  results_anc <- matrix(ncol=8, nrow=length(treesBS_root))
  print(paste("analizing trait=", colnames(traits)[trait_number]), sep=" ")
  
  for (t in 1:length(treesBS_root)) { 
    trait <- traits[,trait_number]
    names(trait) <- traits$Name
    rem <- which(is.na(trait) == TRUE)
    
    trait <- trait[-rem]
    temp <- match.phylo.data(treesBS_root[[t]], trait)    
    trait_BM<-fitContinuous(temp$phy,temp$data,model="BM")
    trait_OU<-fitContinuous(temp$phy,temp$data,model="OU")
    trait_EB<-fitContinuous(temp$phy,temp$data,model="EB")
    
    AIC_mls <- function(LL, K) {
      -2*LL + 2*K	}
    AICc_mls <- function (LL, K,n) {
      x <- AIC_mls(LL, K)
      x + 2*K*(K+1)/(n-K-1) 
    }

    results_anc[t,1] <- names(traits)[trait_number]
    results_anc[t,2] <-  trait_BM$opt$lnL
    results_anc[t,3] <-  trait_BM$opt$aicc
    results_anc[t,4] <-  trait_OU$opt$lnL
    results_anc[t,5] <-  trait_OU$opt$aicc
    results_anc[t,6] <-  trait_EB$opt$lnL
    results_anc[t,7] <-  trait_EB$opt$aicc
    
    AICcs <- as.numeric(results_anc[t,c(3,5,7)])
    names(AICcs) <- c("BM","OU","EB")    
    bestModel <- which (AICcs == min(AICcs))
    bestModel
    AICcs
    
    if (all(bestModel==1)){
      bestModel<-AICcs[1]
    } else {
      diffs<-AICcs-min(AICcs)
      ifelse(diffs[1]<=2, bestModel<-AICcs[1],bestModel<-bestModel)
    }
    results_anc[t,8] <- names(bestModel)
  }
  colnames(results_anc) <- c("Trait", "BM_logL","BM_AICc","OU_logL", "OU_AICc","EB_logL","EB_AICc","Best_model")
  return(results_anc)
}

anc_T2 <- ancestral_states_function(2)

## 6. Perform ancestral state reconstruction
ASR_function_continuous=function(trait_number,bestModel){
  
  # Prepare traits and tree for each specific trait
  trait <- traits[,trait_number]
  names(trait) <- traits$Name
  rem <- which(is.na(trait) == TRUE)
  trait <- trait[-rem]
  temp <- match.phylo.data(pruned_tree, trait)
  
  pal<-c('#461220','#8c2f39','#b23a48','#ef6351','#f38375','#fed0bb','#878787','#4d4d4d','#1a1a1a')
  palrev<-rev(pal)
  pdf(paste("ASR_continuous_newColor",names(traits)[trait_number], ".pdf", sep=""))
  
  fitModel<-anc.ML(temp$phy,temp$data,model=bestModel)
  obj<-contMap(temp$phy,temp$data,anc.states=fitModel$ace,plot = FALSE,res=200,leg.txt=names(traits)[trait_number],lwd=4,outline = FALSE,legend=1)
  obj<-setMap(obj,colors=palrev)
  plot(obj,outline=FALSE,leg.txt=names(traits)[trait_number],legend=1)
  title(main=paste(names(traits)[trait_number]," model=",bestModel),line=-0.6)
  dev.off()
}

ASR_function_continuous(2,"BM")
```

