# Annotation and analysis of repeats in _Bostrychia_ and further red algae

Prior to annotating protein-coding genes in the _Bostrychia_ genome, each assembly was initially soft-masked for repetitive regions and TEs. First, a comprehensive TE library was generated containing consensus sequences of all TE families and repetitive sequences.

## 1. Generating a repeat library 
```bash
## 1. Running RepeatModeler
singularity exec tetools_latest.sif BuildDatabase \
-name Bm3235_M Bm3235_M_assembly_v3.fasta
 
singularity exec tetools_latest.sif RepeatModeler \
-database Bm3235_M -LTRStruct \
-threads 40 > Out_RepeatModeler
```
In the case of the male genome assembly of _Bostrychia_, the initial de novo repeat library was further manually curated to improve TE classification by both the modification of TE annotations and the structural refinement of selected consensus sequences. Each consensus sequence in the repeat library was compared to the eukaryotic repetitive sequences in the RepBase database v29_02. The TE annotation of each consensus sequence was manually reviewed against the assigned RepBase annotation and adjusted in cases of strong e-value support. Selected TE families were also manually curated. In the initial TE library, the most abundant families were classified as Mutator-like DNA transposons (MULEs). To verify and refine this classification, the homology of several of these families was assessed by sequence alignment, which resulted in five distinct curated families. Based on homology to transposase proteins, terminal sequence motifs and target site duplication lengths, these five families were re-classified as EnSpm superfamily TEs of the Plavaka lineage. The original misclassification stemmed from homology between the Plavaka accessory genes found both in _Bostrychia_ Plavaka elements and _Chondrus crispus_ Mutator elements. In addition, the consensus sequence of one abundant family belonging to LTR Gypsy elements was manually curated and added to the repeat library.

## 2. Masking TEs and repeats by using the TE library
```
for i in Bm3235_M_assembly_v3.fasta
 
do
 
NAME=$(basename $i .fasta)
 
echo "2. Running TandemRepeatsFinder"
 
mkdir TRF
cd TRF
 
trf ../$i 2 7 7 80 10 50 2000 -f -d -m -ngs > "$NAME"_TRF.out
perl trf2sats.pl --in "$NAME"_TRF.out --prefix "$NAME"_TRF
 
## Combine satellites and microsatellites
cat "$NAME"_TRF.microsats.bed "$NAME"_TRF.sats.bed | sort -k1,1 -k2n,2n | bedtools merge -i stdin > $NAME.trf.bed
 
echo "3. Running RepeatMasker"

cd .. 
mkdir RepeatMasker_V2
cd RepeatMasker_V2
 
conda activate /ebio/abt5_projects/software/code/miniconda3/envs/RepeatMasker
 
RepeatMasker \
-lib ../RepeatModeler/Bm3235M_MC_lib.fa \
-pa 40 -xsmall -gff -dir Repeatmasking_out \
../$i
 
conda deactivate
 
echo "4. Mask satellites from TRF"
cd ..
mkdir FinalCombinedLibrary
cd FinalCombinedLibrary
 
bedtools maskfasta -fi ../RepeatMasker/Repeatmasking_out/"$i".masked -bed ../TRF/Bm3235_M_assembly_v3.trf.bed -fo "$NAME"_masked.fa -soft
 
echo "DONE"
done
```
In case of the other red algae, the repeatmasking step was directly performed after running RepeatModeler, without the addition of some manual curated consensus sequences of further manual edits.

## 3. Visualization using a TE landscape
To visualize the final TE complement in _Bostrychia_, a landscape of TE divergence was generated based on the percent divergence from the consensus sequence of their assigned family, which was derived from the RepeatMasker output.
The code how to plot the landscape can be found [here](https://github.com/Borg-Lab/Bostrychia_genome/tree/main/scripts/TE_landscape.R)

## 4. Overlap-based approach to detect _Plavaka_ elements across red algae
Since it was infeasible to assign _Plavaka_ elements through manual curation of the TE libraries across all investigated red algal genomes, an overlap-based approach was used as an alternative. Homologs of the Bostrychia PAG genes were first identified in each red algal TE library with tBLASTn (using the evalue 0.001 option). The PAG2 gene hits in each species were then used to determine whether the _Plavaka_ transposase and the other two PAG genes lied within a 15 kb flanking region (total region 30 kb) on the same consensus sequence, to ensure that only complete _Plavaka_ elements were recovered.

To run the analysis, a fasta file containing the PAG2 genes, as well as a fasta file containing all PAG and transposase sequences of as many species as possible is needed. In addition, the genome assembly fasta file and the RepeatMasker output (.out file) file of the species to be analyzed is needed. 
### Run the analysis
During the run, the main script [Overlap_analysis_Plavaka_elements.sh](https://github.com/Borg-Lab/Bostrychia_genome/tree/main/scripts/Overlap_analysis_Plavaka_elements.sh) uses two additional scripts, [add_flanks.sh](https://github.com/Borg-Lab/Bostrychia_genome/tree/main/scripts/add_flanks.sh) and [process_ids.py](https://github.com/Borg-Lab/Bostrychia_genome/tree/main/scripts/process_ids.py). 
```bash
# Overlap_analysis_Plavaka_elements.sh genome_fasta RepeatMasker_out_file 5_letter_abbreviation_of_species_name
Overlap_analysis_Plavaka_elements.sh \
../PYRYE/GCA_009829735.1_ASM982973v1_genomic.fna.unmasked \
../PYRYE/GCA_009829735.1_ASM982973v1_genomic.fna.unmasked.out \
PYRYE
```

## 5. Ancestral state reconstruction of _EnSpm_ elements

To perform the ancestral state reconstruction of _EnSpm_ elements, their presence was treated as a discrete trait. The appropriate model of trait evolution was assessed as described in [Red algal phylogeny and ancestral state reconstruction of genome sizes](https://github.com/Borg-Lab/Bostrychia_genome/blob/main/code/1-Genome_size_expansion.md), with the Early-burst model selected for analysis. 
```
## Required R packages 
library(ape)
library(phytools)
library(corHMM)
library(geiger)

tree<- read.newick("red_algal_phylogeny.tre")

x <- read.table("Plavaka.csv", sep=";", h=T, stringsAsFactors=F)
x$Name
rownames(x)<-x[,1]
trait <- x[,2]
names(trait) <- x$Name
temp <- match.phylo.data(tree, trait)

col_traits <- list(Plavaka=c("#ff7f00", "#984ea3","blue"))
pdf(paste("AncState_recon_",names(x)[2], ".pdf", sep=""), width=6, height=8)
best_model <- "ER"
trait <- x[,2]
names(trait) <- x$Name

temp <- match.phylo.data(tree, trait)
data <- cbind(names(temp$data) , data.frame(temp$data)) 
rownames(data)<-NULL

trait_best <- rayDISC(temp$phy,data,ntraits=1, model=best_model, node.states=c("marginal"), root.p="maddfitz")
col_i <- col_traits [ which(names(col_traits) == colnames(x)[2]) ]
plotRECON(ladderize(temp$phy, right=FALSE),trait_best$states, piecolors=unlist(col_i), pie.cex=0.65, cex=0.6, label.offset=0.015, title=paste(names(x)[2], "best =",best_model, sep=" "),width=6, height=7)
tip_color<- unlist(col_i)
tip_color
names(tip_color)<- names(table(trait))
trait_col <- temp$data
for (l in 1:length(tip_color)){ 
  trait_col [which(trait_col==names(tip_color)[l])]<- tip_color[l]
}
tiplabels(pch=20, col = trait_col, cex = 1.0, adj=0.5)
add.simmap.legend(colors = c("#ff7f00","#984ea3"), x=-0.01, y=23,prompt = FALSE, leg=c("Not present","Present"))
dev.off()
```
