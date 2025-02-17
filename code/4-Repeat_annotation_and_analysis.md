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
