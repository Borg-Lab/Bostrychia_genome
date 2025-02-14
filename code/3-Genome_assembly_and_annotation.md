# Genome assembly and annotation of Bostrychia moritziana

## 1. Genome assembly 
A draft hybrid assembly from ONT and PacBio HiFi reads was generated separately from the male
and female gametophyte sequencing data. Based on the overrepresentation of Nanopore reads compared to the male dataset, the
assembler Flye was used for the female and the assembler hifiasm for the male.
To resolve the assemblies at a chromosome-level and to filter contigs belonging to the associated microbiome, a 3D de novo assembly Hi-C pipeline was applied using juicer and 3D-DNA.

### 1.1 Male assembly
#### 1.1.1 Draft assembly
```bash
singularity exec hifiasm_0.16.1--h2e03b76_0.sif hifiasm \
-o Bm3235_M.asm -l0 --hg-size 1g -t40 \
--ul Bm3235M_reads_total.fastq \
Bm3235_M_PB_DC.fastq
```
Polishing the draft assembly generated with hifiasm using both, PacBio HiFi and Nanopore reads.
```bash
minimap2 \
-ax map-ont Bm3235_M.p_ctg.fa \
Bm3235_M_PB_NP_allReads.fastq \
> Bm3235_M_output.sam
```
```bash
racon \
Bm3235_M_PB_NP_allReads.fastq \
Bm3235_M_output.sam Bm3235_M.p_ctg.fa \
> Bm3235M_hifiasm_racon.fasta
```
#### 1.1.2 Hi-C
Using the generate_site_positions.py script from the Juicer tool to identify the genomic localization of DpnII restriction sites. Subsequently, we employed Juicer to generate Hi-C contact maps.
```bash
generate_site_positions.py DpnII Bm3235M_hifiasm_racon \
Bm3235M_hifiasm_racon.fasta
 
bwa index -a bwtsw Bm3235M_hifiasm_racon.fasta
samtools faidx Bm3235M_hifiasm_racon.fasta
awk '{ print $1 "\t" $2 }' Bm3235M_hifiasm_racon.fasta.fai > Bm3235M_hifiasm_racon.fasta.chrom.sizes
 
#Run the juicer pipeline to generate the contact map
juicer.sh -g Bm3235M_hifiasm_racon -z Bm3235M_hifiasm_racon.fasta -s DpnII -a 'Bm3235M_scaffolding' \
-p Bm3235M_hifiasm_racon.fasta.chrom.sizes -y Bm3235M_hifiasm_racon_DpnII.txt \
-D /ebio/abt5_projects/software/code/miniconda3/envs/Hi-C_scaffolding -t 40
```
Perform 3D-DNA assembly. This step facilitated manual correction of scaffolding a refinement of the order and orientation of each assembly. 
```bash
run-asm-pipeline.sh Bm3235M_hifiasm_racon.fasta merged_nodups.txt
```
The Hi-C contact maps were visualiz using Juicebox and some manual edits were done. Afterwards, we incorporated manual edits and finalized the genome assembly in fasta format.
```bash
run-asm-pipeline-post-review.sh --sort-output -r Bm3235M_hifiasm_racon.0.review.assembly \
Bm3235M_hifiasm_racon.fasta merged_nodups.txt
```
After chromosome 1 (microbial community) was removed, the final gap closing and error correction was performed.
```
tgsgapcloser \
--scaff Bm3235M_hifiasm_racon_HiC.fasta \
--reads Bm3235_M_PB_NP_allReads.fasta \
--output tgs_out \
--racon /ebio/abt5_projects/small_projects/rpetroll/mambaforge/bin/racon \
--tgstype pb \
>pipe.log 2>pipe.err
```

### 1.2 Female assembly
#### 1.2.1 Draft assembly
```bash
flye \
--subassemblies Bm3235_F_NP_total_reads.fastq \
Bm3235_F_PB_DC.fastq \
--out-dir flye_hybrid_female -t 50 --scaffold
```
Polishing the draft assembly generated with flye using both, PacBio HiFi and Nanopore reads.
```bash
minimap2 \
-ax map-ont assembly.fasta \
Bm3235F_NP_PB_allReads.fastq \
> Bm3235_F_output.sam
```
```bash
racon \
Bm3235F_NP_PB_allReads.fastq \
Bm3235_F_output.sam assembly.fasta \
> Bm3235F_flye_racon.fasta
```

#### 1.1.2 Hi-C
Using the generate_site_positions.py script from the Juicer tool to identify the genomic localization of DpnII restriction sites. Subsequently, we employed Juicer to generate Hi-C contact maps.
```bash
generate_site_positions.py DpnII Bm3235F_flye_racon \
Bm3235F_flye_racon.fasta
 
bwa index -a bwtsw Bm3235F_flye_racon.fasta
samtools faidx Bm3235F_flye_racon.fasta
awk '{ print $1 "\t" $2 }' Bm3235F_flye_racon.fasta.fai > Bm3235F_flye_racon.fasta.chrom.sizes
 
#Run the juicer pipeline to generate the contact map
juicer.sh -g Bm3235F_flye_racon -z Bm3235F_flye_racon.fasta -s DpnII -a 'Bm3235F_scaffolding' \
-p Bm3235F_flye_racon.fasta.chrom.sizes -y Bm3235F_flye_racon_DpnII.txt \
-D /ebio/abt5_projects/software/code/miniconda3/envs/Hi-C_scaffolding -t 40
```
Perform 3D-DNA assembly. This step facilitated manual correction of scaffolding a refinement of the order and orientation of each assembly. 
```bash
run-asm-pipeline.sh Bm3235F_flye_racon.fasta merged_nodups.txt
```
The Hi-C contact maps were visualiz using Juicebox and some manual edits were done. Afterwards, we incorporated manual edits and finalized the genome assembly in fasta format.
```bash
run-asm-pipeline-post-review.sh --sort-output -r Bm3235F_flye_racon.0.review.assembly \
Bm3235F_flye_racon.fasta merged_nodups.txt
```
After chromosome 1 (microbial community) was removed, the final gap closing and error correction was performed.
```bash
tgsgapcloser \
--scaff Bm3235F_flye_racon.fasta \
--reads Bm3235_F_PB_NP_allReads.fasta \
--output tgs_out \
--racon /ebio/abt5_projects/small_projects/rpetroll/mambaforge/bin/racon \
--tgstype pb \
>pipe.log 2>pipe.err
```
We standardized chromosome nomenclature of both assemblies using a uniform format (i.e., Chr[number]) by renaming the chromosomes accordingly. Subsequent analyses revealed the integration of a bacterial contig on the female sex chromosome caused by assembly gaps, which was manually removed to ensure assembly accuracy. In addition, the complete mitochondrial genome was found to be integrated into chromosome 10 in both the male and female assembly and was also removed prior to further analyses. Due to the higher contiguity and completeness of the male assembly, as indicated by a higher N50 and fewer gaps, we used it as the reference genome for our analyses. To ensure the representation of both sex-specific regions, we added the female SDR to the male assembly.


## 2. Genome annotation
Prior to annotating protein-coding genes in the Bostrychia genome, each assembly was initially soft-masked for repetitive regions and TEs.
Gene prediction was performed using RNA-seq data alone (BRAKER1) and using both RNA-seq data and a published manually curated orthologous set of red algal protein sequences (BRAKER3).
Since all steps were performed similar in both the male and female, I will just show the code of annotating the male genes here.

### 2.1 Preparing orthologous set of red algal protein sequences
```bash
orthofinder \
-f /ebio/scratch/rpetroll/orthofinder/public_red_algae -M msa -X
```
[OGFilter](https://github.com/pnatsi/OGFilter)
```bash
python OGFilter.py -g Orthogroups.GeneCount.tsv -s Orthogroup_Sequences/ -o ./outdir -min_species 0.8 -max_copies 100
```
### 2.2 Running genome annotation with BRAKER
```bash
echo "1. Running hisat2"
hisat2-build Bm3235_M_assembly_v3_masked.fa Bm3235_M_assembly
 
hisat2 Bm3235_M_assembly \
-1 /tmp/global2/rpetroll/rnaseq/trimming/Bm3235F1_S36_R1_001_val_1.fq.gz,\
/tmp/global2/rpetroll/rnaseq/trimming/Bm3235F2_S37_R1_001_val_1.fq.gz,\
/tmp/global2/rpetroll/rnaseq/trimming/Bm3235F3_S38_R1_001_val_1.fq.gz \
-2 /tmp/global2/rpetroll/rnaseq/trimming/Bm3235F1_S36_R2_001_val_2.fq.gz,\
/tmp/global2/rpetroll/rnaseq/trimming/Bm3235F2_S37_R2_001_val_2.fq.gz,\
/tmp/global2/rpetroll/rnaseq/trimming/Bm3235F3_S38_R2_001_val_2.fq.gz \
--rna-strandness RF -p 50 --dta -S female_RF_trimmed.sam
 
hisat2 Bm3235_M_assembly \
-1 /tmp/global2/rpetroll/rnaseq/trimming/Bm3235M1_S33_R1_001_val_1.fq.gz,\
/tmp/global2/rpetroll/rnaseq/trimming/Bm3235M2_S34_R1_001_val_1.fq.gz,\
/tmp/global2/rpetroll/rnaseq/trimming/Bm3235M3_S35_R1_001_val_1.fq.gz \
-2 /tmp/global2/rpetroll/rnaseq/trimming/Bm3235M1_S33_R2_001_val_2.fq.gz,\
/tmp/global2/rpetroll/rnaseq/trimming/Bm3235M2_S34_R2_001_val_2.fq.gz,\
/tmp/global2/rpetroll/rnaseq/trimming/Bm3235M3_S35_R2_001_val_2.fq.gz \
--rna-strandness RF -p 50 --dta -S male_RF_trimmed.sam
 
hisat2 Bm3235_M_assembly \
-1 /tmp/global2/rpetroll/rnaseq/trimming/Bm3235T1_S39_R1_001_val_1.fq.gz,\
/tmp/global2/rpetroll/rnaseq/trimming/Bm3235T2_S40_R1_001_val_1.fq.gz,\
/tmp/global2/rpetroll/rnaseq/trimming/Bm3235T3_S41_R1_001_val_1.fq.gz \
-2 /tmp/global2/rpetroll/rnaseq/trimming/Bm3235T1_S39_R2_001_val_2.fq.gz,\
/tmp/global2/rpetroll/rnaseq/trimming/Bm3235T2_S40_R2_001_val_2.fq.gz,\
/tmp/global2/rpetroll/rnaseq/trimming/Bm3235T3_S41_R2_001_val_2.fq.gz \
--rna-strandness RF -p 50 --dta -S tetrasporophyte_RF_trimmed.sam
 
echo "2. Converting SAM to BAM"
samtools sort -o tetrasporophyte_RF_trimmed_sort.bam tetrasporophyte_RF_trimmed.sam
rm tetrasporophyte_RF_trimmed.sam
samtools sort -o female_RF_trimmed_sort.bam female_RF_trimmed.sam
rm female_RF_trimmed.sam
samtools sort -o male_RF_trimmed_sort.bam male_RF_trimmed.sam
rm male_RF_trimmed.sam
 
echo "3. Running BRAKER3"
 
### Using both, RNA-seq reads and red algal database
export AUGUSTUS_CONFIG_PATH=/ebio/abt5_projects/small_projects/rpetroll/mambaforge/envs/augustus/config
export AUGUSTUS_BIN_PATH=/ebio/abt5_projects/small_projects/rpetroll/mambaforge/envs/augustus/bin
export SINGULARITY_BIND="/ebio/abt5_projects:/ebio/abt5_projects"
 
singularity exec braker3.sif braker.pl \
--species=Bm3235_M_v3_masked \
--genome=Bm3235_M_assembly_v3_masked.fa \
--prot_seq=RedAlgae_DB.fa \
--bam=female_RF_trimmed_sort.bam,\
male_RF_trimmed_sort.bam,\
tetrasporophyte_RF_trimmed_sort.bam \
--workingdir=braker_out_prot_RNAseq_trimmed --threads=40 --softmasking \
--AUGUSTUS_CONFIG_PATH=/ebio/abt5_projects/small_projects/rpetroll/mambaforge/envs/augustus/config
 
### Using RNA-seq data only
singularity exec braker3.sif braker.pl \
--species=Bm3235_M_v3_masked_RNAonly \
--genome=Bm3235_M_assembly_v3_masked.fa \
--bam=female_RF_trimmed_sort.bam,\
male_RF_trimmed_sort.bam,\
tetrasporophyte_RF_trimmed_sort.bam \
--softmasking --workingdir=braker_out_RNAseq_trimmed --threads=40 \
--AUGUSTUS_CONFIG_PATH=/ebio/abt5_projects/small_projects/rpetroll/mambaforge/envs/augustus/config

echo "done"
```

### 2.3 Combining annotations using TSEBRA
Integrate both rounds of gene prediction.
```bash
singularity exec braker3.sif tsebra.py -g braker_prot_rnaseq.gtf,braker_rnaseq.gtf \
-c default.cfg -e hintsfile_prot_rnaseq.gff,hintsfile_rnaseq.gff -o braker_combined_Male.gtf
```
Using script from [Tsebra](https://github.com/Gaius-Augustus/TSEBRA) to rename merged gtf. 
```bash
rename_gtf.py --gtf braker_combined_Male.gtf --out braker_combined_Male_renamed.gtf
```
Using custom script [rename_gtf.py](https://github.com/Borg-Lab/Bostrychia_genome/tree/main/scripts/rename_gtf.py) to standardize gene IDs to follow a systematic naming convention. The format is structured as Bm01g000010.t1, where Bm refers to the species (Bostrychia moritziana) 01 indicates the chromosome number, g denotes "gene", 000010 represents the unique gene identifier, and .t1 specifies the transcript number.
Afterwards, I am using a script from [Augustus](https://github.com/Gaius-Augustus/Augustus) to get cds and protein fasta out of the gtf file.
```bash
getAnnoFastaFromJoingenes.py -g Bm3235_M_assembly_v3_masked.fa -o braker_Bm3235_M_assembly_v3_combined -f braker_combined_Male_renamed.gtf
```
To simplify downstream analyses, only the longest transcript isoform was used. I used a script from [Genome-Zoo](https://github.com/Rensing-Lab/Genome-Zoo) to filter.

```bash
split_isoform_dot.py -o . -i braker_Bm3235_M_assembly_v3_combined.aa -c "" -p "."
```
The gene set was additionally filtered to exclude TE-related genes, which were defined as sequences that were 100% soft-masked by the repeat annotation. 

Get masked percentage of genes:
```bash
#!/bin/bash

file_path="braker_Bm3235_M_assembly_v3_combined.codingseq"  
output_file="percentage_repeats.txt" 

awk '/^>/ {if (seq) {print header, seq; seq="";} header=$0; next;} {seq = seq $0;} END {print header, seq;}' "$file_path" |
while read -r header sequence; do
    total_length=${#sequence}
    lower_count=$(echo "$sequence" | tr -cd 'a-z' | wc -c)
    lower_percentage=$(awk "BEGIN {printf \"%.2f\", ($lower_count / $total_length) * 100}")
    echo -e "$header\t$lower_percentage%" >> "$output_file"
done
```
