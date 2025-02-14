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
```
tgsgapcloser \
--scaff Bm3235F_flye_racon.fasta \
--reads Bm3235_F_PB_NP_allReads.fasta \
--output tgs_out \
--racon /ebio/abt5_projects/small_projects/rpetroll/mambaforge/bin/racon \
--tgstype pb \
>pipe.log 2>pipe.err
```

## 2. Genome annotation
