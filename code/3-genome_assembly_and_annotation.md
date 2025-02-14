# Genome assembly and annotation of Bostrychia moritziana

## 1. Genome assembly 
A draft hybrid assembly from ONT and PacBio HiFi reads was generated separately from the male
and female gametophyte sequencing data. Based on the overrepresentation of Nanopore reads compared to the male dataset, the
assembler Flye was used for the female and the assembler hifiasm for the male.

### 1.1 Male assembly
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

### 1.2 Female assembly
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
