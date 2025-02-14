# Basecalling of Nanopore (ONT) reads and preparation of PacBio HiFi reads

The isolated HMW DNA was subjected to sequencing using two long-read technologies, 
namely Oxford Nanopore Technologies (ONT) MinION and PacBio SMRTbell sequencing. 

## 1. Basecalling of ONT reads

One male and one female ONT library were prepared and sequenced using FLO-MIN106D flow cells. 
An additional female ONT library was prepared using a FLO-MIN114 flow cell. 
For the two libraries sequenced on FLO-MIN106D flow cells, the dna_r9.4.1_450bps_hac configuration file was used. 
For two runs of the same library using the FLO-MIN114 flow cell, the dna_r10.4.1_e8.2_400bps_sup and the dna_r10.4.1_e8.2_260bps_sup 
configuration files were applied. 

### Basecalling of male reads
```bash
guppy_basecaller  -x "cuda:0 cuda:1" -i /ebio/scratch/rpetroll/Bm3235_M_seq1 \
-s basecall_male_out --recursive --compress_fastq -c dna_r9.4.1_450bps_hac.cfg \
--print_workflows --trim_adapters --trim_primers --num_callers 20 --gpu_runners_per_device 24 --chunks_per_runner 480
```

### Basecalling of female reads
First sequencing run using a FLO-MIN106D flow cell.
```bash
guppy_basecaller  -x "cuda:0 cuda:1" -i /ebio/scratch/rpetroll/Bm3235_F_seq1 \
-s basecalling_female_out --recursive --compress_fastq -c dna_r9.4.1_450bps_hac.cfg \
--print_workflows --trim_adapters --trim_primers --num_callers 20 --gpu_runners_per_device 24 --chunks_per_runner 480
```
Resequencing sequencing run using a FLO-MIN114 flow cell. I used the same flow cell and loaded library twice.
```bash
guppy_basecaller  -x "cuda:0" -i /ebio/scratch/rpetroll/nanopore/Bm3235_F_reseq_run1 \
-s basecalling_female_reseq_out_1 --recursive --compress_fastq -c dna_r10.4.1_e8.2_400bps_sup.cfg \
--print_workflows --trim_adapters --trim_primers \
--num_callers 20 --gpu_runners_per_device 24 --chunks_per_runner 480
```
```bash
guppy_basecaller  -x "cuda:1" -i /ebio/scratch/rpetroll/nanopore/Bm3235_F_reseq_run2 \
-s basecalling_female_reseq_out_2 --recursive --compress_fastq -c dna_r10.4.1_e8.2_260bps_sup.cfg \
--print_workflows --trim_adapters --trim_primers \
--num_callers 20 --gpu_runners_per_device 24 --chunks_per_runner 480
```

## 2. Preparation of PacBio HiFi reads
As an output of sequencing we received multiplexed subreads of the male and the female. 

### 2.1 Process PacBio subreads into HiFi reads 
```bash
ccs --min-rq=0.88 -j 15 m64079_230302_100317.subreads.bam subreads.ccs.bam
```
### 2.2 Demultiplexing of reads 
```bash
lima -j 20 subreads.ccs.bam barcodes.fasta css_reads.demux.bam --split-named
```

### 2.3 Realigning of subreads to the generated HiFi reads 
Male reads:
```bash
actc -j "$(nproc)" m64079_230302_100317.subreads.bam css_reads.demux.bc1002--bc1002.bam /ebio/scratch/rpetroll/actc/male/bc1002_male_subreads_to_css.bam
```
Female reads:
```bash
actc -j "$(nproc)" m64079_230302_100317.subreads.bam css_reads.demux.bc1003--bc1003.bam /ebio/scratch/rpetroll/actc/female/bc1003_female_subreads_to_css.bam
```
### 2.4 Error correction of HiFi reads
Finally, error correction of the HiFi reads was performed using DeepConsensus v1.2.0 with the --batch_size=2048 and --batch_zmws=1000 options114.
Male reads:
```bash
singularity exec --nv ${SINGULARITY_SANDBOX} deepconsensus run \
        --subreads_to_ccs=bc1002_male_subreads_to_css.bam \
        --ccs_bam=css_reads.demux.bc1002--bc1002.bam \
        --checkpoint /ebio/scratch/rpetroll/deepconsensus/male/model/checkpoint \
        --output /ebio/scratch/rpetroll/deepconsensus/male/Bm3235_M_PB_DC.fastq \
        --batch_size=2048 --batch_zmws=1000 \
        --cpus=16 \
```
Female reads:
```bash
singularity exec --nv ${SINGULARITY_SANDBOX} deepconsensus run \
        --subreads_to_ccs=bc1003_female_subreads_to_css.bam \
        --ccs_bam=css_reads.demux.bc1003--bc1003.bam \
        --checkpoint /ebio/scratch/rpetroll/deepconsensus/female/model/checkpoint \
        --output /ebio/scratch/rpetroll/deepconsensus/female/Bm3235_F_PB_DC.fastq \
        --batch_size=2048 --batch_zmws=1000 \
        --cpus=16
```

