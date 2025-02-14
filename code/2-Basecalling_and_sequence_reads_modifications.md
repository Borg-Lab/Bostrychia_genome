# Basecalling of Nanopore (ONT) reads and preparation of PacBio HiFi reads

The isolated HMW DNA was subjected to sequencing using two long-read technologies, 
namely Oxford Nanopore Technologies (ONT) MinION and PacBio SMRTbell sequencing. 

## 1. Basecalling of ONT reads

One male and one female ONT library were prepared and sequenced using FLO-MIN106D flow cells. 
An additional female ONT library was prepared using a FLO-MIN114 flow cell. 
For the two libraries sequenced on FLO-MIN106D flow cells, the dna_r9.4.1_450bps_hac configuration file was used. 
For two runs of the same library using the FLO-MIN114 flow cell, the dna_r10.4.1_e8.2_400bps_sup and the dna_r10.4.1_e8.2_260bps_sup 
configuration files were applied. 



## 2. Preparation of PacBio HiFi reads

### 2.1 Process PacBio subreads into HiFi reads 
using ccs 6.4.0 with the --min-rq=0.88 option111. 

### 2.2 Demultiplexing of reads 
was performed using lima 2.7.1 with the --split-named option112. 

### 2.3 Realigning of subreads to the generated HiFi reads 
separately for the male and female datasets using actc v0.2.0113. 

### 2.4 Error correction of HiFi reads
Finally, error correction of the HiFi reads was performed using DeepConsensus v1.2.0 with the --batch_size=2048 and --batch_zmws=1000 options114.
