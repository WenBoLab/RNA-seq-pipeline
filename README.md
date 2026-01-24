# RNA-seq-pipeline

 This is a reproducible and modular pipeline for bulk RNA-seq data processing and downstream analysis, from raw FASTQ files to biological interpretation.

## ⚠️ Reminders

> **Note 1:** This upstream analysis is designed to run on a Linux system.

> **Note 2:** Make sure Conda is installed and properly configured before running the pipeline, as it is required to manage all dependencies.

> **Note 3:** It is recommended to activate a dedicated Conda environment for reproducibility:

```bash
conda create -n rnaseq_env
conda activate rnaseq_env
conda install bioconda::fastqc
```
## Overview

## This pipeline performs:

### 1. Raw data quality control
This step evaluates the quality of raw RNA-seq FASTQ files before any downstream processing.

Required software：fastqc and multiqc. you can install them using Conda.

```bash
fastqc -t 4 -o data/raw_fastqc/ -f fastq data/raw_data/*.fastq.gz
multiqc data/raw_fastqc/ -o data/raw_multiqc/
```
Export of results to an HTML to quickly assess data
<img width="1266" height="519" alt="image" src="https://github.com/user-attachments/assets/5d3786fa-c349-4638-85e8-832d49be39cd" />

### 2. Read trimming and filtering
This step removes adapter sequences and low-quality bases from raw RNA-seq reads.

Required software：trim-galore.
```bash
trim_galore --paired --cores 4 \
  data/raw_data/sample1_R1.fastq.gz \
  data/raw_data/sample1_R2.fastq.gz \
  -o data/trimmed_fastq/
```
### 3. Read alignment
This step maps trimmed reads to a reference genome to determine the origin of each read.

Required software：bowtie2 and samtools.
```bash
bowtie2 --local --very-sensitive -p 4 \
    -x reference_index_dir/GRCm39.genome.mm \
    -1 data/trimmed_fastq/sample1_R1.fastq.gz \
    -2 data/trimmed_fastq/sample1_R2.fastq.gz \
    -S data/alignment/sam/sample1.sam &> data/alignment/sam/bowtie2_summary/sample1.txt

#Convert SAM files to BAM files.
samtools view -bS -F 0x04 data/alignment/sam/sample1.sam -o data/alignment/bam/sample1.bam  
# Sort BAM file
samtools sort data/alignment/bam/sample1.bam -o data/alignment/bam/sample1.sort.bam
# build index
samtools index data/alignment/bam/sample1.sort.bam
```

### 4. Gene quantification
This step counts how many reads map to each gene.

Required software：subread.
```bash
featureCounts -T 4 -p -a reference_dir/Mus_musculus.GRCm39.114.chr.gtf -g exon_id -f \
	-o data/counts/read_counts.txt \
	data/alignment/bam/sample1.sort.bam
```
### 5. Downstream functional analysis and visualization
This step interprets differential expression results and visualizes patterns in the data, including pathway enrichment, clustering, and principal component analysis.
## ⚠️ Reminders

> **Note:** This step requires R and the necessary Bioconductor/R packages installed before running.

