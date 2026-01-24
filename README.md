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
This will generate an HTML report to quickly assess data
<img width="1266" height="519" alt="image" src="https://github.com/user-attachments/assets/5d3786fa-c349-4638-85e8-832d49be39cd" />

### 2. Read trimming and filtering
This step removes adapter sequences and low-quality bases from raw RNA-seq reads.

Required software：trim-galore or fastp.

the example of trim-galore
```bash
trim_galore --paired --cores 4 \
  data/raw_data/sample1_R1.fastq.gz \
  data/raw_data/sample1_R2.fastq.gz \
  -o data/trimmed_fastq/
```

the example of fastp
```bash
fastp -i data/raw_data/sample1_R1.fastq.gz -o data/clean_data/sample1_R1_clean.fq.gz \
    -I data/raw_data/sample1_R2.fastq.gz -O data/clean_data/sample1_R2_clean.fq.gz \
    -q 20 -u 30 -n 5 -l 30
```	

### 3. Read alignment
This step maps trimmed reads to a reference genome to determine the origin of each read.

Required software：bowtie2 or STAR and samtools.

the example of bowtie2
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
the example of STAR
```bash
STAR --runThreadN 4 \
       --genomeDir index_dir \
       --readFilesIn sample1_R1.fastq.gz sample1_R2.fastq.gz \
       --readFilesCommand zcat \
       --outFileNamePrefix alignment_dir \
       --twopassMode Basic \
       --sjdbGTFfile GTF_dir \
       --quantMode  GeneCounts \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattributes NH HI AS nM XS \
       --outFilterMatchNminOverLread 0.3 \
       --outFilterScoreMinOverLread 0.3
```

### 4. Gene quantification
This step counts how many reads map to each gene.

Required software：subread or HTSeq.

the example of subread
```bash
featureCounts -T 4 -p -a reference_dir/Mus_musculus.GRCm39.114.chr.gtf -g exon_id -f \
	-o data/counts/read_counts.txt \
	data/alignment/bam/sample1.sort.bam
```

the example of HTSeq
```bash
htseq-count \
  -f bam -r name -s yes -t exon -i gene_id -m union -n 4 \
  --with-header \
  data/alignment/bam/sample1.sort.bam \
  GTF_dir > data/counts/read_counts.txt
```

### 5. Downstream functional analysis and visualization
This step interprets differential expression results and visualizes patterns in the data, including pathway enrichment, clustering, and principal component analysis.
## ⚠️ Reminders

> **Note:** This step requires R and the necessary Bioconductor/R packages installed before running.

###### (1) Differentially expressed genes (DEGs)
This step aims to identify differentially expressed genes between conditions and visualize the results.
```R
counts <- read.table("data/counts_matrix.tsv", header=TRUE, row.names=1) # load expression data
meta   <- read.table("data/metadata.tsv", header=TRUE, row.names=1) # sample information
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, meta, design=~condition)
vsd <- vst(dds)
plotPCA(vsd, intgroup="condition")

<img width="570" height="263" alt="image" src="https://github.com/user-attachments/assets/20451504-2091-4f63-96c3-4d3b23788731" />

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","Treatment","Control"))
res.sig <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ] # apply significance threshold
```
<img width="513" height="514" alt="image" src="https://github.com/user-attachments/assets/a4fb3b0f-55f2-4230-af2f-1069c25fbabd" />

###### (2) Functional Enrichment Analysis
This step aims to interpret differentially expressed genes through functional enrichment analysis.
```R
# GO analysis
library(clusterProfiler)
library(org.Mm.eg.db)
genes <- rownames(res.sig)
ego <- enrichGO(gene=genes,
                OrgDb=org.Mm.eg.db,
                keyType="SYMBOL",
                ont="ALL")

<img width="758" height="759" alt="image" src="https://github.com/user-attachments/assets/516407b0-2756-41e3-81fe-d5ede3685e7f" />

# KEGG analysis
library(biomaRt)
gene2symbols <- getBM(attributes = c('external_gene_name','entrezgene_id'),
                      filters = "external_gene_name", values = genes, mart = gene_info)

ekegg <- enrichKEGG(
  gene = gene2symbols$entrezgene_id,
  keyType = 'kegg',
  organism = 'mmu'
)

<img width="547" height="457" alt="image" src="https://github.com/user-attachments/assets/6b2a7a8a-fcf9-4107-8baf-cd8d9acb2d4e" />

```
