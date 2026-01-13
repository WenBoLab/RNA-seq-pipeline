#!/bin/bash
set -euo pipefail

# set the path
data_dir="/data/RNA-seq/01/data/01.RawData"
fastqc_dir="/data/RNA-seq/01/results/fastqc"
multiqc_dir="/data/RNA-seq/01/results/multiqc"

# use conda environment
source /data/software/miniconda3/etc/profile.d/conda.sh
conda activate RNA_env

# create folders
mkdir -p "${fastqc_dir}" "${multiqc_dir}"

echo "start running FastQC..."
# step1. fastqc for raw quality control
for i in {12..19};do
    fastqc -t 4 -o ${fastqc_dir}/ -f fastq ${data_dir}/CDNA${i}/*.fq.gz
done

echo "start running MultiQC..."

# step2. use multiqc to merge fastqc results
multiqc ${fastqc_dir}/ -o ${multiqc_dir}/


# Parameter Description:
# fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] [-c contaminant file] seqfile1 .. seqfileN
# -o: output dir to save results
# -t: threads number
