#!/bin/bash
set -euo pipefail

# set the path
data_dir="/data/RNA-seq/01/data/01.RawData"
fastqc_dir="/data/RNA-seq/01/results/fastqc"
multiqc_dir="/data/RNA-seq/01/results/multiqc"

# ========= 参数 =========
SAMPLES=$1          # 例如：14,18,IgG
CONFIG=$2           # RNA-seq_config.sh

# ========= 读取配置 =========
source ${CONFIG}

# ========= conda =========
set +u
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
set -u

# create folders
mkdir -p "${fastqc_dir}" "${multiqc_dir}"

echo "start running FastQC..."
# step1. fastqc for raw quality control

for i in "${SAMPLE_ARRAY[@]}";do
    fastqc -t ${THREADS} -o ${fastqc_dir}/ -f fastq ${data_dir}/${i}/*.fq.gz
done

# parameters
# {12..19}: the sample name, this can change with different
# -t number of threads to use for parallel processing
# -o directory where FastQC output files will be saved
# -f input file format
# path to the input compressed FASTQ files for each sample

echo "start running MultiQC..."

# step2. use multiqc to merge fastqc results
multiqc ${fastqc_dir}/ -o ${multiqc_dir}/

