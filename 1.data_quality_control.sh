#!/bin/bash
set -euo pipefail

# ========= global parameters =========
SAMPLES=$1          # eg.：SRR12043508,SRR12043509,SRR18909870,SRR18909871
CONFIG=$2           # RNA-seq_config.sh

# ========= config =========
source ${CONFIG}

# ========= conda =========
set +u
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate ${CONDA_ENV}
set -u

# create folders
mkdir -p "${fastqc_dir}" "${multiqc_dir}"

echo "start running FastQC..."

IFS=',' read -ra SAMPLE_ARRAY <<< "${SAMPLES}"

for i in "${SAMPLE_ARRAY[@]}";do
    fastqc -t ${THREADS} -o ${fastqc_dir}/ -f fastq ${data_dir}/${i}/*.fq.gz
done

# parameters
# -t number of threads to use for parallel processing
# -o directory where FastQC output files will be saved
# -f input file format
# path to the input compressed FASTQ files for each sample

echo "start running MultiQC..."
multiqc ${fastqc_dir}/ -o ${multiqc_dir}/
