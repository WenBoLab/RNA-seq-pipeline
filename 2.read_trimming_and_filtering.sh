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
mkdir -p "${clean_dir}"

echo "start running Fastp..."
IFS=',' read -ra SAMPLE_ARRAY <<< "${SAMPLES}"

for i in "${SAMPLE_ARRAY[@]}"; do
    fastp -i ${data_dir}/${i}/"$i"_1.fq.gz -o ${clean_dir}/"$i"_1_clean.fq.gz \
    -I ${data_dir}/${i}/"$i"_2.fq.gz -O ${clean_dir}/"$i"_2_clean.fq.gz \
    -q 20 -u 30 -n 5 -l 30 --thread ${THREADS} --html ${clean_dir}/${i}.fastp.html
done

echo "start running FastQC after clean..."
mkdir -p "${clean_fastqc_dir}" "${clean_multiqc_dir}"
fastqc -t ${THREADS} -o ${clean_fastqc_dir}/ -f fastq ${clean_dir}/*.fq.gz
echo "start running MultiQC after clean..."
multiqc ${clean_fastqc_dir}/ -o ${clean_multiqc_dir}/
