#!/bin/bash
set -eo pipefail

# ========= global parameters =========
SAMPLES=$1          # eg.：SRR12043508,SRR12043509,SRR18909870,SRR18909871
CONFIG=$2           # RNA-seq_config.sh

# ========= config =========
source "${CONFIG}"

# ========= conda =========
set +u
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}"
set -u

echo "Convert SRA Files to FASTQ Files..."
IFS=',' read -ra SAMPLE_ARRAY <<< "${SAMPLES}"

for i in "${SAMPLE_ARRAY[@]}"; do

    if [[ ! -f "${data_dir}/${i}/${i}.sra" ]]; then
        echo "Error: ${data_dir}/${i}/${i}.sra not found!"
        exit 1
    fi

	mkdir -p "${data_dir}/fastq_data/${i}"
	
	${sratoolkit_dir}/fastq-dump --gzip --split-3 "${data_dir}/${i}/${i}.sra" -O ${data_dir}/fastq_data/${i}/
done



