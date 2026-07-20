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

    echo "Processing ${i} ..."

    if [[ "${LibraryLayout}" == "P" ]]; then

        # ==============================
        # Paired-end data
        # ==============================

        if [[ -f "${data_dir}/${i}/${i}_1.fastq.gz" ]]; then

            R1="${data_dir}/${i}/${i}_1.fastq.gz"
            R2="${data_dir}/${i}/${i}_2.fastq.gz"

        elif [[ -f "${data_dir}/${i}/${i}_1.fq.gz" ]]; then

            R1="${data_dir}/${i}/${i}_1.fq.gz"
            R2="${data_dir}/${i}/${i}_2.fq.gz"

        else
            echo "No paired-end FASTQ files found for ${i}"
            continue
        fi


        fastp \
        -i "${R1}" \
        -I "${R2}" \
        -o "${clean_dir}/${i}_1_clean.fq.gz" \
        -O "${clean_dir}/${i}_2_clean.fq.gz" \
        -q 20 \
        -u 30 \
        -n 5 \
        -l 30 \
        --thread "${THREADS}" \
        --html "${clean_dir}/${i}.fastp.html"


    elif [[ "${LibraryLayout}" == "S" ]]; then

        # ==============================
        # Single-end data
        # ==============================

        if [[ -f "${data_dir}/${i}/${i}.fastq.gz" ]]; then

            R1="${data_dir}/${i}/${i}.fastq.gz"

        elif [[ -f "${data_dir}/${i}/${i}.fq.gz" ]]; then

            R1="${data_dir}/${i}/${i}.fq.gz"

        else
            echo "No single-end FASTQ file found for ${i}"
            continue
        fi


        fastp \
        -i "${R1}" \
        -o "${clean_dir}/${i}_clean.fq.gz" \
        -q 20 \
        -u 30 \
        -n 5 \
        -l 30 \
        --thread "${THREADS}" \
        --html "${clean_dir}/${i}.fastp.html"


    else

        echo "Error: LibraryLayout must be P or S"
        exit 1

    fi

done

echo "start running FastQC after clean..."
mkdir -p "${clean_fastqc_dir}" "${clean_multiqc_dir}"
fastqc -t ${THREADS} -o ${clean_fastqc_dir}/ -f fastq ${clean_dir}/*.fq.gz
echo "start running MultiQC after clean..."
multiqc ${clean_fastqc_dir}/ -o ${clean_multiqc_dir}/
