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
mkdir -p "${counts_dir}"

echo "start running featureCounts..."

# Select GTF according to species
if [[ "${SPECIES}" == "mouse" ]]; then
    echo "Species: Mouse"
    annotation="${mouse_gtf_dir}/Mus_musculus.GRCm39.114.chr.gtf"
elif [[ "${SPECIES}" == "human" ]]; then
    echo "Species: Human"
    annotation="${human_gtf_dir}/Homo_sapiens.GRCh38.114.gtf"
else
    echo "ERROR: SPECIES must be mouse or human"
    exit 1
fi

# Library layout
if [[ "${LibraryLayout}" == "P" ]]; then
    paired="-p"
elif [[ "${LibraryLayout}" == "S" ]]; then
    paired=""
else
    echo "ERROR: LibraryLayout must be P or S"
    exit 1
fi

featureCounts -T ${THREADS} ${paired} -a ${annotation} \
-o ${counts_dir}/gene_counts.txt \
-Q 10 \
${alignment_dir}/*.Aligned.sortedByCoord.out.bam

