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

featureCounts -T ${THREADS} -p -a ${GTF_dir}/gencode.v19.annotation.gtf \
-o ${counts_dir}/gene_counts.txt \
-Q 10 \
${alignment_dir}/*.Aligned.sortedByCoord.out.bam

