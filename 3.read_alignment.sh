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
mkdir -p "${alignment_dir}"

echo "start running STAR..."

IFS=',' read -ra SAMPLE_ARRAY <<< "${SAMPLES}"

fa_gz="${reference_dir}/hg19.fa.gz"
fa="${reference_dir}/hg19.fa"

gtf_gz="${GTF_dir}/gencode.v19.annotation.gtf.gz"
gtf="${GTF_dir}/gencode.v19.annotation.gtf"

# parameter -f: if fa_gz file path exists
if [[ -f "$fa_gz" && ! -f "$fa" ]]; then
    echo "Unzipping genome fasta..."
    gunzip -k "$fa_gz"
else
    echo "Genome fasta already unzipped, skip."
fi

if [[ -f "$gtf_gz" && ! -f "$gtf" ]]; then
    echo "Unzipping GTF..."
    gunzip -k "$gtf_gz"
else
    echo "GTF already unzipped, skip."
fi

# create index
index_dir="${reference_dir}/hg19_Index"
# parameter -d: if the index_dir path exists
if [[ ! -d "$index_dir" ]]; then
    echo "create STAR index..."
    mkdir -p "${index_dir}"
    STAR --runMode genomeGenerate \
    --runThreadN ${THREADS} \
    --genomeDir ${index_dir} \
    --genomeFastaFiles ${reference_dir}/hg19.fa \
    --sjdbGTFfile "${GTF_dir}/gencode.v19.annotation.gtf" \
    --sjdbOverhang 149
else
    echo "index already exist, skip."
fi


for i in "${SAMPLE_ARRAY[@]}"; do
  echo ">>> Mapping sample: ${i}.clean.fq.gz"

  STAR --runThreadN ${THREADS} \
       --genomeDir ${index_dir} \
       --readFilesIn "${clean_dir}/${i}_1_clean.fq.gz" "${clean_dir}/${i}_2_clean.fq.gz" \
       --readFilesCommand zcat \
       --outFileNamePrefix "${alignment_dir}/${i}." \
       --twopassMode Basic \
       --sjdbGTFfile "${GTF_dir}/gencode.v19.annotation.gtf" \
       --quantMode  GeneCounts \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattributes NH HI AS nM XS \
       --outFilterMatchNminOverLread 0.3 \
       --outFilterScoreMinOverLread 0.3 \
       --outSAMattrRGline ID:$i SM:$i PL:ILLUMINA

  # create BAM index
  samtools index "${alignment_dir}/${i}.Aligned.sortedByCoord.out.bam"
  mv "${alignment_dir}/${i}.Log.final.out" "${alignment_dir}/${i}.Log.final.out.tsv"
#done
