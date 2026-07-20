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

if [[ "${SPECIES}" == "mouse" ]]; then

    echo "Species: Mouse"
    reference_dir="${mouse_reference_dir}"
    GTF_dir="${mouse_GTF_dir}"
    genome_name="GRCm39.genome"
    fa_gz="${reference_dir}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
    fa="${reference_dir}/Mus_musculus.GRCm39.dna.primary_assembly.fa"
    gtf_gz="${GTF_dir}/Mus_musculus.GRCm39.114.chr.gtf.gz"
    gtf="${GTF_dir}/Mus_musculus.GRCm39.114.chr.gtf"

elif [[ "${SPECIES}" == "human" ]]; then

    echo "Species: Human"
    reference_dir="${human_reference_dir}"
    GTF_dir="${human_GTF_dir}"
    genome_name="GRCh38.genome"
    fa_gz="${reference_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    fa="${reference_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    gtf_gz="${GTF_dir}/Homo_sapiens.GRCh38.114.gtf.gz"
    gtf="${GTF_dir}/Homo_sapiens.GRCh38.114.gtf"

else

    echo "ERROR: SPECIES must be mouse or human"
    exit 1
fi

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
index_dir="${reference_dir}/${genome_name}_index"
# parameter -d: if the index_dir path exists
if [[ ! -d "$index_dir" ]]; then
    echo "create STAR index..."
    mkdir -p "${index_dir}"
    STAR --runMode genomeGenerate \
    --runThreadN "${THREADS}" \
    --genomeDir "${index_dir}" \
    --genomeFastaFiles "${fa}" \
    --sjdbGTFfile "${gtf}" \
    --sjdbOverhang 149
else
    echo "index already exist, skip."
fi

for i in "${SAMPLE_ARRAY[@]}"; do
    echo ">>> Mapping sample: ${i}.clean.fq.gz"
    if [[ "${LibraryLayout}" == "P" ]]; then
        echo "Paired-end data"
        
        STAR runThreadN "${THREADS}" \
        --genomeDir "${index_dir}" \
        --readFilesIn "${clean_dir}/${i}_1_clean.fq.gz" "${clean_dir}/${i}_2_clean.fq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${alignment_dir}/${i}." \
        --twopassMode Basic \
        --sjdbGTFfile "${gtf}" \
        --quantMode GeneCounts \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM XS \
        --outFilterMatchNminOverLread 0.3 \
        --outFilterScoreMinOverLread 0.3 \
        --outSAMattrRGline ID:${i} SM:${i} PL:ILLUMINA

    elif [[ "${LibraryLayout}" == "S" ]]; then
    
        echo "Single-end data"

        STAR --runThreadN "${THREADS}" \
        --genomeDir "${index_dir}" \
        --readFilesIn "${clean_dir}/${i}_clean.fq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${alignment_dir}/${i}." \
        --twopassMode Basic \
        --sjdbGTFfile "${gtf}" \
        --quantMode GeneCounts \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM XS \
        --outFilterMatchNminOverLread 0.3 \
        --outFilterScoreMinOverLread 0.3 \
        --outSAMattrRGline ID:${i} SM:${i} PL:ILLUMINA

    else
        echo "ERROR: LibraryLayout must be P or S"
        exit 1
    fi

  # create BAM index
  samtools index "${alignment_dir}/${i}.Aligned.sortedByCoord.out.bam"
  mv "${alignment_dir}/${i}.Log.final.out" "${alignment_dir}/${i}.Log.final.out.tsv"
done
