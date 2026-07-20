# ========= base directory =========
base_dir="/data/NAFLD/results/mouse/GSE176681"
sra_dir="/data/NAFLD/data/mouse/GSE176681/sra"
data_dir="/data/NAFLD/data/mouse/GSE176681/fastq_data"

# ========= other directory =========
fastqc_dir="${base_dir}/fastqc"
multiqc_dir="${base_dir}/multiqc"

clean_dir="${base_dir}/fastp"
clean_fastqc_dir="${base_dir}/fastp/fastqc"
clean_multiqc_dir="${base_dir}/fastp/multiqc"

mouse_genome_dir="/data/reference/mouse"
mouse_gtf_dir="/data/reference/mouse"

human_genome_dir="/data/reference/human"
human_gtf_dir="/data/reference/human"

alignment_dir="${base_dir}/map"

counts_dir="${base_dir}/counts"

sratoolkit_dir="/data/software/sratoolkit.3.2.1-ubuntu64/bin"

# ========= parameters =========
LibraryLayout="P"  # P: paired-end; S: single-end
SPECIES="mouse"  # mouse or human

THREADS=8
INSERT_MIN=10
INSERT_MAX=700

# ========= conda =========
CONDA_BASE="/data/software/miniconda3"
CONDA_ENV="RNA_env"

