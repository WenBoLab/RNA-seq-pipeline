# ========= base directory =========
data_dir="/data/RNA-seq/01/fastq_data"
fastqc_dir="/data/RNA-seq/01/results/fastqc"
multiqc_dir="/data/RNA-seq/01/results/multiqc"

clean_dir="/data/RNA-seq/01/results/fastp"
clean_fastqc_dir="/data/RNA-seq/01/results/fastp/fastqc"
clean_multiqc_dir="/data/RNA-seq/01/results/fastp/multiqc"

reference_dir="/data/RNA-seq/reference/human/hg19"
GTF_dir="/data/RNA-seq/reference/human/hg19"
alignment_dir="/data/RNA-seq/01/results/map"

counts_dir="/data/RNA-seq/01/results/counts"

# ========= parameters =========
THREADS=8
INSERT_MIN=10
INSERT_MAX=700

# ========= conda =========
CONDA_BASE="/data/software/miniconda3"
CONDA_ENV="RNA_env"

