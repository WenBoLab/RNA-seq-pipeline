# ========= 基础路径 =========
data_dir="/data/xiejing/liuhao/RNA-seq/03/fastq_data"
fastqc_dir="/data/xiejing/liuhao/RNA-seq/03/results/fastqc"
multiqc_dir="/data/xiejing/liuhao/RNA-seq/03/results/multiqc"

clean_dir="/data/xiejing/liuhao/RNA-seq/03/results/fastp"
clean_fastqc_dir="/data/xiejing/liuhao/RNA-seq/03/fastp/fastqc"
clean_multiqc_dir="/data/xiejing/liuhao/RNA-seq/03/fastp/multiqc"

reference_dir="/data/xiejing/liuhao/RNA-seq/reference/human/hg19"
alignment_dir="/data/xiejing/liuhao/RNA-seq/03/results/map"

# ========= 运行参数 =========
THREADS=8
INSERT_MIN=10
INSERT_MAX=700

# ========= conda =========
CONDA_BASE="/data/xiejing/software/miniconda3"
CONDA_ENV="RNA_env"

