#!/bin/bash
set -eo pipefail

source /data/shaolizhen/software/miniconda3/etc/profile.d/conda.sh
conda activate RNA_env

data_dir="/data/shaolizhen/crosstalk/NAFLD/data/mouse/GSE165752"
software_dir="/data/shaolizhen/software/sratoolkit.3.2.1-ubuntu64/bin"
for i in 79 80 81 82 87 88 89 90;do
	${software_dir}/fastq-dump --gzip --split-3 ${data_dir}/SRR135756"$i"/SRR135756"$i".sra -O ${data_dir}/fastq_data/SRR135756"$i"/
done



