bash 1.data_quality_control.sh \
"SRR12043508,SRR12043509,SRR18909870,SRR18909871" \
RNA-seq_config.sh

# Usage:
#   bash <program_file> <sample_ids> <config_file>

# Arguments:
#   program_file   The executable shell script to be run (e.g. 1.data_quality_control.sh, 2.read_trimming_and_filtering.sh)
#   sample_ids   A comma-separated list of SRA accession IDs (e.g. SRR12043508,SRR12043509)
#   config_file  Configuration file specifying paths and parameters for RNA-seq analysis
