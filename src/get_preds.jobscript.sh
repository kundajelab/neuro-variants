#!/bin/bash

module load cuda/11.2
module load cudnn/8.1

# https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
eval "$(conda shell.bash hook)"
conda activate chrombpnet
module load system
module load cairo

ref_fasta=$1
chrom_sizes=$2
peaks=$3
bigwig=$4
bias_model=$5
chrombpnet_model=$6
nobias_model=$7
out_prefix=$8

echo "Live"

chrombpnet pred_bw \
    -bm $bias_model \
    -cm $chrombpnet_model \
    -cmb $nobias_model \
    -r $peaks \
    -bw $bigwig \
    -g $ref_fasta \
    -c $chrom_sizes \
    -op $out_prefix

