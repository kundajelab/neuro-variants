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
model=$4
out_prefix=$5
score_type=$6

echo "Live"

chrombpnet contribs_bw \
            -m $model \
            -r $peaks \
            -g $ref_fasta \
            -c $chrom_sizes \
            -op $out_prefix \
            -pc $score_type

