#!/bin/bash

module load cuda/11.2
module load cudnn/8.1

# https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
eval "$(conda shell.bash hook)"
conda activate chrombpnet
module load system
module load cairo

onehot=$1
mean_shap=$2
num_seqlets=$3
window=$4
output_file=$5

echo "Live"

modisco motifs \
    -s $onehot \
    -a $mean_shap \
    -n $num_seqlets \
    -w $window \
    -o $output_file \
    -v

