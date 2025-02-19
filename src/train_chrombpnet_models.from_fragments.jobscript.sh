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
splits=$3
peaks=$4
negatives=$5
fragments=$6
bias_model=$7
model_dir=$8

echo "Live"

if [[ -d $model_dir/logs ]]
then
    echo "Found logdir. Deleting previous model"
    rm -rf $model_dir/logs
    rm -rf $model_dir/models
    rm -rf $model_dir/evaluation
    rm -rf $model_dir/auxiliary
fi

mkdir -p $model_dir

chrombpnet pipeline \
    -d "ATAC" \
    -g $ref_fasta \
    -c $chrom_sizes \
    -fl $splits \
    -p $peaks \
    -n $negatives \
    -ifrag $fragments \
    -b $bias_model \
    -o $model_dir

