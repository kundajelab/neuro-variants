#!/bin/bash

module load cuda/11.2
module load cudnn/8.1

# https://stackoverflow.com/questions/34534513/calling-conda-source-activate-from-bash-script
eval "$(conda shell.bash hook)"
conda activate finemo
module load system
module load cairo

modisco_file=$1
mean_shap_file=$2
peak_file=$3
alpha=$4
window=$5
conv_tol=$6
out_dir=$7
include_motifs=$8

echo "Live"

finemo extract-regions-chrombpnet-h5 \
    -c $mean_shap_file \
    -o $out_dir/shap.npz

echo Generated .npz files

finemo call-hits \
    -r $out_dir/shap.npz \
    -m $modisco_file \
    -p $peak_file \
    -b 256 \
    -a $alpha \
    -I $include_motifs \
    -N $include_motifs \
    -o $out_dir

echo Finished running hit-caller

finemo report \
    -r $out_dir/shap.npz \
    -m $modisco_file \
    -p $peak_file \
    -W $window \
    -H $out_dir/hits.tsv \
    -I $include_motifs \
    -N $include_motifs \
    -n \
    -s \
    -o $out_dir

echo Generated Report

