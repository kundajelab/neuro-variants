#!/bin/bash

set -e
set -u
set -o pipefail

shap_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/peak_shap/trevino_2021

window=500
num_seqlets=1000000
modisco_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/modisco/trevino_2021
log_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/logs/modisco/trevino_2021

jobscript=/home/groups/akundaje/soumyak/neuro-variants/src/run_modisco.jobscript.sh

time=48:00:00
cpus=24
mem=60G
partitions=akundaje,owners

for peak_type in original; do # specific; do
    echo ${peak_type}_peaks

    for sample in $shap_dir/${peak_type}_peaks/*; do
        sample=$(basename $sample)
        echo $sample

        mkdir -p $modisco_dir/${peak_type}_peaks/$sample
        mkdir -p $log_dir/${peak_type}_peaks/$sample

        for score_type in counts; do
            mean_shap_npz=$shap_dir/${peak_type}_peaks/$sample/mean/$sample.mean.${peak_type}_peaks.${score_type}_scores.${score_type}_scores.npz
            onehot_npz=$shap_dir/${peak_type}_peaks/$sample/mean/$sample.mean.${peak_type}_peaks.${score_type}_scores.onehot.npz

            modisco_h5=$modisco_dir/${peak_type}_peaks/$sample/$sample.${peak_type}_peaks.${score_type}_scores.w${window}.modisco.h5

            [[ -f $modisco_h5 ]] || \
            sbatch -J $sample.$peak_type.$score_type.modisco -t $time -c $cpus --mem=$mem \
                -p $partitions --requeue \
                -o $log_dir/${peak_type}_peaks/$sample/${score_type}_scores.w${window}.modisco.log.txt \
                -e $log_dir/${peak_type}_peaks/$sample/${score_type}_scores.w${window}.modisco.err.txt \
                $jobscript \
                    $onehot_npz \
                    $mean_shap_npz \
                    $num_seqlets \
                    $window \
                    $modisco_h5
        done
    done
done

