#!/bin/bash

set -e
set -u
set -o pipefail

window=500
conv_tol=0.0005
alpha=0.8

shap_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/peak_shap/domcke_2020/specific_peaks
modisco_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/modisco/domcke_2020/original_peaks
finemo_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/finemo_peaks/domcke_2020/specific_peaks/original_modisco
log_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/logs/finemo_peaks/domcke_2020/specific_peaks/original_modisco

jobscript=/home/groups/akundaje/soumyak/neuro-variants/src/run_finemo.original_modisco.peaks.jobscript.sh

time=720
cpus=4
mem=20G
partitions=akundaje,owners,gpu

for sample in $shap_dir/*; do
    sample=$(basename $sample)
    echo $sample

    for score_type in counts; do
        mean_shap_h5=$shap_dir/$sample/mean/$sample.mean.specific_peaks.${score_type}_scores.${score_type}_scores.h5
        modisco_h5=$modisco_dir/$sample/$sample.original_peaks.counts_scores.w500.modisco.h5
        peak_file=$shap_dir/$sample/mean/$sample.mean.specific_peaks.${score_type}_scores.interpreted_regions.bed
        out_dir=$finemo_dir/$sample/$score_type/alpha_$alpha

        mkdir -p $out_dir
        mkdir -p $log_dir/$sample

        if [[ -f $out_dir/report.html ]]; then
            echo Done
        else
            [[ ! -f $out_dir/hits_unique.tsv ]] || rm -r $out_dir
            mkdir -p $out_dir
            sbatch -J $sample.$alpha \
                -t $time -c $cpus --mem=$mem \
                -p $partitions --gpus 1 --requeue \
                -o $log_dir/$sample/$score_type.alpha_$alpha.log.txt \
                -e $log_dir/$sample/$score_type.alpha_$alpha.err.txt \
            $jobscript \
                $modisco_h5 \
                $mean_shap_h5 \
                $peak_file \
                $alpha \
                $window \
                $conv_tol \
                $out_dir
        fi
    done
done

