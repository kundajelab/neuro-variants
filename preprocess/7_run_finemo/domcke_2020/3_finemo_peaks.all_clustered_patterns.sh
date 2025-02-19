#!/bin/bash

set -e
set -u
set -o pipefail

window=500
conv_tol=0.0005
alpha=0.7

modisco_h5=/oak/stanford/groups/akundaje/projects/neuro-variants/motif_compendium/all_data/leiden_96/neuro-variants.all_data.motif_compendium.avg.leiden_96.h5
shap_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/peak_shap/domcke_2020/specific_peaks
finemo_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/finemo_peaks/domcke_2020/specific_peaks/all_data_all_patterns/leiden_96
log_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/logs/finemo_peaks/domcke_2020/specific_peaks/all_data_all_patterns/leiden_96

jobscript=/home/groups/akundaje/soumyak/neuro-variants/src/run_finemo.all_patterns.peaks.jobscript.sh

time=1-0
cpus=4
mem=40G
partitions=akundaje,owners

for sample in $shap_dir/*; do
    sample=$(basename $sample)
    echo $sample

    for score_type in counts; do
        mean_shap_h5=$shap_dir/$sample/mean/$sample.mean.specific_peaks.${score_type}_scores.${score_type}_scores.h5
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

