#!/bin/bash

set -e
set -u
set -o pipefail

window=500
conv_tol=0.0005
alpha=0.8

variants=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_lists/shap_variants.rare.model_inputs.hg38.tsv
shap_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_shap/rare/K562_bias/domcke_2020
modisco_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/modisco/domcke_2020/original_peaks

log_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/logs/finemo_variants/rare/K562_bias/domcke_2020
finemo_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/finemo_variants/rare/K562_bias/domcke_2020

jobscript=/home/groups/akundaje/soumyak/neuro-variants/src/run_finemo.original_modisco.jobscript.sh
finemo_script=/home/groups/akundaje/soumyak/variant-scorer/src/hitcaller_variant.py

time=20
cpus=1
mem=10G
partitions=akundaje,owners

for clust in $shap_dir/*; do
    cluster=$(basename $clust)
    echo
    echo $cluster
    echo

    for score_type in counts; do
        shap_data=$shap_dir/$cluster/mean/$cluster.mean.rare.shap.variant_shap.counts.h5
        modisco_h5=$modisco_dir/$cluster/$cluster.original_peaks.counts_scores.w500.modisco.h5
        out_dir=$finemo_dir/original_modisco/$cluster/mean/$score_type/alpha_$alpha
        hit_calls=$out_dir/variant_hit_calls.tsv

        mkdir -p $log_dir/original_modisco/$cluster/mean
        mkdir -p $out_dir

        [[ -f $hit_calls ]] || sbatch -J $cluster.mean.rare.finemo.variants \
            -t $time -c $cpus --mem=$mem \
            -p $partitions --gpus 1 --requeue \
            -o $log_dir/original_modisco/$cluster/mean/$score_type.alpha_$alpha.log.txt \
            -e $log_dir/original_modisco/$cluster/mean/$score_type.alpha_$alpha.err.txt \
            $jobscript $finemo_script \
                --shap_data $shap_data \
                --input_type h5 \
                --modisco_h5 $modisco_h5 \
                --variant_file $variants \
                --output_dir $out_dir \
                --alpha $alpha
    done
done

