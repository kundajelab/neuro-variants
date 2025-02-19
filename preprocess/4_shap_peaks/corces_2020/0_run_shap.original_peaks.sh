#!/bin/bash

set -e
set -u
set -o pipefail

bias_name=K562
processed_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/corces_2020
overlap_peak_dir=$processed_dir/peaks/overlap
model_input_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/chrombpnet_inputs/corces_2020
negatives_dir=$model_input_dir/negatives
model_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/chrombpnet_models/corces_2020/${bias_name}_bias

shap_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/peak_shap/corces_2020/original_peaks
log_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/logs/peak_shap/corces_2020/original_peaks

hg38_ref_fasta=/oak/stanford/groups/akundaje/soumyak/refs/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
hg38_chrom_sizes=/oak/stanford/groups/akundaje/soumyak/refs/hg38/GRCh38_EBV.chrom.sizes.tsv

jobscript=/home/groups/akundaje/soumyak/neuro-variants/src/run_shap.jobscript.sh

time=24:00:00
cpus=1
mem=60G
partitions=akundaje,owners

for sample in $negatives_dir/*; do
    sample=$(basename $sample)
    echo $sample

    mkdir -p $shap_dir/$sample
    mkdir -p $log_dir/$sample

    for fld in {0..4}; do
        fold=fold_$fld
        echo $fold

        mkdir -p $shap_dir/$sample/$fold

        for score_type in counts; do
            bw_outfile=$shap_dir/$sample/$fold/$sample.$fold.original_peaks.${score_type}_scores.${score_type}_scores.bw

            [[ -f $bw_outfile ]] || \
            sbatch -J $sample.$fold -t $time -c $cpus --mem=$mem \
                -p $partitions --gpus 1 --requeue \
                -C GPU_SKU:A100_SXM4 \
                -o $log_dir/$sample/$fold.log.txt \
                -e $log_dir/$sample/$fold.err.txt \
                $jobscript \
                    $hg38_ref_fasta \
                    $hg38_chrom_sizes \
                    $overlap_peak_dir/$sample.overlap.peaks.bed.gz \
                    $model_dir/$sample/$fold/models/chrombpnet_nobias.h5 \
                    $shap_dir/$sample/$fold/$sample.$fold.original_peaks.${score_type}_scores \
                    $score_type
        done
    done
done

