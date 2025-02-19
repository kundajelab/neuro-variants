#!/bin/bash

set -e
set -u
set -o pipefail

processed_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/domcke_2020
tagalign_dir=$processed_dir/tagaligns
model_input_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/chrombpnet_inputs/domcke_2020
filtered_peak_dir=$model_input_dir/peaks/filtered
negatives_dir=$model_input_dir/negatives

bias_name=K562
bias_model=/oak/stanford/groups/akundaje/soumyak/bias_models/ENCSR868FGK_bias_fold_0.h5
model_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/chrombpnet_models/domcke_2020/${bias_name}_bias
log_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/logs/chrombpnet_models/domcke_2020/${bias_name}_bias

hg19_ref_fasta=/oak/stanford/groups/akundaje/soumyak/refs/hg19/male.hg19.fa
hg19_chrom_sizes=/oak/stanford/groups/akundaje/soumyak/refs/hg19/hg19.chrom.sizes
hg19_splits_dir=/oak/stanford/groups/akundaje/soumyak/refs/hg19/splits

jobscript=/home/groups/akundaje/soumyak/neuro-variants/src/train_chrombpnet_models.from_tagaligns.jobscript.sh

time=24:00:00
cpus=2
mem=60G
partitions=akundaje,owners

for sample in $negatives_dir/*; do
    sample=$(basename $sample)
    echo $sample

    mkdir -p $log_dir/$sample

    for fld in {0..4}; do
        fold=fold_$fld
        echo $fold

        mkdir -p $model_dir/$sample/$fold

        report=$model_dir/$sample/$fold/evaluation/overall_report.html

        [[ -f $report ]] || \
        sbatch -J $sample.$fold -t $time -c $cpus --mem=$mem \
            -p $partitions --gpus 1 --requeue \
            -C GPU_SKU:A100_SXM4 \
            -o $log_dir/$sample/$fold.log.txt \
            -e $log_dir/$sample/$fold.err.txt \
            $jobscript \
                $hg19_ref_fasta \
                $hg19_chrom_sizes \
                $hg19_splits_dir/$fold.json \
                $filtered_peak_dir/$sample.filtered.peaks.bed.gz \
                $negatives_dir/$sample/${fold}_negatives.bed \
                $tagalign_dir/$sample.tagAlign.gz \
                $bias_model \
                $model_dir/$sample/$fold
    done
done

