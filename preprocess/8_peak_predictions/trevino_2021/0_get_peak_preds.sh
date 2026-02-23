#!/bin/bash

set -e
set -u
set -o pipefail

bias_name=K562
dataset=trevino_2021
processed_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed
peak_dir=$processed_dir/$dataset/peaks/overlap
model_input_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/chrombpnet_inputs/$dataset
negatives_dir=$model_input_dir/negatives
model_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/chrombpnet_models/$dataset/${bias_name}_bias

pred_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/peak_preds/$dataset/${bias_name}_bias/overlap_peaks
log_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/logs/peak_preds/$dataset/${bias_name}_bias/overlap_peaks

hg38_ref_fasta=/oak/stanford/groups/akundaje/soumyak/refs/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
hg38_chrom_sizes=/oak/stanford/groups/akundaje/soumyak/refs/hg38/GRCh38_EBV.chrom.sizes.tsv

jobscript=/home/groups/akundaje/soumyak/flare-revision/src/get_preds.jobscript.sh

time=2:00:00
cpus=1
mem=40G
partitions=akundaje,owners

for sample in $negatives_dir/*; do
    sample=$(basename $sample)
    echo $sample

    mkdir -p $log_dir/$sample

    for fld in {0..4}; do
        fold=fold_$fld
        echo $fold

        mkdir -p $pred_dir/$sample/$fold

        bias_bw_file=$pred_dir/$sample/$fold/$sample.$fold.overlap_peaks_bias.bw

        [[ -f $bias_bw_file ]] || sbatch -J $sample.$fold -t $time -c $cpus --mem=$mem \
            -p $partitions --gpus 1 --requeue \
            -o $log_dir/$sample/$fold.log.txt \
            -e $log_dir/$sample/$fold.err.txt \
            $jobscript \
                $hg38_ref_fasta \
                $hg38_chrom_sizes \
                $peak_dir/$sample.overlap.peaks.bed.gz \
                $model_dir/$sample/$fold/auxiliary/data_unstranded.bw \
                $model_dir/$sample/$fold/models/bias_model_scaled.h5 \
                $model_dir/$sample/$fold/models/chrombpnet.h5 \
                $model_dir/$sample/$fold/models/chrombpnet_nobias.h5 \
                $pred_dir/$sample/$fold/$sample.$fold.overlap_peaks
    done
done

