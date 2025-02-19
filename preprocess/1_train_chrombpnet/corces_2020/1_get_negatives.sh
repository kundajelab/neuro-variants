#!/bin/bash

set -e
set -u
set -o pipefail

model_input_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/chrombpnet_inputs/corces_2020
filtered_peak_dir=$model_input_dir/peaks/filtered
negatives_dir=negatives
log_dir=/users/soumyak/transfer_logs/neuro-variants

hg38_exclude=/oak/stanford/groups/akundaje/soumyak/refs/hg38/hg38_exclusion_regions.slop_1057.bed.gz
hg38_ref_fasta=/oak/stanford/groups/akundaje/soumyak/refs/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
hg38_chrom_sizes=/oak/stanford/groups/akundaje/soumyak/refs/hg38/GRCh38_EBV.chrom.sizes.tsv
hg38_splits_dir=/oak/stanford/groups/akundaje/soumyak/refs/hg38/splits

for infile in $filtered_peak_dir/*.peaks.bed.gz; do
    sample=$(basename $infile)
    group=$(echo $sample | cut -d '.' -f1)
    celltype=$(echo $sample | cut -d '.' -f2)
    sample=$group.$celltype
    echo $sample

    # mkdir -p $negatives_dir/$sample

    # for fld in {0..4}
    # do
    #     fold=fold_$fld
    #     chrombpnet prep nonpeaks -g $hg38_ref_fasta \
    #                              -p $infile \
    #                              -c $hg38_chrom_sizes \
    #                              -fl $hg38_splits_dir/$fold.json \
    #                              -br $hg38_exclude \
    #                              -o $negatives_dir/$sample/$fold &
    # done
    # wait

done

echo
echo Copying negatives
echo

scp -r $negatives_dir soumyak@dtn.sherlock.stanford.edu:$model_input_dir/

echo
echo Checking diff after copying
echo

diff -r $negatives_dir $model_input_dir/$negatives_dir

