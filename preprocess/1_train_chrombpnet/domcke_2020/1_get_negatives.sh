#!/bin/bash

set -e
set -u
set -o pipefail

model_input_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/chrombpnet_inputs/domcke_2020
filtered_peak_dir=$model_input_dir/peaks/filtered
negatives_dir=negatives
log_dir=/users/soumyak/transfer_logs/neuro-variants

hg19_exclude=/oak/stanford/groups/akundaje/soumyak/refs/hg19/hg19_exclusion_regions.slop_1057.bed.gz
hg19_ref_fasta=/oak/stanford/groups/akundaje/soumyak/refs/hg19/male.hg19.fa
hg19_chrom_sizes=/oak/stanford/groups/akundaje/soumyak/refs/hg19/hg19.chrom.sizes
hg19_splits_dir=/oak/stanford/groups/akundaje/soumyak/refs/hg19/splits

for infile in $filtered_peak_dir/*.peaks.bed.gz; do
    sample=$(basename $infile)
    group=$(echo $sample | cut -d '.' -f1)
    dataset=$(echo $sample | cut -d '.' -f2)
    celltype=$(echo $sample | cut -d '.' -f3)
    sample=$group.$dataset.$celltype
    echo $sample

    mkdir -p $negatives_dir/$sample

    for fld in {0..4}
    do
        fold=fold_$fld

        negatives_bed=$negatives_dir/$sample/${fold}_negatives.bed

        [[ -f $negatives_bed ]] || rm -rf $negatives_dir/$sample/${fold}_auxiliary

        [[ -f $negatives_bed ]] || chrombpnet prep nonpeaks -g $hg19_ref_fasta \
                                 -p $infile \
                                 -c $hg19_chrom_sizes \
                                 -fl $hg19_splits_dir/$fold.json \
                                 -br $hg19_exclude \
                                 -o $negatives_dir/$sample/$fold &
    done
    wait

done

echo
echo Copying negatives
echo

scp -r $negatives_dir soumyak@dtn.sherlock.stanford.edu:$model_input_dir/

echo
echo Checking diff after copying
echo

diff -r $negatives_dir $model_input_dir/$negatives_dir

