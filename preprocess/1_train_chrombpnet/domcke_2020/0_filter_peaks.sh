#!/bin/bash

set -e
set -u
set -o pipefail

processed_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/domcke_2020
overlap_peak_dir=$processed_dir/peaks/overlap
model_input_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/chrombpnet_inputs/domcke_2020
filtered_peak_dir=$model_input_dir/peaks/filtered

hg19_exclude=/oak/stanford/groups/akundaje/soumyak/refs/hg19/hg19_exclusion_regions.slop_1057.bed.gz

mkdir -p $filtered_peak_dir

for infile in $overlap_peak_dir/*.peaks.bed.gz; do
    sample=$(basename $infile)
    group=$(echo $sample | cut -d '.' -f1)
    dataset=$(echo $sample | cut -d '.' -f2)
    celltype=$(echo $sample | cut -d '.' -f3)
    sample=$group.$dataset.$celltype
    echo $sample

    bedtools intersect -v -a $infile -b $hg19_exclude | bedtools sort -i stdin | bgzip > $filtered_peak_dir/$sample.filtered.peaks.bed.gz
    tabix -p bed $filtered_peak_dir/$sample.filtered.peaks.bed.gz

done

