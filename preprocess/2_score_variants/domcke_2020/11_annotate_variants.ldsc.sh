#!/bin/bash

set -e
set -u
set -o pipefail

variants=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_lists/ldsc.filtered.variants.hg38.tsv
score_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_scores/ldsc/domcke_2020/K562_bias
merged_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/merged_variant_scores/ldsc/domcke_2020/K562_bias
summary_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_summary/ldsc/domcke_2020/K562_bias

processed_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/domcke_2020
overlap_peak_dir=$processed_dir/peaks/hg38

annotated_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_annotations/ldsc/domcke_2020/K562_bias
log_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/logs/variant_annotations/ldsc/domcke_2020/K562_bias

jobscript=/home/groups/akundaje/soumyak/neuro-variants/src/annotate_variants.jobscript.sh

time=120
cpus=1
mem=30G
partitions=akundaje,owners

for clust in $score_dir/*; do
    cluster=$(basename $clust)
    echo
    echo $cluster
    echo

    mkdir -p $annotated_dir
    mkdir -p $log_dir

    summary_file=$summary_dir/$cluster.mean.variant_scores.tsv

    expected_lines=$(wc -l < $variants)
    if [[ -f $summary_file ]]; then
        observed_lines=$(cat $summary_file | grep -v variant_id | wc -l)
    else
        observed_lines=0
    fi

    echo Expected Lines: $expected_lines
    echo Observed Lines: $observed_lines

    annotated_file=$annotated_dir/$cluster.annotations.tsv

    if [[ $expected_lines -eq $observed_lines ]]; then
        [[ -f $annotated_file ]] || \
        sbatch -J $cluster.ldsc -t $time -c $cpus --mem=$mem \
            -p $partitions --requeue \
            -o $log_dir/$cluster.log \
            -e $log_dir/$cluster.err \
            $jobscript \
                $summary_dir \
                $annotated_dir \
                $overlap_peak_dir/$cluster.overlap.peaks.hg38.bed \
                $cluster
    fi
done

