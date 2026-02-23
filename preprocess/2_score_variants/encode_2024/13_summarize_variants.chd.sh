#!/bin/bash

set -e
set -u
set -o pipefail

variants=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_lists/chd_snv_list.tsv
score_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_scores/chd/encode_2024/K562_bias
merged_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/merged_variant_scores/chd/encode_2024/K562_bias
summary_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_summary/chd/encode_2024/K562_bias
log_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/logs/variant_summary/chd/encode_2024/K562_bias

jobscript=/home/groups/akundaje/soumyak/neuro-variants/src/summarize_variants.jobscript.sh

time=60
cpus=1
mem=30G
partitions=akundaje,owners
num_folds=5
latest_fold=fold_$((num_folds - 1))

for clust in $score_dir/*; do
    cluster=$(basename $clust)
    echo
    echo $cluster
    echo

    expected_lines=$(wc -l < $variants)
    if [[ -n $(ls -A $score_dir/$cluster/$latest_fold) ]]; then
        observed_lines=$(cat $score_dir/$cluster/$latest_fold/*.variant_scores.tsv | grep -v variant_id | wc -l)
    else
        observed_lines=0
    fi

    echo Expected Lines: $expected_lines
    echo Observed Lines: $observed_lines

    summary_file=$summary_dir/$cluster.mean.variant_scores.tsv

    if [[ $expected_lines -eq $observed_lines ]]; then

        mkdir -p $merged_dir/$cluster
        mkdir -p $summary_dir
        mkdir -p $log_dir

        [[ -f $summary_file ]] || \
        sbatch -J $cluster.rare -t $time -c $cpus --mem=$mem \
            -p $partitions --requeue \
            -o $log_dir/$cluster.log \
            -e $log_dir/$cluster.err \
            $jobscript \
                $score_dir \
                $merged_dir \
                $summary_dir \
                $cluster \
                $num_folds
    fi
done

