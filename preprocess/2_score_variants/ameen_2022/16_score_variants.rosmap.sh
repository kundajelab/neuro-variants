#!/bin/bash

set -e
set -u
set -o pipefail

variants=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_lists/rosmap_variants.tsv
model_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/chrombpnet_models/ameen_2022/K562_bias
peaks_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/chrombpnet_inputs/ameen_2022/peaks/filtered
score_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_scores/rosmap/ameen_2022/K562_bias
log_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/logs/variant_scores/rosmap/ameen_2022/K562_bias

hg38_ref_fasta=/oak/stanford/groups/akundaje/soumyak/refs/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
hg38_chrom_sizes=/oak/stanford/groups/akundaje/soumyak/refs/hg38/GRCh38_EBV.chrom.sizes.tsv

jobscript=/home/groups/akundaje/soumyak/neuro-variants/src/score_variants.jobscript.sh
scoring_script=/home/groups/akundaje/soumyak/variant-scorer/src/variant_scoring.per_chrom.py

time=24:00:00
cpus=1
mem=60G
partitions=akundaje,owners

for clust in $model_dir/*; do
    cluster=$(basename $clust)
    echo
    echo $cluster
    echo

    mkdir -p $log_dir/$cluster

    for fld in {0..4}; do
        fold=fold_$fld
        echo
        echo $fold
        echo

        mkdir -p $score_dir/$cluster/$fold

        expected_lines=$(wc -l < $variants)
        if [[ -n $(ls -A $score_dir/$cluster/$fold) ]]; then
            observed_lines=$(cat $score_dir/$cluster/$fold/*.variant_scores.tsv | grep -v variant_id | wc -l)
        else
            observed_lines=0
        fi

        echo Expected Lines: $expected_lines
        echo Observed Lines: $observed_lines

        [[ $expected_lines -eq $observed_lines ]] || \
        sbatch -J $cluster.$fold.rosmap -t $time -c $cpus --mem=$mem \
            -p $partitions --gpus 1 --requeue \
            -o $log_dir/$cluster/$fold.log.txt \
            -e $log_dir/$cluster/$fold.err.txt \
            $jobscript $scoring_script \
                -l $variants \
                -g $hg38_ref_fasta \
                -s $hg38_chrom_sizes \
                -m $model_dir/$cluster/$fold/models/chrombpnet_nobias.h5 \
                -p $peaks_dir/$cluster.filtered.peaks.bed.gz \
                -pg $hg38_ref_fasta \
                -ps $hg38_chrom_sizes \
                -o $score_dir/$cluster/$fold/$cluster.$fold.rosmap \
                -t 1000000 \
                -sc chrombpnet \
                --no_hdf5
    done
done

