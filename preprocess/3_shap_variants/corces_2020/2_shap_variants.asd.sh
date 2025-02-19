#!/bin/bash

set -e
set -u
set -o pipefail

variants=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_lists/neuro.asd.variants.shap.tsv
model_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/chrombpnet_models/corces_2020/K562_bias
shap_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_shap/asd/K562_bias/corces_2020
log_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/logs/variant_shap/asd/K562_bias/corces_2020

hg38_ref_fasta=/oak/stanford/groups/akundaje/soumyak/refs/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
hg38_chrom_sizes=/oak/stanford/groups/akundaje/soumyak/refs/hg38/GRCh38_EBV.chrom.sizes.tsv

jobscript=/home/groups/akundaje/soumyak/neuro-variants/src/score_variants.jobscript.sh
scoring_script=/home/groups/akundaje/soumyak/variant-scorer/src/variant_scoring.py
shap_script=/home/groups/akundaje/soumyak/variant-scorer/src/variant_shap.py

time=30
cpus=1
mem=10G
partitions=akundaje,owners

# for clust in $model_dir/*; do
#     cluster=$(basename $clust)

for cluster in corces_2020.Cluster1 corces_2020.Cluster8 corces_2020.Cluster11 corces_2020.Cluster13 corces_2020.Cluster19 corces_2020.Cluster24; do
    echo
    echo $cluster
    echo

    mkdir -p $log_dir/$cluster

    for fld in {0..4}; do
        fold=fold_$fld
        echo
        echo $fold
        echo

        mkdir -p $shap_dir/$cluster/$fold

        score_output_file=$shap_dir/$cluster/$fold/$cluster.$fold.asd.shap.variant_scores.tsv
        shap_output_file=$shap_dir/$cluster/$fold/$cluster.$fold.asd.shap.variant_shap.counts.h5

        expected_lines=$(wc -l < $variants)
        if [[ -n $(ls -A $shap_dir/$cluster/$fold) ]]; then
            observed_lines=$(cat $shap_dir/$cluster/$fold/*.variant_scores.tsv | grep -v variant_id | wc -l)
        else
            observed_lines=0
        fi

        echo Expected Lines: $expected_lines
        echo Observed Lines: $observed_lines

        [[ $expected_lines -eq $observed_lines ]] || \
        sbatch -J $cluster.$fold.asd.shap.variant_scores -t $time -c $cpus --mem=$mem \
            -p $partitions --gpus 1 --requeue \
            -o $log_dir/$cluster/$fold.variant_scores.log.txt \
            -e $log_dir/$cluster/$fold.variant_scores.err.txt \
            $jobscript $scoring_script \
                -l $variants \
                -g $hg38_ref_fasta \
                -s $hg38_chrom_sizes \
                -m $model_dir/$cluster/$fold/models/chrombpnet_nobias.h5 \
                -o $shap_dir/$cluster/$fold/$cluster.$fold.asd.shap \
                -t 0 \
                -sc chrombpnet

        [[ -f $shap_output_file ]] || \
        sbatch -J $cluster.$fold.asd.shap.variant_shap -t $time -c $cpus --mem=$mem \
            -p $partitions --gpus 1 --requeue \
            -o $log_dir/$cluster/$fold.variant_shap.log.txt \
            -e $log_dir/$cluster/$fold.variant_shap.err.txt \
            $jobscript $shap_script \
                -l $variants \
                -g $hg38_ref_fasta \
                -s $hg38_chrom_sizes \
                -m $model_dir/$cluster/$fold/models/chrombpnet_nobias.h5 \
                -o $shap_dir/$cluster/$fold/$cluster.$fold.asd.shap \
                -sc chrombpnet

    done
done

