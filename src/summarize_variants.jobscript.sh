#!/bin/bash

summary_script=/home/groups/akundaje/soumyak/variant-scorer/src/variant_summary_across_folds.py

score_dir=$1
merged_dir=$2
summary_dir=$3
sample=$4
num_folds=$5

mkdir -p $merged_dir
mkdir -p $summary_dir

out_prefix=$summary_dir/$sample

for fld in $(seq 0 $((num_folds - 1)))
do
    fold=fold_$fld
    echo $fold
    cat $score_dir/$sample/$fold/*.variant_scores.tsv | head -n 1 > $merged_dir/$sample/$fold.variant_scores.tsv
    cat $score_dir/$sample/$fold/*.variant_scores.tsv | grep -v "allele1" >> $merged_dir/$sample/$fold.variant_scores.tsv
    rm -f $merged_dir/$sample/$fold.variant_scores.tsv.gz
    bgzip $merged_dir/$sample/$fold.variant_scores.tsv
done

echo "Done Concatenating"

sl_args=$(seq -f "fold_%g.variant_scores.tsv.gz" 0 $((num_folds - 1)) | xargs)

python $summary_script \
    -sd $merged_dir/$sample \
    -sl $sl_args \
    -o $out_prefix \
    -sc chrombpnet

echo "Done Summarizing"
echo

# bgzip $out_prefix.mean.variant_scores.tsv

