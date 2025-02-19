#!/bin/bash

annotate_script=/home/groups/akundaje/soumyak/variant-scorer/src/variant_annotation.py
# tss_gtex_file=/oak/stanford/groups/akundaje/soumyak/refs/gtex/hg38.gtex.protein_coding.tss.bed
# tss_gencode_file=/oak/stanford/groups/akundaje/soumyak/refs/gencode/hg38/hg38.gencode.all.tss.bed
tss_gencode_coding_file=/oak/stanford/groups/akundaje/soumyak/refs/gencode/hg38/hg38.gencode.protein_coding.tss.bed

summary_dir=$1
annotated_dir=$2
peaks=$3
sample=$4

mkdir -p $annotated_dir

out_prefix=$annotated_dir/$sample

python $annotate_script \
    -l $summary_dir/$sample.mean.variant_scores.tsv \
    -o $out_prefix \
    -p $peaks \
    -g $tss_gencode_coding_file \
    -sc chrombpnet

echo "Done Annotating"
echo

# bgzip $out_prefix.annotations.tsv

