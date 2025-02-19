#!/bin/bash

set -e
set -u
set -o pipefail

processed_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/corces_2020
croo_dir=$processed_dir/croo
overlap_peak_dir=$processed_dir/peaks/overlap
idr_peak_dir=$processed_dir/peaks/idr
bigwig_dir=$processed_dir/bigwigs

mkdir -p $overlap_peak_dir
mkdir -p $idr_peak_dir
mkdir -p $bigwig_dir

for sample in $croo_dir/*; do
    sample=$(basename $sample)
    echo
    echo $sample
    echo

    zcat $croo_dir/$sample/peak/overlap_reproducibility/overlap.optimal_peak.narrowPeak.gz | bedtools sort -i stdin | bgzip > $overlap_peak_dir/$sample.overlap.peaks.bed.gz
    tabix -p bed $overlap_peak_dir/$sample.overlap.peaks.bed.gz

    zcat $croo_dir/$sample/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz | bedtools sort -i stdin | bgzip > $idr_peak_dir/$sample.idr.peaks.bed.gz
    tabix -p bed $idr_peak_dir/$sample.idr.peaks.bed.gz

    fc_bigwig_file=$(find $croo_dir/$sample/signal -maxdepth 2 -type f -name *.fc.signal.bigwig)
    echo $fc_bigwig_file

    pval_bigwig_file=$(find $croo_dir/$sample/signal -maxdepth 2 -type f -name *.pval.signal.bigwig)
    echo $pval_bigwig_file

    ln -s $fc_bigwig_file $bigwig_dir/$sample.fc.bigwig
    ln -s $pval_bigwig_file $bigwig_dir/$sample.pval.bigwig

done

