#!/bin/bash

set -e
set -u
set -o pipefail

processed_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/trevino_2021
croo_dir=$processed_dir/croo
qc_dir=$processed_dir/qc_reports

mkdir -p $qc_dir

for sample in $croo_dir/*; do
    sample=$(basename $sample)
    echo
    echo $sample
    echo

    qc_html_file=$croo_dir/$sample/qc/qc.html
    echo $qc_html_file

    ln -s $qc_html_file $qc_dir/$sample.qc.html

done

