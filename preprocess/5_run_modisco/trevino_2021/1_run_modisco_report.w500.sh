#!/bin/bash

set -e
set -u
set -o pipefail

window=500
num_seqlets=1000000
modisco_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/modisco/trevino_2021
modisco_report_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/modisco_reports/trevino_2021
meme_file=/oak/stanford/groups/akundaje/soumyak/motifs/latest/hocomoco_v12/H12CORE_meme_format.meme

for peak_type in original; do # specific; do
    echo ${peak_type}_peaks

    for sample in $modisco_dir/${peak_type}_peaks/*; do
        sample=$(basename $sample)
        echo $sample

        mkdir -p $modisco_report_dir/${peak_type}_peaks/$sample

        for score_type in counts; do
            modisco_h5=$modisco_dir/${peak_type}_peaks/$sample/$sample.${peak_type}_peaks.${score_type}_scores.w${window}.modisco.h5
            out_dir=$modisco_report_dir/${peak_type}_peaks/$sample/$score_type
            modisco_html=$modisco_report_dir/${peak_type}_peaks/$sample/$score_type/motifs.html

            [[ -f $modisco_html ]] || \
            if [[ -f $modisco_h5 ]]; then
                modisco report \
                    -i $modisco_h5 \
                    -o $out_dir \
                    -n 10 \
                    -m $meme_file
            fi
        done
    done
done

