#!/bin/bash

set -e
set -u
set -o pipefail

orig_dir=/mnt/lab_data3/salil512/adult_heart/2_chrombpnet/2_celltype_fragments

raw_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/raw/encode_2024/fragments
log_dir=/users/soumyak/transfer_logs/neuro-variants

# rsync -vahP --log-file $log_dir/encode_2024.transfer.log.txt $orig_dir/fragments/*_LV_*_sorted.tsv $raw_dir/

stdbuf -oL cp --verbose $orig_dir/fragments/*_LV_*_sorted.tsv $raw_dir/ > $log_dir/encode_2024.transfer.log.txt 2>&1

