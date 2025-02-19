#!/bin/bash

set -e
set -u
set -o pipefail

orig_dir=/mnt/lab_data3/salil512/cardiogenesis_reprocess/3_celltype_fragments

raw_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/raw/ameen_2022/fragments
log_dir=/users/soumyak/transfer_logs/neuro-variants

# rsync -vahP --log-file $log_dir/ameen_2022.transfer.log.txt $orig_dir/fragments/*_sorted.tsv $raw_dir/

stdbuf -oL cp --verbose $orig_dir/fragments/*_sorted.tsv $raw_dir/ > $log_dir/ameen_2022.transfer.log.txt 2>&1

