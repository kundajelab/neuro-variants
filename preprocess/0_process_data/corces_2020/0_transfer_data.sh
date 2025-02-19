#!/bin/bash

set -e
set -u
set -o pipefail

orig_dir=/oak/stanford/groups/akundaje/projects/alzheimers_parkinsons/pseudobulk/pseudobulk_rds

raw_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/raw/corces_2020
log_dir=/home/groups/akundaje/soumyak/transfer_logs/neuro-variants

rsync -varhP --log-file $log_dir/corces_2020.transfer.log.txt $orig_dir $raw_dir/

