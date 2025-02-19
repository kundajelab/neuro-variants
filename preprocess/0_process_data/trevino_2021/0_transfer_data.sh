#!/bin/bash

set -e
set -u
set -o pipefail

orig_dir=/oak/stanford/groups/smontgom/amarder/data/fetal_brain_fragments/pseudobulk_all_data

raw_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/raw/trevino_2021/fragments
log_dir=/home/groups/akundaje/soumyak/transfer_logs/neuro-variants

rsync -varhP --log-file $log_dir/trevino_2021.transfer.log.txt $orig_dir/*.txt.gz $raw_dir/

