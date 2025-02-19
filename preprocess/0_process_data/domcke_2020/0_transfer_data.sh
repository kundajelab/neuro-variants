#!/bin/bash

set -e
set -u
set -o pipefail

orig_brain_dir=/oak/stanford/groups/smontgom/amarder/data/Domcke_2020_Nature/pseudobulk/cluster
orig_heart_dir=/oak/stanford/groups/smontgom/amarder/data/Domcke_2020_Nature_HEART/pseudobulk/cluster

raw_brain_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/raw/domcke_2020/fetal_brain
raw_heart_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/raw/domcke_2020/fetal_heart
log_dir=/home/groups/akundaje/soumyak/transfer_logs/neuro-variants

mkdir -p $raw_brain_dir
mkdir -p $raw_heart_dir

rsync -varhP --log-file $log_dir/domcke_2020.fetal_brain.transfer.log.txt $orig_brain_dir/* $raw_brain_dir/
rsync -varhP --log-file $log_dir/domcke_2020.fetal_heart.transfer.log.txt $orig_heart_dir/* $raw_heart_dir/

