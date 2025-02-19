#!/bin/bash

set -e
set -u
set -o pipefail

processed_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/domcke_2020/peaks
orig_peak_dir=/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/peaks/domcke_2020/hg38

echo Copying Peaks

cp -r $orig_peak_dir $processed_dir

echo Checking Diff

diff -r $processed_dir/hg38 $orig_peak_dir

