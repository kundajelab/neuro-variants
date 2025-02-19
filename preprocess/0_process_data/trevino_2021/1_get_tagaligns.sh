#!/bin/bash

set -e
set -u
set -o pipefail

fragments_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/raw/trevino_2021/fragments
tagalign_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/trevino_2021/tagaligns

mkdir -p $tagalign_dir

for fragment_file in $fragments_dir/*.txt.gz
do
    sample=$(basename $fragment_file .all_data.pb.txt.gz)
    echo $sample

    output_file=$tagalign_dir/trevino_2021.$sample.tagAlign.gz
    zcat $fragment_file | grep ^chr | awk -v OFS='\t' '{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+"; print $1,mid,$3,"N",1000,"-"}' | bgzip > $output_file &
done
wait

