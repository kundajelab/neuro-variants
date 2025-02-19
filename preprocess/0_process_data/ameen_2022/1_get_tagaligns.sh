#!/bin/bash

set -e
set -u
set -o pipefail

fragments_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/raw/ameen_2022/fragments
tagalign_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/ameen_2022/tagaligns

mkdir -p $tagalign_dir

for fragment_file in $fragments_dir/*.tsv
do
    sample=$(basename $fragment_file _sorted.tsv)
    echo $sample

    output_file=$tagalign_dir/ameen_2022.$sample.tagAlign.gz
    cat $fragment_file | grep ^chr | awk -v OFS='\t' '{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+"; print $1,mid,$3,"N",1000,"-"}' | bgzip > $output_file &
done
wait

