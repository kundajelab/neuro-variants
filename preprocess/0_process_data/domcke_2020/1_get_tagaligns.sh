#!/bin/bash

set -e
set -u
set -o pipefail

raw_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/raw/domcke_2020
tagalign_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/domcke_2020/tagaligns

src_dir=/users/soumyak/neuro-variants/src

mkdir -p $tagalign_dir

for dataset in fetal_brain fetal_heart
do
    echo
    echo $dataset
    echo
    fragments_dir=$raw_dir/$dataset/frags

    for fragment_file in $fragments_dir/*
    do
        sample=$(basename "$fragment_file" .frags.txt.gz | sed 's/ /_/g' | sed 's/\?//g' | sed 's/\./_/g')
        echo $sample

        output_file=$tagalign_dir/domcke_2020.$dataset.$sample.tagAlign.gz
        zcat "$fragment_file" | grep ^chr | awk -v OFS='\t' '{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+"; print $1,mid,$3,"N",1000,"-"}' | bgzip > $output_file &

    done
done
wait

