#!/bin/bash

set -e
set -u
set -o pipefail

fragments_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/raw/encode_2024/fragments
tagalign_dir=/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/encode_2024/tagaligns

mkdir -p $tagalign_dir

for fragment_file in $fragments_dir/*.tsv
do
    base_name=$(basename $fragment_file _sorted.tsv)
    region=$(echo $base_name | rev | cut -d '_' -f2 | rev)
    condition=$(echo $base_name | rev | cut -d '_' -f1 | rev)
    celltype=$(basename $base_name _${region}_$condition)
    sample=$region.$celltype.$condition
    echo $sample

    output_file=$tagalign_dir/encode_2024.$sample.tagAlign.gz
    cat $fragment_file | grep ^chr | awk -v OFS='\t' '{mid=int(($2+$3)/2); print $1,$2,mid,"N",1000,"+"; print $1,mid,$3,"N",1000,"-"}' | bgzip > $output_file &
done
wait

