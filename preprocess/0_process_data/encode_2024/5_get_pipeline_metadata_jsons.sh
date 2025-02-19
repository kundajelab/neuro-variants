#!/bin/bash

set -e
set -u
set -o pipefail

pipeline_output_dir=gs://neuro-variants/data/processed/atac_pipeline_outputs/atac
metadata_json_dir=pipeline_metadata_jsons
hide_result_before=2024-06-09T10:04:27.714Z

mkdir -p $metadata_json_dir

caper list --hide-result-before $hide_result_before > caper_list.tsv

while IFS=$'\t' read -r id status name str_label _
do
    if [ "$id" == "id" ]; then
        continue
    fi

    if [[ $status == "Succeeded" ]]; then
        prefix=$(echo "$str_label" | cut -d '.' -f1)
        if [[ $prefix == "encode_2024" ]]; then
            src=$pipeline_output_dir/$id/metadata.json
            dest=$metadata_json_dir/$str_label.metadata.json

            echo $id : $str_label
            echo $src
            echo $dest
            echo

            gsutil cp $src $dest
        fi
    fi
done < caper_list.tsv

