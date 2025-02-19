#!/bin/bash

# This script processes variant annotation files for different variant sets and datasets.
# It extracts specific columns from each file and saves the output to a designated directory.

# Loop over different variant sets: 'common', 'rare', 'asd'
for variantSet in common rare asd; do
# for variantSet in asd; do # need to asd! not ran yet

# Extract chr, pos, and gene columns once and store in a file
input_file=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_annotations/$variantSet/corces_2020/K562_bias/corces_2020.Cluster1.annotations.tsv
metainfofile=/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/$variantSet.metainfo.tsv
cut -f1-5,28-33 $input_file > "$metainfofile"

# Loop over each dataset
for dataset in ameen_2022 domcke_2020 encode_2024 trevino_2021 corces_2020; do
echo $dataset

# Set the 'bias' directory
bias=K562_bias

# Define the input directory path where the variant annotations are stored
dir=/oak/stanford/groups/akundaje/projects/neuro-variants/variant_annotations/$variantSet/$dataset/$bias

# Define the output directory path for the processed data
outdir=/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/$variantSet/$dataset/$bias

# Create the output directory if it doesn't already exist
mkdir -p $outdir

# Loop over each file in the input directory
for file in $dir/*; do

# Print the full path of the current file being processed
echo $file

# Extract the filename from the full path (without the directory prefix)
filename=$(basename $file)

# Define the output file path, appending '.sub' to the filename
out=$outdir/${filename}.sub

# Extract cell type identifier
prefix="${dataset}."
suffix=".annotations.tsv"
cell=$(basename "$file" | sed -n "s/${prefix}\(.*\)${suffix}/\1/p")

# Use awk to extract specific columns from the input file 
# Rename them with the cell type suffix
# And with the dataset id
# Then write to the output file:
# BEGIN {OFS="\t"} sets the output field separator to a tab character
# {print $8, $9, $34} extracts columns 8, 9, and 34 from each line of the input file
awk 'BEGIN {OFS="\t"} {print $8, $9, $34}' $file | awk -v cell="$cell" -v dataset="$dataset" 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=3;i++)$i=$i"."cell"."dataset} NR>1{for(i=1;i<=3;i++)$i=$i} {print $1, $2, $3}' > $out

# Print the full path of the output file that is now saved
echo $out

done

# Combine all the columns into a single file
output_file=/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/$variantSet.$dataset.$bias.txt
paste "$metainfofile" "$outdir"/*.sub > "$output_file"

echo "Output written to $output_file"

done

# Combine all datasets into a single "all_dataset" file
dir=/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize
paste $dir/$variantSet.ameen_2022.K562_bias.txt \
    <(cut -f12- $dir/$variantSet.corces_2020.K562_bias.txt) \
    <(cut -f12- $dir/$variantSet.domcke_2020.K562_bias.txt) \
    <(cut -f12- $dir/$variantSet.encode_2024.K562_bias.txt) \
    <(cut -f12- $dir/$variantSet.trevino_2021.K562_bias.txt) > $dir/$variantSet.all_dataset.K562_bias.txt

echo "Final output merged and written to $dir/$variantSet.all_dataset.K562_bias.txt"

done




