#!/bin/bash

# Rare vs common analysis - SLURM submission script

# cd /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/scripts/slurm
# sbatch --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/scripts/4_Rare_vs_Common/slurm_submit_rare_vs_common.sh

# Scripts to run:

dir=/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/scripts/4_Rare_vs_Common

source ~/.bashrc
conda activate r

echo "Script 1..."
Rscript $dir/1a_rare_vs_common.R

echo "Script 2..."
Rscript $dir/2a_rare_vs_common_celltype_specificity.R
