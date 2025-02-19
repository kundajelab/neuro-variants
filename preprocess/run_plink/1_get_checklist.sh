#!/bin/bash

vcf_dir=/oak/stanford/groups/akundaje/soumyak/refs/1kg_30x_hg38_snvs_indels
manifest=$vcf_dir/20201028_genotyped-manifest.tsv
checklist=$vcf_dir/1kg_30x_hg38_snvs_indels.vcfs.checklist.chk

cd $vcf_dir

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_genotyped-manifest.tsv

awk '{print $3 "  " $1}' $manifest | grep recalibrated_variants.vcf.gz$ | grep -v chrX | grep -v chrY | grep -v others > $checklist

