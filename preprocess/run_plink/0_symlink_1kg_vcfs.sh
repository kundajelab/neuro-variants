#!/bin/bash

orig_vcf_dir=/oak/stanford/groups/smontgom/amarder/bin/1kg/30x
new_vcf_dir=/oak/stanford/groups/akundaje/soumyak/refs/1kg_30x_hg38_snvs_indels

cd $new_vcf_dir

for i in {1..22}
do
    echo chr$i
    ln -s $orig_vcf_dir/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.recalibrated_variants.vcf.gz
done

