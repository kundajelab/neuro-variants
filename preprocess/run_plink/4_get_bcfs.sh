#!/bin/bash

vcf_dir=/oak/stanford/groups/akundaje/soumyak/refs/1kg_30x_hg38_snvs_indels
bcf_dir=/oak/stanford/groups/akundaje/soumyak/refs/bcf_1kg_30x_hg38_snvs_indels
bcf_eur_dir=/oak/stanford/groups/akundaje/soumyak/refs/bcf_eur_1kg_30x_hg38_snvs_indels
eur_ids=/oak/stanford/groups/akundaje/soumyak/refs/1kg_30x_hg38_snvs_indels/eur_ids.txt
hg38_fasta=/oak/stanford/groups/akundaje/soumyak/refs/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
chrom=$1

echo chr$chrom

bcftools norm -Ou -m -any $vcf_dir/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrom}.recalibrated_variants.vcf.gz |
bcftools norm -Ou -f $hg38_fasta |
bcftools annotate -Ob -x ID \
    -I +'%CHROM:%POS:%REF:%ALT' > $bcf_dir/chr${chrom}.all.bcf

#------------------------------#

bcftools view -Ob --samples-file $eur_ids \
    $bcf_dir/chr${chrom}.all.bcf > $bcf_eur_dir/chr${chrom}.eur.bcf

