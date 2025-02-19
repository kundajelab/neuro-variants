#!/bin/bash

bcf_dir=/oak/stanford/groups/akundaje/soumyak/refs/bcf_1kg_30x_hg38_snvs_indels
bcf_eur_dir=/oak/stanford/groups/akundaje/soumyak/refs/bcf_eur_1kg_30x_hg38_snvs_indels
plink_dir=/oak/stanford/groups/akundaje/soumyak/refs/plink_1kg_30x_hg38_snvs_indels
plink_eur_dir=/oak/stanford/groups/akundaje/soumyak/refs/plink_eur_1kg_30x_hg38_snvs_indels
chrom=$1

echo chr$chrom

plink --bcf $bcf_dir/chr${chrom}.all.bcf \
    --memory 40000 \
    --keep-allele-order \
    --allow-extra-chr 0 \
    --chr 1-22 \
    --mac 1 \
    --make-bed \
    --out $plink_dir/chr${chrom}.all

cut -f 2 $plink_dir/chr${chrom}.all.bim | grep \* | sort | uniq > $plink_dir/chr${chrom}.all.exclude.txt

plink --bfile $plink_dir/chr${chrom}.all \
    --memory 40000 \
    --exclude $plink_dir/chr${chrom}.all.exclude.txt \
    --keep-allele-order \
    --mac 1 \
    --make-bed \
    --out $plink_dir/chr${chrom}.all.filtered

#-----------------------------------------------#

plink --bcf $bcf_eur_dir/chr${chrom}.eur.bcf \
    --memory 40000 \
    --keep-allele-order \
    --allow-extra-chr 0 \
    --chr 1-22 \
    --mac 1 \
    --make-bed \
    --out $plink_eur_dir/chr${chrom}.eur

cut -f 2 $plink_eur_dir/chr${chrom}.eur.bim | grep \* | sort | uniq > $plink_eur_dir/chr${chrom}.eur.exclude.txt

plink --bfile $plink_eur_dir/chr${chrom}.eur \
    --memory 40000 \
    --exclude $plink_eur_dir/chr${chrom}.eur.exclude.txt \
    --keep-allele-order \
    --mac 1 \
    --make-bed \
    --out $plink_eur_dir/chr${chrom}.eur.filtered

