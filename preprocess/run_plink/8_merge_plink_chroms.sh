#!/bin/bash

plink_eur_dir=/oak/stanford/groups/akundaje/soumyak/refs/plink_eur_1kg_30x_hg38_snvs_indels
mergelist=plink.eur.mergelist.txt
outdir=plink_outputs

echo chr1
echo $plink_eur_dir/chr1.eur.filtered > $mergelist

for chrom in {2..22}
do
    echo chr$chrom
    echo $plink_eur_dir/chr${chrom}.eur.filtered >> $mergelist
done

plink --merge-list $mergelist \
    --memory 400000 \
    --keep-allele-order \
    --mac 1 \
    --make-bed \
    --out $outdir/1kg_30x_hg38_snvs_indels.eur.filtered.merged

