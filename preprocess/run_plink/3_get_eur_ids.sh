#!/bin/bash

vcf_dir=/oak/stanford/groups/akundaje/soumyak/refs/1kg_30x_hg38_snvs_indels
samples=$vcf_dir/20130606_g1k_3202_samples_ped_population.txt
eur_ids=$vcf_dir/eur_ids.txt

cd $vcf_dir

wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt

grep -w EUR $samples | cut -f 2 -d ' ' | sort | uniq > $eur_ids

