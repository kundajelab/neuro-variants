#!/bin/bash

vcf_dir=/oak/stanford/groups/akundaje/soumyak/refs/1kg_30x_hg38_snvs_indels
checklist=$vcf_dir/1kg_30x_hg38_snvs_indels.vcfs.checklist.chk

cd $vcf_dir

md5sum -c $checklist > $vcf_dir/md5sum_check.txt

