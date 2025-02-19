#!/bin/bash

hg19_exclude=/oak/stanford/groups/akundaje/soumyak/refs/hg19/hg19_exclusion_regions.bed.gz
hg19_chrom_sizes=/oak/stanford/groups/akundaje/soumyak/refs/hg19/hg19.chrom.sizes
hg19_refs=/oak/stanford/groups/akundaje/soumyak/refs/hg19

bedtools slop -i $hg19_exclude -g $hg19_chrom_sizes -b 1057 | bgzip > $hg19_refs/hg19_exclusion_regions.slop_1057.bed.gz

