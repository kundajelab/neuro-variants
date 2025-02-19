#!/bin/bash

logdir=/oak/stanford/groups/akundaje/soumyak/refs/logs/1kg_bcfs
jobscript=jobscript.sh
bcf_script=4_get_bcfs.sh

for chrom in {1..22}
do
    echo chr$chrom

    sbatch -J bcfs.chr$chrom \
           -t 1-0 -c 4 --mem=60G -p akundaje,owners --requeue \
           -o $logdir/chr$chrom.log.txt \
           -e $logdir/chr$chrom.err.txt  \
           $jobscript bash $bcf_script $chrom
done

