#!/bin/bash

logdir=/oak/stanford/groups/akundaje/soumyak/refs/logs/1kg_plink
jobscript=jobscript.sh
plink_script=6_run_plink.sh

for chrom in {1..22}
do
    echo chr$chrom

    sbatch -J plink.chr$chrom \
           -t 360 -c 4 --mem=50G -p akundaje,owners --requeue \
           -o $logdir/chr$chrom.log.txt \
           -e $logdir/chr$chrom.err.txt  \
           $jobscript bash $plink_script $chrom
done

