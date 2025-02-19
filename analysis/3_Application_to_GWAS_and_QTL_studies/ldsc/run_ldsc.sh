#!/bin/bash

# initialize
source ~/micromamba/etc/profile.d/conda.sh

path_to_ldsc=/oak/stanford/groups/smontgom/amarder/bin/ldsc

####################################################################################

# command line arguments
trait=$1
ANNOTNAME=$2

echo "Running trait ${trait} using annotation set $ANNOTNAME ..."

HEADDIR=/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects

# N_participants=`awk -F, -v pat=$trait '$1~pat {print $5}' /oak/stanford/groups/smontgom/amarder/neuro-variants/scripts/snps/helper_func/colocal.csv`
 
samplesizeFile=/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/scripts/sample_size.csv
 
ldcts=$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/$ANNOTNAME.ldcts

####################################################################################

ldsc=$path_to_ldsc/ldsc.py
baselinedir=$path_to_ldsc/baseline
outpath=$HEADDIR/output/ldsc/results/$ANNOTNAME
TMPDIR=/oak/stanford/groups/smontgom/amarder/tmp/$trait
hapmap=$path_to_ldsc/munging/w_hm3.snplist
munge="python2 $path_to_ldsc/munge_sumstats.py"
sumstatpath=/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/out/gwas/munge/hg38/$trait/$trait.v2.txt.gz

mkdir -p $outpath
mkdir -p $TMPDIR

# N_CAS=100204
# N_CON=154587

if [ ! -e $TMPDIR/gwas_munge.sumstats.gz ]; then

conda activate r
zcat $sumstatpath | head -10000 > $TMPDIR/tmpfile
Rscript /oak/stanford/groups/smontgom/amarder/LDSC_pipeline/scripts/etc/effect_detect.R $TMPDIR
signedsumstat=`cat $TMPDIR/signedsumstat`
conda deactivate

echo "Pre-munging..."
conda activate ldsc
zcat $sumstatpath | awk 'NR==1 {print "chr\tsnp_pos\tA0\tA1\trsid\tbeta\tse\tpvalue"; for (i=1; i<=NF; i++) col[$i] = i} NR>1 {print $col["chr"], $col["snp_pos"], $col["A0"], $col["A1"], $col["rsid"], $col["beta"], $col["se"], $col["pvalue"]}' OFS="\t" | gzip > $TMPDIR/pre_munge.txt.gz
# zcat $sumstatpath | cut -f1,2,13,14,16,5,6,7 | gzip > $TMPDIR/pre_munge.txt.gz

echo "Munging..."
N_INDIV=`awk -F, -v pat=$trait '$1~pat {print $2}' $samplesizeFile`
N_CAS=`awk -F, -v pat=$trait '$1~pat {print $3}' $samplesizeFile`
N_CON=`awk -F, -v pat=$trait '$1~pat {print $4}' $samplesizeFile`
if [ "$N_CAS" = "NA" ]; then
echo "No case-control information. Running with N_participants instead..."
$munge --N $N_INDIV --sumstats $TMPDIR/pre_munge.txt.gz --out $TMPDIR/gwas_munge --merge-alleles $hapmap --chunksize 500000 --signed-sumstats $signedsumstat --a1 A0 --a2 A1
else
$munge --N-cas $N_CAS --N-con $N_CON --sumstats $TMPDIR/pre_munge.txt.gz --out $TMPDIR/gwas_munge --merge-alleles $hapmap --chunksize 500000 --signed-sumstats $signedsumstat --a1 A0 --a2 A1
fi

else

echo "Skipping munging since already complete!"

fi

echo "Running LDSC..."

newsumstatpath=$TMPDIR/gwas_munge.sumstats.gz # output from munging
# joint_set=$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/joint_set/joint_set.chr_
joint_set=$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/union.chr
origldcts=$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/$ANNOTNAME.ldcts

mkdir -p $outpath/single-annot
python2 $ldsc \
    --h2-cts $newsumstatpath \
    --ref-ld-chr $HEADDIR/output/ldsc/annot/hg19/baseline/baseline. \
    --out $outpath/single-annot/$trait \
    --ref-ld-chr-cts $origldcts \
    --w-ld-chr $baselinedir/weights_hm3_no_hla/weights.

mkdir -p $outpath/union
python2 $ldsc \
    --h2-cts $newsumstatpath \
    --ref-ld-chr $HEADDIR/output/ldsc/annot/hg19/baseline/baseline.,$joint_set \
    --out $outpath/union/$trait \
    --ref-ld-chr-cts $origldcts \
    --w-ld-chr $baselinedir/weights_hm3_no_hla/weights.


