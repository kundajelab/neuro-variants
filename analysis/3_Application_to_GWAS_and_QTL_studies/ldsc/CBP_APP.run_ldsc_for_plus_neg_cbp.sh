source ~/micromamba/etc/profile.d/conda.sh

path_to_ldsc=/oak/stanford/groups/smontgom/amarder/bin/ldsc


conda activate ldsc

# try to run with ldsc_loop.sh with:

annotation="extended"
# annotation="extended2"

# and
trait="Alzheimers_Bellenguez_2022"
# trait="Schizophrenia_PGCWave3_2022"

for annotation in "extended" "extended2"; do
for trait in "Alzheimers_Bellenguez_2022" "Schizophrenia_PGCWave3_2022"; do

ANNOTNAME=$annotation

#######################

echo "Running trait ${trait} using annotation set $ANNOTNAME ..."

HEADDIR=/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects
 
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


echo "Running LDSC..."

newsumstatpath=$TMPDIR/gwas_munge.sumstats.gz # output from munging
# joint_set=$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/joint_set/joint_set.chr_
joint_set=$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/union.chr
origldcts=$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/$ANNOTNAME.ldcts

# mkdir -p $outpath/single-annot
# python2 $ldsc \
#     --h2-cts $newsumstatpath \
#     --ref-ld-chr $HEADDIR/output/ldsc/annot/hg19/baseline/baseline. \
#     --out $outpath/single-annot/$trait \
#     --ref-ld-chr-cts $origldcts \
#     --w-ld-chr $baselinedir/weights_hm3_no_hla/weights.

mkdir -p $outpath/joint-annot
python2 $ldsc \
  --h2 $newsumstatpath \
  --ref-ld-chr \
$HEADDIR/output/ldsc/annot/hg19/baseline/baseline.,\
$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/plus_cbp.chr \
  --w-ld-chr $baselinedir/weights_hm3_no_hla/weights. \
  --overlap-annot \
  --frqfile-chr $baselinedir/1000G_Phase3_frq/1000G.EUR.QC. \
  --print-coefficients \
  --out $outpath/joint-annot/$trait.plus

python2 $ldsc \
  --h2 $newsumstatpath \
  --ref-ld-chr \
$HEADDIR/output/ldsc/annot/hg19/baseline/baseline.,\
$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/neg_cbp.chr \
  --w-ld-chr $baselinedir/weights_hm3_no_hla/weights. \
  --overlap-annot \
  --frqfile-chr $baselinedir/1000G_Phase3_frq/1000G.EUR.QC. \
  --print-coefficients \
  --out $outpath/joint-annot/$trait.neg

python2 $ldsc \
  --h2 $newsumstatpath \
  --ref-ld-chr \
$HEADDIR/output/ldsc/annot/hg19/baseline/baseline.,\
$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/peak.chr \
  --w-ld-chr $baselinedir/weights_hm3_no_hla/weights. \
  --overlap-annot \
  --frqfile-chr $baselinedir/1000G_Phase3_frq/1000G.EUR.QC. \
  --print-coefficients \
  --out $outpath/joint-annot/$trait.peak

done
done

# 
# python2 $ldsc \
#   --h2 $newsumstatpath \
#   --ref-ld-chr \
# $HEADDIR/output/ldsc/annot/hg19/baseline/baseline.,\
# $HEADDIR/output/ldsc/annot/hg19/extended/plus_cbp.chr,\
# $HEADDIR/output/ldsc/annot/hg19/extended/neg_cbp.chr \
#   --w-ld-chr $baselinedir/weights_hm3_no_hla/weights. \
#   --overlap-annot \
#   --frqfile-chr $baselinedir/1000G_Phase3_frq/1000G.EUR.QC. \
#   --out $outpath/joint-annot/$trait.plus
# 
