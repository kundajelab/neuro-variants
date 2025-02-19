ANNOTNAME=$1
HEADDIR=/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects
ldcts=$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/$ANNOTNAME.ldcts

rm $ldcts # start from scratch
celltypes=$(ls $HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/*.chr22.l2.ldscore.gz | sed 's/\.chr22.l2.ldscore.*//g' | uniq | xargs -n 1 basename | grep -v '^union$')
for CELLTYPE in $celltypes; do
echo -e "${CELLTYPE}\t$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/${CELLTYPE}.chr@" >> $ldcts
done
wc -l $ldcts
