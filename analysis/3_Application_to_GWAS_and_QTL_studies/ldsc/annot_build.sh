#!/bin/bash

source ~/micromamba/etc/profile.d/conda.sh

################################################################################
################################################################################

# Initialize ## These are dataset-specific! ##
for i in {5..5}; do
  case $i in
    1)
      ANNOTNAME="CRC"
      INDIR="/oak/stanford/groups/akundaje/projects/CRC_finemap/data/processed/peaks/overlap"
      build="hg38"
      ;;
    2)
      ANNOTNAME="ameen_2022"
      INDIR="/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/$ANNOTNAME/peaks/overlap"
      build="hg38"
      ;;
    3)
      ANNOTNAME="corces_2020"
      INDIR="/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/$ANNOTNAME/peaks/overlap"
      build="hg38"
      ;;
    4)
      ANNOTNAME="domcke_2020"
      INDIR="/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/$ANNOTNAME/peaks/overlap"
      build="hg19"
      ;;
    5)
      ANNOTNAME="encode_2024"
      INDIR="/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/$ANNOTNAME/peaks/overlap"
      build="hg38"
      ;;
    6)
      ANNOTNAME="roussos_2024"
      INDIR="/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/$ANNOTNAME/peaks/overlap"
      build="hg38"
      ;;
    7)
      ANNOTNAME="trevino_2021"
      INDIR="/oak/stanford/groups/akundaje/projects/neuro-variants/data/processed/$ANNOTNAME/peaks/overlap"
      build="hg38"
      ;;
    8)
      ANNOTNAME="IGVFCoronaryArteriesMultiome_2024"
      INDIR="/oak/stanford/groups/akundaje/projects/cad/data/processed/human_multiome_washu/peaks/overlap"
      build="hg38"
      ;;
    9)
      ANNOTNAME="turner_2022"
      INDIR="/oak/stanford/groups/akundaje/projects/cad/data/processed/clint/peaks/overlap"
      build="hg38"
      ;;
    *)
      echo "Invalid value for i: $i"
      exit 1
      ;;
  esac

  # Here you can add the commands to process each set of parameters
  echo "Processing $ANNOTNAME with build $build at $INDIR"


# Initialize ## These are cluster-specific! ##
hg38ToHg19chain=/oak/stanford/groups/smontgom/amarder/bin/chain/hg38ToHg19.over.chain.gz
path_to_ldsc=/oak/stanford/groups/smontgom/amarder/bin/ldsc
HEADDIR=/oak/stanford/groups/smontgom/amarder/neuro-variants

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# Preprocessing
mkdir -p $HEADDIR/output/ldsc/annot/hg19
mkdir -p $HEADDIR/output/ldsc/annot/hg38
mkdir -p $HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME
mkdir -p $HEADDIR/output/ldsc/annot/hg38/$ANNOTNAME

# ldsc files
ldsc=$path_to_ldsc/ldsc.py
make_annot=$path_to_ldsc/make_annot.py

################################################################################
################################################################################
# 
# # LiftOver!
# 
# echo "LiftOver for $ANNOTNAME..."
# 
# celltypes=$(ls $INDIR | sed 's/\.overlap.*//g' | uniq)
# 
# # Cycle through diff cell types:
# for cell in $celltypes; do
# 
# # Set parameters:
# echo $cell
# infile=$INDIR/$cell.overlap.peaks.bed.gz
# hg38peaks=$HEADDIR/output/ldsc/annot/hg38/$ANNOTNAME/$cell.overlap.peaks.bed
# hg19peaks=$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/$cell.overlap.peaks.bed
# 
# # Need to liftover if in hg38!
# if [ "$build" = "hg38" ]; then
# 
# zcat "$infile" | cut -f1,2,3,4 > $hg38peaks
# conda activate kent-tools
# liftOver $hg38peaks "$hg38ToHg19chain" $hg19peaks $hg19peaks.unmapp
# conda deactivate
# rm $hg19peaks.unmapp
# 
# # Do not need to liftover hg19:
# elif [ "$build" = "hg19" ]; then
# 
# zcat "$infile" | cut -f1,2,3,4 > $hg19peaks
# 
# fi
# done
# 
# 
# # # Create union set: 
# # for some reason, sometimes throws error: 
# # Error: Type checker found wrong number of fields while tokenizing data line.
# # Perhaps you have extra TAB at the end of your line? Check with "cat -t"
# cat $HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/*.overlap.peaks.bed > $HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/union.bed
# sort -k1,1 -k2,2n $HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/union.bed | sed 's/\t*$//' > $HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/union_sorted.bed
# conda activate bedtools
# bedtools merge -i $HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/union_sorted.bed > $HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/union.overlap.peaks.bed
# conda deactivate
# rm $HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/union.bed
# rm $HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/union_sorted.bed
# 
################################################################################
################################################################################

# Build LDSC annotation:

conda activate ldsc
celltypes=$(ls $INDIR | sed 's/\.overlap.*//g' | uniq)
# for cell in $celltypes; do
# for cell in union; do
for cell in $celltypes union; do

for chrnum in {1..22}; do

# set paths:
hg19peaks=$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/$cell.overlap.peaks.bed
annotation_output=$HEADDIR/output/ldsc/annot/hg19/$ANNOTNAME/$cell.chr$chrnum
hapmap=$path_to_ldsc/hapmap/hapmap3_snps/hm.${chrnum}.snp

if [ -e $annotation_output.l2.ldscore.gz ]; then
    echo "$ANNOTNAME exists."
else

# Genotype file - you want this to be on the same chr build as the BEDFILE
# (but: don't need GWAS and BEDFILE to be the same build, since SNP ids are used, not positions)
# geno=/oak/stanford/groups/smontgom/amarder/bin/1kg/1kg_eur/plink/1kg.all.b37.chr${chrnum}
geno=/oak/stanford/groups/smontgom/amarder/bin/ldsc/1kg/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrnum}

# the annotation is basically just 0 or 1 determining whether the 1000G SNPs are in the BEDFILE
# $annotation_output.gz has the same number of rows as $geno.bim
echo "running make_annot.."
$make_annot --bed-file $hg19peaks --bimfile $geno.bim --annot-file $annotation_output.annot.gz

echo "running ldsc..."
# Note: having --out be the same as --annot allows the same file prefix
# It will NOT overwrite the annot.gz file.
# --print-snps $hapmap writes only hapmap snps to the final file, which allows a smaller annotation file. the l2.ldscore.gz is what is used for ld score regression
$ldsc --l2 --bfile $geno --ld-wind-cm 1 --annot $annotation_output.annot.gz --thin-annot --out $annotation_output --print-snps $hapmap
echo "job complete."

fi

done
done

done

conda deactivate

