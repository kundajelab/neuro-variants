# this script converts a bed into a vcf

library(bedr)
library(dplyr)
library(data.table)

inp="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/chd_snv_list.tsv"
# inp="~/Downloads/chd_snv_list.tsv"
r<-fread(inp,data.table = F,stringsAsFactors = F)
bed = r
bed[,6] = bed$V2 - 1
bed = bed[,c(1,6,2,3,4,5)]
colnames(bed) = paste0("V",1:6)

# Assuming 'bed' is your data frame
vcf <- data.frame(
  CHROM = substring(bed$V1,4),
  POS = bed$V2 + 1, # VCF uses 1-based position
  ID = ".",
  REF = bed$V4,
  ALT = bed$V5,
  QUAL = ".",
  FILTER = ".",
  INFO = "."
)

# Write VCF to file
vcf_file <- sprintf("%s.vcf", fs::path_ext_remove(inp))
writeLines("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", vcf_file)
fwrite(vcf, file = vcf_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE,na = "NA")

################################################

conda deactivate
cd /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep
conda activate /oak/stanford/groups/smontgom/amarder/micromamba/envs/ws
# head -100 /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/chd_snv_list.vcf | bgzip > /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep/data/vcf/chd_snv_list.vcf.gz
# cat /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/chd_snv_list.vcf | bcftools sort -Oz -o /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/chd_snv_list.sorted.vcf.gz /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/chd_snv_list.vcf

# split header and body
dir=/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists
grep '^#' $dir/chd_snv_list.vcf > $dir/header.vcf
grep -v '^#' $dir/chd_snv_list.vcf \
| sort -k1,1V -k2,2n \
> $dir/body.sorted.vcf

# recombine, compress, index
cat $dir/header.vcf $dir/body.sorted.vcf \
| bgzip -c > $dir/chd_snv_list.sorted.vcf.gz
tabix -p vcf $dir/chd_snv_list.sorted.vcf.gz
cp $dir/chd_snv_list.sorted.vcf.gz /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep/data/vcf/chd_snv_list.sort.rare.vcf.gz

# cat /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/chd_snv_list.vcf | sort -k1,1 -k2,2 | bgzip > /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep/data/vcf/chd_snv_list.vcf.gz
tabix /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep/data/vcf/chd_snv_list.sort.rare.vcf.gz
VEP_PLUGINS_DIR=$CONDA_PREFIX/share/ensembl-vep-112.0-0
snakemake --cores 1 data/vcf/chd_snv_list.sort.rare.CADD.vcf.gz --force --rerun-incomplete
# # didn't run for some reason:
# snakemake --cores 1 data/vcf/chd_snv_list.VEP.gnomad.vcf.gz --force
# updated paths so this is the old one:
# vep --verbose --vcf -i data/vcf/chd_snv_list.vcf.gz -o data/vcf/chd_snv_list.VEP.gnomad.vcf --distance 10000 --no_stats --force_overwrite --offline --dir_cache data/vep --regulatory --dir_plugins ${VEP_PLUGINS_DIR} --custom file=data/gnomad/gnomad.genomes.v4.0.0.sites.afonly.vcf.gz,short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AF_joint_afr%AF_joint_amr%AF_joint_asj%AF_joint_eas%AF_joint_sas%AF_joint_fin%AF_joint_nfe --plugin LoF,human_ancestor_fa:data/vep/hg38/human_ancestor.fa.gz,loftee_path:${VEP_PLUGINS_DIR},conservation_file:data/vep/hg38/loftee.sql,gerp_bigwig:data/vep/hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw --cache_version 110
vep --verbose --vcf -i data/vcf/chd_snv_list.sort.rare.vcf.gz -o data/vcf/chd_snv_list.sort.rare.VEP.gnomad.vcf --distance 10000 --no_stats --force_overwrite --offline --dir_cache data/annotations/vep --regulatory --dir_plugins ${VEP_PLUGINS_DIR} --custom file=data/annotations/gnomad/gnomad.genomes.v4.0.0.sites.afonly.vcf.gz,short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AF_joint_afr%AF_joint_amr%AF_joint_asj%AF_joint_eas%AF_joint_sas%AF_joint_fin%AF_joint_nfe --plugin LoF,human_ancestor_fa:data/annotations/vep/hg38/human_ancestor.fa.gz,loftee_path:${VEP_PLUGINS_DIR},conservation_file:data/annotations/vep/hg38/loftee.sql,gerp_bigwig:data/annotations/vep/hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw --cache_version 110
bgzip --keep data/vcf/chd_snv_list.sort.rare.VEP.gnomad.vcf
tabix data/vcf/chd_snv_list.sort.rare.VEP.gnomad.vcf.gz
snakemake --cores 1 data/vcf/chd_snv_list.sort.rare.CADD.VEP.gnomad.vcf.gz
snakemake --cores 1 data/vcf/chd_snv_list.sort.rare.CADD.VEP.gnomad.split.vcf.gz
sh scripts/format_v2.sh data/vcf/chd_snv_list.sort.rare.CADD.VEP.gnomad.split.vcf.gz config/format_v2 > data/watershed/chd_snv_list.sort.all.tsv

