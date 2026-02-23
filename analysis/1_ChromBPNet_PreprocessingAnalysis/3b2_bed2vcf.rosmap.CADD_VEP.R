# this script converts a bed into a vcf

# you really need to make sure ref/alt and position is correct, or cadd wont return results

library(bedr)
library(dplyr)
library(data.table)

# inp="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/variant_lists/Trevino_et_al_AllMutations.sort.bed"
inp=paste0("/oak/stanford/groups/akundaje/projects/neuro-variants/variant_lists/rosmap_variants.tsv")
r<-fread(inp,data.table = F,stringsAsFactors = F)
# bed <- r %>% mutate(V6=V2+1) %>% relocate(V6, .after=V2)
# r$i = 1:nrow(r)
bed = r
# bed = head(bed)
# Assuming 'bed' is your data frame
vcf <- data.frame(
  CHROM = substring(bed$V1,4),
  POS = bed$V2, # VCF uses 1-based position
  ID = ".",
  REF = bed$V3,
  ALT = bed$V4,
  QUAL = ".",
  FILTER = ".",
  INFO = "."
)

# Write VCF to file
vcf_file <- paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/rosmap.vcf")
writeLines("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", vcf_file)
options(scipen = 999)
fwrite(vcf, file = vcf_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE,na = "NA")
options(scipen = 0)
q()

################################################

conda deactivate
cd /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep
conda activate /oak/stanford/groups/smontgom/amarder/micromamba/envs/ws
# head -100 /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/rosmap.vcf | bgzip > /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep/data/vcf/rosmap.vcf.gz
# cat /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/rosmap.vcf | bcftools sort -Oz -o /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/rosmap.sorted.vcf.gz /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/rosmap.vcf

# split header and body
dir=/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists
grep '^#' $dir/rosmap.vcf > $dir/header.vcf
grep -v '^#' $dir/rosmap.vcf \
| sort -k1,1V -k2,2n \
> $dir/body.sorted.vcf

# recombine, compress, index
cat $dir/header.vcf $dir/body.sorted.vcf \
| bgzip -c > $dir/rosmap.sorted.vcf.gz
# tabix -p vcf $dir/rosmap.sorted.vcf.gz
cp $dir/rosmap.sorted.vcf.gz /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep/data/vcf/rosmap.sort.rare.vcf.gz

# cat /oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/rosmap.vcf | sort -k1,1 -k2,2 | bgzip > /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep/data/vcf/rosmap.vcf.gz
tabix /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep/data/vcf/rosmap.sort.rare.vcf.gz
VEP_PLUGINS_DIR=$CONDA_PREFIX/share/ensembl-vep-112.0-0
snakemake --cores 1 data/vcf/rosmap.sort.rare.CADD.vcf.gz --force --rerun-incomplete
# # didn't run for some reason:
# snakemake --cores 1 data/vcf/rosmap.VEP.gnomad.vcf.gz --force
# vep --verbose --vcf -i data/vcf/rosmap.vcf.gz -o data/vcf/rosmap.VEP.gnomad.vcf --distance 10000 --no_stats --force_overwrite --offline --dir_cache data/vep --regulatory --dir_plugins ${VEP_PLUGINS_DIR} --custom file=data/gnomad/gnomad.genomes.v4.0.0.sites.afonly.vcf.gz,short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AF_joint_afr%AF_joint_amr%AF_joint_asj%AF_joint_eas%AF_joint_sas%AF_joint_fin%AF_joint_nfe --plugin LoF,human_ancestor_fa:data/vep/hg38/human_ancestor.fa.gz,loftee_path:${VEP_PLUGINS_DIR},conservation_file:data/vep/hg38/loftee.sql,gerp_bigwig:data/vep/hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw --cache_version 110
vep --verbose --vcf -i data/vcf/rosmap.sort.rare.vcf.gz -o data/vcf/rosmap.sort.rare.VEP.gnomad.vcf --distance 10000 --no_stats --force_overwrite --offline --dir_cache data/annotations/vep --regulatory --dir_plugins ${VEP_PLUGINS_DIR} --custom file=data/annotations/gnomad/gnomad.genomes.v4.0.0.sites.afonly.vcf.gz,short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AF_joint_afr%AF_joint_amr%AF_joint_asj%AF_joint_eas%AF_joint_sas%AF_joint_fin%AF_joint_nfe --plugin LoF,human_ancestor_fa:data/annotations/vep/hg38/human_ancestor.fa.gz,loftee_path:${VEP_PLUGINS_DIR},conservation_file:data/annotations/vep/hg38/loftee.sql,gerp_bigwig:data/annotations/vep/hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw --cache_version 110

bgzip --keep data/vcf/rosmap.sort.rare.VEP.gnomad.vcf
tabix data/vcf/rosmap.sort.rare.VEP.gnomad.vcf.gz
snakemake --cores 1 data/vcf/rosmap.sort.rare.CADD.VEP.gnomad.vcf.gz
snakemake --cores 1 data/vcf/rosmap.sort.rare.CADD.VEP.gnomad.split.vcf.gz
sh scripts/format_v2.sh data/vcf/rosmap.sort.rare.CADD.VEP.gnomad.split.vcf.gz config/format_v2 > data/watershed/rosmap.sort.all.tsv

