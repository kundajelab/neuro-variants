# this script converts a bed into a vcf

file=/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/gt_0.05/chrALL.filter.score.bed

module load bedtools
cat $file | bedtools sort -i - > /oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/gt_0.05/chrALL.filter.score.sort.bed

##################

library(bedr)
library(dplyr)
library(data.table)

# inp="/oak/stanford/groups/smontgom/amarder/neuro-variants/output/data/variant_lists/Trevino_et_al_AllMutations.sort.bed"
inp=paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/gt_0.05/chrALL.filter.score.sort.bed")
r<-fread(inp,data.table = F,stringsAsFactors = F)
# bed <- r %>% mutate(V6=V2+1) %>% relocate(V6, .after=V2)
# r$i = 1:nrow(r)
bed = r

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
# writeLines("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", vcf_file, append=TRUE)
options(scipen = 999)
fwrite(vcf, file = vcf_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE,na = "NA")
options(scipen = 0)
q()


conda deactivate
cd /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep
conda activate /oak/stanford/groups/smontgom/amarder/micromamba/envs/ws
cat /oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/1kg_variants/gt_0.05/chrALL.filter.score.sort.vcf | bgzip > /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep/data/vcf/1kg.common.gt_0.05.sort.rare.vcf.gz
tabix /oak/stanford/groups/smontgom/amarder/bin/watershed_dataprep/data/vcf/1kg.common.gt_0.05.sort.rare.vcf.gz
snakemake --cores 1 data/vcf/1kg.common.gt_0.05.sort.rare.CADD.vcf.gz --force --rerun-incomplete
VEP_PLUGINS_DIR=$CONDA_PREFIX/share/ensembl-vep-112.0-0
# # didn't run for some reason:
# snakemake --cores 1 data/vcf/Trevino_et_al_AllMutations.sort.rare.VEP.gnomad.vcf.gz --force
vep --verbose --vcf -i data/vcf/1kg.common.gt_0.05.sort.rare.vcf.gz -o data/vcf/1kg.common.gt_0.05.sort.rare.VEP.gnomad.vcf --distance 10000 --no_stats --force_overwrite --offline --dir_cache data/vep --regulatory --dir_plugins ${VEP_PLUGINS_DIR} --custom file=data/gnomad/gnomad.genomes.v4.0.0.sites.afonly.vcf.gz,short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AF_joint_afr%AF_joint_amr%AF_joint_asj%AF_joint_eas%AF_joint_sas%AF_joint_fin%AF_joint_nfe --plugin LoF,human_ancestor_fa:data/vep/hg38/human_ancestor.fa.gz,loftee_path:${VEP_PLUGINS_DIR},conservation_file:data/vep/hg38/loftee.sql,gerp_bigwig:data/vep/hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw --cache_version 110
bgzip --keep data/vcf/1kg.common.gt_0.05.sort.rare.VEP.gnomad.vcf
tabix data/vcf/1kg.common.gt_0.05.sort.rare.VEP.gnomad.vcf.gz
snakemake --cores 1 data/vcf/1kg.common.gt_0.05.sort.rare.CADD.VEP.gnomad.vcf.gz --force --rerun-incomplete
snakemake --cores 1 data/vcf/1kg.common.gt_0.05.sort.rare.CADD.VEP.gnomad.split.vcf.gz --force --rerun-incomplete
sh scripts/format_v2.sh data/vcf/1kg.common.gt_0.05.sort.rare.CADD.VEP.gnomad.split.vcf.gz config/format_v2 > data/watershed/1kg.common.gt_0.05.sort.all.tsv