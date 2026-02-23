
# need 256 to merge common and rare:
# srun --account=smontgom --partition=batch --time=24:00:00 --mem=256G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash 
# conda activate r

library(data.table)
library(stringr)
library(fgsea)

common = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/common.all_dataset.K562_bias.annot2.txt",data.table = F,stringsAsFactors = F)
rare = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/rare.all_dataset.K562_bias.annot2.txt",data.table = F,stringsAsFactors = F)

# Make column names match
common <- common[, colnames(rare), drop = FALSE]
rare <- rare[, colnames(common), drop = FALSE]

common$set = "common"
rare$set = "rare"
# colnames(rare)[colnames(rare)=="cbp_min_pval_a"] = "cbp_min_pval_fb"
df = as.data.frame(rbind(common,rare))
df$gene_distance_1_log10 = log10(df$gene_distance_1 + 1)

# extract unique cols (dataset, cell)
string = grep("abs_logfc.mean.pval.",colnames(df),value = TRUE)
dataset_lst <- str_extract(string, "(?<=\\.)[^.]+$")
cell_lst <- str_extract(string, "(?<=abs_logfc\\.mean\\.pval\\.).*(?=\\.[^.]+$)")

t.test(gene_distance_1_log10~set,df)
t.test(gene_distance_1_log10~set,df)$p.value

cor.test(ct.mg2$s_het_1,ct.mg2$y)

# no rare exclusive genes (every gene has at least 1 common variant)
# however: found 30 genes where there is no rare variant but there are common variants
# these are under less constraint
i = df$set=="rare"
gene_set1 = sort(unique(df$closest_gene_1[i]))
gene_set2 = sort(unique(df$closest_gene_1[!i]))
rare_exclusive = gene_set1[!(gene_set1 %in% gene_set2)]
common_rare_shared = gene_set1[(gene_set1 %in% gene_set2)]
common_exclusive = gene_set2[!(gene_set2 %in% gene_set1)]
length(rare_exclusive)
length(common_rare_shared)
length(common_exclusive)
# y = unique(df[df$closest_gene_1%in%common_exclusive,c("closest_gene_1","s_het_1")])
# y2 = unique(df[,c("closest_gene_1","s_het_1")])
# y2 = subset(y2,!(closest_gene_1 %in% common_exclusive))
y = unique(df[df$closest_gene_1%in%rare_exclusive,c("closest_gene_1","s_het_1")])
# y2 = unique(df[,c("closest_gene_1","s_het_1")])
# y2 = subset(y2,!(closest_gene_1 %in% common_exclusive))
t.test(y[,2],y2[,2])

# run GSEA:
library(enrichR)
# Check available databases
# dbs <- enrichR::listEnrichrDbs()
dbs_of_interest <- c("GO_Biological_Process_2023", "GO_Molecular_Function_2023","KEGG_2021_Human")

# Run enrichment for foreground genes
enrichment_results <- enrichr(rare_exclusive, dbs_of_interest)
lapply(enrichment_results,head)
enrichment_results <- enrichr(common_exclusive, dbs_of_interest)
lapply(enrichment_results,head)

i = df$set=="rare"
ct1 = data.frame(table(df$closest_gene_1[i]))
ct2 = data.frame(table(df$closest_gene_1[!i]))
ct.mg = merge(ct1,ct2,all=T,by="Var1")
ct.mg[is.na(ct.mg)] = 0
ct.mg$y = log2((ct.mg[,2] + 1)/(ct.mg[,3] + 1))
ct.mg[1]
y2 = unique(df[,c("closest_gene_1","s_het_1")])
y2 = y2[order(y2$s_het_1,decreasing = T),]
y2 = subset(y2,!duplicated(closest_gene_1))
ct.mg2 = merge(ct.mg,y2,by.x="Var1",by.y="closest_gene_1",all=T)
cor.test(ct.mg2$s_het_1,ct.mg2$y)
ct.mg2 = ct.mg2[order(ct.mg2$y,decreasing = T),]
fwrite(ct.mg2,"/oak/stanford/groups/smontgom/amarder/tmp/fgsea_input_common_v_rare.txt",quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

########################################################################################
########################################################################################
########################################################################################

df = fread("~/Downloads/fgsea_input_common_v_rare.txt",data.table = F,stringsAsFactors = F)

# Create ranked vector
ranking <- setNames(df$y, df$Var1)
library(fgsea)
library(msigdbr)

# Get Hallmark gene sets for Homo sapiens
msigdb_sets <- msigdbr(species = "Homo sapiens", category = "H") 
pathways <- split(msigdb_sets$gene_symbol, msigdb_sets$gs_name)

fgsea_res <- fgsea(
  pathways = pathways,
  stats = ranking,
  minSize = 15,
  maxSize = 500,
  nperm = 10000
)

# Order by p-value
fgsea_res <- fgsea_res[order(fgsea_res$pval), ]

# Print top results
head(fgsea_res)

# Save to file
fwrite(fgsea_res, "~/Downloads/fgsea_results.tsv", sep = "\t")
