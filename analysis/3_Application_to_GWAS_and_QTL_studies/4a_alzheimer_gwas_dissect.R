# srun --account=smontgom --partition=batch --time=24:00:00 --mem=64G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# srun --account=default --partition=interactive --time=24:00:00 --mem=64G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash
# conda activate r

library(data.table)
library(stringr)
library(enrichR)
run_enrichr <- function(genes,max_fdr = 0.2,min_gene = 3) {
  dbs = c("GO_Biological_Process_2021","GO_Molecular_Function_2021")
  enrichment_results <- enrichr(genes, databases = dbs)

  # Filter significant results with adjusted p-value < 0.05 and at least 3 overlapping genes
  filtered_results.lst = list()
  for (db in dbs) {
    go_results <- enrichment_results[[db]]
    go_results$GeneCount <- as.numeric(sapply(strsplit(go_results$Overlap, "/"), `[`, 1))
    filtered_results <- go_results[go_results$Adjusted.P.value < max_fdr & go_results$GeneCount >= min_gene, ]
    if (nrow(filtered_results)==0) {next}
    filtered_results$db = db
    filtered_results.lst[[db]] = filtered_results
  }
  go_res = as.data.frame(do.call(rbind,filtered_results.lst)[,c("db","Term","Adjusted.P.value","Genes","GeneCount")])
  rownames(go_res) = NULL
  
  go_res = go_res[order(go_res$Adjusted.P.value),]
  return(go_res)
}

################

gwas = fread("/oak/stanford/groups/smontgom/amarder/HarmonizeGWAS/out/gwas/Alzheimers_Bellenguez_2022/Alzheimers_Bellenguez_2022.txt.gz",data.table = F,stringsAsFactors = F)

################################################################################

# Read in common variant data:

common = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/common.all_dataset.K562_bias.annot2.txt",data.table = F,stringsAsFactors = F)

################################################################################


traitName = "Alzheimers_Bellenguez_2022"

fm = fread(paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/all/",traitName,".susie_all.txt"),data.table = F,stringsAsFactors = F)



common$gwas = NA
fm.sub = subset(fm,PIP.L_10 > 0.2 | PIP.L_1 > 0.2)
idx = common$snp_id %in% fm.sub$snp
common$gwas[idx] = TRUE
fm.sub = subset(fm,PIP.L_10 < 0.01 & PIP.L_1 < 0.01)
idx = common$snp_id %in% fm.sub$snp
common$gwas[idx] = FALSE

# tmp = subset(common,gwas & abs_logfc.mean.pval.Cluster24.corces_2020 < 0.1 & peak_overlap.Cluster24.corces_2020)
# nrow(tmp)
# tmp$snp_id
# tmp$abs_logfc.mean.pval.Cluster24.corces_2020
# 
# y = subset(common,snp_id %in% subset(fm,sentinel=="chr5:180201150:G:A" & CS.L_1 > -1)$snp)

# fm.sub = subset(fm,(PIP.L_10 > 0.001 | PIP.L_1 > 0.001) & sentinel=="chr11:86156833:A:G")
# # picalm_tmp = subset(common,(closest_gene_3=="PICALM" | closest_gene_2=="PICALM" | closest_gene_1=="PICALM") & abs_logfc.mean.pval.Cluster24.corces_2020 < 0.05 & peak_overlap.Cluster24.corces_2020)
# picalm_tmp = subset(common,snp_id %in% fm.sub$snp & abs_logfc.mean.pval.Cluster24.corces_2020 < 0.05 & peak_overlap.Cluster24.corces_2020)
# length(picalm_tmp$snp_id)
# subset(fm.sub,snp %in% picalm_tmp$snp_id)
# 
# fm.sub = subset(fm,(PIP.L_10 > 0.005 | PIP.L_1 > 0.005))
# picalm_tmp = subset(common,snp_id %in% fm.sub$snp & abs_logfc.mean.pval.Cluster24.corces_2020 < 0.05 & peak_overlap.Cluster24.corces_2020)
# y = merge(picalm_tmp,fm.sub,by.x="snp_id",by.y="snp")
# table(y$sentinel)
# 
# 
# length(picalm_tmp$snp_id)
# picalm_tmp = subset(common,snp_id %in% fm.sub$snp)
# length(picalm_tmp$snp_id)
# 
# y = subset(common,pos==86151686)
# 
# subset(fm,snp %in% picalm_tmp$snp_id)
# 
# fm.sub = subset(fm,sentinel=="chr11:86156833:A:G")
# gwas.sub = subset(gwas,variant_alternate_id %in% fm.sub$snp)
# nrow(fm.sub)
# nrow(gwas.sub)
# sum(gwas.sub$p_value < 5e-8)
# 



# suppressMessages(suppressWarnings(library(clusterProfiler)))
# suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
# go_result <- enrichGO(gene = genes, 
#                       OrgDb = org.Hs.eg.db, 
#                       keyType = "SYMBOL",  # Specify gene identifier type
#                       ont = "BP",  # Options: "BP" (Biological Process), "MF", "CC"
#                       pAdjustMethod = "BH", 
#                       pvalueCutoff = 0.05, 
#                       qvalueCutoff = 0.2)
# go_result = as.data.frame(go_result)


library(enrichR)

run_enrichr <- function(genes,max_fdr = 0.2,min_gene = 3) {
  dbs = c("GO_Biological_Process_2021","GO_Molecular_Function_2021")
  enrichment_results <- enrichr(genes, databases = dbs)
  go_results <- enrichment_results[["GO_Biological_Process_2021"]]
  
  # Filter significant results with adjusted p-value < 0.05 and at least 3 overlapping genes
  filtered_results.lst = list()
  for (db in dbs) {
    go_results <- enrichment_results[[db]]
    go_results$GeneCount <- as.numeric(sapply(strsplit(go_results$Overlap, "/"), `[`, 1))
    filtered_results <- go_results[go_results$Adjusted.P.value < max_fdr & go_results$GeneCount >= min_gene, ]
    if (nrow(filtered_results)==0) {next}
    filtered_results$db = db
    filtered_results.lst[[db]] = filtered_results
  }
  go_res = as.data.frame(do.call(rbind,filtered_results.lst)[,c("db","Term","Adjusted.P.value","Genes","GeneCount")])
  rownames(go_res) = NULL
  
  go_res = go_res[order(go_res$Adjusted.P.value),]
  return(go_res)
}

cell_types = colnames(common)[grep("peak_overlap.",colnames(common))]
cell_types <- gsub("peak_overlap\\.", "", cell_types)
ind1 <- grepl("corces_2020", cell_types)
adult_brain <- cell_types[ind1]
adult_brain_no_microglia = adult_brain[adult_brain!="Cluster24.corces_2020"]
common$min_pval_ab = apply(common[,paste0("abs_logfc.mean.pval.",adult_brain_no_microglia)],1,min)

tmp = subset(common,(snp_id %in% subset(fm,PIP.L_10 > 0.005 | PIP.L_1 > 0.005)$snp) & 
               (abs_logfc.mean.pval.Cluster24.corces_2020 < 0.05) &
               (min_pval_ab > 0.05) &
               (chr != "chr6"))
gene_lst = sort(unique(c(tmp$closest_gene_1,tmp$closest_gene_2)))
enrichr_results = run_enrichr(gene_lst,min_gene = 0,max_fdr = Inf)


tmp2 = subset(common,(snp_id %in% subset(fm,PIP.L_10 > 0.005 | PIP.L_1 > 0.005)$snp) & 
               (chr != "chr6"))
gene_lst = sort(unique(c(tmp2$closest_gene_1,tmp2$closest_gene_2)))
enrichr_results_base = run_enrichr(gene_lst,min_gene = 0,max_fdr = Inf)

df.mg = merge(enrichr_results_base,enrichr_results,by=c("db","Term"),all=TRUE)
df.mg$Adjusted.P.value.y[is.na(df.mg$Adjusted.P.value.y)] = 1
df.mg$Adjusted.P.value.x[is.na(df.mg$Adjusted.P.value.x)] = 1

cor.test(-log10(df.mg$Adjusted.P.value.x),-log10(df.mg$Adjusted.P.value.y))

df.mg[order(df.mg$Adjusted.P.value.y)[1:10],-1]
subset(enrichr_results,grepl("amyloid",Term))

f.out = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/FinalAnalysis/Alzheimers_GSEA.txt"
fwrite(df.mg,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)



