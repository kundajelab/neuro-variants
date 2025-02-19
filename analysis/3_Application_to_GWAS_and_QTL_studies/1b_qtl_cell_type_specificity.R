# srun --account=smontgom --partition=batch --time=24:00:00 --mem=64G --nodes=1 --ntasks=1 --cpus-per-task=1 --pty bash 
# conda activate r

library(data.table)
library(stringr)
# 
common = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/common.all_dataset.K562_bias.annot2.txt",data.table = F,stringsAsFactors = F)

common$chr_pos = paste0(common$chr,"_",common$pos)

cell_types = colnames(common)[grep("peak_overlap.",colnames(common))]
cell_types <- gsub("peak_overlap\\.", "", cell_types)
ind1 <- grepl("domcke_2020", cell_types) & grepl("fetal_brain", cell_types)
ind2 <- grepl("trevino_2021", cell_types)
fetal_brain <- cell_types[ind1 | ind2]

ind1 <- grepl("corces_2020", cell_types)
adult_brain <- cell_types[ind1]

ind1 <- grepl("domcke_2020", cell_types) & grepl("fetal_brain", cell_types)
ind2 <- grepl("ameen_2022", cell_types)
fetal_heart <- cell_types[ind1 | ind2]

ind1 <- grepl("encode_2024", cell_types)
adult_heart <- cell_types[ind1]

# organ = "brain"
organ = "heart_artery"
for (organ in c("brain","heart_artery")) { 
  finemap = fread(paste0("/oak/stanford/groups/smontgom/amarder/data/finemap_data/gtex/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.",organ,".txt.gz"),data.table = F,stringsAsFactors = F)
  finemap$chr_pos <- sub("(_[A-Z]+_[A-Z]+_b38)", "", finemap$V5)
  finemap = finemap[order(finemap$V6,decreasing = T),]
  finemap = finemap[!duplicated(finemap$chr_pos),]
  
  common$gwas = NA
  finemap.sub = subset(finemap,V6 > 0.9)
  idx = common$chr_pos %in% finemap.sub$chr_pos
  common$gwas[idx] = TRUE
  finemap.sub = subset(finemap,V6 < 0.01)
  idx = common$chr_pos %in% finemap.sub$chr_pos
  common$gwas[idx] = FALSE
  table(common$gwas)
  
  # df = common
  # rm(common)
  
  df = subset(common,!is.na(gwas))
  
  keep = adult_brain
  df$num_peaks_ab = apply(df[,paste0("peak_overlap.",keep)],1,sum)
  tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
  for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
  df$num_peakscbp_ab = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})
  
  keep = fetal_brain
  df$num_peaks_fb = apply(df[,paste0("peak_overlap.",keep)],1,sum)
  tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
  for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
  df$num_peakscbp_fb = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})
  
  keep = fetal_heart
  df$num_peaks_fh = apply(df[,paste0("peak_overlap.",keep)],1,sum)
  tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
  for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
  df$num_peakscbp_fh = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})
  
  keep = adult_heart
  df$num_peaks_ah = apply(df[,paste0("peak_overlap.",keep)],1,sum)
  tmp = df[,paste0("abs_logfc.mean.pval.",keep)]
  for (i in 1:length(keep)) {tmp[,i][!df[,paste0("peak_overlap.",keep[i])]] = 1}
  df$num_peakscbp_ah = apply(tmp,1,function(x) {sum(x<0.01,na.rm = T)})
  
  df$gene_distance_1_log10 = log10(df$gene_distance_1 + 1)
  
  # res = list()
  # mod = (lm(scale(num_peakscbp_fb)~s_het_1 + gene_distance_1_log10 + gwas,subset(df,num_peakscbp_fb > 0)))
  # res[[1]] = summary(mod)$coef["gwasTRUE",];res
  # mod = (lm(scale(num_peakscbp_ab)~s_het_1 + gene_distance_1_log10 + gwas,subset(df,num_peakscbp_ab > 0)))
  # res[[2]] = summary(mod)$coef["gwasTRUE",];res
  # mod = (lm(scale(num_peakscbp_fh)~s_het_1 + gene_distance_1_log10 + gwas,subset(df,num_peakscbp_fh > 0)))
  # res[[3]] = summary(mod)$coef["gwasTRUE",];res
  # mod = (lm(scale(num_peakscbp_ah)~s_het_1 + gene_distance_1_log10 + gwas,subset(df,num_peakscbp_ah > 0)))
  # res[[4]] = summary(mod)$coef["gwasTRUE",];res
  
  # res = list()
  # mod = (lm(scale(num_peakscbp_fb)~s_het_1 + gene_distance_1_log10 + gwas,df))
  # res[[1]] = summary(mod)$coef["gwasTRUE",];res
  # mod = (lm(scale(num_peakscbp_ab)~s_het_1 + gene_distance_1_log10 + gwas,df))
  # res[[2]] = summary(mod)$coef["gwasTRUE",];res
  # mod = (lm(scale(num_peakscbp_fh)~s_het_1 + gene_distance_1_log10 + gwas,df))
  # res[[3]] = summary(mod)$coef["gwasTRUE",];res
  # mod = (lm(scale(num_peakscbp_ah)~s_het_1 + gene_distance_1_log10 + gwas,df))
  # res[[4]] = summary(mod)$coef["gwasTRUE",];res
  
  res = list()
  mod = (lm(scale(num_peakscbp_fb)~s_het_1 + gene_distance_1_log10 + gwas,subset(df,num_peaks_fb > 0)))
  res[[1]] = summary(mod)$coef["gwasTRUE",];res
  mod = (lm(scale(num_peakscbp_ab)~s_het_1 + gene_distance_1_log10 + gwas,subset(df,num_peaks_ab > 0)))
  res[[2]] = summary(mod)$coef["gwasTRUE",];res
  mod = (lm(scale(num_peakscbp_fh)~s_het_1 + gene_distance_1_log10 + gwas,subset(df,num_peaks_fh > 0)))
  res[[3]] = summary(mod)$coef["gwasTRUE",];res
  mod = (lm(scale(num_peakscbp_ah)~s_het_1 + gene_distance_1_log10 + gwas,subset(df,num_peaks_ah > 0)))
  res[[4]] = summary(mod)$coef["gwasTRUE",];res
  
  res.df = as.data.frame(do.call(rbind,res))
  res.df$context = c("fb","ab","fh","ah")
  res.df
  
  # 
  # mod = (lm(scale(num_peakscbp_ab)~s_het_1 + gene_distance_1_log10 + gwas,subset(df,num_peaks_ab > 0)))
  # summary(mod)$coef["gwasTRUE",]
  # mod = (lm(scale(num_peakscbp_ah)~s_het_1 + gene_distance_1_log10 + gwas,subset(df,num_peaks_ah > 0)))
  # summary(mod)$coef["gwasTRUE",]
  # 
  
  
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/FinalAnalysis/qtl.",organ,".cell_type_specificity.txt")
  fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}


