library(data.table)
library(stringr)

################################################################################

# Read in common variant data:

common = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/common.all_dataset.K562_bias.annot2.txt",data.table = F,stringsAsFactors = F)

################################################################################

# This section analyzes fine-mapped GWAS meta-analysis results that have been fine-mapped using an out-of-sample 1000G reference panel.

trait_lst = c("Alzheimers_Bellenguez_2022","CAD_Tcheandjieu_2022","CAD_Aragam_2022")
# traitName = "Alzheimers_Bellenguez_2022"

for (traitName in trait_lst) {
  print(traitName)
  fm = fread(paste0("/oak/stanford/groups/smontgom/amarder/VariantPrioritization/out/finemap/all/",traitName,".susie_all.txt"),data.table = F,stringsAsFactors = F)
  
  common$gwas = NA
  fm.sub = subset(fm,PIP.L_10 > 0.2 | PIP.L_1 > 0.2)
  idx = common$snp_id %in% fm.sub$snp
  common$gwas[idx] = TRUE
  fm.sub = subset(fm,PIP.L_10 < 0.01 & PIP.L_1 < 0.01)
  idx = common$snp_id %in% fm.sub$snp
  common$gwas[idx] = FALSE
  
  # contexts
  string = grep("abs_logfc.mean.pval.",colnames(common),value = TRUE)
  dataset_lst <- str_extract(string, "(?<=\\.)[^.]+$")
  cell_lst <- str_extract(string, "(?<=abs_logfc\\.mean\\.pval\\.).*(?=\\.[^.]+$)")
  #
  res.lst=list();j=1;
  for (i in 1:length(cell_lst)) {
    # cell type to test:
    cell=cell_lst[i]
    dataset=dataset_lst[i]
    print(paste0(i,"/",length(cell_lst),": ",cell," (",dataset,")"))
    idx = common[,paste0("peak_overlap.",cell,".",dataset)]
    # test predicting fine map variants
    if (sum(common$gwas[idx],na.rm=T)>1) {
      mod = lm((common[idx,paste0("abs_logfc.mean.",cell,".",dataset)])~common$gwas[idx] + log10(common$gene_distance_1[idx]+1) + common$s_het_1[idx])
      res = summary(mod)$coef[2,]
      j = j + 1
      res.lst[[j]] = data.frame(trait=traitName,dataset,cell,Estimate=res["Estimate"],se = res["Std. Error"],tvalue=res["t value"],pval=res["Pr(>|t|)"])
    }
  }
  
  res.df = as.data.frame(do.call(rbind,res.lst))
  if (nrow(res.df) > 0) {
    print(res.df[order(res.df$pval,decreasing = F),][1:10,])
  }
  
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_gwas_binarizedpip_enrichments.",traitName,".txt")
  fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}


###################################

# This section analyzes GTEx QTL results that have been fine-mapped using in-sample UKB LD.

# QTL

for (organ in c("heart_artery","brain")) {
  
  finemap = fread(paste0("/oak/stanford/groups/smontgom/amarder/data/finemap_data/gtex/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.",organ,".txt.gz"),data.table = F,stringsAsFactors = F)
  finemap$chr_pos <- sub("(_[A-Z]+_[A-Z]+_b38)", "", finemap$V5)
  finemap = finemap[order(finemap$V6,decreasing = T),]
  finemap = finemap[!duplicated(finemap$chr_pos),]
  common$chr_pos = paste0(common$chr,"_",common$pos)
  
  common$gwas = NA
  finemap.sub = subset(finemap,V6 > 0.9)
  idx = common$chr_pos %in% finemap.sub$chr_pos
  common$gwas[idx] = TRUE
  finemap.sub = subset(finemap,V6 < 0.01)
  idx = common$chr_pos %in% finemap.sub$chr_pos
  common$gwas[idx] = FALSE
  table(common$gwas)
  
  # contexts
  string = grep("abs_logfc.mean.pval.",colnames(common),value = TRUE)
  dataset_lst <- str_extract(string, "(?<=\\.)[^.]+$")
  cell_lst <- str_extract(string, "(?<=abs_logfc\\.mean\\.pval\\.).*(?=\\.[^.]+$)")
  #
  res.lst=list();j=1;
  for (i in 1:length(cell_lst)) {
    # cell type to test:
    cell=cell_lst[i]
    dataset=dataset_lst[i]
    print(paste0(i,"/",length(cell_lst),": ",cell," (",dataset,")"))
    idx = common[,paste0("peak_overlap.",cell,".",dataset)]
    # test predicting fine map variants
    if (sum(common$gwas[idx],na.rm=T)>1) {
      mod = lm((common[idx,paste0("abs_logfc.mean.",cell,".",dataset)])~common$gwas[idx] + log10(common$gene_distance_1[idx]+1) + common$s_het_1[idx])
      res = summary(mod)$coef[2,]
      j = j + 1
      res.lst[[j]] = data.frame(trait=organ,dataset,cell,Estimate=res["Estimate"],se = res["Std. Error"],tvalue=res["t value"],pval=res["Pr(>|t|)"])
    }
  }
  
  res.df = as.data.frame(do.call(rbind,res.lst))
  if (nrow(res.df) > 0) {
    print(res.df[order(res.df$pval,decreasing = F),][1:10,])
  }
  
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_gwas_binarizedpip_enrichments.qtl.",organ,".txt")
  fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
  
}


#########################################

# This section analyzes UKB GWAS results that have been fine-mapped using in-sample UKB LD.

f="/oak/stanford/groups/smontgom/amarder/data/finemap_data/ukbb/release1.1/UKBB_94traits_release1.bed.gz"
finemap.ukb = fread(f,data.table = F,stringsAsFactors = F)
colnamesUse=fread("/oak/stanford/groups/smontgom/amarder/data/finemap_data/ukbb/release1.1/UKBB_94traits_release1.cols",data.table = F,stringsAsFactors = F,header = F)
colnames(finemap.ukb) = colnamesUse[,1]
coords = fread("/oak/stanford/groups/smontgom/amarder/data/finemap_data/ukbb/release1.1/ukb.hg38.bed",data.table = F,stringsAsFactors = F)
colnames(coords) = colnames(finemap.ukb)[1:4]
colnames(coords)[3] = paste0(colnames(coords)[3] ,"_hg38")
finemap.ukb = merge(finemap.ukb,coords[,c(3,4)],by="variant")

trait_lst = c("BMI","Neuroticism","AFib")
for (traitName in trait_lst) {
  print(traitName)
  # traitName = "Alzheimers_Bellenguez_2022"
  finemap.ukb.sub = subset(finemap.ukb,trait==traitName)
  finemap.ukb.sub$chr_pos = paste0(finemap.ukb.sub$chromosome,"_",finemap.ukb.sub$end_hg38)
  
  common$gwas = NA
  fm.sub = subset(finemap.ukb.sub,pip > 0.2)
  idx = common$chr_pos %in% fm.sub$chr_pos
  common$gwas[idx] = TRUE
  fm.sub = subset(finemap.ukb.sub,pip < 0.01)
  idx = common$chr_pos %in% fm.sub$chr_pos
  common$gwas[idx] = FALSE
  
  # contexts
  string = grep("abs_logfc.mean.pval.",colnames(common),value = TRUE)
  dataset_lst <- str_extract(string, "(?<=\\.)[^.]+$")
  cell_lst <- str_extract(string, "(?<=abs_logfc\\.mean\\.pval\\.).*(?=\\.[^.]+$)")
  #
  res.lst=list();j=1;
  for (i in 1:length(cell_lst)) {
    # cell type to test:
    cell=cell_lst[i]
    dataset=dataset_lst[i]
    print(paste0(i,"/",length(cell_lst),": ",cell," (",dataset,")"))
    idx = common[,paste0("peak_overlap.",cell,".",dataset)]
    # test predicting fine map variants
    if (sum(common$gwas[idx],na.rm=T)>1) {
      mod = lm((common[idx,paste0("abs_logfc.mean.",cell,".",dataset)])~common$gwas[idx] + log10(common$gene_distance_1[idx]+1) + common$s_het_1[idx])
      res = summary(mod)$coef[2,]
      j = j + 1
      res.lst[[j]] = data.frame(trait=traitName,dataset,cell,Estimate=res["Estimate"],se = res["Std. Error"],tvalue=res["t value"],pval=res["Pr(>|t|)"])
    }
  }
  
  res.df = as.data.frame(do.call(rbind,res.lst))
  if (nrow(res.df) > 0) {
    print(res.df[order(res.df$pval,decreasing = F),][1:10,])
  }
  
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_gwas_binarizedpip_enrichments.",traitName,".txt")
  fwrite(res.df,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}



