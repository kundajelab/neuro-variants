# conda activate r
library(data.table)
trait_lst = c("Alzheimers_Bellenguez_2022","CAD_Tcheandjieu_2022","CAD_Aragam_2022","BMI","Neuroticism","AFib")
traitName = "Alzheimers_Bellenguez_2022"

j = 0
df.mg.lst = list()
for (traitName in trait_lst) {
  print(traitName)
  f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_gwas_binarizedpip_enrichments.",traitName,".txt")
  df1 = fread(f,data.table = F,stringsAsFactors = F)
  df1$fdr = p.adjust(df1$pval,method = 'fdr')
  
  f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_gwas_binarizedpip_enrichments.peaks.",traitName,".txt")
  df2 = fread(f,data.table = F,stringsAsFactors = F)
  df2$fdr = p.adjust(df2$pval,method = 'fdr')
  
  df.mg = merge(df1,df2,by=c('trait','dataset','cell'))
  j = j + 1; df.mg.lst[[j]] = df.mg
  print(cor.test(df.mg$Estimate.x,df.mg$Estimate.y))
  print(cor.test(-log10(df.mg$pval.x),-log10(df.mg$pval.y)))
  print(table(df.mg$fdr.x < 0.1,df.mg$fdr.y < 0.1))
  subset(df.mg,fdr.y < 0.1)
}

traitName="brain"
trait_lst = c("heart_artery","brain")
for (traitName in trait_lst) {
  print(traitName)
  f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_gwas_binarizedpip_enrichments.qtl.",traitName,".txt")
  df1 = fread(f,data.table = F,stringsAsFactors = F)
  df1$fdr = p.adjust(df1$pval,method = 'fdr')
  
  f = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/ukb_fm_gwas_binarizedpip_enrichments.peaks.qtl.",traitName,".txt")
  df2 = fread(f,data.table = F,stringsAsFactors = F)
  df2$fdr = p.adjust(df2$pval,method = 'fdr')
  
  df.mg = merge(df1,df2,by=c('trait','dataset','cell'))
  j = j + 1; df.mg.lst[[j]] = df.mg
  # print(cor.test(df.mg$Estimate.x,df.mg$Estimate.y))
  # print(cor.test(-log10(df.mg$pval.x),-log10(df.mg$pval.y)))
  # print(table(df.mg$fdr.x < 0.1,df.mg$fdr.y < 0.1))
  # subset(df.mg,fdr.y < 0.1)
}

  
df.mg.all = as.data.frame(do.call(rbind,df.mg.lst))
df.mg.all = subset(df.mg.all,!(trait %in% c("brain","heart_artery")))
cor.test(df.mg.all$Estimate.x,df.mg.all$Estimate.y)
cor.test(-log10(df.mg.all$pval.x),-log10(df.mg.all$pval.y))
table(df.mg.all$fdr.x < 0.1,df.mg.all$fdr.y < 0.1)

df.mg.all = as.data.frame(do.call(rbind,df.mg.lst))
df.mg.all = subset(df.mg.all,(trait %in% c("brain","heart_artery")))
cor.test(df.mg.all$Estimate.x,df.mg.all$Estimate.y)
cor.test(-log10(df.mg.all$pval.x),-log10(df.mg.all$pval.y))
table(df.mg.all$fdr.x < 0.1,df.mg.all$fdr.y < 0.1)

