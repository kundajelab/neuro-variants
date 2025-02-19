# Use 1_FLARE_Predict.R to get 'df' for ASD variants as input to this script.
# can run this script locally...

# ASD mutation info
asd_info = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/Trevino_et_al_info.txt",data.table = F,stringsAsFactors = F)
asd_info$GeneSymbol = "."
asd_info$GeneSymbol[asd_info$SYMBOL!="."] = asd_info$SYMBOL[asd_info$SYMBOL!="."]
asd_info$GeneSymbol[asd_info$NEAREST!="."] = asd_info$NEAREST[asd_info$NEAREST!="."]
asd_info$Pheno = factor(asd_info$Pheno,c("control","case"))
df.mg = merge(df,asd_info[,c("snp_id","Pheno","GeneSymbol","Gene","SampleID")],by="snp_id")

# SFARI gene info
sfari = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/gene_lists/SFARI-Gene_genes_08-19-2024release_08-22-2024export.csv",data.table = F,stringsAsFactors = F)
sfari.sub = subset(sfari,syndromic==1 & (`number-of-reports` >= 10))
# sfari.sub = subset(sfari,syndromic==1 & `number-of-reports` >= 10 & `gene-score` %in% c(1,2))
# sfari gene needs to be closest TSS or closest gene body
df.mg$sfari = df.mg$closest_gene_1 %in% sfari.sub[,"gene-symbol"] | #|
  df.mg$GeneSymbol %in% sfari.sub[,"gene-symbol"]

# subset to near SFARI genes
tmp = subset(df.mg,sfari)

# compare predictive value
t.test(FLARE_heart~Pheno,tmp)
t.test(FLARE_fb~Pheno,tmp)
wilcox.test(FLARE_fb~Pheno,tmp) # shift in the mean but not the median
fisher.test(tmp$FLARE_fb > quantile(tmp$FLARE_fb,0.99),tmp$Pheno)

# other metrics
t.test(phylop~Pheno,tmp)
t.test(PHRED~Pheno,tmp)
t.test(abs_logfc.mean.c11.trevino_2021~Pheno,tmp)
fisher.test(tmp$peak_overlap.c11.trevino_2021,tmp$Pheno)
fisher.test(tmp$peak_overlap.c8.trevino_2021,tmp$Pheno)
tmp.sub = subset(tmp,peak_overlap.c8.trevino_2021)
t.test(abs_logfc.mean.c8.trevino_2021~Pheno,tmp.sub)

#########

outlier_case_imbalance = function(score_of_interest) {
  tmp.sub = tmp[order(tmp[,score_of_interest],decreasing = T),] # sort by score ranking
  rng = c(seq(3,1000,by=1)) 
  save_y = list()
  for (k in 1:length(rng)) {
    j = rng[k]
    df.mg.sub = tmp.sub[1:j,]
    tab = (table(df.mg.sub$Pheno))
    save_y[[j]] = tab
  }
  
  save_y.all = as.data.frame(do.call(rbind,save_y))
  save_y.all$n = rng
  save_y.all$prop = save_y.all$case/save_y.all$n
  
  total_prob = mean(tmp.sub$Pheno=="case")
  for (k in 1:length(rng)) {
    res = binom.test(save_y.all$case[k],save_y.all$n[k],total_prob)
    save_y.all$l[k] = res$conf.int[1]
    save_y.all$h[k] = res$conf.int[2]
    save_y.all$pval[k] = res$p.value
  }
  return(save_y.all)
}

predictor = "FLARE_fb"
res = outlier_case_imbalance(predictor)
head(res,40)

predictor = "FLARE_fb_peaks"
res = outlier_case_imbalance(predictor)
head(res,40)

# compare predictive value
t.test(FLARE_fb~Pheno,tmp)
t.test(FLARE_fb_peaks~Pheno,tmp)

for (predictor in c("FLARE_fb","FLARE_fb_peaks","FLARE_heart","FLARE_ab","PHRED","phylop","abs_logfc.mean.c11.trevino_2021","s_het_1")) {
  res = outlier_case_imbalance(predictor)
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/",variantSet,".",predictor,".extrema.txt")
  fwrite(res,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}

predictor = "FLARE_fb"
tmp.sub = tmp[order(tmp[,predictor],decreasing = T),] # sort by score ranking
tmp2 = subset(tmp.sub,GeneSymbol=="CNTNAP2")[,c("snp_id","SampleID","Pheno","closest_gene_1","GeneSymbol","Consequence",predictor)]
cat("Identified ",nrow(tmp)," variants in CNTNAP2.\n",sep = '')
wilcox.test(tmp2[,predictor] ~ tmp2$Pheno)

f.out = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/scores_sfari.txt"
fwrite(tmp[,c("snp_id","Pheno","FLARE_fb","FLARE_heart","phylop","PHRED","abs_logfc.mean.c11.trevino_2021")],f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
tmp2 = subset(tmp.sub,GeneSymbol=="CNTNAP2")
f.out = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/scores_cntnap2.txt"
fwrite(tmp2[,c("snp_id","Pheno","FLARE_fb","FLARE_heart","phylop","PHRED","abs_logfc.mean.c11.trevino_2021")],f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

# # exploratory commands below:
# subset(sfari,sfari[,2]=="CNTNAP2")
# subset(sfari,sfari[,2]=="NFIB")
# fisher.test(tmp$pred > quantile(tmp$pred,0.99),tmp$Pheno)
# tmp.sub[1:50,c("snp_id","Pheno","closest_gene_1","GeneSymbol","Consequence","pred")]
# df.mg$gene_interest = df.mg$closest_gene_1
# idx = grepl("intron",df.mg$Consequence); df.mg$gene_interest[idx] = df.mg$GeneSymbol[idx]
# subset(df.mg,snp_id=="chr14:27886425:C:T")[,c("GeneSymbol","Gene","GeneName","closest_gene_1","closest_gene_2","closest_gene_3")]
# subset(df.mg,snp_id=="chr7:148230432:G:T")[,c("GeneSymbol","Gene","GeneName","closest_gene_1","closest_gene_2","closest_gene_3","sfari")]

# f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/",i,".",model,".",variantSet,".",modelType,".extrema.txt")
# fwrite(save_y.all,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
