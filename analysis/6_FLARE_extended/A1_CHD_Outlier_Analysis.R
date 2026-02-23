library(data.table)
f="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/summarize/chd.all_dataset.K562_bias.annot2.txt"
df=fread(f,data.table = F,stringsAsFactors = F)

f="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/variant_lists/chd_snv_case_info.tsv"
info = fread(f,data.table = F,stringsAsFactors = F)

genes = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/gene_lists/chdgene_table.csv",data.table = F,stringsAsFactors = F)
df = merge(df,info,by.x="snp_id",by.y='id')

flare = fread("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/chd.FLARE.2.txt",data.table = F,stringsAsFactors = F)
flare = flare[,c(1,3:(ncol(flare)-1))]
df = merge(df,flare,by='snp_id')

df$case=ifelse(df$case,"case","control")
df$gene_distance_1.log10 = log10(df$gene_distance_1 + 1)

outlier_case_imbalance = function(score_of_interest,tmp,nrank=1000) {
  tmp.sub = tmp[order(tmp[,score_of_interest],decreasing = T),] # sort by score ranking
  rng = c(seq(3,nrank,by=1)) 
  save_y = list()
  for (k in 1:length(rng)) {
    j = rng[k]
    df.mg.sub = tmp.sub[1:j,]
    tab <- table(factor(df.mg.sub$case, levels = c("case","control")))
    save_y[[k]] = tab
  }
  
  save_y.all = as.data.frame(do.call(rbind,save_y))
  save_y.all$n = rng
  save_y.all$prop = save_y.all$case/save_y.all$n
  
  total_prob = mean(tmp.sub$case=="case")
  for (k in 1:length(rng)) {
    res = binom.test(save_y.all$case[k],save_y.all$n[k],total_prob)
    save_y.all$l[k] = res$conf.int[1]
    save_y.all$h[k] = res$conf.int[2]
    save_y.all$pval[k] = res$p.value
  }
  return(save_y.all)
}

genes.sub = genes[genes$`Supporting References`>=3,]
df.sub = subset(df,
                (closest_gene_1 %in% genes.sub$Gene |
                   closest_gene_2 %in% genes.sub$Gene |
                   closest_gene_3 %in% genes.sub$Gene) & is.na(`Known cause of CHD (cases)`)) # ignore GeneName

predictor = "FLARE_heart"
res = outlier_case_imbalance(predictor,df.sub,nrank=100)
print(head(res,30))

exact2x2::exact2x2(df.sub$FLARE_h > quantile(df.sub$FLARE_h,0.99),df.sub$case=="case")
subset(df.sub,FLARE_heart > quantile(df.sub$FLARE_h,0.99) & case=="case")

for (predictor in c("FLARE_fb","FLARE_heart","PHRED","phylop","abs_logfc.mean.fetal_heart.Endocardial_cells.domcke_2020","abs_logfc.mean.fetal_heart.Cardiomyocytes.domcke_2020","gene_distance_1.log10")) {
  print(predictor)
  res = outlier_case_imbalance(predictor,df.sub,nrank=1000)
  # print(head(res,30))
  f.out = paste0("/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/chd_results/",predictor,".extrema.txt")
  fwrite(res,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)
}

f.out = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/pred/results/scores_chd.txt"
fwrite(df.sub[,c("snp_id","case","FLARE_fb","FLARE_heart","phylop","PHRED","abs_logfc.mean.fetal_heart.Endocardial_cells.domcke_2020","abs_logfc.mean.fetal_heart.Cardiomyocytes.domcke_2020","gene_distance_1.log10")],f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)

score_of_interest="FLARE_heart"
tmp = df.sub[order(df.sub[,score_of_interest],decreasing = T),] # sort by score ranking
tmp2 <- head(subset(tmp,case=="case"),9)[,c("snp_id","closest_gene_1","closest_gene_2","closest_gene_3","phylop","gene_distance_1","abs_logfc.mean.pval.fetal_heart.Endocardial_cells.domcke_2020","abs_logfc.mean.pval.fetal_heart.Cardiomyocytes.domcke_2020")]

# head(tmp,12)$case
head(subset(tmp,case=="case"),9)$abs_logfc.mean.fetal_heart.Endocardial_cells.domcke_2020
head(subset(tmp,case=="case"),9)[,c("abs_logfc.mean.fetal_heart.Endocardial_cells.domcke_2020","abs_logfc.mean.fetal_heart.Cardiomyocytes.domcke_2020")]
subset(tmp2,closest_gene_1 %in% genes$Gene)

head(tmp$abs_logfc.mean.fetal_heart.Cardiomyocytes.domcke_2020)


score_of_interest="PHRED"
tmp=df.sub
nrank=30

df.sub = subset(df,
                closest_gene_1 %in% genes$Gene |
                  closest_gene_2 %in% genes$Gene |
                  closest_gene_3 %in% genes$Gene) # ignore GeneName
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.99))
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.99))

df.sub = subset(df,
                closest_gene_1 %in% genes$Gene |
                  closest_gene_2 %in% genes$Gene) # ignore GeneName
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.99))
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.999))

df.sub = subset(df,
                closest_gene_1 %in% genes$Gene) # ignore GeneName
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.99))
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.999))

df.sub[order(df.sub$abs_logfc.mean.pval.fetal_heart.Endocardial_cells.domcke_2020,decreasing = F),][1:10,'case']
df.sub[order(df.sub$abs_logfc.mean.pval.fetal_heart.Endocardial_cells.domcke_2020,decreasing = F),][1:10,'closest_gene_1']

df.sub = df
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.99))
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.999))


df.sub = subset(df,is.na(`Known cause of CHD (cases)`))
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.99))
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.999))

genes.sub = genes[genes$`Supporting References`>=3,]
df.sub = subset(df,
                (closest_gene_1 %in% genes.sub$Gene |
                   closest_gene_2 %in% genes.sub$Gene |
                   closest_gene_3 %in% genes.sub$Gene) & is.na(`Known cause of CHD (cases)`)) # ignore GeneName
# table(df$`Known cause of CHD (cases)`)
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.98))
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.99))
fisher.test(df.sub$case,df.sub$FLARE_fb > quantile(df.sub$FLARE_fb,0.99))
fisher.test(df.sub$case,df.sub$FLARE_heart > quantile(df.sub$FLARE_heart,0.999))
df.sub[order(df.sub$FLARE_heart,decreasing = T),][1:30,'case']
df.sub[order(df.sub$phylop.x,decreasing = T),][1:30,'case']
df.sub[order(df.sub$PHRED,decreasing = T),][1:30,'case']

fisher.test(df.sub$case,df.sub$PHRED > quantile(df.sub$PHRED,0.99))

t.test(df.sub$abs_logfc.mean.fetal_heart.Endocardial_cells.domcke_2020~df.sub$case)
t.test(df.sub$abs_logfc.mean.fetal_heart.Cardiomyocytes.domcke_2020~df.sub$case)
t.test(df.sub$FLARE_heart~df.sub$case)

t.test(df.sub$abs_logfc.mean.fetal_heart.Cardiomyocytes.domcke_2020~df.sub$case)
fisher.test(df.sub$abs_logfc.mean.pval.fetal_heart.Cardiomyocytes.domcke_2020<0.01,df.sub$case)
fisher.test(df.sub$peak_overlap.fetal_heart.Cardiomyocytes.domcke_2020,df.sub$case)
fisher.test(df.sub$abs_logfc.mean.pval.fetal_heart.Vascular_endothelial_cells.domcke_2020<0.005,df.sub$case)


t.test(df.sub$abs_logfc.mean.fetal_heart.Endocardial_cells.domcke_2020~df.sub$case)
t.test(df.sub$abs_logfc.mean.fetal_heart.Cardiomyocytes.domcke_2020~df.sub$case)
fisher.test(df.sub$abs_logfc.mean.pval.fetal_heart.Cardiomyocytes.domcke_2020<0.05,df.sub$case)
fisher.test(df.sub$abs_logfc.mean.pval.fetal_heart.Endocardial_cells.domcke_2020<0.1,df.sub$case)
fisher.test(df.sub$abs_logfc.mean.pval.fetal_heart.Vascular_endothelial_cells.domcke_2020<0.005,df.sub$case)




