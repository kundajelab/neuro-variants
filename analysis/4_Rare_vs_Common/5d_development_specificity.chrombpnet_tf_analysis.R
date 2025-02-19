library(data.table)

f = "/oak/stanford/groups/akundaje/projects/neuro-variants/variant_motif_matrices/rare/K562_bias/alpha_0.8/domcke_2020/domcke_2020.fetal_brain.Excitatory_neurons.tsv"
fb = fread(f,data.table = F,stringsAsFactors = F)
f = "/oak/stanford/groups/akundaje/projects/neuro-variants/variant_motif_matrices/rare/K562_bias/alpha_0.8/corces_2020/corces_2020.Cluster1.tsv"
ab = fread(f,data.table = F,stringsAsFactors = F)

f="/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/cbp/analysis/excitatory_neurons_dev.dataset.txt"
df2 = fread(f,data.table = F,stringsAsFactors = F)

df2 = df2[match(ab$V1,df2$snp_id),]

rownames(ab)=ab[,1]
ab = ab[,-1]
rownames(fb)=fb[,1]
fb = fb[,-1]

df2[(ab[,"NEUROD1-ATOH_1"]==1) & (fb[,"NEUROD1-ATOH_1"]==1) & df2$fetal_specific & df2$abs_logfc.mean.fetal_brain.Excitatory_neurons.domcke_2020 > 0.65,]
df2[(ab[,"NEUROD1-ATOH_1"]==1) & (fb[,"NEUROD1-ATOH_1"]==1) & df2$fetal_specific & df2$abs_logfc.mean.Cluster1.corces_2020 > 0.12,]

subset(df2,fetal_specific & abs_logfc.mean.fetal_brain.Excitatory_neurons.domcke_2020 > 0.4 & abs_logfc.mean.Cluster1.corces_2020 > 0.14)
tmp = fb["chr7:131508050:C:T",]; colnames(tmp)[tmp==1]
tmp = fb["chr2:15944130:T:C",]; colnames(tmp)[tmp==1]
tmp = ab["chr2:15944130:T:C",]; colnames(tmp)[tmp==1]

tmp = fb["chr2:191020486:C:T",]; colnames(tmp)[tmp==1]
tmp = fb["chr10:27242206:C:T",]; colnames(tmp)[tmp==1]
tmp = fb["chr12:2255516:C:G",]; colnames(tmp)[tmp==1]
tmp = fb["chr14:24216122:A:T",]; colnames(tmp)[tmp==1]


res = list()
for (tfNum_use in 1:ncol(ab)) {
  tf_adult_specific_ab = ab[df2$adult_specific,tfNum_use]
  tf_fetal_specific_ab = ab[df2$fetal_specific,tfNum_use]
  tf_adult_specific_fb = fb[df2$adult_specific,tfNum_use]
  tf_fetal_specific_fb = fb[df2$fetal_specific,tfNum_use]
  res[[tfNum_use]] = data.frame(tf = colnames(ab)[tfNum_use],
                                tf_adult_specific_ab = sum(tf_adult_specific_ab),
                                tf_fetal_specific_ab = sum(tf_fetal_specific_ab),
                                tf_adult_specific_fb = sum(tf_adult_specific_fb),
                                tf_fetal_specific_fb = sum(tf_fetal_specific_fb))
}
res.df = as.data.frame(do.call(rbind,res))
sum(df2$adult_specific)
sum(df2$fetal_specific)
res.df[order(res.df$tf_adult_specific_ab,decreasing = T),][1:10,]
res.df[order(res.df$tf_adult_specific_fb,decreasing = T),][1:10,]
res.df[order(res.df$tf_fetal_specific_ab,decreasing = T),][1:10,]
res.df[order(res.df$tf_fetal_specific_fb,decreasing = T),][1:10,]

res.df$tf_adult_specific_ab = res.df$tf_adult_specific_ab / sum(df2$adult_specific)
res.df$tf_adult_specific_fb = res.df$tf_adult_specific_fb / sum(df2$adult_specific)
res.df$tf_fetal_specific_ab = res.df$tf_fetal_specific_ab / sum(df2$fetal_specific)
res.df$tf_fetal_specific_fb = res.df$tf_fetal_specific_fb / sum(df2$fetal_specific)

tmp = subset(res.df,tf_adult_specific_ab > 0.05 | tf_fetal_specific_fb > 0.05)
f.out = "/oak/stanford/groups/smontgom/amarder/chrombpnet_variant_effects/output/data/FinalAnalysis/dev_specific_chrombpnet_motif_prop.txt"
fwrite(tmp,f.out,quote = F,na = "NA",sep = '\t',row.names = F,col.names = T)






